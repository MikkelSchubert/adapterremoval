// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2021 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "demultiplexing.hpp"
#include "barcode_table.hpp" // for barcode_table
#include "debug.hpp"         // for AR_REQUIRE, AR_REQUIRE_SINGLE_THREAD
#include "fastq_io.hpp"      // for chunk_ptr, fastq...
#include "output.hpp"        // for output_files
#include "sequence_sets.hpp" // for adapter_set
#include "userconfig.hpp"    // for userconfig, ar_command, ar_command::demul...
#include <cstddef>           // for size_t
#include <memory>            // for make_unique, unique_ptr
#include <utility>           // for move

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Implementations for `demultiplex_reads`

demultiplex_reads::demultiplex_reads(const userconfig& config,
                                     const post_demux_steps& steps,
                                     demux_stats_ptr stats)
  : analytical_step(processing_order::ordered, "demultiplex_reads")
  , m_samples(config.samples)
  , m_barcode_table(m_samples,
                    config.barcode_mm,
                    config.barcode_mm_r1,
                    config.barcode_mm_r2)
  , m_config(config)
  , m_steps(steps)
  , m_cache(steps)
  , m_statistics(std::move(stats))
{
  AR_REQUIRE(m_samples.size());
  AR_REQUIRE(m_samples.size() == m_steps.samples.size());
  AR_REQUIRE(m_statistics);

  AR_REQUIRE(m_statistics->samples.empty());
  m_statistics->samples.resize(m_samples.size());

  // Map global barcode offsets to sample and relative barcode offsets
  for (size_t i = 0; i < m_samples.size(); ++i) {
    const size_t barcodes = m_samples.at(i).size();

    m_statistics->samples.at(i).resize_up_to(barcodes);
    for (size_t j = 0; j < barcodes; ++j) {
      m_barcodes.emplace_back(i, j);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////

demultiplex_se_reads::demultiplex_se_reads(const userconfig& config,
                                           const post_demux_steps& steps,
                                           demux_stats_ptr stats)
  : demultiplex_reads(config, steps, std::move(stats))
{
}

chunk_vec
demultiplex_se_reads::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  for (auto& read : chunk->reads_1) {
    const auto [sample, barcode] = m_barcode_table.identify(read);
    // TODO: We should keep reads even if we cannot identify the exact barcode
    if (sample < 0 || barcode < 0) {
      switch (sample) {
        case barcode_key::unidentified:
          m_statistics->unidentified += 1;
          break;
        case barcode_key::ambiguous:
          m_statistics->ambiguous += 1;
          break;
        default:
          AR_FAIL("invalid barcode match sample");
      }

      m_cache.add_unidentified_1(std::move(read));
    } else {
      read.truncate(m_barcode_table.length_1());

      m_statistics->samples.at(sample).inc(barcode);
      m_cache.add_read_1(std::move(read), sample, barcode);
    }
  }

  return m_cache.flush(chunk->eof, chunk->mate_separator);
}

///////////////////////////////////////////////////////////////////////////////

demultiplex_pe_reads::demultiplex_pe_reads(const userconfig& config,
                                           const post_demux_steps& steps,
                                           demux_stats_ptr stats)
  : demultiplex_reads(config, steps, std::move(stats))
{
}

chunk_vec
demultiplex_pe_reads::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size());

  auto it_1 = chunk->reads_1.begin();
  auto it_2 = chunk->reads_2.begin();
  for (; it_1 != chunk->reads_1.end(); ++it_1, ++it_2) {
    const auto [sample, barcode] = m_barcode_table.identify(*it_1, *it_2);

    // TODO: We should keep reads even if we cannot identify the exact barcode
    if (sample < 0 || barcode < 0) {
      switch (sample) {
        case barcode_key::unidentified:
          m_statistics->unidentified += 2;
          break;
        case barcode_key::ambiguous:
          m_statistics->ambiguous += 2;
          break;
        default:
          AR_FAIL("invalid barcode match sample");
      }

      m_cache.add_unidentified_1(std::move(*it_1));
      m_cache.add_unidentified_2(std::move(*it_2));
    } else {
      it_1->truncate(m_barcode_table.length_1());
      m_cache.add_read_1(std::move(*it_1), sample, barcode);

      it_2->truncate(m_barcode_table.length_2());
      m_cache.add_read_2(std::move(*it_2), sample);

      m_statistics->samples.at(sample).inc(barcode, 2);
    }
  }

  return m_cache.flush(chunk->eof, chunk->mate_separator);
}

///////////////////////////////////////////////////////////////////////////////

process_demultiplexed::process_demultiplexed(const userconfig& config,
                                             const sample_output_files& output,
                                             const size_t sample,
                                             trim_stats_ptr sink)
  : analytical_step(processing_order::unordered, "process_demultiplexed")
  , m_config(config)
  , m_output(output)
  , m_samples(config.samples)
  , m_sample(sample)
  , m_stats_sink(std::move(sink))
{
  m_stats.emplace_back_n(m_config.max_threads, m_config.report_sample_rate);
}

chunk_vec
process_demultiplexed::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  processed_reads chunks{ m_output };
  chunks.set_sample(m_samples.at(m_sample));
  chunks.set_mate_separator(chunk->mate_separator);
  chunks.set_demultiplexing_only(true);

  if (chunk->first) {
    chunks.write_headers(m_config.args);
  }

  auto stats = m_stats.acquire();

  AR_REQUIRE(chunk->reads_1.size() == chunk->barcodes.size());
  auto barcode = chunk->barcodes.begin();

  if (m_config.paired_ended_mode) {
    AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size());

    auto it_1 = chunk->reads_1.begin();
    auto it_2 = chunk->reads_2.begin();
    for (; it_1 != chunk->reads_1.end(); ++it_1, ++it_2, ++barcode) {
      it_1->add_prefix_to_name(m_config.prefix_read_1);
      stats->read_1->process(*it_1);
      chunks.add(*it_1, read_type::pe_1, *barcode);

      it_2->add_prefix_to_name(m_config.prefix_read_2);
      stats->read_2->process(*it_2);
      chunks.add(*it_2, read_type::pe_2, *barcode);
    }
  } else {
    for (auto& read : chunk->reads_1) {
      read.add_prefix_to_name(m_config.prefix_read_1);

      stats->read_1->process(read);
      chunks.add(read, read_type::se, *barcode++);
    }
  }

  m_stats.release(stats);

  return chunks.finalize(chunk->eof);
}

void
process_demultiplexed::finalize()
{
  m_stats.merge_into(*m_stats_sink);
}

///////////////////////////////////////////////////////////////////////////////

processes_unidentified::processes_unidentified(const userconfig& config,
                                               const output_files& output,
                                               demux_stats_ptr stats)
  : analytical_step(processing_order::unordered, "processes_unidentified")
  , m_config(config)
  , m_statistics(std::move(stats))
{
  m_output.set_file(read_file::mate_1, output.unidentified_1);
  m_output.set_file(read_file::mate_2, output.unidentified_2);

  if (output.unidentified_1_step != output_files::disabled) {
    m_output.set_step(read_file::mate_1, output.unidentified_1_step);
  }

  if (output.unidentified_1_step != output.unidentified_2_step &&
      output.unidentified_2_step != output_files::disabled) {
    m_output.set_step(read_file::mate_2, output.unidentified_2_step);
  }

  m_stats_1.emplace_back_n(m_config.max_threads, config.report_sample_rate);
  m_stats_2.emplace_back_n(m_config.max_threads, config.report_sample_rate);
}

chunk_vec
processes_unidentified::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  processed_reads chunks{ m_output };
  chunks.set_sample(m_config.samples.unidentified());
  chunks.set_mate_separator(chunk->mate_separator);

  if (chunk->first) {
    chunks.write_headers(m_config.args);
  }

  auto stats_1 = m_stats_1.acquire();
  auto stats_2 = m_stats_2.acquire();

  if (m_config.paired_ended_mode) {
    AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size());

    auto it_1 = chunk->reads_1.begin();
    auto it_2 = chunk->reads_2.begin();

    for (; it_1 != chunk->reads_1.end(); ++it_1, ++it_2) {
      it_1->add_prefix_to_name(m_config.prefix_read_1);
      stats_1->process(*it_1);
      chunks.add(*it_1, read_type::pe_1);

      it_2->add_prefix_to_name(m_config.prefix_read_2);
      stats_2->process(*it_2);
      chunks.add(*it_2, read_type::pe_2);
    }
  } else {
    for (auto& read : chunk->reads_1) {
      read.add_prefix_to_name(m_config.prefix_read_1);

      stats_1->process(read);
      chunks.add(read, read_type::se);
    }
  }

  m_stats_1.release(stats_1);
  m_stats_2.release(stats_2);

  return chunks.finalize(chunk->eof);
}

void
processes_unidentified::finalize()
{
  m_stats_1.merge_into(*m_statistics->unidentified_stats_1);
  m_stats_2.merge_into(*m_statistics->unidentified_stats_2);
}

} // namespace adapterremoval
