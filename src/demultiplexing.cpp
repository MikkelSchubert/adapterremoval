/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2021 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
\*************************************************************************/
#include "demultiplexing.hpp"
#include "adapterset.hpp" // for adapter_set
#include "debug.hpp"      // for AR_REQUIRE, AR_REQUIRE_SINGLE_THREAD
#include "fastq_io.hpp"   // for chunk_ptr, fastq...
#include "output.hpp"     // for output_files
#include "serializer.hpp" // for fastq_flags
#include "userconfig.hpp" // for userconfig, ar_command, ar_command::demul...
#include <cstddef>        // for size_t
#include <memory>         // for make_unique, unique_ptr
#include <utility>        // for move

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Implementations for `demultiplex_reads`

demultiplex_reads::demultiplex_reads(const userconfig& config,
                                     const post_demux_steps& steps,
                                     demux_stats_ptr stats)
  : analytical_step(processing_order::ordered, "demultiplex_reads")
  , m_barcodes(config.adapters.get_barcodes())
  , m_barcode_table(m_barcodes,
                    config.barcode_mm,
                    config.barcode_mm_r1,
                    config.barcode_mm_r2)
  , m_config(config)
  , m_steps(steps)
  , m_cache(steps)
  , m_statistics(std::move(stats))
{
  AR_REQUIRE(!m_barcodes.empty());
  AR_REQUIRE(m_barcodes.size() == m_steps.samples.size());
  AR_REQUIRE(m_statistics);
  AR_REQUIRE(m_statistics->barcodes.size() == m_barcodes.size());
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
    const int best_barcode = m_barcode_table.identify(read);

    if (best_barcode < 0) {
      // Always add prefix since unidentified reads are not processed further
      read.add_prefix_to_name(m_config.prefix_read_1);

      if (best_barcode == -1) {
        m_statistics->unidentified += 1;
      } else {
        m_statistics->ambiguous += 1;
      }

      m_cache.add_unidentified_1(std::move(read));
    } else {
      read.truncate(m_barcodes.at(best_barcode).first.length());

      // Prefixing with user supplied prefixes is also done during trimming
      if (m_config.run_type == ar_command::demultiplex_only) {
        read.add_prefix_to_name(m_config.prefix_read_1);
      }

      m_statistics->barcodes.at(best_barcode) += 1;
      m_cache.add_read_1(std::move(read), best_barcode);
    }
  }

  return m_cache.flush(chunk->eof);
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
    const int best_barcode = m_barcode_table.identify(*it_1, *it_2);

    if (best_barcode < 0) {
      // Always add prefix since unidentified reads are not processed further
      it_1->add_prefix_to_name(m_config.prefix_read_1);
      it_2->add_prefix_to_name(m_config.prefix_read_2);

      if (best_barcode == -1) {
        m_statistics->unidentified += 2;
      } else {
        m_statistics->ambiguous += 2;
      }

      m_cache.add_unidentified_1(std::move(*it_1));
      m_cache.add_unidentified_2(std::move(*it_2));
    } else {
      // Prefixing with user supplied prefixes is also done during trimming
      if (m_config.run_type == ar_command::demultiplex_only) {
        it_1->add_prefix_to_name(m_config.prefix_read_1);
        it_2->add_prefix_to_name(m_config.prefix_read_2);
      }

      it_1->truncate(m_barcodes.at(best_barcode).first.length());
      m_cache.add_read_1(std::move(*it_1), best_barcode);

      it_2->truncate(m_barcodes.at(best_barcode).second.length());
      m_cache.add_read_2(std::move(*it_2), best_barcode);

      m_statistics->barcodes.at(best_barcode) += 2;
    }
  }

  return m_cache.flush(chunk->eof);
}

processes_unidentified::processes_unidentified(const userconfig& config,
                                               const output_files& output,
                                               demux_stats_ptr stats)
  : analytical_step(processing_order::unordered, "processes_unidentified")
  , m_config(config)
  , m_statistics(std::move(stats))
{
  m_output.set_file(read_type::mate_1, output.unidentified_1);
  m_output.set_file(read_type::mate_2, output.unidentified_2);

  m_output.set_step(read_type::mate_1, output.unidentified_1_step);
  if (output.unidentified_1_step != output.unidentified_2_step &&
      output.unidentified_2_step != output_files::disabled) {
    m_output.set_step(read_type::mate_2, output.unidentified_2_step);
  }

  m_stats_1.emplace_back_n(m_config.max_threads);
  m_stats_2.emplace_back_n(m_config.max_threads);
}

chunk_vec
processes_unidentified::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  processed_reads chunks{ m_output, chunk->first };

  auto stats_1 = m_stats_1.acquire();
  auto stats_2 = m_stats_2.acquire();

  if (chunk->reads_2.empty()) {
    for (const auto& read : chunk->reads_1) {
      stats_1->process(read);
      chunks.add(read, read_type::mate_1, fastq_flags::se);
    }
  } else {
    AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size());
    auto it_1 = chunk->reads_1.begin();
    auto it_2 = chunk->reads_2.begin();

    for (; it_1 != chunk->reads_1.end(); ++it_1, ++it_2) {
      stats_1->process(*it_1);
      chunks.add(*it_1, read_type::mate_1, fastq_flags::pe_1);

      stats_2->process(*it_2);
      chunks.add(*it_2, read_type::mate_2, fastq_flags::pe_2);
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
