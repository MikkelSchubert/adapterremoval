/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
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
#include "simd.hpp"       // for size_t
#include "userconfig.hpp" // for userconfig, ar_command, ar_command::demul...
#include <cstddef>        // for size_t
#include <memory>         // for make_unique, unique_ptr
#include <utility>        // for move

namespace adapterremoval {

namespace {

template<typename T>
void
flush_chunk(chunk_vec& output, std::unique_ptr<T>& ptr, size_t step, bool eof)
{
  if (eof || ptr->nucleotides >= INPUT_BLOCK_SIZE) {
    ptr->eof = eof;
    output.push_back(chunk_pair(step, std::move(ptr)));
    ptr = std::make_unique<T>();
  }
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for `post_demux_steps`

const size_t post_demux_steps::disabled = output_files::disabled;

post_demux_steps::post_demux_steps(const output_files& output)
  : unidentified_1(output.unidentified_1_step)
  , unidentified_2(output.unidentified_2_step)
{
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `demultiplex_reads`

demultiplex_reads::demultiplex_reads(const userconfig& config,
                                     post_demux_steps steps,
                                     demux_stats_ptr stats)
  : analytical_step(processing_order::ordered, "demultiplex_reads")
  , m_barcodes(config.adapters.get_barcodes())
  , m_barcode_table(m_barcodes,
                    config.barcode_mm,
                    config.barcode_mm_r1,
                    config.barcode_mm_r2)
  , m_config(config)
  , m_steps(std::move(steps))
  , m_statistics(std::move(stats))
{
  AR_REQUIRE(!m_barcodes.empty());
  AR_REQUIRE(m_barcodes.size() == m_steps.samples.size());
  AR_REQUIRE(m_statistics);
  AR_REQUIRE(m_statistics->barcodes.size() == m_barcodes.size());

  if (m_steps.unidentified_1 != post_demux_steps::disabled) {
    m_unidentified_1 = std::make_unique<analytical_chunk>();
  }

  if (m_steps.unidentified_2 != post_demux_steps::disabled) {
    m_unidentified_2 = std::make_unique<analytical_chunk>();
  }

  for (const auto next_step : m_steps.samples) {
    AR_REQUIRE(next_step != post_demux_steps::disabled);

    m_cache.push_back(std::make_unique<analytical_chunk>());
  }
}

chunk_vec
demultiplex_reads::flush_cache(bool eof)
{
  chunk_vec output;

  if (m_unidentified_1) {
    flush_chunk(output, m_unidentified_1, m_steps.unidentified_1, eof);
  }

  if (m_unidentified_2) {
    flush_chunk(output, m_unidentified_2, m_steps.unidentified_2, eof);
  }

  for (size_t nth = 0; nth < m_cache.size(); ++nth) {
    flush_chunk(output, m_cache.at(nth), m_steps.samples.at(nth), eof);
  }

  return output;
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

      if (m_unidentified_1) {
        m_unidentified_1->add(read);
      }

      if (best_barcode == -1) {
        m_statistics->unidentified += 1;
      } else {
        m_statistics->ambiguous += 1;
      }

      m_statistics->unidentified_stats_1->process(read);
    } else {
      // Prefixing with user supplied prefixes is also done during trimming
      if (m_config.run_type == ar_command::demultiplex_only) {
        read.add_prefix_to_name(m_config.prefix_read_1);
      }

      auto& dst = *m_cache.at(best_barcode);
      read.truncate(m_barcodes.at(best_barcode).first.length());
      dst.nucleotides += read.length();
      dst.reads_1.push_back(read);

      m_statistics->barcodes.at(best_barcode) += 1;
    }
  }

  return flush_cache(chunk->eof);
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

      if (m_unidentified_1) {
        m_unidentified_1->add(*it_1);
      }

      if (m_config.interleaved_output) {
        if (m_unidentified_1) {
          m_unidentified_1->add(*it_2);
        }
      } else if (m_unidentified_2) {
        m_unidentified_2->add(*it_2);
      }

      if (best_barcode == -1) {
        m_statistics->unidentified += 2;
      } else {
        m_statistics->ambiguous += 2;
      }

      m_statistics->unidentified_stats_1->process(*it_1);
      m_statistics->unidentified_stats_2->process(*it_2);
    } else {
      // Prefixing with user supplied prefixes is also done during trimming
      if (m_config.run_type == ar_command::demultiplex_only) {
        it_1->add_prefix_to_name(m_config.prefix_read_1);
        it_2->add_prefix_to_name(m_config.prefix_read_2);
      }

      const chunk_ptr& dst = m_cache.at(best_barcode);

      it_1->truncate(m_barcodes.at(best_barcode).first.length());
      dst->nucleotides += it_1->length();
      dst->reads_1.push_back(*it_1);
      it_2->truncate(m_barcodes.at(best_barcode).second.length());
      dst->nucleotides += it_2->length();
      dst->reads_2.push_back(*it_2);

      m_statistics->barcodes.at(best_barcode) += 2;
    }
  }

  return flush_cache(chunk->eof);
}

} // namespace adapterremoval
