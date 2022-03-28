/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <memory>  // for unique_ptr
#include <utility> // for move

#include "adapterset.hpp"  // for adapter_set
#include "commontypes.hpp" // for fastq_vec
#include "debug.hpp"       // for AR_DEBUG_ASSERT, AR_DEBUG_LOCK
#include "demultiplexing.hpp"
#include "fastq_io.hpp"   // for fastq_read_chunk, fastq_output_chunk, rea...
#include "userconfig.hpp" // for userconfig, fastq_encoding_ptr

///////////////////////////////////////////////////////////////////////////////

template<typename T>
void
flush_chunk(chunk_vec& output, std::unique_ptr<T>& ptr, size_t step, bool eof)
{
  if (eof || ptr->nucleotides >= INPUT_BLOCK_SIZE) {
    ptr->eof = eof;
    output.push_back(chunk_pair(step, std::move(ptr)));
    ptr.reset(new T());
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `post_demux_steps`

post_demux_steps::post_demux_steps()
  : unidentified_1(post_demux_steps::disabled)
  , unidentified_2(post_demux_steps::disabled)
  , samples()
{}

const size_t post_demux_steps::disabled = static_cast<size_t>(-1);

///////////////////////////////////////////////////////////////////////////////
// Implementations for `demultiplex_reads`

demultiplex_reads::demultiplex_reads(const userconfig& config,
                                     const post_demux_steps& steps,
                                     demux_stats_ptr stats)
  : analytical_step(processing_order::ordered)
  , m_barcodes(config.adapters.get_barcodes())
  , m_barcode_table(m_barcodes,
                    config.barcode_mm,
                    config.barcode_mm_r1,
                    config.barcode_mm_r2)
  , m_config(config)
  , m_cache()
  , m_unidentified_1()
  , m_unidentified_2()
  , m_steps(steps)
  , m_statistics(stats)
  , m_lock()
{
  AR_DEBUG_ASSERT(!m_barcodes.empty());
  AR_DEBUG_ASSERT(m_barcodes.size() == m_steps.samples.size());
  AR_DEBUG_ASSERT(m_statistics);
  AR_DEBUG_ASSERT(m_statistics->barcodes.size() == m_barcodes.size());

  AR_DEBUG_ASSERT(m_steps.unidentified_1 != post_demux_steps::disabled);
  m_unidentified_1.reset(new fastq_output_chunk());

  if (m_steps.unidentified_2 != post_demux_steps::disabled &&
      m_steps.unidentified_1 != m_steps.unidentified_2) {
    m_unidentified_2.reset(new fastq_output_chunk());
  }

  for (const auto next_step : m_steps.samples) {
    AR_DEBUG_ASSERT(next_step != post_demux_steps::disabled);

    m_cache.push_back(read_chunk_ptr(new fastq_read_chunk()));
  }
}

demultiplex_reads::~demultiplex_reads() {}

chunk_vec
demultiplex_reads::flush_cache(bool eof)
{
  chunk_vec output;

  flush_chunk(output, m_unidentified_1, m_steps.unidentified_1, eof);

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
  : demultiplex_reads(config, steps, stats)
{}

chunk_vec
demultiplex_se_reads::process(chunk_ptr chunk)
{
  AR_DEBUG_LOCK(m_lock);
  auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);

  for (auto& read : read_chunk.reads_1) {
    const int best_barcode = m_barcode_table.identify(read);

    if (best_barcode < 0) {
      m_unidentified_1->add(read);

      if (best_barcode == -1) {
        m_statistics->unidentified += 1;
      } else {
        m_statistics->ambiguous += 1;
      }

      m_statistics->unidentified_stats_1->process(read);
    } else {
      read_chunk_ptr& dst = m_cache.at(best_barcode);
      read.truncate(m_barcodes.at(best_barcode).first.length());
      dst->nucleotides += read.length();
      dst->reads_1.push_back(read);

      m_statistics->barcodes.at(best_barcode) += 1;
    }
  }

  return flush_cache(read_chunk.eof);
}

///////////////////////////////////////////////////////////////////////////////

demultiplex_pe_reads::demultiplex_pe_reads(const userconfig& config,
                                           const post_demux_steps& steps,
                                           demux_stats_ptr stats)
  : demultiplex_reads(config, steps, stats)
{}

chunk_vec
demultiplex_pe_reads::process(chunk_ptr chunk)
{
  AR_DEBUG_LOCK(m_lock);
  auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  AR_DEBUG_ASSERT(read_chunk.reads_1.size() == read_chunk.reads_2.size());

  fastq_vec::iterator it_1 = read_chunk.reads_1.begin();
  fastq_vec::iterator it_2 = read_chunk.reads_2.begin();
  for (; it_1 != read_chunk.reads_1.end(); ++it_1, ++it_2) {
    const int best_barcode = m_barcode_table.identify(*it_1, *it_2);

    if (best_barcode < 0) {
      m_unidentified_1->add(*it_1);
      if (m_config.interleaved_output) {
        m_unidentified_1->add(*it_2);
      } else {
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
      read_chunk_ptr& dst = m_cache.at(best_barcode);

      it_1->truncate(m_barcodes.at(best_barcode).first.length());
      dst->nucleotides += it_1->length();
      dst->reads_1.push_back(*it_1);
      it_2->truncate(m_barcodes.at(best_barcode).second.length());
      dst->nucleotides += it_2->length();
      dst->reads_2.push_back(*it_2);

      m_statistics->barcodes.at(best_barcode) += 2;
    }
  }

  return flush_cache(read_chunk.eof);
}
