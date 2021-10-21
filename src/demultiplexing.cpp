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
#include <algorithm>
#include <iostream>

#include "commontypes.hpp"
#include "debug.hpp"
#include "demultiplexing.hpp"
#include "fastq_io.hpp"
#include "strutils.hpp"
#include "userconfig.hpp"

///////////////////////////////////////////////////////////////////////////////

demultiplex_reads::demultiplex_reads(const userconfig* config,
                                     demultiplexing_statistics* statistics)
  : analytical_step(analytical_step::ordering::ordered)
  , m_barcodes(config->adapters.get_barcodes())
  , m_barcode_table(m_barcodes,
                    config->barcode_mm,
                    config->barcode_mm_r1,
                    config->barcode_mm_r2)
  , m_config(config)
  , m_cache()
  , m_unidentified_1(new fastq_output_chunk())
  , m_unidentified_2()
  , m_statistics(statistics)
  , m_lock()
{
  AR_DEBUG_ASSERT(!m_barcodes.empty());
  AR_DEBUG_ASSERT(m_statistics);
  AR_DEBUG_ASSERT(m_statistics->empty());
  m_statistics->resize(m_barcodes.size());

  if (!config->interleaved_output) {
    m_unidentified_2.reset(new fastq_output_chunk());
  }

  for (size_t i = 0; i < m_barcodes.size(); ++i) {
    m_cache.push_back(read_chunk_ptr(new fastq_read_chunk()));
  }
}

demultiplex_reads::~demultiplex_reads() {}

chunk_vec
demultiplex_reads::flush_cache(bool eof)
{
  chunk_vec output;

  if (eof || m_unidentified_1->count >= FASTQ_CHUNK_SIZE) {
    m_unidentified_1->eof = eof;
    output.push_back(
      chunk_pair(ai_write_unidentified_1, std::move(m_unidentified_1)));
    m_unidentified_1 = output_chunk_ptr(new fastq_output_chunk());
  }

  if (m_config->paired_ended_mode && !m_config->interleaved_output &&
      (eof || m_unidentified_2->count >= FASTQ_CHUNK_SIZE)) {
    m_unidentified_2->eof = eof;
    output.push_back(
      chunk_pair(ai_write_unidentified_2, std::move(m_unidentified_2)));
    m_unidentified_2 = output_chunk_ptr(new fastq_output_chunk());
  }

  for (size_t nth = 0; nth < m_cache.size(); ++nth) {
    read_chunk_ptr& chunk = m_cache.at(nth);
    if (eof || chunk->reads_1.size() >= FASTQ_CHUNK_SIZE) {
      chunk->eof = eof;

      const size_t step_id = (nth + 1) * ai_analyses_offset;
      output.push_back(chunk_pair(step_id, std::move(chunk)));
      chunk = read_chunk_ptr(new fastq_read_chunk());
    }
  }

  return output;
}

///////////////////////////////////////////////////////////////////////////////

demultiplex_se_reads::demultiplex_se_reads(
  const userconfig* config,
  demultiplexing_statistics* statistics)
  : demultiplex_reads(config, statistics)
{}

chunk_vec
demultiplex_se_reads::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
  read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));

  for (const auto& read : read_chunk->reads_1) {
    const int best_barcode = m_barcode_table.identify(read);

    if (best_barcode < 0) {
      m_unidentified_1->add(*m_config->quality_output_fmt, read);

      if (best_barcode == -1) {
        m_statistics->unidentified += 1;
      } else {
        m_statistics->ambiguous += 1;
      }

      m_statistics->unidentified_stats.process(read);
    } else {
      read_chunk_ptr& dst = m_cache.at(best_barcode);
      dst->reads_1.push_back(read);
      dst->reads_1.back().truncate(m_barcodes.at(best_barcode).first.length());

      m_statistics->barcodes.at(best_barcode) += 1;
    }
  }

  return flush_cache(read_chunk->eof);
}

///////////////////////////////////////////////////////////////////////////////

demultiplex_pe_reads::demultiplex_pe_reads(
  const userconfig* config,
  demultiplexing_statistics* statistics)
  : demultiplex_reads(config, statistics)
{}

chunk_vec
demultiplex_pe_reads::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
  read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
  AR_DEBUG_ASSERT(read_chunk->reads_1.size() == read_chunk->reads_2.size());

  fastq_vec::iterator it_1 = read_chunk->reads_1.begin();
  fastq_vec::iterator it_2 = read_chunk->reads_2.begin();
  for (; it_1 != read_chunk->reads_1.end(); ++it_1, ++it_2) {
    const int best_barcode = m_barcode_table.identify(*it_1, *it_2);

    if (best_barcode < 0) {
      m_unidentified_1->add(*m_config->quality_output_fmt, *it_1);
      if (m_config->interleaved_output) {
        m_unidentified_1->add(*m_config->quality_output_fmt, *it_2);
      } else {
        m_unidentified_2->add(*m_config->quality_output_fmt, *it_2);
      }

      if (best_barcode == -1) {
        m_statistics->unidentified += 2;
      } else {
        m_statistics->ambiguous += 2;
      }

      m_statistics->unidentified_stats.process(*it_1);
      m_statistics->unidentified_stats.process(*it_2);
    } else {
      read_chunk_ptr& dst = m_cache.at(best_barcode);

      it_1->truncate(m_barcodes.at(best_barcode).first.length());
      dst->reads_1.push_back(*it_1);
      it_2->truncate(m_barcodes.at(best_barcode).second.length());
      dst->reads_2.push_back(*it_2);

      m_statistics->barcodes.at(best_barcode) += 2;
    }
  }

  return flush_cache(read_chunk->eof);
}
