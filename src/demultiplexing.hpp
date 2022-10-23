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
#pragma once

#include <algorithm> // for max, copy
#include <mutex>     // for mutex
#include <stddef.h>  // for size_t
#include <vector>    // for vector

#include "barcode_table.hpp" // for barcode_table
#include "fastq.hpp"         // for fastq_pair_vec
#include "fastq_io.hpp"      // for read_chunk_ptr, output_chunk_ptr
#include "scheduler.hpp"     // for chunk_vec, analytical_step
#include "statistics.hpp"    // for demux_statistics, fastq_statistics

namespace adapterremoval {

class userconfig;
class demux_statistics;

/** Map of samples to downstream FASTQ processing/writing steps. */
class post_demux_steps
{
public:
  post_demux_steps();

  /* Step used to write unidentified mate 1 reads. */
  size_t unidentified_1;
  /* Step used to write unidentified mate 2 reads; may be disabled. */
  size_t unidentified_2;

  /* Processing step for each sample. */
  std::vector<size_t> samples;

  /** Constant indicating that a step has been disabled. */
  static const size_t disabled;
};

/**
 * Baseclass for demultiplexing of reads; responsible for building the quad-tree
 * representing the set of adapter sequences, and for maintaining the cache of
 * demultiplexed reads.
 */
class demultiplex_reads : public analytical_step
{
public:
  /** Setup demultiplexer; keeps reference to config object. */
  demultiplex_reads(const userconfig& config,
                    const post_demux_steps& steps,
                    demux_stats_ptr stats);

  /** Frees any unflushed caches. */
  virtual ~demultiplex_reads();

  //! Copy construction not supported
  demultiplex_reads(const demultiplex_reads&) = delete;
  //! Assignment not supported
  demultiplex_reads& operator=(const demultiplex_reads&) = delete;

protected:
  //! List of barcode (pairs) supplied by caller
  const fastq_pair_vec& m_barcodes;
  //! Quad-tree representing all mate 1 adapters; for search with n mismatches
  const barcode_table m_barcode_table;
  //! Pointer to user settings used for output format for unidentified reads
  const userconfig& m_config;

  //! Returns a chunk-list with any set of reads exceeding the max cache size
  //! If 'eof' is true, all chunks are returned, and the 'eof' values in the
  //! chunks are set to true.
  chunk_vec flush_cache(bool eof = false);

  typedef std::vector<read_chunk_ptr> demultiplexed_cache;

  //! Cache of demultiplex reads; used to reduce the number of output chunks
  //! generated from each processed chunk, which would otherwise increase
  //! linearly with the number of barcodes.
  demultiplexed_cache m_cache;
  //! Cache of unidentified mate 1 reads
  output_chunk_ptr m_unidentified_1;
  //! Cache of unidentified mate 2 reads
  output_chunk_ptr m_unidentified_2;

  //! Map of steps for output chunks;
  post_demux_steps m_steps;

  //! Sink for demultiplexing statistics; used by subclasses.
  demux_stats_ptr m_statistics;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
};

/** Demultiplexer for single-end reads. */
class demultiplex_se_reads : public demultiplex_reads
{
public:
  /** See demultiplex_reads::demultiplex_reads. */
  demultiplex_se_reads(const userconfig& config,
                       const post_demux_steps& steps,
                       demux_stats_ptr stats);

  /**
   * Processes a read chunk, and forwards chunks to downstream steps, with
   * the IDs corresponding to ai_analyses_offset * (nth + 1) for the nth
   * barcode (pair). Unidentified reads are sent to ai_write_unidentified_1.
   */
  chunk_vec process(chunk_ptr chunk) override;
};

/** Demultiplexer for paired-end reads. */
class demultiplex_pe_reads : public demultiplex_reads
{
public:
  /** See demultiplex_reads::demultiplex_reads. */
  demultiplex_pe_reads(const userconfig& config,
                       const post_demux_steps& steps,
                       demux_stats_ptr stats);

  /**
   * Processes a read chunk, and forwards chunks to downstream steps, with
   * the IDs corresponding to ai_analyses_offset * (nth + 1) for the nth
   * barcode (pair). Unidentified reads are sent to ai_write_unidentified_1
   * and ai_write_unidentified_2.
   */
  chunk_vec process(chunk_ptr chunk) override;
};

} // namespace adapterremoval
