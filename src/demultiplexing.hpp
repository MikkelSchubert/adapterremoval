// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "barcode_table.hpp" // for barcode_table
#include "output.hpp"        // for demultiplexed_reads
#include "scheduler.hpp"     // for analytical_step
#include "statistics.hpp"    // for fastq_statistics, trimming_statistics
#include "threading.hpp"     // for threadstate
#include <cstddef>           // for size_t
#include <memory>            // for shared_ptr
#include <mutex>             // for mutex
#include <utility>           // for pair
#include <vector>            // for vector

namespace adapterremoval {

class demux_statistics;
class output_files;
class post_demux_steps;
class sample_set;
class trimming_statistics;
class userconfig;
class fastq_statistics;

using demux_stats_ptr = std::shared_ptr<demux_statistics>;
using trim_stats_ptr = std::shared_ptr<trimming_statistics>;

/**
 * Base-class for demultiplexing of reads; responsible for building the
 * quad-tree representing the set of adapter sequences, and for maintaining the
 * cache of demultiplexed reads.
 */
class demultiplex_reads : public analytical_step
{
public:
  /** Setup demultiplexer; keeps reference to config object. */
  demultiplex_reads(const userconfig& config,
                    const post_demux_steps& steps,
                    demux_stats_ptr stats);

protected:
  //! Set of samples to be demultiplexed
  const threadsafe_data<sample_set> m_samples;
  //! Quad-tree representing all mate 1 adapters; for search with n mismatches
  const barcode_table m_barcode_table;
  //! Pointer to user settings used for output format for unidentified reads
  const userconfig& m_config;
  //! Map of steps for output chunks;
  const post_demux_steps& m_steps;

  //! Mapping of barcode global offsets to sample and sample barcode offsets
  std::vector<std::pair<size_t, size_t>> m_barcodes{};
  //! Cache of reads used to buffer chunks for downstream processing
  demultiplexed_reads m_cache;

  //! Sink for demultiplexing statistics; used by subclasses.
  demux_stats_ptr m_statistics{};

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock{};
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
  chunk_vec process(chunk_ptr data) override;
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
  chunk_vec process(chunk_ptr data) override;
};

class process_demultiplexed : public analytical_step
{
public:
  process_demultiplexed(const userconfig& config,
                        const sample_output_files& output,
                        size_t sample,
                        trim_stats_ptr sink);

  chunk_vec process(chunk_ptr data) override;

  void finalize() override;

private:
  const userconfig& m_config;
  const sample_output_files& m_output;
  //! Set of samples to be demultiplexed
  const threadsafe_data<sample_set> m_samples;
  //! The sample processed by this step
  const size_t m_sample;

  threadstate<trimming_statistics> m_stats{};
  trim_stats_ptr m_stats_sink{};
};

/** Collects statistics for and serializes unidentified reads */
class processes_unidentified : public analytical_step
{
public:
  /** Setup step; keeps reference to config object. */
  processes_unidentified(const userconfig& config,
                         const output_files& output,
                         demux_stats_ptr stats);

  /** Collect statistics and serialize (if kept) unidentified reads */
  chunk_vec process(chunk_ptr data) override;

  void finalize() override;

protected:
  //! Pointer to user settings used for output format for unidentified reads
  const userconfig& m_config;
  //! Mapping of output files; unidentified reads are treated as a pseudo-sample
  sample_output_files m_output{};
  //! Sink for demultiplexing statistics
  demux_stats_ptr m_statistics{};

  //! Per thread statistics collected from unidentified mate 1 reads
  threadstate<fastq_statistics> m_stats_1{};
  //! Per thread statistics collected from unidentified mate 2 reads
  threadstate<fastq_statistics> m_stats_2{};
};
} // namespace adapterremoval
