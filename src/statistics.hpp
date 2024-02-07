/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 * Copyright (C) 2010-17 by Simon Andrews
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

#include "commontypes.hpp" // for string_vec
#include "counts.hpp"      // for counts, indexed_count, indexed_counts, rates
#include "fastq.hpp"       // for ACGTN, ACGT
#include "robin_hood.hpp"  // for unordered_flat_map
#include <cstdlib>         // for size_t
#include <memory>          // for shared_ptr
#include <random>          // for mt19937
#include <stdint.h>        // for uint64_t
#include <string>          // for string
#include <vector>          // for vector

namespace adapterremoval {

class demux_statistics;
class duplication_statistics;
class fastq_statistics;
class statistics;
class trimming_statistics;

using demux_stats_ptr = std::shared_ptr<demux_statistics>;
using duplication_stats_ptr = std::shared_ptr<duplication_statistics>;
using fastq_stats_ptr = std::shared_ptr<fastq_statistics>;
using stats_ptr = std::shared_ptr<statistics>;
using trim_stats_ptr = std::shared_ptr<trimming_statistics>;

/** Simple counter class for read/bases statistics */
class reads_and_bases
{
public:
  reads_and_bases(uint64_t reads = 0, uint64_t bases = 0);
  reads_and_bases& operator+=(const reads_and_bases& other);
  reads_and_bases operator+(const reads_and_bases& other);

  /** Increments the number of reads by one the number of bases by an amount */
  void inc(uint64_t bases = 1);

  /** Increments the number of reads */
  void inc_reads(uint64_t reads = 1);

  /** Increments the number of bases */
  void inc_bases(uint64_t bases = 1);

  /** Returns the total number of reads (number of increments) */
  uint64_t reads() const;
  /** Returns the total number of bases */
  uint64_t bases() const;

private:
  //! Number of reads/times inc was called
  uint64_t m_reads;
  //! Sum of the number bases inc'd
  uint64_t m_bases;
};

/**
 * Estimation of the fraction of duplicated sequences.
 *
 * Adapted from FastQC v0.11.9 by Simon Andrews under GPLv3+.
 */
class duplication_statistics
{
public:
  /** Represents the results of a duplication estimate. **/
  struct summary
  {
    summary();

    //! X-axis labels; the number of duplicates : 1, 2, .., >10, >50, >100, ..
    string_vec labels;
    //! Fraction of sequences for X-axis label
    rates total_sequences;
    //! Fraction of sequences for X-axis label after de-duplication
    rates unique_sequences;
    //! Estimated fraction of unique sequences in input
    double unique_frac;
  };

  explicit duplication_statistics(size_t max_unique_sequences);

  /** **/
  void process(const fastq& read);

  summary summarize() const;

  size_t max_unique() const;

private:
  void insert(const std::string& key);

  /** Attempts to correct the number of observations for a given bin */
  double correct_count(size_t bin, size_t count) const;

  using string_counts = robin_hood::unordered_flat_map<std::string, size_t>;

  /** Maximum size of m_sequence_counts. */
  size_t m_max_unique_sequences;
  /** Map of truncated sequences to sequence counts. */
  string_counts m_sequence_counts;

  size_t m_sequences_counted;
  size_t m_sequences_counted_at_max;
};

/** Class used to collect statistics about pre/post-processed FASTQ reads. */
class fastq_statistics
{
public:
  explicit fastq_statistics(double sample_rate = 1.0);

  void process(const fastq& read, size_t num_input_reads = 1);

  inline size_t number_of_input_reads() const
  {
    return m_number_of_input_reads;
  }

  inline size_t number_of_output_reads() const
  {
    return m_number_of_output_reads;
  }

  inline size_t number_of_sampled_reads() const
  {
    return m_number_of_sampled_reads;
  }

  /** Distribution of read lengths */
  inline const counts& length_dist() const { return m_length_dist; }

  /** Distribution of Phred quality scores (offset 0) */
  inline const counts& quality_dist() const { return m_quality_dist; }

  /** Smoothed distribution of GC content */
  inline const rates& gc_content() const { return m_gc_content_dist; }

  /** Counts of ACGTN nucleotides by position */
  inline counts nucleotides_pos(char nuc) const
  {
    return m_nucleotide_pos.to_counts(nuc);
  }

  /** Sum of nucleotide counts by position */
  counts nucleotides_pos() const;

  /** Sum of nucleotide counts by position for GC only */
  inline counts nucleotides_gc_pos() const
  {
    return nucleotides_pos('G') + nucleotides_pos('C');
  }

  /** Sum of base qualities for each nucleotide (ACGTN) by position */
  inline counts qualities_pos(char nuc) const
  {
    return m_quality_pos.to_counts(nuc);
  }

  /** Sum of base qualities for ACGTN by position */
  counts qualities_pos() const;

  /** Count duplications with the given max number of unique sequences. */
  void init_duplication_stats(size_t max_unique);
  /** Return the duplication statistics, if any. */
  const duplication_stats_ptr& duplication() const;

  /** Sum statistics, e.g. those used by different threads. */
  fastq_statistics& operator+=(const fastq_statistics& other);

private:
  //! Sample every nth read for statistics.
  double m_sample_rate;
  //! RNG used to downsample sample reads for curves
  std::mt19937 m_rng;

  //! Number of input reads used to produce these statistics
  size_t m_number_of_input_reads;
  //! The actual number of reads written, e.g. input / 2 for merged.
  size_t m_number_of_output_reads;
  //! Number of reads sampled for computationally expensive stats
  size_t m_number_of_sampled_reads;

  /** Length distribution. */
  counts m_length_dist;
  /** Quality distribution. */
  counts m_quality_dist;
  /** GC content distribution. */
  rates m_gc_content_dist;

  /** Per position nucleotide counts indexed using ACGTN_TO_IDX. */
  indexed_counts<ACGTN> m_nucleotide_pos;
  /** Per position quality score sums indexed using ACGTN_TO_IDX. */
  indexed_counts<ACGTN> m_quality_pos;

  //! Maximum size of read processed; used to resize counters as needed
  size_t m_max_sequence_len;

  //! Optional duplication statistics
  duplication_stats_ptr m_duplication;

  //! Copy construction not supported
  fastq_statistics(const fastq_statistics&) = delete;
  //! Assignment not supported
  fastq_statistics& operator=(const fastq_statistics&) = delete;
};

/** Object used to collect summary statistics for trimming. */
class trimming_statistics
{
public:
  explicit trimming_statistics(double sample_rate = 1.0);

  //! Statistics for first reads
  fastq_stats_ptr read_1;
  //! Statistics for second reads
  fastq_stats_ptr read_2;
  //! Statistics for discarded reads
  fastq_stats_ptr singleton;
  //! Statistics for merged reads
  fastq_stats_ptr merged;
  //! Statistics for discarded reads
  fastq_stats_ptr discarded;

  /** Insert size distribution. */
  counts insert_sizes;

  //! Number of reads with adapters trimmed
  counts adapter_trimmed_reads;
  //! Number of bases trimmed for a given adapter (pair)
  counts adapter_trimmed_bases;

  //! Number of reads that overlap/can be merged
  size_t overlapping_reads;

  //! Total number of reads/bases merged
  reads_and_bases reads_merged;

  //! Number of reads/bases trimmed with --pre-trim5p/3p
  reads_and_bases terminal_pre_trimmed;
  //! Number of reads/bases trimmed with --post-trim5p/3p
  reads_and_bases terminal_post_trimmed;

  //! Number of reads trimmed with --pre-trim-polyx
  indexed_count<ACGT> poly_x_pre_trimmed_reads;
  //! Number of 3' bases trimmed with --pre-trim-polyx
  indexed_count<ACGT> poly_x_pre_trimmed_bases;

  //! Number of reads trimmed with --post-trim-polyx
  indexed_count<ACGT> poly_x_post_trimmed_reads;
  //! Number of 3' bases trimmed with --post-trim-polyx
  indexed_count<ACGT> poly_x_post_trimmed_bases;

  //! Number of reads/bases trimmed for low quality bases
  reads_and_bases low_quality_trimmed;
  //! Total number of reads/bases trimmed
  reads_and_bases total_trimmed;

  //! Number of reads/bases filtered due to length (min)
  reads_and_bases filtered_min_length;
  reads_and_bases filtered_max_length;
  //! Number of reads filtered due to too many Ns
  reads_and_bases filtered_ambiguous;
  //! Number of reads filtered due to low complexity
  reads_and_bases filtered_low_complexity;

  /** Combine statistics objects, e.g. those used by different threads. */
  trimming_statistics& operator+=(const trimming_statistics& other);

  //! Copy construction not supported
  trimming_statistics(const trimming_statistics&) = delete;
  //! Assignment not supported
  trimming_statistics& operator=(const trimming_statistics&) = delete;
};

/** Object used to collect summary statistics for demultiplexing. */
class demux_statistics
{
public:
  explicit demux_statistics(double sample_rate = 1.0);

  size_t total() const;

  //! Number of reads identified for for each barcode (pair)
  std::vector<size_t> barcodes;
  //! Number of reads with no hits
  size_t unidentified;
  //! Number of reads with no single best hit
  size_t ambiguous;

  //! Statistics for unidentified mate 1 reads
  fastq_stats_ptr unidentified_stats_1;
  //! Statistics for unidentified mate 2 reads
  fastq_stats_ptr unidentified_stats_2;

  //! Copy construction not supported
  demux_statistics(const demux_statistics&) = delete;
  //! Assignment not supported
  demux_statistics& operator=(const demux_statistics&) = delete;

  friend class statistics_builder;
};

class statistics
{
public:
  fastq_stats_ptr input_1;
  fastq_stats_ptr input_2;

  demux_stats_ptr demultiplexing;
  std::vector<trim_stats_ptr> trimming;

private:
  explicit statistics(double sample_rate);

  friend class statistics_builder;
};

class statistics_builder
{
public:
  statistics_builder();

  /** The number of barcodes used for demultiplexing (0 if disabled). */
  statistics_builder& demultiplexing(size_t barcodes);
  /** Sampling rate of sequences used for sequence composition statistics. */
  statistics_builder& sample_rate(double rate);
  /** Estimate duplication rate using the algorithm implemented in FastQC. */
  statistics_builder& estimate_duplication(size_t max_unique);

  statistics initialize() const;

private:
  //! Number of demultiplexing barcodes
  size_t m_barcode_count;
  //! Fraction of reads sampled for costly statistics (quality distrib., etc.).
  double m_sample_rate;
  //! The max number of unique sequences counted when estimating duplication.
  size_t m_max_unique;
};

////////////////////////////////////////////////////////////////////////////////
inline reads_and_bases::reads_and_bases(uint64_t reads, uint64_t bases)
  : m_reads(reads)
  , m_bases(bases)
{
}

inline reads_and_bases&
reads_and_bases::operator+=(const reads_and_bases& other)
{
  m_reads += other.m_reads;
  m_bases += other.m_bases;

  return *this;
}

inline reads_and_bases
reads_and_bases::operator+(const reads_and_bases& other)
{
  return { m_reads + other.m_reads, m_bases + other.m_bases };
}

inline void
reads_and_bases::inc(uint64_t bases)
{
  m_reads += 1;
  m_bases += bases;
}

inline void
reads_and_bases::inc_reads(uint64_t reads)
{
  m_reads += reads;
}

inline void
reads_and_bases::inc_bases(uint64_t bases)
{
  m_bases += bases;
}

inline uint64_t
reads_and_bases::reads() const
{
  return m_reads;
}

inline uint64_t
reads_and_bases::bases() const
{
  return m_bases;
}

} // namespace adapterremoval
