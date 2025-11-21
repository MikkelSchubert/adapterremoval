// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
// SPDX-FileCopyrightText: 2010-17 Simon Andrews
#pragma once

#include "counts.hpp"    // for counts, indexed_count, indexed_counts, rates
#include "fastq_enc.hpp" // for ACGTN, ACGT
#include <cstdint>       // for uint64_t
#include <cstdlib>       // for size_t
#include <memory>        // for shared_ptr
#include <random>        // for mt19937
#include <string>        // for string
#include <vector>        // for vector

namespace adapterremoval {

class adapter_id_statistics;
class demux_statistics;
class duplication_statistics;
class fastq_statistics;
class fastq;
class statistics;
class string_counts;
class trimming_statistics;

using demux_stats_ptr = std::shared_ptr<demux_statistics>;
using duplication_stats_ptr = std::shared_ptr<duplication_statistics>;
using fastq_stats_ptr = std::shared_ptr<fastq_statistics>;
using stats_ptr = std::shared_ptr<statistics>;
using trim_stats_ptr = std::shared_ptr<trimming_statistics>;
using adapter_id_stats_ptr = std::shared_ptr<adapter_id_statistics>;

/** Simple counter class for read/bases statistics */
class reads_and_bases
{
public:
  explicit reads_and_bases(uint64_t reads = 0, uint64_t bases = 0);
  reads_and_bases& operator+=(const reads_and_bases& other);
  reads_and_bases operator+(const reads_and_bases& other) const;

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
    std::vector<std::string> labels;
    //! Fraction of sequences for X-axis label
    rates total_sequences;
    //! Fraction of sequences for X-axis label after de-duplication
    rates unique_sequences;
    //! Estimated fraction of unique sequences in input
    double unique_frac;
  };

  explicit duplication_statistics(size_t max_unique_sequences);
  ~duplication_statistics();

  duplication_statistics(const duplication_statistics&) = delete;
  duplication_statistics(duplication_statistics&&) = delete;
  duplication_statistics& operator=(const duplication_statistics&) = delete;
  duplication_statistics& operator=(duplication_statistics&&) = delete;

  /** **/
  void process(const fastq& read);

  summary summarize() const;

private:
  /** Attempts to correct the number of observations for a given bin */
  double correct_count(size_t bin, size_t count) const;

  /** Map of truncated sequences to sequence counts. */
  std::unique_ptr<string_counts> m_sequence_counts;
};

/** Class used to collect statistics about pre/post-processed FASTQ reads. */
class fastq_statistics
{
public:
  explicit fastq_statistics(double sample_rate = 1.0);
  fastq_statistics(double sample_rate, unsigned int seed);
  ~fastq_statistics() = default;

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
  counts nucleotides_pos() const { return m_nucleotide_pos.merge(); }

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
  counts qualities_pos() const { return m_quality_pos.merge(); }

  /** Sum statistics, e.g. those used by different threads. */
  fastq_statistics& operator+=(const fastq_statistics& other);

  fastq_statistics(const fastq_statistics&) = delete;
  fastq_statistics(fastq_statistics&&) = delete;
  fastq_statistics& operator=(const fastq_statistics&) = delete;
  fastq_statistics& operator=(fastq_statistics&&) = delete;

private:
  //! Sample every nth read for statistics.
  double m_sample_rate = 1.0;
  //! RNG used to downsample sample reads for curves
  std::mt19937 m_rng{};

  //! Number of input reads used to produce these statistics
  size_t m_number_of_input_reads = 0;
  //! The actual number of reads written, e.g. input / 2 for merged.
  size_t m_number_of_output_reads = 0;
  //! Number of reads sampled for computationally expensive stats
  size_t m_number_of_sampled_reads = 0;

  /** Length distribution. */
  counts m_length_dist{};
  /** Quality distribution. */
  counts m_quality_dist{};
  /** GC content distribution. */
  rates m_gc_content_dist{};

  /** Per position nucleotide counts indexed using ACGTN::to_index. */
  indexed_counts<ACGTN> m_nucleotide_pos{};
  /** Per position quality score sums indexed using ACGTN::to_index. */
  indexed_counts<ACGTN> m_quality_pos{};
};

/** Object used to collect summary statistics for trimming. */
class trimming_statistics
{
public:
  explicit trimming_statistics(double sample_rate = 1.0);
  ~trimming_statistics() = default;

  //! Statistics for first reads
  fastq_stats_ptr read_1{};
  //! Statistics for second reads
  fastq_stats_ptr read_2{};
  //! Statistics for discarded reads
  fastq_stats_ptr singleton{};
  //! Statistics for merged reads
  fastq_stats_ptr merged{};
  //! Statistics for discarded reads
  fastq_stats_ptr discarded{};

  /** Insert size distribution. */
  counts insert_sizes{};

  //! Number of reads with adapters trimmed
  counts adapter_trimmed_reads{};
  //! Number of bases trimmed for a given adapter (pair)
  counts adapter_trimmed_bases{};

  //! Total number of reads/bases merged
  reads_and_bases reads_merged{};
  //! Number of mismatches where a higher quality base was picked during merging
  size_t mismatches_resolved = 0;
  //! Number of mismatches not resolved during merging; bases set to 'N'
  size_t mismatches_unresolved = 0;
  //! Number of 'N's that could be resolved using the other base
  size_t ns_resolved = 0;

  //! Number of reads/bases trimmed with --pre-trim5p/3p
  reads_and_bases terminal_pre_trimmed{};
  //! Number of reads/bases trimmed with --post-trim5p/3p
  reads_and_bases terminal_post_trimmed{};

  //! Number of reads trimmed with --pre-trim-polyx
  indexed_count<ACGT> poly_x_pre_trimmed_reads{};
  //! Number of 3' bases trimmed with --pre-trim-polyx
  indexed_count<ACGT> poly_x_pre_trimmed_bases{};

  //! Number of reads trimmed with --post-trim-polyx
  indexed_count<ACGT> poly_x_post_trimmed_reads{};
  //! Number of 3' bases trimmed with --post-trim-polyx
  indexed_count<ACGT> poly_x_post_trimmed_bases{};

  //! Number of reads/bases trimmed for low quality bases
  reads_and_bases low_quality_trimmed{};
  //! Total number of reads/bases trimmed
  reads_and_bases total_trimmed{};

  //! Number of reads/bases filtered due to length (min)
  reads_and_bases filtered_min_length{};
  reads_and_bases filtered_max_length{};
  //! Number of reads filtered due to too many Ns
  reads_and_bases filtered_ambiguous{};
  //! Number of reads filtered due to low mean quality
  reads_and_bases filtered_mean_quality{};
  //! Number of reads filtered due to low complexity
  reads_and_bases filtered_low_complexity{};

  /** Combine statistics objects, e.g. those used by different threads. */
  trimming_statistics& operator+=(const trimming_statistics& other);

  trimming_statistics(const trimming_statistics&) = delete;
  trimming_statistics(trimming_statistics&&) = delete;
  trimming_statistics& operator=(const trimming_statistics&) = delete;
  trimming_statistics& operator=(trimming_statistics&&) = delete;
};

/** Object used to collect summary statistics for demultiplexing. */
class demux_statistics
{
public:
  explicit demux_statistics(double sample_rate = 1.0);
  ~demux_statistics() = default;

  size_t total() const;

  //! Number of reads identified for each barcode (pair) for each sample
  std::vector<counts> samples{};
  //! Number of reads with no hits
  size_t unidentified = 0;
  //! Number of reads with no single best hit
  size_t ambiguous = 0;

  //! Statistics for unidentified mate 1 reads
  fastq_stats_ptr unidentified_stats_1{};
  //! Statistics for unidentified mate 2 reads
  fastq_stats_ptr unidentified_stats_2{};

  demux_statistics(const demux_statistics&) = delete;
  demux_statistics(demux_statistics&&) = delete;
  demux_statistics& operator=(const demux_statistics&) = delete;
  demux_statistics& operator=(demux_statistics&&) = delete;

  friend class statistics_builder;
};

class statistics
{
public:
  fastq_stats_ptr input_1{};
  fastq_stats_ptr input_2{};

  demux_stats_ptr demultiplexing{};
  std::vector<trim_stats_ptr> trimming{};

  //! Optional duplication statistics
  duplication_stats_ptr duplication_1{};
  //! Optional duplication statistics
  duplication_stats_ptr duplication_2{};

  //! Optional adapter identification statistics
  adapter_id_stats_ptr adapter_id{};

private:
  explicit statistics(double sample_rate);

  friend class statistics_builder;
};

class statistics_builder
{
public:
  statistics_builder();

  /** The number of barcodes used for demultiplexing (0 if disabled). */
  statistics_builder& demultiplexing(size_t samples);
  /** Sampling rate of sequences used for sequence composition statistics. */
  statistics_builder& sample_rate(double rate);
  /** Estimate duplication rate using the algorithm implemented in FastQC. */
  statistics_builder& estimate_duplication(size_t max_unique);
  /** Enable adapter identification if set to > 0 */
  statistics_builder& adapter_identification(size_t max_length);

  statistics initialize() const;

private:
  //! Number of demultiplexing samples
  size_t m_sample_count;
  //! Fraction of reads sampled for costly statistics (quality distrib., etc.).
  double m_sample_rate;
  //! The max number of unique sequences counted when estimating duplication.
  size_t m_max_unique;
  //! The number of bases to infer for adapter sequences
  size_t m_adapter_id;
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
reads_and_bases::operator+(const reads_and_bases& other) const
{
  return reads_and_bases{ m_reads + other.m_reads, m_bases + other.m_bases };
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
