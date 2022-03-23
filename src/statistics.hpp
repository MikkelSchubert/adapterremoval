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
#pragma once

#include <cstdlib> // for size_t
#include <memory>  // for unique_ptr
#include <random>  // for mt19937
#include <vector>  // for vector

#include "counts.hpp" // for counts
#include "fastq.hpp"  // for ACGT_TO_IDX

class fastq_statistics;
class demux_statistics;
class trimming_statistics;
class statistics;

typedef std::shared_ptr<fastq_statistics> fastq_stats_ptr;
typedef std::shared_ptr<demux_statistics> demux_stats_ptr;
typedef std::shared_ptr<trimming_statistics> trim_stats_ptr;
typedef std::shared_ptr<statistics> stats_ptr;

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

  inline const counts& length_dist() const { return m_length_dist; }
  inline const counts& quality_dist() const { return m_quality_dist; }
  inline const counts& gc_content() const { return m_gc_content_dist; }
  inline const counts& uncalled_pos() const { return m_uncalled_pos; }
  inline const counts& uncalled_quality_pos() const
  {
    return m_uncalled_quality_pos;
  }
  inline const counts& nucleotides_pos(char nuc) const
  {
    return m_called_pos.at(ACGT_TO_IDX(nuc));
  }
  inline const counts& qualities_pos(char nuc) const
  {
    return m_quality_pos.at(ACGT_TO_IDX(nuc));
  }

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
  counts m_gc_content_dist;
  /** Count of uncalled bases per position. */
  counts m_uncalled_pos;
  /** Sum of qualities of Ns per position. */
  counts m_uncalled_quality_pos;
  /** Count of A/C/G/T per position; indexed using ACGT_TO_IDX. */
  std::vector<counts> m_called_pos;
  /** Sum of qualities of A/C/G/Ts per position; indexed using ACGT_TO_IDX. */
  std::vector<counts> m_quality_pos;

  //! Maximum size of read processed; used to resize counters as needed
  size_t m_max_sequence_len;

  //! Copy construction not supported
  fastq_statistics(const fastq_statistics&) = delete;
  //! Assignment not supported
  fastq_statistics& operator=(const fastq_statistics&) = delete;
};

/** Object used to collect summary statistics for trimming. */
class trimming_statistics
{
public:
  trimming_statistics(double sample_rate = 1.0);

  //! Statistics for first reads
  fastq_stats_ptr read_1;
  //! Statistics for second reads
  fastq_stats_ptr read_2;
  //! Statistics for merged reads
  fastq_stats_ptr merged;
  //! Statistics for discarded reads
  fastq_stats_ptr discarded;

  //! Number of reads with adapters trimmed
  counts adapter_trimmed_reads;
  //! Number of bases trimmed for a given adapter (pair)
  counts adapter_trimmed_bases;

  //! Number of paired reads merged
  size_t overlapping_reads_merged;

  //! Number of bases 5p/3p bases trimmed with --trim5p/3p
  size_t terminal_bases_trimmed;

  //! Number of reads/bases trimmed for low quality bases
  size_t low_quality_trimmed_reads;
  size_t low_quality_trimmed_bases;

  //! Number of reads/bases filtered due to length (min)
  size_t filtered_min_length_reads;
  size_t filtered_min_length_bases;
  size_t filtered_max_length_reads;
  size_t filtered_max_length_bases;

  size_t filtered_ambiguous_reads;
  size_t filtered_ambiguous_bases;

  size_t filtered_low_complexity_reads;
  size_t filtered_low_complexity_bases;

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
  demux_statistics(double sample_rate = 1.0);

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
  statistics(double sample_rate);

  friend class statistics_builder;
};

struct statistics_builder
{
public:
  statistics_builder();

  statistics_builder& demultiplexing(size_t barcodes);
  statistics_builder& sample_rate(double rate);

  statistics initialize() const;

private:
  //! Number of demultiplexing barcodes
  size_t m_barcode_count;
  //! Fraction of reads sampled for costly statistics (quality distrib., etc.).
  double m_sample_rate;
};
