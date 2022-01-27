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
#include <cstdlib> // for size_t
#include <string>  // for string

#include "fastq.hpp"     // for ACGT_TO_IDX, fastq
#include "fastq_enc.hpp" // for PHRED_OFFSET_33
#include "statistics.hpp"
#include "utilities.hpp" // for prng_seed

fastq_statistics::fastq_statistics(double sample_rate)
  : m_sample_rate(sample_rate)
  , m_rng(prng_seed())
  , m_number_of_input_reads()
  , m_number_of_output_reads()
  , m_number_of_sampled_reads()
  , m_length_dist()
  , m_quality_dist(MAX_PHRED_SCORE + 1)
  , m_gc_content_dist(101)
  , m_uncalled_pos()
  , m_uncalled_quality_pos()
  , m_called_pos(4)
  , m_quality_pos(4)
  , m_max_sequence_len()
{}

void
fastq_statistics::process(const fastq& read, size_t num_input_reads)
{
  m_number_of_input_reads += num_input_reads;
  m_number_of_output_reads++;

  if (read.length() >= m_max_sequence_len) {
    m_max_sequence_len = read.length();

    m_length_dist.resize_up_to(m_max_sequence_len + 1);

    m_uncalled_pos.resize_up_to(m_max_sequence_len);
    m_uncalled_quality_pos.resize_up_to(m_max_sequence_len);

    for (size_t nuc_i = 0; nuc_i < 4; ++nuc_i) {
      m_called_pos.at(nuc_i).resize_up_to(m_max_sequence_len);
      m_quality_pos.at(nuc_i).resize_up_to(m_max_sequence_len);
    }
  }

  m_length_dist.inc(read.length());

  if (std::generate_canonical<float, 32>(m_rng) <= m_sample_rate) {
    m_number_of_sampled_reads += num_input_reads;

    const std::string& sequence = read.sequence();
    const std::string& qualities = read.qualities();

    size_t acgt_counts[] = { 0, 0, 0, 0 };
    for (size_t i = 0; i < sequence.length(); ++i) {
      const auto nuc = sequence.at(i);
      if (nuc == 'N') {
        m_uncalled_pos.inc(i);
        m_uncalled_quality_pos.inc(i, qualities.at(i));
      } else {
        const auto nuc_i = ACGT_TO_IDX(nuc);

        acgt_counts[nuc_i]++;
        m_called_pos.at(nuc_i).inc(i);
        m_quality_pos.at(nuc_i).inc(i, qualities.at(i) - PHRED_OFFSET_33);
      }

      m_quality_dist.inc(qualities.at(i) - PHRED_OFFSET_33);
    }

    auto n_at = acgt_counts[ACGT_TO_IDX('A')] + acgt_counts[ACGT_TO_IDX('T')];
    auto n_gc = acgt_counts[ACGT_TO_IDX('G')] + acgt_counts[ACGT_TO_IDX('C')];
    if (n_at || n_gc) {
      m_gc_content_dist.inc((100.0 * n_gc) / (n_at + n_gc) + 0.5);
    }
  }
}

fastq_statistics&
fastq_statistics::operator+=(const fastq_statistics& other)
{
  m_number_of_input_reads += other.m_number_of_input_reads;
  m_number_of_output_reads += other.m_number_of_output_reads;
  m_number_of_sampled_reads += other.m_number_of_sampled_reads;
  m_length_dist += other.m_length_dist;
  m_quality_dist += other.m_quality_dist;
  m_uncalled_pos += other.m_uncalled_pos;
  m_uncalled_quality_pos += other.m_uncalled_quality_pos;

  for (size_t i = 0; i < 4; ++i) {
    m_called_pos.at(i) += other.m_called_pos.at(i);
    m_quality_pos.at(i) += other.m_quality_pos.at(i);
  }

  return *this;
}

trimming_statistics::trimming_statistics(double sample_rate)
  : read_1(sample_rate)
  , read_2(sample_rate)
  , merged(sample_rate)
  , discarded(sample_rate)
  , adapter_trimmed_reads()
  , adapter_trimmed_bases()
  , overlapping_reads_merged()
  , terminal_bases_trimmed()
  , low_quality_trimmed_reads()
  , low_quality_trimmed_bases()
  , filtered_min_length_reads()
  , filtered_min_length_bases()
  , filtered_max_length_reads()
  , filtered_max_length_bases()
  , filtered_ambiguous_reads()
  , filtered_ambiguous_bases()
  , filtered_low_complexity_reads()
  , filtered_low_complexity_bases()
{}

trimming_statistics&
trimming_statistics::operator+=(const trimming_statistics& other)
{
  read_1 += other.read_1;
  read_2 += other.read_2;
  merged += other.merged;
  discarded += other.discarded;

  adapter_trimmed_reads += other.adapter_trimmed_reads;
  adapter_trimmed_bases += other.adapter_trimmed_bases;
  overlapping_reads_merged += other.overlapping_reads_merged;
  terminal_bases_trimmed += other.terminal_bases_trimmed;
  low_quality_trimmed_reads += other.low_quality_trimmed_reads;
  low_quality_trimmed_bases += other.low_quality_trimmed_bases;
  filtered_min_length_reads += other.filtered_min_length_reads;
  filtered_min_length_bases += other.filtered_min_length_bases;
  filtered_max_length_reads += other.filtered_max_length_reads;
  filtered_max_length_bases += other.filtered_max_length_bases;
  filtered_ambiguous_reads += other.filtered_ambiguous_reads;
  filtered_ambiguous_bases += other.filtered_ambiguous_bases;
  filtered_low_complexity_reads += other.filtered_low_complexity_reads;
  filtered_low_complexity_bases += other.filtered_low_complexity_bases;

  return *this;
}

demultiplexing_statistics::demultiplexing_statistics(double sample_rate)
  : barcodes()
  , unidentified(0)
  , ambiguous(0)
  , unidentified_stats_1(sample_rate)
  , unidentified_stats_2(sample_rate)
{}

bool
demultiplexing_statistics::empty() const
{
  return barcodes.empty();
}

void
demultiplexing_statistics::resize(size_t n)
{
  barcodes.resize(n);
}

size_t
demultiplexing_statistics::total() const
{
  size_t total = unidentified + ambiguous;
  for (size_t i = 0; i < barcodes.size(); ++i) {
    total += barcodes.at(i);
  }

  return total;
}
