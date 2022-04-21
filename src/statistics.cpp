/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 * Copyright (C) 2010-17 by Simon Andrews
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
#include <iostream>
#include <string> // for string

#include "debug.hpp"     // for AR_DEBUG_ASSERT
#include "fastq.hpp"     // for ACGT_TO_IDX, fastq
#include "fastq_enc.hpp" // for PHRED_OFFSET_33
#include "statistics.hpp"
#include "userconfig.hpp" // for userconfig
#include "utilities.hpp"  // for prng_seed

////////////////////////////////////////////////////////////////////////////////

const size_t DUPLICATION_LEVELS = 16;
const string_vec DUPLICATION_LABELS = { "1",    "2",   "3",   "4",
                                        "5",    "6",   "7",   "8",
                                        "9",    ">10", ">50", ">100",
                                        ">500", ">1k", ">5k", ">10k" };

size_t
duplication_level(size_t count)
{
  if (count > 10000) {
    return 15;
  } else if (count > 50000) {
    return 14;
  } else if (count > 10000) {
    return 13;
  } else if (count > 500) {
    return 12;
  } else if (count > 100) {
    return 11;
  } else if (count > 50) {
    return 10;
  } else if (count > 10) {
    return 9;
  } else {
    return count - 1;
  }
}

duplication_statistics::summary::summary()
  : labels(DUPLICATION_LABELS)
  , total_sequences(DUPLICATION_LEVELS)
  , unique_sequences(DUPLICATION_LEVELS)
  , unique_frac()
{
}

duplication_statistics::duplication_statistics(size_t max_unique_sequences)
  : m_max_unique_sequences(max_unique_sequences)
  , m_sequence_counts()
  , m_sequences_counted()
  , m_sequences_counted_at_max()
{
  m_sequence_counts.reserve(max_unique_sequences);
}

void
duplication_statistics::process(const fastq& read)
{
  if (m_max_unique_sequences) {
    if (read.length() > 75) {
      insert(read.sequence().substr(0, 50));
    } else {
      insert(read.sequence());
    }
  }
}

duplication_statistics::summary
duplication_statistics::summarize() const
{
  // Histogram of sequence duplication counts
  robin_hood::unordered_flat_map<size_t, size_t> histogram;
  for (const auto& it : m_sequence_counts) {
    histogram[it.second]++;
  }

  duplication_statistics::summary summary;

  summary.total_sequences.resize_up_to(DUPLICATION_LEVELS);
  summary.unique_sequences.resize_up_to(DUPLICATION_LEVELS);

  for (const auto& it : histogram) {
    const auto level = it.first;
    const auto count = correct_count(level, it.second);

    const auto slot = duplication_level(level);
    summary.total_sequences.inc(slot, count * level);
    summary.unique_sequences.inc(slot, count);
  }

  const auto unique_count = summary.unique_sequences.sum();
  const auto total_count = summary.total_sequences.sum();

  if (total_count) {
    summary.unique_frac = unique_count / static_cast<double>(total_count);
    summary.total_sequences = summary.total_sequences / total_count;
    summary.unique_sequences = summary.unique_sequences / unique_count;
  } else {
    summary.unique_frac = 1.0;
  }

  return summary;
}

size_t
duplication_statistics::max_unique() const
{
  return m_max_unique_sequences;
}

double
duplication_statistics::correct_count(size_t bin, size_t count) const
{
  // Adapted from
  //   https://github.com/s-andrews/FastQC/blob/v0.11.9/uk/ac/babraham/FastQC/Modules/DuplicationLevel.java#L158

  // See if we can bail out early
  if (m_sequences_counted_at_max == m_sequences_counted) {
    return count;
  }

  // If there aren't enough sequences left to hide another sequence with this
  // count then we can also skip the calculation
  if (m_sequences_counted - count < m_sequences_counted_at_max) {
    return count;
  }

  // If not then we need to see what the likelihood is that we had another
  // sequence with this number of observations which we would have missed.

  // We'll start by working out the probability of NOT seeing a sequence with
  // this duplication level within the first countAtLimit sequences of
  // count.  This is easier than calculating the probability of
  // seeing it.
  double p_not_seeing_at_limit = 1.0;

  // To save doing long calculations which are never going to produce anything
  // meaningful we'll set a limit to our p-value calculation.  This is the
  // probability below which we won't increase our count by 0.01 of an
  // observation.  Once we're below this we stop caring about the corrected
  // value since it's going to be so close to the observed value that we can
  // just return that instead.
  const double limit_of_caring = 1.0 - (count / (count + 0.01));

  for (size_t i = 0; i < m_sequences_counted_at_max; i++) {
    p_not_seeing_at_limit *= ((m_sequences_counted - i) - bin) /
                             static_cast<double>(m_sequences_counted - i);

    if (p_not_seeing_at_limit < limit_of_caring) {
      p_not_seeing_at_limit = 0;
      break;
    }
  }

  // Scale by the chance of seeing
  return count / (1 - p_not_seeing_at_limit);
}

void
duplication_statistics::insert(const std::string& key)
{
  m_sequences_counted++;

  const auto unique_seqs = m_sequence_counts.size();
  if (unique_seqs >= m_max_unique_sequences) {
    auto it = m_sequence_counts.find(key);
    if (it != m_sequence_counts.end()) {
      it->second++;
    }
  } else {
    m_sequence_counts[key]++;
    m_sequences_counted_at_max++;
  }
}

////////////////////////////////////////////////////////////////////////////////

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
  , m_duplication()
{
}

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

  if (m_duplication) {
    m_duplication->process(read);
  }
}

void
fastq_statistics::init_duplication_stats(size_t max_unique_sequences)
{
  AR_DEBUG_ASSERT(!m_duplication);
  m_duplication = make_shared<duplication_statistics>(max_unique_sequences);
}

const duplication_stats_ptr&
fastq_statistics::duplication() const
{
  return m_duplication;
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
  : read_1(make_shared<fastq_statistics>(sample_rate))
  , read_2(make_shared<fastq_statistics>(sample_rate))
  , singleton(make_shared<fastq_statistics>(sample_rate))
  , merged(make_shared<fastq_statistics>(sample_rate))
  , discarded(make_shared<fastq_statistics>(sample_rate))
  , insert_sizes()
  , adapter_trimmed_reads()
  , adapter_trimmed_bases()
  , overlapping_reads()
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
{
}

trimming_statistics&
trimming_statistics::operator+=(const trimming_statistics& other)
{
  *read_1 += *other.read_1;
  *read_2 += *other.read_2;
  *singleton += *other.singleton;
  *merged += *other.merged;
  *discarded += *other.discarded;

  insert_sizes += other.insert_sizes;
  adapter_trimmed_reads += other.adapter_trimmed_reads;
  adapter_trimmed_bases += other.adapter_trimmed_bases;
  overlapping_reads += other.overlapping_reads;
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

demux_statistics::demux_statistics(double sample_rate)
  : barcodes()
  , unidentified(0)
  , ambiguous(0)
  , unidentified_stats_1(make_shared<fastq_statistics>(sample_rate))
  , unidentified_stats_2(make_shared<fastq_statistics>(sample_rate))
{
}

size_t
demux_statistics::total() const
{
  size_t total = unidentified + ambiguous;
  for (size_t i = 0; i < barcodes.size(); ++i) {
    total += barcodes.at(i);
  }

  return total;
}

statistics::statistics(double sample_rate)
  : input_1(make_shared<fastq_statistics>(sample_rate))
  , input_2(make_shared<fastq_statistics>(sample_rate))
  , demultiplexing(make_shared<demux_statistics>(sample_rate))
  , trimming()
{
}

statistics_builder::statistics_builder()
  : m_barcode_count(0)
  , m_sample_rate(1.0)
  , m_max_unique(0)
{
}

statistics_builder&
statistics_builder::sample_rate(double value)
{
  m_sample_rate = value;

  return *this;
}

statistics_builder&
statistics_builder::demultiplexing(size_t barcodes)
{
  m_barcode_count = barcodes;

  return *this;
}

statistics_builder&
statistics_builder::estimate_duplication(size_t max_unique)
{
  m_max_unique = max_unique;

  return *this;
}

statistics
statistics_builder::initialize() const
{
  statistics stats(m_sample_rate);
  if (m_barcode_count) {
    stats.demultiplexing->barcodes.resize(m_barcode_count);
  }

  stats.input_1->init_duplication_stats(m_max_unique);
  stats.input_2->init_duplication_stats(m_max_unique);

  return stats;
}
