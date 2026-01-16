// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
// SPDX-FileCopyrightText: 2010-17 Simon Andrews
#include "statistics.hpp" // declarations
#include "adapter_id.hpp" // for adapter_id_statistics
#include "fastq.hpp"      // for ACGTN, fastq, ACGT
#include "fastq_enc.hpp"  // for PHRED_OFFSET_MIN, PHRED_SCORE_MAX
#include "robin_hood.hpp" // for unordered_flat_map
#include "utilities.hpp"  // for prng_seed
#include <algorithm>      // for max, min
#include <array>          // for array
#include <cstdlib>        // for size_t
#include <functional>     // for equal_to
#include <memory>         // for make_shared, __shared_ptr_access, __shared_...
#include <string>         // for basic_string, string
#include <string_view>    // for string_view

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////

const size_t DUPLICATION_LEVELS = 16;
const std::array<std::string_view, DUPLICATION_LEVELS> DUPLICATION_LABELS = {
  "1", "2",   "3",   "4",    "5",    "6",   "7",   "8",
  "9", ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k"
};

namespace {

size_t
duplication_level(size_t count)
{
  if (count > 10000) {
    return 15;
  } else if (count > 5000) {
    return 14;
  } else if (count > 1000) {
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

} // namespace

class string_counts
{
  using map = robin_hood::unordered_flat_map<std::string, size_t>;

public:
  explicit string_counts(size_t max_unique_sequences)
    : m_max_unique_sequences(max_unique_sequences)
  {
    m_counts.reserve(max_unique_sequences);
  }

  void add(const std::string& key)
  {
    m_counted++;

    const auto unique_seqs = m_counts.size();
    if (unique_seqs >= m_max_unique_sequences) {
      auto it = m_counts.find(key);
      if (it != m_counts.end()) {
        it->second++;
      } else if (!m_counted_at_max_set) {
        // Small difference from FastQC: Continue counting until the max unique
        // sequence count is exceeded, rather than until the max unique sequence
        // count is reached. For unique counts close to the max, this should
        // give slightly more accurate results
        m_counted_at_max_set = true;
        m_counted_at_max = m_counted;
      }
    } else {
      m_counts[key]++;
    }
  }

  [[nodiscard]] auto counted() const { return m_counted; }

  [[nodiscard]] auto counted_at_max() const
  {
    return m_counted_at_max_set ? m_counted_at_max : m_counted;
  }

  [[nodiscard]] auto begin() const { return m_counts.begin(); }

  [[nodiscard]] auto end() const { return m_counts.end(); }

private:
  /** Maximum size of m_sequence_counts. */
  size_t m_max_unique_sequences{};
  /** Total number of sequences counted */
  size_t m_counted{};
  /** Sequences counted before max number of unique sequences was exceeded */
  size_t m_counted_at_max{};
  /** If m_counted_at_max has been set, when capacity is first exceeded */
  bool m_counted_at_max_set = false;

  /** Map of truncated sequences to sequence counts. */
  map m_counts{};
};

duplication_statistics::summary::summary()
  : labels(DUPLICATION_LABELS.begin(), DUPLICATION_LABELS.end())
  , total_sequences(DUPLICATION_LEVELS)
  , unique_sequences(DUPLICATION_LEVELS)
  , unique_frac()
{
}

duplication_statistics::duplication_statistics(size_t max_unique_sequences)
  : m_sequence_counts()
{
  if (max_unique_sequences) {
    m_sequence_counts = std::make_unique<string_counts>(max_unique_sequences);
  }
}

duplication_statistics::~duplication_statistics()
{
  m_sequence_counts.reset();
}

void
duplication_statistics::process(const fastq& read)
{
  if (m_sequence_counts) {
    if (read.length() > 75) {
      m_sequence_counts->add(read.sequence().substr(0, 50));
    } else {
      m_sequence_counts->add(read.sequence());
    }
  }
}

duplication_statistics::summary
duplication_statistics::summarize() const
{
  // Histogram of sequence duplication counts
  robin_hood::unordered_flat_map<size_t, size_t> histogram;
  for (const auto& it : *m_sequence_counts) {
    histogram[it.second]++;
  }

  duplication_statistics::summary result;

  result.total_sequences.resize_up_to(DUPLICATION_LEVELS);
  result.unique_sequences.resize_up_to(DUPLICATION_LEVELS);

  for (const auto& it : histogram) {
    const auto level = it.first;
    const auto count = correct_count(level, it.second);

    const auto slot = duplication_level(level);
    result.total_sequences.inc(slot, count * level);
    result.unique_sequences.inc(slot, count);
  }

  const auto unique_count = result.unique_sequences.sum();
  const auto total_count = result.total_sequences.sum();

  if (total_count) {
    result.unique_frac = unique_count / static_cast<double>(total_count);
    result.total_sequences = result.total_sequences / total_count;
    result.unique_sequences = result.unique_sequences / unique_count;
  } else {
    result.unique_frac = 1.0;
  }

  return result;
}

double
duplication_statistics::correct_count(size_t bin, size_t count) const
{
  // Adapted from
  //   https://github.com/s-andrews/FastQC/blob/v0.11.9/uk/ac/babraham/FastQC/Modules/DuplicationLevel.java#L158

  const auto sequences_counted = m_sequence_counts->counted();
  const auto sequences_counted_at_max = m_sequence_counts->counted_at_max();

  // See if we can bail out early
  if (sequences_counted_at_max == sequences_counted) {
    return count;
  }

  // If there aren't enough sequences left to hide another sequence with this
  // count then we can also skip the calculation
  if (sequences_counted - count < sequences_counted_at_max) {
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

  for (size_t i = 0; i < sequences_counted_at_max; i++) {
    p_not_seeing_at_limit *= ((sequences_counted - i) - bin) /
                             static_cast<double>(sequences_counted - i);

    if (p_not_seeing_at_limit < limit_of_caring) {
      p_not_seeing_at_limit = 0;
      break;
    }
  }

  // Scale by the chance of seeing
  return count / (1 - p_not_seeing_at_limit);
}

////////////////////////////////////////////////////////////////////////////////

namespace {

void
smoothed_gc_count(rates& distribution, size_t count, size_t length)
{
  // counts are smoothed across adjacent (percentage) bins
  const auto lower = std::max<double>(0, 100 * (count - 0.5) / length);
  const auto upper = std::min<double>(100, 100 * (count + 0.5) / length);

  const size_t lower_i = lower;
  const size_t upper_i = upper;

  const auto lower_f = lower - lower_i;
  const auto upper_f = upper - upper_i;

  if (lower_i == upper_i) {
    // Count falls in a single bin
    distribution.inc(lower_i, (upper - lower));
  } else {
    // Increment first bin/partially overlapped adjacent bins
    distribution.inc(lower_i, (1 - lower_f));
    distribution.inc(upper_i, upper_f);

    // Ranges are half open, except for the final bin
    const auto final_i = count == length ? 101 : upper_i;
    for (size_t i = lower + 1; i < final_i; ++i) {
      distribution.inc(i);
    }
  }
}

} // namespace

fastq_statistics::fastq_statistics(double sample_rate)
  : fastq_statistics(sample_rate, prng_seed())
{
}

fastq_statistics::fastq_statistics(double sample_rate, uint32_t seed)
  : m_sample_rate(sample_rate)
  , m_rng(seed)
  , m_quality_dist(PHRED_SCORE_MAX + 1)
  , m_gc_content_dist(101)
{
}

void
fastq_statistics::process(const fastq& read, size_t num_input_reads)
{
  m_number_of_input_reads += num_input_reads;
  m_number_of_output_reads++;

  const auto sequence_len = read.length();
  m_length_dist.resize_up_to(sequence_len + 1);
  m_length_dist.inc(sequence_len);

  if (std::generate_canonical<float, 32>(m_rng) <= m_sample_rate) {
    m_nucleotide_pos.resize_up_to(sequence_len);
    m_quality_pos.resize_up_to(sequence_len);

    m_number_of_sampled_reads++;

    const std::string_view sequence = read.sequence();
    const std::string_view qualities = read.qualities();

    indexed_count<ACGTN> nucls;
    for (size_t i = 0; i < sequence_len; ++i) {
      const auto nuc = sequence[i];

      nucls.inc(nuc);
      m_nucleotide_pos.inc_unsafe(nuc, i);

      const auto quality = qualities[i] - PHRED_OFFSET_MIN;
      m_quality_pos.inc_unsafe(nuc, i, quality);
      m_quality_dist.inc(quality);
    }

    auto n_at = nucls.get('A') + nucls.get('T');
    auto n_gc = nucls.get('G') + nucls.get('C');
    if (n_at || n_gc) {
      // uncalled bases are not considered for the length of the reads
      smoothed_gc_count(m_gc_content_dist, n_gc, n_gc + n_at);
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
  m_gc_content_dist += other.m_gc_content_dist;
  m_nucleotide_pos += other.m_nucleotide_pos;
  m_quality_pos += other.m_quality_pos;

  return *this;
}

trimming_statistics::trimming_statistics(double sample_rate)
  : singleton(std::make_shared<fastq_statistics>(sample_rate))
  , merged(std::make_shared<fastq_statistics>(sample_rate))
  , discarded(std::make_shared<fastq_statistics>(sample_rate))
{
  // Synchronize sampling of mate 1 and mate 2 reads
  const auto seed = prng_seed();
  read_1 = std::make_shared<fastq_statistics>(sample_rate, seed);
  read_2 = std::make_shared<fastq_statistics>(sample_rate, seed);
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
  reads_merged += other.reads_merged;
  mismatches_resolved += other.mismatches_resolved;
  mismatches_unresolved += other.mismatches_unresolved;
  ns_resolved += other.ns_resolved;
  terminal_pre_trimmed += other.terminal_pre_trimmed;
  terminal_post_trimmed += other.terminal_post_trimmed;
  poly_x_pre_trimmed_reads += other.poly_x_pre_trimmed_reads;
  poly_x_pre_trimmed_bases += other.poly_x_pre_trimmed_bases;
  poly_x_post_trimmed_reads += other.poly_x_post_trimmed_reads;
  poly_x_post_trimmed_bases += other.poly_x_post_trimmed_bases;
  low_quality_trimmed += other.low_quality_trimmed;
  total_trimmed += other.total_trimmed;
  filtered_min_length += other.filtered_min_length;
  filtered_max_length += other.filtered_max_length;
  filtered_ambiguous += other.filtered_ambiguous;
  filtered_mean_quality += other.filtered_mean_quality;
  filtered_low_complexity += other.filtered_low_complexity;

  return *this;
}

demux_statistics::demux_statistics(double sample_rate)
  : unidentified_stats_1(std::make_shared<fastq_statistics>(sample_rate))
  , unidentified_stats_2(std::make_shared<fastq_statistics>(sample_rate))
{
}

size_t
demux_statistics::total() const
{
  size_t total = unidentified + ambiguous;
  for (const auto& count : samples) {
    total += count.sum();
  }

  return total;
}

statistics::statistics(double sample_rate)
  : input_1(std::make_shared<fastq_statistics>(sample_rate))
  , input_2(std::make_shared<fastq_statistics>(sample_rate))
  , demultiplexing(std::make_shared<demux_statistics>(sample_rate))
{
}

statistics_builder&
statistics_builder::sample_rate(double rate)
{
  m_sample_rate = rate;

  return *this;
}

statistics_builder&
statistics_builder::demultiplexing(size_t samples)
{
  m_sample_count = samples;

  return *this;
}

statistics_builder&
statistics_builder::estimate_duplication(size_t max_unique)
{
  m_max_unique = max_unique;

  return *this;
}

/** Enable adapter identification */
statistics_builder&
statistics_builder::adapter_identification(size_t max_length)
{
  m_adapter_id = max_length;

  return *this;
}

statistics
statistics_builder::initialize() const
{
  statistics stats(m_sample_rate);

  stats.duplication_1 = std::make_shared<duplication_statistics>(m_max_unique);
  stats.duplication_2 = std::make_shared<duplication_statistics>(m_max_unique);

  if (m_adapter_id) {
    stats.adapter_id = std::make_shared<adapter_id_statistics>(m_adapter_id);
  }

  for (size_t i = 0; i < m_sample_count; ++i) {
    stats.trimming.push_back(std::make_shared<trimming_statistics>());
  }

  return stats;
}

} // namespace adapterremoval
