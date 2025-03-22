// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_id.hpp" // declarations
#include "utilities.hpp"  // for merge
#include <algorithm>      // for max, copy, min, reverse
#include <queue>          // for priority_queue
#include <sstream>        // for basic_ostringstream

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Utility functions

namespace {

using adapter_kmer = std::pair<size_t, uint32_t>;
using adapter_kmer_vec = std::vector<adapter_kmer>;

/**
 * Hashing function for string consisting of the chars "ACGT" (uppercase only).
 * Will return a unique number in the range 0 to 4^N - 1 for a given nucleotide
 * sequence. Passing characters other than "ACGT" (uppercase only) will result
 * in hash collisions.
 */
inline size_t
kmer_to_size_t(const std::string& kmer)
{
  size_t index = 0;
  for (const char c : kmer) {
    index = (index << 2) | ACGT::to_index(c);
  }

  return index;
}

/** Translates a hash generated using kmer_to_size_t into a NT sequence. */
inline std::string
size_t_to_kmer(size_t kmer)
{
  std::string kmer_s(consensus_adapter_stats::kmer_length, 'N');
  for (size_t i = 1; i <= consensus_adapter_stats::kmer_length; ++i) {
    kmer_s.at(consensus_adapter_stats::kmer_length - i) =
      ACGT::to_value(kmer & 0x3);
    kmer = kmer >> 2;
  }

  return kmer_s;
}

/** Functor for sorting kmers by frequency. */
struct cmp_nt_count
{
  bool operator()(const adapter_kmer& a, const adapter_kmer& b) const
  {
    return (a.second > b.second);
  }
};

using kmer_queue =
  std::priority_queue<adapter_kmer, adapter_kmer_vec, cmp_nt_count>;

/**
 * Takes an indexed_counts object, and returns a fastq sequence containing the
 * majority nucleotide at each position.The quality score of each base is
 * calculated as the proportion of the bases which match the majority nucleotide
 * (p = m / (N + 1)). If no majority nucleotide can be found 'N' is used
 * instead.
 */
fastq
build_consensus_sequence(const indexed_counts<ACGTN>& consensus)
{
  std::ostringstream sequence;
  std::ostringstream qualities;

  for (size_t i = 0; i < consensus.size(); ++i) {
    // Always assume one non-consensus observation; this is more reasonable
    // than allowing an error-rate of 0, especially for few observations.
    size_t total_count = 1;

    char best_nuc = 'N';
    size_t best_count = 0;
    for (const char nuc : ACGT::values) {
      const size_t cur_count = consensus.get(nuc, i);
      total_count += cur_count;

      if (cur_count > best_count) {
        best_nuc = nuc;
        best_count = cur_count;
      } else if (cur_count == best_count) {
        best_nuc = 'N';
      }
    }

    sequence << best_nuc;

    if (best_nuc == 'N') {
      qualities << static_cast<char>(PHRED_OFFSET_MIN);
    } else {
      const double pvalue = 1.0 - best_count / static_cast<double>(total_count);
      qualities << fastq::p_to_phred_33(pvalue);
    }
  }

  return { "consensus", sequence.str(), qualities.str() };
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// consensus_adapter

consensus_adapter::consensus_adapter(const indexed_counts<ACGTN>& consensus,
                                     const kmer_map& kmers,
                                     const size_t n_kmers)
  : m_adapter(build_consensus_sequence(consensus))
{
  kmer_queue queue;
  for (size_t i = 0; i < kmers.size(); ++i) {
    adapter_kmer value(i, kmers.at(i));
    m_total_kmers += value.second;

    if (queue.size() >= n_kmers) {
      // The top value will be the currently lowest value in the queue
      if (queue.top().second < value.second) {
        queue.pop();
        queue.push(value);
      }
    } else if (value.second) {
      queue.push(value);
    }
  }

  while (!queue.empty()) {
    const auto pair = queue.top();

    m_top_kmers.emplace_back(size_t_to_kmer(pair.first), pair.second);
    queue.pop();
  }

  std::reverse(m_top_kmers.begin(), m_top_kmers.end());
}

std::string
consensus_adapter::compare_with(const std::string& other) const
{
  const auto& adapter = m_adapter.sequence();
  const auto size = std::min(adapter.size(), other.size());

  std::ostringstream identity;
  for (size_t i = 0; i < size; ++i) {
    if (adapter.at(i) == 'N' || other.at(i) == 'N') {
      identity << '*';
    } else {
      identity << (adapter.at(i) == other.at(i) ? '|' : ' ');
    }
  }

  return identity.str();
}

///////////////////////////////////////////////////////////////////////////////
// consensus_adapter_stats

consensus_adapter_stats::consensus_adapter_stats(size_t max_length)
  : m_max_length(max_length)
  , m_kmers(kmer_count, 0)
{
}

consensus_adapter_stats&
consensus_adapter_stats::operator+=(const consensus_adapter_stats& other)
{
  m_consensus += other.m_consensus;
  merge(m_kmers, other.m_kmers);

  return *this;
}

void
consensus_adapter_stats::process(const std::string& sequence)
{
  const auto length = std::min(m_max_length, sequence.length());

  m_consensus.resize_up_to(length);
  for (size_t i = 0; i < length; ++i) {
    m_consensus.inc(sequence.at(i), i);
  }

  if (sequence.length() >= consensus_adapter_stats::kmer_length) {
    const std::string kmer =
      sequence.substr(0, consensus_adapter_stats::kmer_length);
    if (kmer.find('N') == std::string::npos) {
      m_kmers.at(kmer_to_size_t(kmer)) += 1;
    }
  }
}

consensus_adapter
consensus_adapter_stats::summarize(size_t n_kmers) const
{
  return { m_consensus, m_kmers, n_kmers };
}

adapter_id_statistics::adapter_id_statistics(size_t max_length)
  : adapter1(max_length)
  , adapter2(max_length)
  , aligned_pairs(0)
  , pairs_with_adapters(0)
{
}

/** Merge overall trimming_statistics, consensus, and k-mer counts. */
adapter_id_statistics&
adapter_id_statistics::operator+=(const adapter_id_statistics& other)
{
  adapter1 += other.adapter1;
  adapter2 += other.adapter2;
  aligned_pairs += other.aligned_pairs;
  pairs_with_adapters += other.pairs_with_adapters;

  return *this;
}

} // namespace adapterremoval
