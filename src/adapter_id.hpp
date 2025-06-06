// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "counts.hpp"    // for indexed_counts
#include "fastq.hpp"     // for fastq
#include "fastq_enc.hpp" // for ACGTN
#include <cstddef>       // for size_t
#include <cstdint>       // for uint32_t
#include <string>        // for string
#include <utility>       // for pair
#include <vector>        // for vector

namespace adapterremoval {

using kmer_map = std::vector<uint32_t>;

/** Structure containing results of consensus adapter inference */
class consensus_adapter
{
public:
  using kmer = std::pair<std::string, size_t>;
  using kmer_vec = std::vector<kmer>;

  /** Constructs consensus adapter sequence and selects the top N kmers */
  consensus_adapter(const indexed_counts<ACGTN>& consensus,
                    const kmer_map& kmers,
                    size_t n_kmers);

  /**
   * Build representation of identity between two sequences.
   *
   * The resulting string represents N with wildcards ('*'), matching bases with
   * pipes ('|') and mismatches with spaces (' '). Only overlapping bases are
   * compared.
   */
  std::string compare_with(const std::string& other) const;

  /** Returns the consensus adapter sequence */
  inline const fastq& adapter() const { return m_adapter; }

  /** Returns vector containing the top N KMers */
  const kmer_vec& top_kmers() const { return m_top_kmers; }

  /** Returns the total number of KMers recorded */
  size_t total_kmers() const { return m_total_kmers; }

private:
  //! Consensus adapter sequence
  fastq m_adapter{};
  //! Vector of the top N k-mer sequences and the number of observations
  kmer_vec m_top_kmers{};
  //! Total number of kmers observed
  size_t m_total_kmers = 0;
};

/** Raw statistics for consensus adapter */
class consensus_adapter_stats
{
public:
  //! Length of kmers to collect to find common kmers
  static const size_t kmer_length = 9;
  //! Size of vector needed for k-mer counts
  static const size_t kmer_count = 2LLU << (2 * kmer_length);
  //! The N most common kmers to print
  static const size_t top_n_kmers = 5;

  explicit consensus_adapter_stats(size_t max_length);
  /** Merge overall trimming_statistics, consensus, and k-mer counts. */
  consensus_adapter_stats& operator+=(const consensus_adapter_stats& other);

  /** Returns the max size of inferred consensus adapter sequences */
  size_t max_length() const { return m_max_length; }

  /** Process an adapter fragment, assumed to only contain bases ACGTN */
  void process(const std::string& sequence);
  /** Constructs consensus adapter sequence and selects the top N kmers */
  consensus_adapter summarize(size_t n_kmers = top_n_kmers) const;

private:
  //! Maximum length of consensus adapter sequence
  size_t m_max_length = 0;
  //! Nucleotide frequencies of putative adapter fragments
  indexed_counts<ACGTN> m_consensus{};
  //! 5' k-mer frequencies of putative adapter fragments
  kmer_map m_kmers{};
};

/** Struct for collecting adapter fragment statistics */
class adapter_id_statistics
{
public:
  explicit adapter_id_statistics(size_t max_length);
  ~adapter_id_statistics() = default;

  /** Merge overall trimming_statistics, consensus, and k-mer counts. */
  adapter_id_statistics& operator+=(const adapter_id_statistics& other);

  //! Statistics based on putative adapter 1 fragments
  consensus_adapter_stats adapter1;
  //! Statistics based on putative adapter 2 fragments
  consensus_adapter_stats adapter2;
  //! Number of properly aligned reads
  size_t aligned_pairs;
  //! Number of reads with adapter fragments
  size_t pairs_with_adapters;

  adapter_id_statistics(const adapter_id_statistics&) = delete;
  adapter_id_statistics(adapter_id_statistics&&) = delete;
  adapter_id_statistics& operator=(const adapter_id_statistics&) = delete;
  adapter_id_statistics& operator=(adapter_id_statistics&&) = delete;
};

} // namespace adapterremoval
