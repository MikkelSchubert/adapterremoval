/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include "counts.hpp" // for indexed_counts
#include "fastq.hpp"  // for ACGTN
#include <cstddef>    // for size_t
#include <utility>    // for pair
#include <vector>     // for vector

namespace adapterremoval {

using kmer_map = std::vector<uint32_t>;

/** Structure containing results of consensus adapter inference */
struct consensus_adapter
{
public:
  using kmer = std::pair<std::string, size_t>;
  using kmer_vec = std::vector<kmer>;

  /** Constructs consensus adapter sequence and selects the top N kmers */
  consensus_adapter(const indexed_counts<ACGTN>& consensus,
                    const kmer_map& kmers,
                    const size_t n_kmers);

  //! Consensus adapter sequence
  fastq adapter;
  //! Vector of the top N kmer sequences and the number of observations
  kmer_vec top_kmers;
  //! Total number of kmers observed
  size_t total_kmers;
};

/** Raw statistics for consensus adapter */
struct consensus_adapter_stats
{
public:
  //! Length of kmers to collect to find common kmers
  static const size_t kmer_length = 9;
  //! Size of vector needed for kmer counts
  static const size_t kmer_count = 2 << (2 * kmer_length);

  consensus_adapter_stats();

  /** Merge overall trimming_statistics, consensus, and k-mer counts. */
  consensus_adapter_stats& operator+=(const consensus_adapter_stats& other);

  /** Process an adapter fragment, assumed to only contain bases ACGTN */
  void process(const std::string& sequence);

  /** Constructs consensus adapter sequence and selects the top N kmers */
  consensus_adapter summarize(size_t n_kmers) const;

private:
  //! Nucleotide frequencies of putative adapter fragments
  indexed_counts<ACGTN> consensus;
  //! 5' KMer frequencies of putative adapter fragments
  kmer_map kmers;
};

/** Struct for collecting adapter fragments, kmer frequencies, and read stats */
struct adapter_id_stats
{
public:
  adapter_id_stats();

  /** Merge overall trimming_statistics, consensus, and k-mer counts. */
  adapter_id_stats& operator+=(const adapter_id_stats& other);

  //! Statistics based on putative adapter 1 fragments
  consensus_adapter_stats adapter1;
  //! Statistics based on putative adapter 2 fragments
  consensus_adapter_stats adapter2;
  //! Number of properly aligned reads
  size_t aligned_pairs;
  //! Number of reads that could not be aligned
  size_t unaligned_pairs;
  //! Number of reads with adapter fragments
  size_t pairs_with_adapters;

  //! Copy construction not supported
  adapter_id_stats(const adapter_id_stats&) = delete;
  //! Assignment not supported
  adapter_id_stats& operator=(const adapter_id_stats&) = delete;
};

} // namespace adapterremoval
