/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <algorithm> // for max, copy, min, count, reverse
#include <iomanip>   // for operator<<, setw
#include <iostream>  // for operator<<, basic_ostream, char_traits
#include <memory>    // for unique_ptr
#include <queue>     // for priority_queue
#include <sstream>   // for stringstream
#include <stddef.h>  // for size_t
#include <string>    // for string, basic_string, operator<<
#include <utility>   // for pair
#include <vector>    // for vector, vector<>::iterator, vector<>::con...

#include "adapterset.hpp"  // for adapter_set
#include "alignment.hpp"   // for align_paired_ended_sequences, extract_ada...
#include "commontypes.hpp" // for fastq_vec
#include "debug.hpp"       // for AR_DEBUG_ASSERT
#include "fastq.hpp"       // for fastq, ACGT_TO_IDX, fastq_pair_vec, IDX_T...
#include "fastq_io.hpp"    // for fastq_read_chunk, read_fastq, read_chunk_ptr
#include "scheduler.hpp"   // for threadstate, scheduler, analytical_step
#include "timer.hpp"       // for progress_timer
#include "userconfig.hpp"  // for userconfig, fastq_encoding_ptr
#include "vecutils.hpp"    // for merge_vectors

///////////////////////////////////////////////////////////////////////////////
// KMer related functions and constants

//! Length of kmers to collect to find common kmers
const size_t KMER_LENGTH = 9;
//! Size of vector needed for kmer counts
const size_t N_KMERS = 2 << (2 * KMER_LENGTH);
//! The N most common kmers to print
const size_t TOP_N_KMERS = 5;

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
  for (size_t i = 0; i < kmer.length(); ++i) {
    index = (index << 2) | ACGT_TO_IDX(kmer.at(i));
  }

  return index;
}

/** Translates a hash generated using kmer_to_size_t into a NT sequence. */
inline std::string
size_t_to_kmer(size_t kmer)
{
  std::string kmer_s(KMER_LENGTH, 'N');
  for (size_t i = 1; i <= KMER_LENGTH; ++i) {
    kmer_s.at(KMER_LENGTH - i) = "ACTG"[kmer & 0x3];
    kmer = kmer >> 2;
  }

  return kmer_s;
}

/** Simple structure for counting the frequency of A, C, G, and Ts. */
struct nt_counts
{
  nt_counts()
    : counts(4, 0)
  {}

  /** Increment count of a nucleotide A, C, G, or T (uppercase only). */
  void increment(char nt)
  {
    if (nt != 'N') {
      ++counts.at(ACGT_TO_IDX(nt));
    }
  }

  /** Merge count objects. */
  nt_counts& operator+=(const nt_counts& other)
  {
    merge_vectors(counts, other.counts);
    return *this;
  }

  //! Fixed sized vector (4)
  std::vector<size_t> counts;
};

typedef std::vector<nt_counts> nt_count_vec;
typedef std::vector<unsigned> kmer_map;
typedef std::pair<size_t, unsigned> nt_count;

/** Functor for sorting kmers by frequency. */
struct cmp_nt_count
{
  bool operator()(const nt_count& a, const nt_count& b) const
  {
    return (a.second > b.second);
  }
};

typedef std::priority_queue<nt_count, std::vector<nt_count>, cmp_nt_count>
  kmer_queue;
typedef std::vector<nt_count> kmer_vector;

/** Prints the N top kmers in a kmer_map, including sequence and frequency. */
void
print_most_common_kmers(const kmer_map& kmers, size_t print_n = TOP_N_KMERS)
{
  size_t total = 0;
  kmer_queue queue;
  for (size_t i = 0; i < kmers.size(); ++i) {
    nt_count value(i, kmers.at(i));
    total += value.second;

    if (queue.size() >= print_n) {
      // The top value will be the currently lowest value in the queue
      if (queue.top().second < value.second) {
        queue.pop();
        queue.push(value);
      }
    } else if (value.second) {
      queue.push(value);
    }
  }

  kmer_vector top_n_kmers;
  while (!queue.empty()) {
    top_n_kmers.push_back(queue.top());
    queue.pop();
  }

  std::cout.precision(2);
  std::cout << std::fixed;
  std::cout << "    Top 5 most common " << KMER_LENGTH << "-bp 5'-kmers:\n";

  std::reverse(top_n_kmers.begin(), top_n_kmers.end());
  for (size_t i = 0; i < top_n_kmers.size(); ++i) {
    const nt_count& count = top_n_kmers.at(i);
    std::string kmer_s = size_t_to_kmer(count.first);

    std::cout << std::string(12, ' ');
    std::cout << i + 1 << ": " << kmer_s << " = " << std::right << std::setw(5)
              << (100.0 * count.second) / total << "% (" << count.second << ")"
              << "\n";
  }
}

///////////////////////////////////////////////////////////////////////////////
// Consensus adapter related functions and constants

/**
 * Build representation of identity between an adapter and a consensus sequence.
 *
 * The resulting string represents N with wildcards ('*'), matching bases with
 * pipes ('|') and mismatches with spaces (' '). Only overlapping bases are
 * compared.
 */
std::string
compare_consensus_with_ref(const std::string& ref, const std::string& consensus)
{
  std::stringstream identity;
  for (size_t i = 0, size = std::min(consensus.size(), ref.size()); i < size;
       ++i) {
    if (ref.at(i) == 'N' || consensus.at(i) == 'N') {
      identity << '*';
    } else {
      identity << (ref.at(i) == consensus.at(i) ? '|' : ' ');
    }
  }

  return identity.str();
}

/**
 * Takes a nt_counts object, and returns a pair containing the majority nt and
 * the Phred encoded quality score of the consensus, defined as the proportion
 * of the bases which match the majority nucleotide (p = m / (N + 1)). If no
 * majority nucleotide can be found, 'N' is returned instead.
 */
std::pair<char, char>
get_consensus_nt(const nt_counts& nts)
{
  // Always assume one non-consensus observation; this is more reasonable
  // than allowing an error-rate of 0, especially for few observations.
  size_t total_count = 1;

  char best_nt_i = -1;
  size_t best_count = 0;
  for (char nt_i = 0; nt_i < 4; ++nt_i) {
    const size_t cur_count = nts.counts.at(nt_i);
    total_count += cur_count;

    if (cur_count > best_count) {
      best_nt_i = nt_i;
      best_count = cur_count;
    } else if (cur_count == best_count) {
      best_nt_i = -1;
    }
  }

  const double pvalue = 1.0 - best_count / static_cast<double>(total_count);
  const char phred = fastq::p_to_phred_33(pvalue);
  const char best_nt = (best_nt_i == -1) ? 'N' : IDX_TO_ACGT(best_nt_i);

  return std::pair<char, char>(best_nt, phred);
}

/**
 * Prints description of consensus adapter sequence.
 *
 * @param counts Observed nucleotide frequencies.
 * @param kmers Observed kmer frequencies.
 * @param name Argument name for adapter (--adapter1 / adapter2)
 * @param ref Default sequence for the inferred adapter
 */
void
print_consensus_adapter(const nt_count_vec& counts,
                        const kmer_map& kmers,
                        const std::string& name,
                        const std::string& ref)
{
  std::stringstream sequence;
  std::stringstream qualities;

  for (nt_count_vec::const_iterator it = counts.begin(); it != counts.end();
       ++it) {
    const std::pair<char, char> consensus = get_consensus_nt(*it);

    sequence << consensus.first;
    qualities << consensus.second;
  }

  const std::string consensus = sequence.str();
  const std::string identity = compare_consensus_with_ref(ref, consensus);

  std::cout << "  " << name << ":  " << ref << "\n"
            << "               " << identity << "\n"
            << "   Consensus:  " << consensus << "\n"
            << "     Quality:  " << qualities.str() << "\n\n";

  print_most_common_kmers(kmers);
}

///////////////////////////////////////////////////////////////////////////////
//

/** Struct for collecting adapter fragments, kmer frequencies, and read stats.
 */
struct adapter_stats
{
public:
  adapter_stats()
    : pcr1_counts()
    , pcr2_counts()
    , pcr1_kmers(N_KMERS, 0)
    , pcr2_kmers(N_KMERS, 0)
    , aligned_pairs(0)
    , unaligned_pairs(0)
    , pairs_with_adapters(0)
  {}

  /** Merge overall trimming_statistics, consensus, and k-mer counts. */
  adapter_stats& operator+=(const adapter_stats& other)
  {
    merge_vectors(pcr1_counts, other.pcr1_counts);
    merge_vectors(pcr2_counts, other.pcr2_counts);
    merge_vectors(pcr1_kmers, other.pcr1_kmers);
    merge_vectors(pcr2_kmers, other.pcr2_kmers);

    aligned_pairs += other.aligned_pairs;
    unaligned_pairs += other.unaligned_pairs;
    pairs_with_adapters += other.pairs_with_adapters;

    return *this;
  }

  //! Nucleotide frequencies of putative adapter 1 fragments
  nt_count_vec pcr1_counts;
  //! Nucleotide frequencies of putative adapter 2 fragments
  nt_count_vec pcr2_counts;
  //! 5' KMer frequencies of putative adapter 1 fragments
  kmer_map pcr1_kmers;
  //! 5' KMer frequencies of putative adapter 2 fragments
  kmer_map pcr2_kmers;
  //! Number of properly aligned reads
  size_t aligned_pairs;
  //! Number of reads that could not be aligned
  size_t unaligned_pairs;
  //! Number of reads with adapter fragments
  size_t pairs_with_adapters;

  //! Copy construction not supported
  adapter_stats(const adapter_stats&) = delete;
  //! Assignment not supported
  adapter_stats& operator=(const adapter_stats&) = delete;
};

///////////////////////////////////////////////////////////////////////////////
// Threaded adapter identification step

class adapter_identification : public analytical_step
{
public:
  adapter_identification(const userconfig& config)
    : analytical_step(analytical_step::ordering::unordered)
    , m_config(config)
    , m_timer("reads")
    , m_stats()
  {
    for (size_t i = 0; i < m_config.max_threads; ++i) {
      m_stats.emplace_back();
    }
  }

  chunk_vec process(analytical_chunk* chunk)
  {
    AR_DEBUG_ASSERT(chunk);
    read_chunk_ptr file_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
    AR_DEBUG_ASSERT(file_chunk);

    const fastq empty_adapter("dummy", "", "");
    fastq_pair_vec adapters;
    adapters.push_back(fastq_pair(empty_adapter, empty_adapter));

    auto aligner = sequence_aligner(adapters);
    aligner.set_mismatch_threshold(m_config.mismatch_threshold);

    auto stats = m_stats.acquire();

    AR_DEBUG_ASSERT(file_chunk->reads_1.size() == file_chunk->reads_2.size());
    fastq_vec::iterator read_1 = file_chunk->reads_1.begin();
    fastq_vec::iterator read_2 = file_chunk->reads_2.begin();

    while (read_1 != file_chunk->reads_1.end()) {
      process_reads(aligner, *stats, *read_1++, *read_2++);
    }

    m_stats.release(stats);
    m_timer.increment(file_chunk->reads_1.size() * 2);

    return chunk_vec();
  }

  /** Prints summary of inferred consensus sequences. */
  void finalize()
  {
    m_timer.finalize();

    auto stats = m_stats.acquire();
    while (!m_stats.empty()) {
      *stats += *m_stats.acquire();
    }

    std::cout << "   Found " << stats->aligned_pairs << " overlapping pairs\n"
              << "   Of which " << stats->pairs_with_adapters
              << " contained adapter sequence(s)\n\n"
              << "Printing adapter sequences, including poly-A tails:"
              << std::endl;

    print_consensus_adapter(
      stats->pcr1_counts,
      stats->pcr1_kmers,
      "--adapter1",
      m_config.adapters.get_raw_adapters().front().first.sequence());
    std::cout << "\n\n";

    fastq adapter2 = m_config.adapters.get_raw_adapters().front().second;
    adapter2.reverse_complement();
    print_consensus_adapter(
      stats->pcr2_counts, stats->pcr2_kmers, "--adapter2", adapter2.sequence());
  }

private:
  void process_reads(const sequence_aligner& aligner,
                     adapter_stats& stats,
                     fastq& read1,
                     fastq& read2)
  {
    // Throws if read-names or mate numbering does not match
    fastq::validate_paired_reads(read1, read2, m_config.mate_separator);

    // Reverse complement to match the orientation of read1
    read2.reverse_complement();

    const auto alignment =
      aligner.align_paired_end(read1, read2, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      stats.aligned_pairs++;
      if (m_config.can_merge_alignment(alignment)) {
        if (extract_adapter_sequences(alignment, read1, read2)) {
          stats.pairs_with_adapters++;

          process_adapter(
            read1.sequence(), stats.pcr1_counts, stats.pcr1_kmers);

          read2.reverse_complement();
          process_adapter(
            read2.sequence(), stats.pcr2_counts, stats.pcr2_kmers);
        }
      }
    } else {
      stats.unaligned_pairs++;
    }
  }

  void process_adapter(const std::string& sequence,
                       nt_count_vec& counts,
                       kmer_map& kmers)
  {
    if (counts.size() < sequence.length()) {
      counts.resize(sequence.length());
    }

    for (size_t i = 0; i < std::min(counts.size(), sequence.length()); ++i) {
      counts.at(i).increment(sequence.at(i));
    }

    if (sequence.length() >= KMER_LENGTH) {
      const std::string kmer = sequence.substr(0, KMER_LENGTH);
      if (!std::count(kmer.begin(), kmer.end(), 'N')) {
        kmers.at(kmer_to_size_t(kmer)) += 1;
      }
    }
  }

  const userconfig& m_config;

  progress_timer m_timer;
  threadstate<adapter_stats> m_stats;
};

int
identify_adapter_sequences(const userconfig& config)
{
  std::cout << "Attempting to identify adapter sequences" << std::endl;

  scheduler sch;

  // Step 3: Attempt to identify adapters through pair-wise alignments
  const size_t identification_step =
    sch.add_step("identify_adapters", new adapter_identification(config));

  // Step 2: Post-process and validate FASTQ reads
  const size_t postproc_step = sch.add_step(
    "post_process_fastq",
    new post_process_fastq(config.io_encoding, identification_step));

  // Step 1: Read input file(s)
  sch.add_step("read_fastq", new read_fastq(config, postproc_step));

  return !sch.run(config.max_threads);
}
