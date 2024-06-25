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
#include "adapter_id.hpp" // for adapter_id_stats
#include "alignment.hpp"  // for extract_adapter_sequences, sequence_aligner
#include "counts.hpp"     // for indexed_counts
#include "debug.hpp"      // for AR_REQUIRE
#include "fastq.hpp"      // for ACGTN, fastq, fastq_pair_vec, ACGT, ACGT:...
#include "fastq_io.hpp"   // for read_fastq, fastq_read_chunk
#include "scheduler.hpp"  // for threadstate, scheduler, analytical_step
#include "userconfig.hpp" // for userconfig
#include <cstddef>        // for size_t
#include <iomanip>        // for operator<<, setw
#include <iostream>       // for operator<<, basic_ostream, ostringstream
#include <string>         // for string, operator<<, char_traits

namespace adapterremoval {

//! The N most common kmers to print
const size_t TOP_N_KMERS = 5;

/** Prints the N top kmers in a kmer_map, including sequence and frequency. */
void
print_most_common_kmers(const consensus_adapter& adapter)
{
  std::cout.precision(2);
  std::cout << std::fixed;
  std::cout << "    Top 5 most common " << adapter.top_kmers.size()
            << "-bp 5'-kmers:\n";

  for (size_t i = 0; i < adapter.top_kmers.size(); ++i) {
    const auto& count = adapter.top_kmers.at(i);

    std::cout << std::string(12, ' ');
    std::cout << i + 1 << ": " << count.first << " = " << std::right
              << std::setw(5) << (100.0 * count.second) / adapter.total_kmers
              << "% (" << count.second << ")\n";
  }
}

/**
 * Build representation of identity between two sequences.
 *
 * The resulting string represents N with wildcards ('*'), matching bases with
 * pipes ('|') and mismatches with spaces (' '). Only overlapping bases are
 * compared.
 */
std::string
compare_consensus_with_ref(const std::string& a, const std::string& b)
{
  std::ostringstream identity;
  for (size_t i = 0, size = std::min(b.size(), a.size()); i < size; ++i) {
    if (a.at(i) == 'N' || b.at(i) == 'N') {
      identity << '*';
    } else {
      identity << (a.at(i) == b.at(i) ? '|' : ' ');
    }
  }

  return identity.str();
}

/**
 * Prints description of consensus adapter sequence.
 *
 * @param nt_counts Observed nucleotide frequencies.
 * @param kmers Observed kmer frequencies.
 * @param name Argument name for adapter (--adapter1 / adapter2)
 * @param ref Default sequence for the inferred adapter
 */
void
print_consensus_adapter(const consensus_adapter_stats& stats,
                        const std::string& name,
                        const std::string& ref)
{
  const auto consensus = stats.summarize(TOP_N_KMERS);
  const auto& adapter = consensus.adapter;

  const std::string identity =
    compare_consensus_with_ref(ref, adapter.sequence());

  std::cout << "  " << name << ":  " << ref << "\n"
            << "               " << identity << "\n"
            << "   Consensus:  " << adapter.sequence() << "\n"
            << "     Quality:  " << adapter.qualities() << "\n\n";

  print_most_common_kmers(consensus);
}

///////////////////////////////////////////////////////////////////////////////
// Threaded adapter identification step

class adapter_identification : public analytical_step
{
public:
  explicit adapter_identification(const userconfig& config)
    : analytical_step(processing_order::unordered, "adapter_identification")
    , m_config(config)
    , m_stats()
  {
    for (size_t i = 0; i < m_config.max_threads; ++i) {
      m_stats.emplace_back();
    }
  }

  ~adapter_identification() override = default;

  chunk_vec process(chunk_ptr chunk) override
  {
    AR_REQUIRE(chunk);
    auto& file_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);

    const fastq empty_adapter("dummy", "", "");
    fastq_pair_vec adapters;
    adapters.emplace_back(empty_adapter, empty_adapter);

    const auto aligner = sequence_aligner(adapters, m_config.simd);
    auto stats = m_stats.acquire();

    AR_REQUIRE(file_chunk.reads_1.size() == file_chunk.reads_2.size());
    auto read_1 = file_chunk.reads_1.begin();
    auto read_2 = file_chunk.reads_2.begin();

    while (read_1 != file_chunk.reads_1.end()) {
      process_reads(aligner, *stats, *read_1++, *read_2++);
    }

    m_stats.release(stats);

    return chunk_vec();
  }

  /** Prints summary of inferred consensus sequences. */
  void finalize() override
  {
    auto stats = m_stats.acquire();
    while (auto next = m_stats.try_acquire()) {
      *stats += *next;
    }

    std::cout << "   Found " << stats->aligned_pairs << " overlapping pairs\n"
              << "   Of which " << stats->pairs_with_adapters
              << " contained adapter sequence(s)\n\n"
              << "Printing adapter sequences, including poly-A tails:"
              << std::endl;

    auto reference = m_config.adapters.get_raw_adapters().front();
    // Revert to user-supplied orientation
    reference.second.reverse_complement();

    print_consensus_adapter(
      stats->adapter1, "--adapter1", reference.first.sequence());
    std::cout << "\n\n";
    print_consensus_adapter(
      stats->adapter2, "--adapter2", reference.second.sequence());
  }

private:
  void process_reads(const sequence_aligner& aligner,
                     adapter_id_stats& stats,
                     fastq& read1,
                     fastq& read2) const
  {
    read1.post_process(m_config.io_encoding);
    read2.post_process(m_config.io_encoding);

    // Reverse complement to match the orientation of read1
    read2.reverse_complement();

    const auto alignment =
      aligner.align_paired_end(read1, read2, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      stats.aligned_pairs++;
      if (m_config.can_merge_alignment(alignment) &&
          extract_adapter_sequences(alignment, read1, read2)) {
        stats.pairs_with_adapters++;

        stats.adapter1.process(read1.sequence());
        read2.reverse_complement();
        stats.adapter2.process(read2.sequence());
      }
    } else {
      stats.unaligned_pairs++;
    }
  }

  const userconfig& m_config;

  threadstate<adapter_id_stats> m_stats;
};

int
identify_adapter_sequences(const userconfig& config)
{
  scheduler sch;

  // Step 2: Attempt to identify adapters through pair-wise alignments
  const size_t id_step = sch.add<adapter_identification>(config);

  // Step 1: Read input file(s)
  sch.add<read_fastq>(config, id_step);

  return !sch.run(config.max_threads);
}

} // namespace adapterremoval
