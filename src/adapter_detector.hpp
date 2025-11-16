// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "adapter_database.hpp" // for adapter_database
#include "alignment.hpp"        // for sequence_aligner
#include "sequence.hpp"         // for dna_sequence
#include <vector>               // for vector

namespace adapterremoval {

class fastq;
class userconfig;
class adapter_database;

//! Minimum overlap required for adapter fragments
const size_t ADAPTER_DETECT_MIN_OVERLAP = 8;

/**
 * Class used to record adapter detection statistics; these statistics are
 * external to the detector to enable parallelization
 **/
class adapter_detection_stats
{
public:
  struct hits
  {
    //! The number of times an adapter was among the set that aligned the best
    size_t hits = 0;
    //! The total number of aligned bases: Excluding Ns but including mismatches
    size_t aligned = 0;
    //! The total number of mismatches
    size_t mismatches = 0;
  };

  /** Creates zero'd stats object */
  adapter_detection_stats() = default;
  /** Create stats object with the specified values */
  adapter_detection_stats(size_t reads,
                          std::vector<hits> mate_1,
                          std::vector<hits> mate_2 = {});

  /** Merge results in `other` into this object */
  void merge(const adapter_detection_stats& other);

  //! Per adapter candidate statistics
  using values = std::vector<hits>;

  /** Returns the number of reads 1 sequences processed */
  [[nodiscard]] size_t reads_1() const noexcept { return m_reads_1; }

  /** Returns the number of reads 2 sequences processed */
  [[nodiscard]] size_t reads_2() const noexcept { return m_reads_2; }

  /** Returns overlapping mate 1 hits for all sequences for a detector */
  [[nodiscard]] const values& mate_1() const noexcept { return m_mate_1; }

  /** Returns overlapping mate 2 hits for all sequences for a detector */
  [[nodiscard]] const values& mate_2() const noexcept { return m_mate_2; }

private:
  /** Only adapter_detector needs write access to the statistics */
  friend class adapter_detector;

  /** Helper function for use with utilities.hpp:merge */
  friend void merge(hits&, const hits&);

  /** Returns true if all fields match */
  friend bool operator==(const hits&, const hits&);

  /** Stream operator for debugging output */
  friend std::ostream& operator<<(std::ostream& os, const hits&);

  //! The number of read 1 sequences processed
  size_t m_reads_1 = 0;
  //! The number of read 2 sequences processed
  size_t m_reads_2 = 0;
  //! Match statistics for each adapter candidate for mate 1 reads
  values m_mate_1{};
  //! Match statistics for each adapter candidate for mate 2 reads
  values m_mate_2{};
};

/**
 * Class for identifying possible adapter sequences in untrimmed reads. No
 * assumptions are made about which adapters are found in which reads, in case
 * of swapped mate 1 / 2 files or swapped sequences.
 */
class adapter_detector
{
public:
  /** Create selector based on known and user-specified adapters */
  adapter_detector(adapter_database database,
                   simd::instruction_set is,
                   double mismatch_threshold);

  /** Detect adapters based only on mate 1 reads, writing results to `stats` */
  void detect_se(adapter_detection_stats& stats, const fastq& read);
  /** Detect adapters based paired reads, writing results to `stats` */
  void detect_pe(adapter_detection_stats& stats,
                 const fastq& read_1,
                 const fastq& read_2);

  /** Selects the best candidate(s) for the provided stats */
  [[nodiscard]] identified_adapter_pair select_best(
    const adapter_detection_stats& stats) const;

  /** Returns the number of unique adapter sequences */
  [[nodiscard]] size_t size() const noexcept { return m_adapters.size(); }

  /** Returns unique adapters sequences, all stored as read 1 adapters */
  [[nodiscard]] const sequence_vec& sequences() const noexcept
  {
    return m_adapters;
  }

private:
  //! Detect adapter and record the statistics into the provided object
  void detect_adapters(const fastq& read,
                       adapter_detection_stats::values& stats);

  //! Database of known and/or user-provided adapter sequences
  adapter_database m_database;
  //! Unique known and user-provided adapters, all stored as read 1 adapters
  sequence_vec m_adapters;
  //! Aligner for full set of adapter sequences, all stored as read 1 adapters
  sequence_aligner m_aligner;
  //! Adapter sequence overlap with adjacent sequences; each sub-vector contains
  // counts of sequences before/after, with lengths - min adapter overlap >= 0
  std::vector<std::vector<std::pair<size_t, size_t>>> m_common_prefixes{};
};

} // namespace adapterremoval
