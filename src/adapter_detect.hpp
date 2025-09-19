// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "alignment.hpp"     // for sequence_aligner
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for adapter_set
#include <vector>            // for vector

namespace adapterremoval {

class fastq;
class userconfig;

/**
 * Class used to contain adapter detection statistics. These statistics are
 * external to the detector to enable parallelization
 **/
struct adapter_detection_stats
{
  adapter_detection_stats() = default;

  /** Merge results in `other` into this object */
  void merge(const adapter_detection_stats& other);

  struct stats
  {
    size_t hits = 0;
    size_t length = 0;
    size_t mismatches = 0;
  };

  using values = std::vector<stats>;

private:
  friend class adapter_detector;

  //! The number of reads processed
  size_t reads = 0;
  //! Matches for each adapter candidate for mate 1 reads
  values mate_1{};
  //! Matches for each adapter candidate for mate 2 reads
  values mate_2{};
};

class adapter_detector
{
public:
  /** Create selector based on known adapters */
  adapter_detector(simd::instruction_set is, double mismatch_threshold);
  /** Create selector based on known and user-specified adapters */
  explicit adapter_detector(const adapter_set& user_sequences,
                            simd::instruction_set is,
                            double mismatch_threshold);

  /** Detect adapters based only on mate 1 reads, writing results to `stats` */
  void detect_se(adapter_detection_stats& stats, const fastq& read);
  /** Detect adapters based paired reads, writing results to `stats` */
  void detect_pe(adapter_detection_stats& stats,
                 const fastq& read1,
                 const fastq& read2);

  /** Selects the best candidate for the provided stats */
  sequence_pair select_best(adapter_detection_stats& stats) const;

private:
  //! Detect adapter and record the statistics into the provided object
  void detect_adapters(const fastq& read,
                       adapter_detection_stats::values& stats);

  //! Full set of known and user-provided adapter sequences
  adapter_set m_adapters;
  //! Aligner for full set of forward and reverse adapter sequences
  sequence_aligner m_aligner;
  //! Overview of adapter sequence overlap with sequences before/after
  std::vector<std::vector<std::pair<size_t, size_t>>> m_common_prefixes{};
};

} // namespace adapterremoval
