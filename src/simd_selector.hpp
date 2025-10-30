// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "scheduler.hpp"     // for analytical_step, chunk_ptr, chunk_vec
#include "sequence_sets.hpp" // for adapter_set
#include "simd.hpp"          // for instruction_set
#include "threading.hpp"     // for threadsafe_data
#include <chrono>            // for
#include <fastq.hpp>         // for fastq

namespace adapterremoval {

class userconfig;
class sequence_aligner;

class simd_selector : public analytical_step
{
public:
  explicit simd_selector(threadsafe_data<sample_set> samples,
                         threadsafe_data<simd::instruction_set> is,
                         double mismatch_threshold,
                         uint32_t shift,
                         size_t next_step);

  chunk_vec process(chunk_ptr data) override;

protected:
  using clock = std::chrono::steady_clock;

  [[nodiscard]] clock::duration process_chunk(sequence_aligner& aligner,
                                              const fastq_chunk& chunk) const;

  //! Next step in the pipeline
  size_t m_next_step;
  //! Sample set. Required to obtain selected adapters
  threadsafe_data<sample_set> m_samples;
  //! Default or selected SIMD instruction set
  threadsafe_data<simd::instruction_set> m_simd;
  //! Mismatch threshold for alignments
  double m_mismatch_threshold;
  //! Maximum shift allowed in alignments
  uint32_t m_shift;

  /** Reads processed in total nano-seconds */
  using stats = std::pair<size_t, clock::duration>;
  //! Per candidate statistics
  std::vector<stats> m_candidates{};
  //! Number of reads processed so far
  size_t m_processed_reads = 0;
};

} // namespace adapterremoval
