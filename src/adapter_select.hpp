// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "alignment.hpp"     // for sequence_aligner
#include "scheduler.hpp"     // for analytical_step, threadstate
#include "sequence_sets.hpp" // for adapter_set
#include <atomic>            // for atomic_size_t
#include <vector>

namespace adapterremoval {

class userconfig;

class adapter_preselector : public analytical_step
{
public:
  explicit adapter_preselector(size_t next_step);
  ~adapter_preselector() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr chunk) override;

  adapter_preselector(const adapter_preselector&) = delete;
  adapter_preselector(adapter_preselector&&) = delete;
  adapter_preselector& operator=(const adapter_preselector&) = delete;
  adapter_preselector& operator=(adapter_preselector&&) = delete;

private:
  //! Next step in the pipeline, expected to be adapter_detector
  const size_t m_next_step;
  //! Total number of cached reads
  size_t m_cached_reads = 0;
  //! Has selection been performed and should the step be bypassed
  bool m_done = false;
};

class adapter_detector : public analytical_step

{
public:
  adapter_detector(const userconfig& config, size_t next_step);
  ~adapter_detector() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr chunk) override;

  adapter_detector(const adapter_detector&) = delete;
  adapter_detector(adapter_detector&&) = delete;
  adapter_detector& operator=(const adapter_detector&) = delete;
  adapter_detector& operator=(adapter_detector&&) = delete;

private:
  void process_reads(const std::vector<fastq>& reads,
                     sequence_aligner& aligner,
                     adapter_selection_values& stats) const;

  const userconfig& m_config;
  //! Next step in the pipeline, either demultiplexing or trimming
  const size_t m_next_step;
  //! Full set of forward and reverse adapter sequences
  adapter_set m_adapters{};
  //! Cached aligners to take advantage of frequency optimizations
  threadstate<sequence_aligner> m_aligners{};
  //! Overview of adapter sequence overlap with sequences before/after
  std::vector<std::vector<std::pair<size_t, size_t>>> m_common_prefixes{};
};

class adapter_selector : public analytical_step
{
public:
  adapter_selector(const userconfig& config,
                   size_t next_step,
                   std::shared_ptr<std::atomic_size_t> entrypoint);
  ~adapter_selector() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr chunk) override;

  adapter_selector(const adapter_selector&) = delete;
  adapter_selector(adapter_selector&&) = delete;
  adapter_selector& operator=(const adapter_selector&) = delete;
  adapter_selector& operator=(adapter_selector&&) = delete;

private:
  const userconfig& m_config;
  //! Next step in the pipeline, either demultiplexing or trimming
  const size_t m_next_step;
  //! Step number containing the entry-point of the detector mini-pipeline
  std::shared_ptr<std::atomic_size_t> m_entrypoint{};
  //! Cached reads to be used for adapter selection
  chunk_vec m_cache{};

  //! Total number of matches for each adapter candidate for mate 1 reads
  adapter_selection_values m_mate_1{};
  //! Total number of matches for each adapter candidate for mate 2 reads
  adapter_selection_values m_mate_2{};
  //! Has selection been performed and should the step be bypassed
  bool m_done = false;
};

} // namespace adapterremoval
