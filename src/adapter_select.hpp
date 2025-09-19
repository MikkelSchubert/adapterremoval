// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "adapter_detect.hpp" // for adapter_detector
#include "scheduler.hpp"      // for analytical_step, threadstate
#include "sequence_sets.hpp"
#include <atomic> // for atomic_size_t

namespace adapterremoval {

class fastq;
class userconfig;

/**
 * This class forwards chunks to the until a pre-determined of reads have been
 * processed, after which forwarded chunks are flagged to prevent further
 * processing. This is required due to the processing step running in parallel.
 */
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

class adapter_selector : public analytical_step
{
public:
  adapter_selector(const userconfig& config, size_t next_step);
  ~adapter_selector() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr chunk) override;

  adapter_selector(const adapter_selector&) = delete;
  adapter_selector(adapter_selector&&) = delete;
  adapter_selector& operator=(const adapter_selector&) = delete;
  adapter_selector& operator=(adapter_selector&&) = delete;

private:
  //! Next step in the pipeline, either demultiplexing or trimming
  const size_t m_next_step;
  //! Cached selectors to take advantage of frequency optimizations
  threadstate<adapter_detector> m_selectors{};
};

class adapter_finalizer : public analytical_step
{
public:
  adapter_finalizer(const userconfig& config,
                    size_t next_step,
                    std::shared_ptr<std::atomic_size_t> entrypoint);
  ~adapter_finalizer() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr chunk) override;

  adapter_finalizer(const adapter_finalizer&) = delete;
  adapter_finalizer(adapter_finalizer&&) = delete;
  adapter_finalizer& operator=(const adapter_finalizer&) = delete;
  adapter_finalizer& operator=(adapter_finalizer&&) = delete;

private:
  //! Next step in the pipeline, either demultiplexing or trimming
  const size_t m_next_step;
  //! Step number containing the entry-point of the detector mini-pipeline
  std::shared_ptr<std::atomic_size_t> m_entrypoint{};
  //! Cached reads to be used for adapter selection
  chunk_vec m_cache{};
  //! User provided adapter sequences for detection
  adapter_set m_adapters;
  //! Detector used to pick the best match once results have been merged
  adapter_detector m_detector;
  //! Total number of matches for each adapter candidate for mate 1/2 reads
  adapter_detection_stats m_stats{};
  //! Is adapterremoval running in PE mode
  bool m_paired_ended_mode = false;
  //! Has selection been performed and should the step be bypassed
  bool m_done = false;
};

} // namespace adapterremoval
