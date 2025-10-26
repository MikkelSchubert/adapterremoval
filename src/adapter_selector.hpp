// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "adapter_database.hpp" // for adapter_database
#include "adapter_detector.hpp" // for adapter_detector_stats
#include "commontypes.hpp"      // for adapter_fallback
#include "scheduler.hpp"        // for analytical_step
#include "sequence_sets.hpp"    // for adapter_set
#include "simd.hpp"             // for instruction_set
#include "threading.hpp"        // for threadstate, threadsafe_data

namespace adapterremoval {

class adapter_database;
class fastq;

/** Chunk containing adapter detection stats in addition to fastq reads */
class adapter_chunk : public analytical_chunk
{
public:
  adapter_chunk() = default;
  ~adapter_chunk() override = default;

  //! The chunk of FASTQ reads being processed
  fastq_chunk_ptr data{};
  //! (Optional) adapter detection stats
  adapter_detection_stats adapters{};
  //! Indicates if this is the the last chunk to run detection on
  bool last_adapter_selection_block = false;

  adapter_chunk(const adapter_chunk&) = delete;
  adapter_chunk(adapter_chunk&&) = delete;
  adapter_chunk& operator=(const adapter_chunk&) = delete;
  adapter_chunk& operator=(adapter_chunk&&) = delete;
};

/**
 * This class wraps fastq_chunks in adapter_chunks until a pre-determined
 * number of reads have been selected, after which chunks are forwarded as is.
 * This step is sequential to ensure that the selection of fastq_chunks is
 * deterministic, as the selection of chunks would otherwise depend on the order
 * in which pre-processing completed.
 */
class adapter_preselector : public analytical_step
{
public:
  explicit adapter_preselector(size_t next_step);
  ~adapter_preselector() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr data) override;

  /** Finalizer; checks that pre-selection is done */
  void finalize() override;

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

/**
 * Central worker in the adapter selection workflow; performs adapter detection
 * on adapter_chunks using the specified adapter database. fastq_chunks are
 * forwarded as is.
 */
class adapter_selector : public analytical_step
{
public:
  adapter_selector(const adapter_database& database,
                   simd::instruction_set is,
                   double mismatch_threshold,
                   size_t next_step,
                   size_t max_threads = 1);
  ~adapter_selector() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr data) override;

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

/**
 * Aggregates the results from adapter_selector, and caches the associated
 * fastq_chunks, until the final adapter_chunk has been received. At this point,
 * the best matching adapter (pair) is selected and assigned to the sample_set.
 * Subsequent fastq_chunks are forwarded as is.
 */
class adapter_finalizer : public analytical_step
{
public:
  adapter_finalizer(adapter_database database,
                    threadsafe_data<sample_set> sample,
                    adapter_fallback fallback,
                    size_t next_step);
  ~adapter_finalizer() override = default;

  /** Cache chunks, select adapters, and/or forward the chunks as is */
  chunk_vec process(chunk_ptr data) override;

  /** Finalizer; checks that selection is done */
  void finalize() override;

  adapter_finalizer(const adapter_finalizer&) = delete;
  adapter_finalizer(adapter_finalizer&&) = delete;
  adapter_finalizer& operator=(const adapter_finalizer&) = delete;
  adapter_finalizer& operator=(adapter_finalizer&&) = delete;

private:
  //! Next step in the pipeline, either demultiplexing or trimming
  const size_t m_next_step;
  //! Fallback strategy if no adapters could be identified
  const adapter_fallback m_fallback;
  //! Set of adapters to potentially detect
  adapter_database m_database;
  //! Sampleset for which adapters should be configured
  threadsafe_data<sample_set> m_samples;
  //! Cached reads held back until adapter selection is complete
  chunk_vec m_cache{};
  //! Total number of matches for each adapter candidate for mate 1/2 reads
  adapter_detection_stats m_stats{};
  //! Has selection been performed and should the step be bypassed
  bool m_done = false;
};

} // namespace adapterremoval
