// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>

#include "adapter_select.hpp" // declarations
#include "adapter_detect.hpp" // for adapter_detection_stats, adapter_detector
#include "adapters.hpp"       // for adapter_database, identified_adapter
#include "commontypes.hpp"    // for read_mate
#include "debug.hpp"          // AR_REQUIRE
#include "logging.hpp"        // for log
#include "scheduler.hpp"      // for analytical_step, chunk_ptr, chunk_vec, ...
#include "sequence_sets.hpp"  // for adapter_set
#include "userconfig.hpp"     // for userconfig
#include <cstddef>            // for size_t
#include <memory>             // for make_unique, shared_ptr
#include <string_view>        // for string_view
#include <utility>            // for move

namespace adapterremoval {

namespace {

//! Max reads/read pairs to evaluate when selecting adapters
const size_t ADAPTER_SELECT_N = 100'000;

} // namespace

////////////////////////////////////////////////////////////////////////////////
// adapter_preselector

adapter_preselector::adapter_preselector(size_t next_step)
  : analytical_step(processing_order::ordered, "adapter_preselector")
  , m_next_step(next_step)
{
}

/** Cache chunks, select adapters, and/or forward the chunks as is */
chunk_vec
adapter_preselector::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  AR_REQUIRE(!chunk->adapters);

  if (!m_done) {
    m_cached_reads += chunk->reads_1.size();
    m_done = (m_cached_reads >= ADAPTER_SELECT_N) || chunk->eof;

    chunk->adapters = std::make_unique<adapter_detection_stats>();
    chunk->last_adapter_selection_block = m_done;
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

////////////////////////////////////////////////////////////////////////////////

adapter_selector::adapter_selector(const userconfig& config, size_t next_step)
  : analytical_step(processing_order::unordered, "adapter_selector")
  , m_next_step(next_step)
{

  for (size_t i = 0; i < config.max_threads; ++i) {
    m_selectors.release(
      std::make_unique<adapter_detector>(config.samples->adapters(),
                                         config.simd,
                                         config.mismatch_threshold));
  }
}

chunk_vec
adapter_selector::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  if (chunk->adapters) {
    auto& stats = *chunk->adapters;
    auto selector = m_selectors.acquire();

    if (chunk->reads_2.empty()) {
      for (const auto& read : chunk->reads_1) {
        selector->detect_se(stats, read);
      }
    } else {
      AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size());
      auto it_1 = chunk->reads_1.begin();
      auto it_2 = chunk->reads_2.begin();
      auto end_1 = chunk->reads_1.end();

      for (; it_1 != end_1; ++it_1, ++it_2) {
        selector->detect_pe(stats, *it_1, *it_2);
      }
    }

    m_selectors.release(std::move(selector));
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

////////////////////////////////////////////////////////////////////////////////
// adapter_selector

namespace {

void
describe_best_match(const identified_adapter& match,
                    std::string_view key,
                    bool required)
{
  if (match.sequence.empty()) {
    if (required) {
      log::warn() << "Could not identify adapter sequence for " << key;
    }
  } else {
    log::info() << "Using '" << match.source << "' " << match.mate
                << " adapter for " << key << ": " << match.sequence.as_string();
  }
}

} // namespace

adapter_finalizer::adapter_finalizer(
  const userconfig& config,
  size_t next_step,
  std::shared_ptr<std::atomic_size_t> entrypoint)
  : analytical_step(processing_order::ordered, "adapter_finalizer")
  , m_next_step(next_step)
  , m_entrypoint(std::move(entrypoint))
  , m_adapters(config.samples->adapters())
  , m_detector(m_adapters, simd::instruction_set::none, 0.0)
  , m_paired_ended_mode(config.paired_ended_mode)
{
  AR_REQUIRE(m_entrypoint);
}

chunk_vec
adapter_finalizer::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  if (!chunk->adapters) {
    chunk_vec chunks;
    chunks.emplace_back(m_next_step, std::move(chunk));

    return chunks;
  }

  AR_REQUIRE(!m_done);
  AR_REQUIRE(chunk->adapters);
  m_done = chunk->last_adapter_selection_block;
  m_stats.merge(*chunk->adapters);
  m_cache.emplace_back(m_next_step, std::move(chunk));

  if (m_done) {
    // Let subsequent chunks skip this mini-pipeline
    *m_entrypoint = m_next_step;

    size_t total_reads_1 = 0;
    size_t total_reads_2 = 0;
    for (const auto& it : m_cache) {
      total_reads_1 += it.second->reads_1.size();
      total_reads_2 += it.second->reads_2.size();
    }

    AR_REQUIRE((total_reads_1 == total_reads_2) || (total_reads_2 == 0));
    auto [seq1, seq2] = m_detector.select_best(m_stats);

    const adapter_database candidates{ m_adapters };
    const auto selection = candidates.identify(seq1, seq2);

    if (selection.first.sequence.empty() && selection.second.sequence.empty()) {
      log::info() << "Could not identify adapter sequences automatically. "
                     "Trimming will be performed without pre-defined "
                     "adapter sequences, equivalent to --adapter1 '' and "
                     "--adapter2 ''";
    } else {
      describe_best_match(selection.first, "--adapter1", true);
      describe_best_match(selection.second, "--adapter2", m_paired_ended_mode);

      if (selection.first.mate != read_mate::_1 ||
          (m_paired_ended_mode && selection.second.mate != read_mate::_2)) {
        log::warn()
          << "The orientation of identified adapters does not match "
             "expectations: Expected to find forward adapter in --in-file1 and "
             "reverse adapter in --in-file2 files (if any). Input files may be "
             "swapped";
      }
    }

    return std::move(m_cache);
  } else {
    return {};
  }
}

} // namespace adapterremoval
