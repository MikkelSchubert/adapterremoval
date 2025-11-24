// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_selector.hpp" // declarations
#include "adapter_database.hpp" // for adapter_database, identified_adapter
#include "adapter_detector.hpp" // for adapter_detection_stats, adapter_detector
#include "alignment.hpp"        // for sequence_aligner
#include "debug.hpp"            // AR_REQUIRE
#include "errors.hpp"           // for fatal_error
#include "logging.hpp"          // for log
#include "scheduler.hpp"        // for analytical_step, chunk_ptr, ...
#include "sequence.hpp"         // for dna_sequence
#include "sequence_sets.hpp"    // for adapter_set
#include "simd.hpp"             // for instruction_set
#include "threading.hpp"        // for threadsafe_data
#include "userconfig.hpp"       // for userconfig
#include "utilities.hpp"        // for dynamic_cast_unique
#include <cstddef>              // for size_t
#include <memory>               // for make_unique, shared_ptr
#include <string>               // for string
#include <string_view>          // for string_view
#include <utility>              // for move
#include <vector>               // for vector

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
adapter_preselector::process(chunk_ptr data)
{
  chunk_vec chunks;
  if (m_done) {
    chunks.emplace_back(m_next_step, std::move(data));
  } else {
    auto chunk = std::make_unique<adapter_chunk>();
    chunk->data = dynamic_cast_unique<fastq_chunk>(data);
    AR_REQUIRE(chunk->data);

    m_cached_reads += chunk->data->reads_1.size();
    m_done = (m_cached_reads >= ADAPTER_SELECT_N) || chunk->data->eof;

    chunk->last_adapter_selection_block = m_done;
    chunks.emplace_back(m_next_step, std::move(chunk));
  }

  return chunks;
}

void
adapter_preselector::finalize()
{
  AR_REQUIRE(m_done);
}

////////////////////////////////////////////////////////////////////////////////

adapter_selector::adapter_selector(const adapter_database& database,
                                   simd::instruction_set is,
                                   double mismatch_threshold,
                                   size_t next_step,
                                   size_t max_threads)
  : analytical_step(processing_order::unordered, "adapter_selector")
  , m_next_step(next_step)
{
  m_selectors.emplace_back_n(max_threads, database, is, mismatch_threshold);
}

chunk_vec
adapter_selector::process(chunk_ptr data)
{
  chunk_vec chunks;
  if (auto chunk = dynamic_cast_unique<adapter_chunk>(data)) {
    AR_REQUIRE(chunk->data);

    auto selector = m_selectors.acquire();
    auto& stats = chunk->adapters;
    if (chunk->data->reads_2.empty()) {
      for (const auto& read : chunk->data->reads_1) {
        selector->detect_se(stats, read);
      }
    } else {
      AR_REQUIRE(chunk->data->reads_1.size() == chunk->data->reads_2.size());
      auto it_1 = chunk->data->reads_1.begin();
      auto it_2 = chunk->data->reads_2.begin();
      auto end_1 = chunk->data->reads_1.end();

      for (; it_1 != end_1; ++it_1, ++it_2) {
        selector->detect_pe(stats, *it_1, *it_2);
      }
    }

    m_selectors.release(std::move(selector));

    chunks.emplace_back(m_next_step, std::move(chunk));
  } else if (auto fq_chunk = dynamic_cast_unique<fastq_chunk>(data)) {
    chunks.emplace_back(m_next_step, std::move(fq_chunk));
  } else {
    AR_FAIL("invalid data chunk");
  }

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
      log::warn() << "Could not identify " << key << " sequence";
    }
  } else {
    log::info() << "Using '" << match.name << "' " << match.mate
                << " adapter for " << key << ": " << match.sequence.as_string();
  }
}

} // namespace

adapter_finalizer::adapter_finalizer(adapter_database database,
                                     threadsafe_data<sample_set> samples,
                                     const adapter_fallback fallback,
                                     const size_t next_step)
  : analytical_step(processing_order::ordered, "adapter_finalizer")
  , m_next_step(next_step)
  , m_fallback(fallback)
  , m_database(std::move(database))
  , m_samples(std::move(samples))
{
}

chunk_vec
adapter_finalizer::process(chunk_ptr data)
{
  chunk_vec chunks;
  if (auto chunk = dynamic_cast_unique<adapter_chunk>(data)) {
    AR_REQUIRE(!m_done);
    AR_REQUIRE(chunk->data);

    m_stats.merge(chunk->adapters);
    m_cache.emplace_back(m_next_step, std::move(chunk->data));

    if (chunk->last_adapter_selection_block) {
      m_done = true;

      if (m_stats.reads_1() == 0 && m_stats.reads_2() == 0) {
        // Adapters must be set to *something*, so later steps can use them, but
        // there's no point in printing diagnostics if there is no data
        m_samples.get_writer()->set_adapters(adapter_set{});
      } else {
        const adapter_detector detector{ m_database,
                                         simd::instruction_set::none,
                                         0.0 };
        const auto [seq_1, seq_2] = detector.select_best(m_stats);

        if (seq_1.sequence.empty() && seq_2.sequence.empty()) {
          switch (m_fallback) {
            case adapter_fallback::abort: {
              throw fatal_error(
                "Could not select adapter sequences automatically; aborting "
                "since '--adapter-fallback' is set to 'abort'. To proceed, "
                "either select a different fallback strategy, or manually "
                "specify adapter sequences using --adapter1 / --adapter2");
            }
            case adapter_fallback::undefined: {
              log::warn()
                << "Could not select adapter sequences automatically; falling "
                   "back to trimming putative adapter sequences based on "
                   "overlap between read pairs, equivalent to `--adapter1 ''` "
                   "and `--adapter2 ''";

              m_samples.get_writer()->set_adapters(adapter_set{ { "", "" } });
              break;
            }
            case adapter_fallback::none: {
              log::warn()
                << "Could not select adapter sequences automatically; falling "
                   "back to assuming that reads do not contain adapter "
                   "sequences";

              m_samples.get_writer()->set_adapters(adapter_set{});
              break;
            }
            default:                                     // GCOVR_EXCL_LINE
              AR_FAIL("invalid adapter_fallback value"); // GCOVR_EXCL_LINE
          }
        } else {
          describe_best_match(seq_1, "--adapter1", true);
          describe_best_match(seq_2, "--adapter2", m_stats.reads_2());

          if (seq_1.mate != read_mate::_1 ||
              (m_stats.reads_2() && seq_2.mate != read_mate::_2)) {
            log::warn() << "The orientation of identified adapters does not "
                           "match expectations: Expected to find forward "
                           "adapter in --in-file1 and reverse adapter in "
                           "--in-file2 files (if any). Input files may be "
                           "swapped";
          }

          adapter_set adapters{ { seq_1.sequence.as_string(),
                                  seq_2.sequence.as_string() } };
          m_samples.get_writer()->set_adapters(std::move(adapters));
        }
      }

      std::swap(m_cache, chunks);
    }
  } else if (auto fq_chunk = dynamic_cast_unique<fastq_chunk>(data)) {
    chunks.emplace_back(m_next_step, std::move(fq_chunk));
  } else {
    AR_FAIL("invalid data chunk");
  }

  return chunks;
}

void
adapter_finalizer::finalize()
{
  AR_REQUIRE(m_done);
  AR_REQUIRE(m_cache.empty());
}

} // namespace adapterremoval
