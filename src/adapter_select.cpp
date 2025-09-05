// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>

#include "adapter_select.hpp" // declarations
#include "adapters.hpp"       // for known_adapters
#include "alignment.hpp"      // for sequence_aligner
#include "commontypes.hpp"
#include "debug.hpp" // AR_REQUIRE
#include "logging.hpp"
#include "scheduler.hpp"
#include "sequence.hpp" // for dna_sequence
#include "sequence_sets.hpp"
#include "strutils.hpp"
#include "userconfig.hpp" // for userconfig
#include "utilities.hpp"  // for merge
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility> // for sort
#include <vector>

namespace adapterremoval {

namespace {

//! Max reads/read pairs to evaluate when selecting adapters
const size_t ADAPTER_SELECT_N = 100'000;
//! Minimum overlap required for adapter fragments
const size_t ADAPTER_SELECT_MIN_OVERLAP = 8;

const double SELECTION_THRESHOLD = 10;
const double SELECTION_THRESHOLD_LENIENT = 1 / 100.0;
const double SELECTION_THRESHOLD_STRICT = 1 / 1000.0;
const double SELECTION_THRESHOLD_STRICT_MISMATCHES = 1 / 10.0;

void
make_unique(std::vector<dna_sequence>& adapters)
{
  std::sort(adapters.begin(), adapters.end());

  size_t last_unique = 0;
  for (size_t i = last_unique + 1; i < adapters.size(); ++i) {
    if (!(adapters.at(last_unique) == adapters.at(i))) {
      last_unique++;

      if (last_unique != i) {
        std::swap(adapters.at(last_unique), adapters.at(i));
      }
    }
  }

  if (adapters.size() > last_unique + 1) {
    adapters.resize(last_unique + 1);
  }
}

adapter_set
collect_adapters(const userconfig& config)
{
  // Database of known and user user-provided sequences
  adapter_database adapters{ config.samples->uninitialized_adapters() };

  std::vector<dna_sequence> sequences;
  for (const auto& it : adapters) {
    sequences.insert(sequences.end(),
                     it.adapter_1().begin(),
                     it.adapter_1().end());
    sequences.insert(sequences.end(),
                     it.adapter_2().begin(),
                     it.adapter_2().end());
  }

  make_unique(sequences);

  adapter_set selection;
  for (const auto& it : sequences) {
    if (it.length() >= ADAPTER_SELECT_MIN_OVERLAP) {
      selection.add(it, {});
    } else {
      log::warn() << "Adapter " << log_escape(it.as_string())
                  << " is too short for auto-detection; adapters must be at "
                  << ADAPTER_SELECT_MIN_OVERLAP << " bp long";
    }
  }

  return selection;
}

std::vector<std::vector<std::pair<size_t, size_t>>>
collect_adapter_prefixes(adapter_set adapters, size_t min_overlap)
{
  std::vector<std::vector<std::pair<size_t, size_t>>> prefixes;
  prefixes.reserve(adapters.size());

  for (size_t i = 0; i < adapters.size(); ++i) {
    std::string_view adapter = adapters.at(i).first.as_string();
    size_t max_overlapping_before = 0;
    size_t max_overlapping_after = 0;

    std::vector<std::pair<size_t, size_t>> overlaps;
    overlaps.reserve(adapter.length() -
                     std::min(min_overlap, adapter.length()));

    for (size_t len = min_overlap; len <= adapter.length(); ++len) {
      auto prefix = adapter.substr(0, len);

      size_t overlapping_before = 0;
      for (size_t j = i; j; --j) {
        if (starts_with(adapters.at(j - 1).first.as_string(), prefix)) {
          overlapping_before++;
        } else {
          break;
        }
      }

      size_t overlapping_after = 0;
      for (size_t j = i + 1; j < adapters.size(); ++j) {
        if (starts_with(adapters.at(j).first.as_string(), prefix)) {
          overlapping_after++;
        } else {
          break;
        }
      }

      max_overlapping_before =
        std::max(max_overlapping_before, overlapping_before);
      max_overlapping_after =
        std::max(max_overlapping_after, overlapping_after);

      overlaps.emplace_back(overlapping_before, overlapping_after);
    }

    prefixes.emplace_back(std::move(overlaps));

    log::debug() << "Adapter " << adapter << " overlaps "
                 << max_overlapping_before << " sequences before, "
                 << max_overlapping_after << " after";
  }

  return prefixes;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

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

    chunk->adapters =
      std::make_unique<analytical_chunk::adapter_selection_stats>();
    chunk->adapters->is_last = m_done;
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

////////////////////////////////////////////////////////////////////////////////

namespace {

std::unique_ptr<sequence_aligner>
create_aligner(const userconfig& config, const adapter_set& adapters)
{
  auto aligner = std::make_unique<sequence_aligner>(adapters,
                                                    config.simd,
                                                    config.mismatch_threshold);
  aligner->set_min_se_overlap(ADAPTER_SELECT_MIN_OVERLAP);
  return aligner;
}

} // namespace

adapter_detector::adapter_detector(const userconfig& config, size_t next_step)
  : analytical_step(processing_order::unordered, "adapter_detector")
  , m_config(config)
  , m_next_step(next_step)
  , m_adapters(collect_adapters(config))
{
  m_common_prefixes =
    collect_adapter_prefixes(m_adapters, ADAPTER_SELECT_MIN_OVERLAP);

  for (size_t i = 0; i < config.max_threads; ++i) {
    m_aligners.release(create_aligner(config, m_adapters));
  }
}

chunk_vec
adapter_detector::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  if (chunk->adapters) {
    auto& adapters = *chunk->adapters;
    auto aligner = m_aligners.acquire();

    process_reads(chunk->reads_1, *aligner, adapters.mate_1);
    process_reads(chunk->reads_2, *aligner, adapters.mate_2);

    m_aligners.release(std::move(aligner));
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

void
adapter_detector::process_reads(const std::vector<fastq>& reads,
                                sequence_aligner& aligner,
                                adapter_selection_values& stats) const
{
  stats.resize(aligner.size());
  for (const auto& read : reads) {
    const auto alignment = aligner.align_single_end(read, m_config.shift);

    if (alignment.type() >= alignment_type::good) {
      auto idx = alignment.adapter_id();
      auto idx_end = idx;

      const auto adapter_len = m_adapters.at(idx).first.length();
      if (alignment.offset() + adapter_len >= read.length()) {
        auto length = alignment.length() - ADAPTER_SELECT_MIN_OVERLAP;

        // idx -= m_common_prefixes.at(idx).at(length).first;
        idx_end += m_common_prefixes.at(idx).at(length).second;
      }

      for (; idx <= idx_end; ++idx) {
        stats.at(idx).inc(alignment.length(), alignment.n_mismatches());
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// adapter_selector

namespace {

dna_sequence
select_best_match(const read_mate mate,
                  const adapter_set& adapters,
                  const adapter_selection_values& values,
                  const size_t total_reads)
{
  auto best_idx = std::numeric_limits<size_t>::max();
  adapter_selection_values::value_type best;

  if (total_reads) {
    for (size_t i = 0; i < adapters.size(); ++i) {
      const auto& current = values.at(i);

      if (current.hits > best.hits ||
          (current.hits == best.hits && (current.length - current.mismatches) >
                                          (best.length - best.mismatches))) {
        best = current;
        best_idx = i;
      }

      if (current.hits) {
        log::debug() << mate << " adapter has " << current.hits << " hits, "
                     << current.mismatches << "/" << current.length
                     << " mismatches: " << adapters.at(i).first.as_string();
      }
    }

    if (best.hits >= SELECTION_THRESHOLD &&
        (best.hits >= total_reads * SELECTION_THRESHOLD_LENIENT ||
         ((best.hits >= total_reads * SELECTION_THRESHOLD_STRICT) &&
          (best.length * SELECTION_THRESHOLD_STRICT_MISMATCHES >=
           best.mismatches)))) {
      log::debug() << "Selected " << adapters.at(best_idx).first;
      return adapters.at(best_idx).first;
    }
  }

  return {};
}

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

void
merge(adapter_selection_values::value_type& dst,
      const adapter_selection_values::value_type& src)
{
  dst.hits += src.hits;
  dst.length += src.length;
  dst.mismatches += src.mismatches;
}

adapter_selector::adapter_selector(
  const userconfig& config,
  size_t next_step,
  std::shared_ptr<std::atomic_size_t> entrypoint)
  : analytical_step(processing_order::ordered, "adapter_selector")
  , m_config(config)
  , m_next_step(next_step)
  , m_entrypoint(std::move(entrypoint))
{
  AR_REQUIRE(m_entrypoint);
}

chunk_vec
adapter_selector::process(chunk_ptr chunk)
{
  AR_REQUIRE(chunk);
  if (!chunk->adapters) {
    chunk_vec chunks;
    chunks.emplace_back(m_next_step, std::move(chunk));

    return chunks;
  }

  AR_REQUIRE(!m_done);
  AR_REQUIRE(chunk->adapters);
  m_done = chunk->adapters->is_last;

  merge(m_mate_1, chunk->adapters->mate_1);
  merge(m_mate_2, chunk->adapters->mate_2);

  m_cache.emplace_back(m_next_step, std::move(chunk));

  if (m_done) {
    // Let subsequent chunks skip this mini-pipeline
    *m_entrypoint = m_next_step;

    const auto adapters = collect_adapters(m_config);

    size_t total_reads_1 = 0;
    size_t total_reads_2 = 0;
    for (const auto& it : m_cache) {
      total_reads_1 += it.second->reads_1.size();
      total_reads_2 += it.second->reads_2.size();
    }

    AR_REQUIRE((total_reads_1 == total_reads_2) || (total_reads_2 == 0));
    auto seq1 =
      select_best_match(read_mate::_1, adapters, m_mate_1, total_reads_1);
    auto seq2 =
      select_best_match(read_mate::_2, adapters, m_mate_2, total_reads_2);

    const auto candidates =
      adapter_database{ m_config.samples->uninitialized_adapters() };
    const auto selection = candidates.identify(seq1, seq2);

    if (selection.first.sequence.empty() && selection.second.sequence.empty()) {
      switch (m_config.adapter_fallback_strategy) {
        case adapter_fallback::abort:
          throw std::runtime_error(
            "Could not identify adapter sequences automatically, and "
            "--adapter-fallback is set to 'abort'. Either specify adapter "
            "sequences manually or set change the --adapter-fallback option. "
            "See the documentation for more information.");

        case adapter_fallback::unknown:
          if (m_config.paired_ended_mode) {
            log::warn()
              << "Could not identify adapter sequences automatically. "
                 "Trimming will be performed without pre-defined "
                 "adapter sequences, equivalent to --adapter1 '' and "
                 "--adapter2 ''";
            break;
          } else {
            [[fallthrough]];
          }

        case adapter_fallback::none:
          log::warn() << "Could not identify adapter sequences automatically. "
                         "Trimming will be performed under the assumption that "
                         "the reads do not contain adapter sequences";

          m_config.adapter_selection_strategy = adapter_selection::none;
          break;

        default:
          AR_FAIL("invalid adapter_fallback value");
      }
    } else {
      describe_best_match(selection.first, "--adapter1", true);
      describe_best_match(selection.second,
                          "--adapter2",
                          m_config.paired_ended_mode);

      if (selection.first.mate != read_mate::_1 ||
          (m_config.paired_ended_mode &&
           selection.second.mate != read_mate::_2)) {
        log::warn()
          << "The orientation of identified adapters does not match "
             "expectations: Expected to find forward adapter in --in-file1 and "
             "reverse adapter in --in-file2 files (if any). Input files may be "
             "swapped";
      }
    }

    // FIXME: No other thread should be touching samples at this point, but
    //       that may change
    m_config.samples->set_adapters(
      adapter_set{ { selection.first.sequence.as_string(),
                     selection.second.sequence.as_string() } });
    m_config.samples->set_uninitialized_adapters(false);

    return std::move(m_cache);
  } else {
    return {};
  }
}

} // namespace adapterremoval
