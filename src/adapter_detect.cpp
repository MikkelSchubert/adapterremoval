// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>

#include "adapter_detect.hpp" // declarations
#include "adapters.hpp"       // for known_adapters
#include "alignment.hpp"      // for sequence_aligner
#include "debug.hpp"          // AR_REQUIRE
#include "fastq.hpp"          // for fastq
#include "logging.hpp"        // for log
#include "sequence.hpp"       // for dna_sequence
#include "sequence_sets.hpp"  // adapter_set
#include "strutils.hpp"       // for log_escape, starts_with
#include "utilities.hpp"      // for merge
#include <algorithm>          // for sort
#include <cstddef>            // for size_t
#include <limits>             // for numeric_limits
#include <string_view>        // for string_view
#include <vector>             // for vector

namespace adapterremoval {

namespace {

//! Minimum overlap required for adapter fragments
const size_t ADAPTER_DETECT_MIN_OVERLAP = 8;

//! The minimum number of hits required for a successful detection
const double DETECTION_THRESHOLD = 10;
//! The fraction of hits/reads required for detection of common sequences
const double DETECTION_THRESHOLD_COMMON = 1 / 100.0;
//! The minimum fraction of hits/reads required for detection of rare sequences
const double DETECTION_THRESHOLD_RARE = 1 / 1000.0;
//! The maximum error rate allowed for detection of rare sequences
const double DETECTION_THRESHOLD_RARE_MISMATCHES = 1 / 20.0;

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
collect_adapters(const adapter_set& user_adapters)
{
  // Database of known and user user-provided sequences
  adapter_database adapters{ user_adapters };

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
    if (it.length() >= ADAPTER_DETECT_MIN_OVERLAP) {
      selection.add(it, {});
    } else {
      log::warn() << "Adapter " << log_escape(it.as_string())
                  << " is too short for auto-detection; adapters must be at "
                  << ADAPTER_DETECT_MIN_OVERLAP << " bp long";
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

      overlaps.emplace_back(overlapping_before, overlapping_after);
    }

    prefixes.emplace_back(std::move(overlaps));
  }

  return prefixes;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// adapter_detection_stats

void
merge(adapter_detection_stats::stats& dst,
      const adapter_detection_stats::stats& src)
{
  dst.hits += src.hits;
  dst.length += src.length;
  dst.mismatches += src.mismatches;
}

void
adapter_detection_stats::merge(const adapter_detection_stats& other)
{
  AR_REQUIRE(mate_1.size() == 0 || mate_1.size() == other.mate_1.size());
  AR_REQUIRE(mate_2.size() == 0 || mate_2.size() == other.mate_2.size());

  reads += other.reads;
  ::adapterremoval::merge(mate_1, other.mate_1);
  ::adapterremoval::merge(mate_2, other.mate_2);
}

////////////////////////////////////////////////////////////////////////////////
// adapter_detector

namespace {

dna_sequence
select_best_match(const read_mate mate,
                  const adapter_set& adapters,
                  const adapter_detection_stats::values& values,
                  const size_t total_reads)
{
  if (values.empty()) {
    // no reads were processed, due to empty input or SE mode
    return {};
  }
  auto best_idx = std::numeric_limits<size_t>::max();
  adapter_detection_stats::stats best;

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

  if (best.hits >= DETECTION_THRESHOLD &&
      (best.hits >= total_reads * DETECTION_THRESHOLD_COMMON ||
       ((best.hits >= total_reads * DETECTION_THRESHOLD_RARE) &&
        (best.length * DETECTION_THRESHOLD_RARE_MISMATCHES >=
         best.mismatches)))) {
    log::debug() << "Selected " << adapters.at(best_idx).first.as_string();
    return adapters.at(best_idx).first;
  } else {
    return {};
  }
}

} // namespace

adapter_detector::adapter_detector(simd::instruction_set is,
                                   double mismatch_threshold)
  : adapter_detector({}, is, mismatch_threshold)
{
}

adapter_detector::adapter_detector(const adapter_set& user_sequences,
                                   simd::instruction_set is,
                                   double mismatch_threshold)
  : m_adapters(collect_adapters(user_sequences))
  , m_aligner(m_adapters, is, mismatch_threshold)
  , m_common_prefixes(
      collect_adapter_prefixes(m_adapters, ADAPTER_DETECT_MIN_OVERLAP))
{
  m_aligner.set_min_se_overlap(ADAPTER_DETECT_MIN_OVERLAP);
}

void
adapter_detector::detect_se(adapter_detection_stats& stats, const fastq& read)
{
  stats.reads++;
  detect_adapters(read, stats.mate_1);
}

void
adapter_detector::detect_pe(adapter_detection_stats& stats,
                            const fastq& read1,
                            const fastq& read2)
{
  // For now PE reads are processed like SE reads
  stats.reads++;
  detect_adapters(read1, stats.mate_1);
  detect_adapters(read2, stats.mate_2);
}

sequence_pair
adapter_detector::select_best(adapter_detection_stats& stats) const
{
  return {
    select_best_match(read_mate::_1, m_adapters, stats.mate_1, stats.reads),
    select_best_match(read_mate::_2, m_adapters, stats.mate_2, stats.reads),
  };
}

void
adapter_detector::detect_adapters(const fastq& read,
                                  adapter_detection_stats::values& stats)
{
  AR_REQUIRE(stats.empty() || stats.size() == m_adapters.size());
  stats.resize(m_adapters.size());

  // No shift, to avoid misdetection of adapters with common postfixes
  const auto alignment = m_aligner.align_single_end(read, 0);
  if (alignment.type() >= alignment_type::good) {
    auto idx = alignment.adapter_id();
    auto idx_start = idx;
    auto idx_end = idx;

    const auto adapter_len = m_adapters.at(idx).first.length();
    if (alignment.offset() + adapter_len >= read.length()) {
      auto length = alignment.length() - ADAPTER_DETECT_MIN_OVERLAP;

      idx_start -= m_common_prefixes.at(idx).at(length).first;
      idx_end += m_common_prefixes.at(idx).at(length).second;
    }

    for (idx = idx_start; idx <= idx_end; ++idx) {
      auto& it = stats.at(idx);

      it.hits++;
      it.length += alignment.length();
      it.mismatches += alignment.n_mismatches();
    }
  }
}

} // namespace adapterremoval
