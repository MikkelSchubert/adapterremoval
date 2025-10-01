// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_detector.hpp" // declarations
#include "adapter_database.hpp" // for known_adapters
#include "alignment.hpp"        // for sequence_aligner
#include "debug.hpp"            // AR_REQUIRE
#include "fastq.hpp"            // for fastq
#include "logging.hpp"          // for log
#include "sequence.hpp"         // for dna_sequence
#include "sequence_sets.hpp"    // adapter_set
#include "strutils.hpp"         // for log_escape, starts_with
#include "utilities.hpp"        // for merge
#include <algorithm>            // for sort
#include <cstddef>              // for size_t
#include <limits>               // for numeric_limits
#include <string_view>          // for string_view

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

/** Remove duplicate sequences from a vector; values are sorted */
sequence_vec
collect_unique(sequence_vec adapters)
{
  std::sort(adapters.begin(), adapters.end());
  sequence_vec unique;

  for (auto&& it : adapters) {
    if (it.length() >= ADAPTER_DETECT_MIN_OVERLAP) {
      if (unique.empty() || !(unique.back() == it)) {
        unique.emplace_back(std::move(it));
      }
    } else {
      // Check is here to avoid repeat messages for copies of the same sequence
      log::warn() << "Adapter sequence " << log_escape(it.as_string())
                  << " is too short for auto-detection; adapters must be at "
                  << ADAPTER_DETECT_MIN_OVERLAP << " bp long";
    }
  }

  return unique;
}

/** Collect all unique adapter sequence as read 1 adapters */
sequence_vec
collect_adapters(const adapter_database& database)
{
  sequence_vec sequences;
  for (const auto& it : database) {
    sequences.insert(sequences.end(),
                     it.adapter_1().begin(),
                     it.adapter_1().end());
    sequences.insert(sequences.end(),
                     it.adapter_2().begin(),
                     it.adapter_2().end());
  }

  return collect_unique(std::move(sequences));
}

adapter_set
adapters_to_set(sequence_vec sequences)
{
  adapter_set adapters;
  for (auto&& it : sequences) {
    adapters.add(std::move(it), {});
  }

  return adapters;
}

std::vector<std::vector<std::pair<size_t, size_t>>>
common_prefixes(const sequence_vec& adapters, size_t min_overlap)
{
  std::vector<std::vector<std::pair<size_t, size_t>>> prefixes;
  prefixes.reserve(adapters.size());

  for (size_t i = 0; i < adapters.size(); ++i) {
    std::string_view adapter = adapters.at(i).as_string();

    std::vector<std::pair<size_t, size_t>> overlaps;
    overlaps.reserve(adapter.length() -
                     std::min(min_overlap, adapter.length()));

    for (size_t len = min_overlap; len <= adapter.length(); ++len) {
      auto prefix = adapter.substr(0, len);

      size_t overlapping_before = 0;
      for (size_t j = i; j; --j) {
        if (starts_with(adapters.at(j - 1).as_string(), prefix)) {
          overlapping_before++;
        } else {
          break;
        }
      }

      size_t overlapping_after = 0;
      for (size_t j = i + 1; j < adapters.size(); ++j) {
        if (starts_with(adapters.at(j).as_string(), prefix)) {
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

dna_sequence
select_best_match(const read_mate mate,
                  const sequence_vec& adapters,
                  const adapter_detection_stats::values& values,
                  const size_t total_reads)
{
  AR_REQUIRE(values.empty() || values.size() == adapters.size());
  if (values.empty()) {
    // no reads were processed, due to empty input or SE mode
    return {};
  }

  adapter_detection_stats::hits best;
  auto best_idx = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < adapters.size(); ++i) {
    const auto& current = values.at(i);

    if (current.hits > best.hits ||
        (current.hits == best.hits && (current.aligned - current.mismatches) >
                                        (best.aligned - best.mismatches))) {
      best = current;
      best_idx = i;
    }

    if (current.hits) {
      log::debug() << mate << " adapter has " << current.hits << " hits, "
                   << current.aligned << " bases, " << current.mismatches
                   << " mismatches: " << adapters.at(i).as_string();
    }
  }

  if (best.hits >= DETECTION_THRESHOLD &&
      (best.hits >= total_reads * DETECTION_THRESHOLD_COMMON ||
       ((best.hits >= total_reads * DETECTION_THRESHOLD_RARE) &&
        (best.aligned * DETECTION_THRESHOLD_RARE_MISMATCHES >=
         best.mismatches)))) {
    log::debug() << "Selected " << adapters.at(best_idx).as_string();
    return adapters.at(best_idx);
  } else {
    return {};
  }
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// adapter_detection_stats

void
merge(adapter_detection_stats::hits& dst,
      const adapter_detection_stats::hits& src)
{
  dst.hits += src.hits;
  dst.aligned += src.aligned;
  dst.mismatches += src.mismatches;
}

adapter_detection_stats::adapter_detection_stats(size_t reads,
                                                 std::vector<hits> mate_1,
                                                 std::vector<hits> mate_2)
  : m_reads_1(reads)
  , m_reads_2(mate_2.size() ? reads : 0)
  , m_mate_1{ std::move(mate_1) }
  , m_mate_2{ std::move(mate_2) }
{
  AR_REQUIRE(mate_2.empty() || (mate_1.size() == mate_2.size()));
}

void
adapter_detection_stats::merge(const adapter_detection_stats& other)
{
  AR_REQUIRE(m_mate_1.size() == 0 || m_mate_1.size() == other.m_mate_1.size());
  AR_REQUIRE(m_mate_2.size() == 0 || m_mate_2.size() == other.m_mate_2.size());

  m_reads_1 += other.m_reads_1;
  m_reads_2 += other.m_reads_2;
  ::adapterremoval::merge(m_mate_1, other.m_mate_1);
  ::adapterremoval::merge(m_mate_2, other.m_mate_2);
}

bool
operator==(const adapter_detection_stats::hits& a,
           const adapter_detection_stats::hits& b)
{
  return (a.hits == b.hits) && (a.aligned == b.aligned) &&
         (a.mismatches == b.mismatches);
}

/** Stream operator for debugging output */
std::ostream&
operator<<(std::ostream& os, const adapter_detection_stats::hits& value)
{
  return os << "adapter_detection_stats::hits{hits=" << value.hits
            << ", aligned=" << value.aligned
            << ", mismatches=" << value.mismatches << "}";
}

////////////////////////////////////////////////////////////////////////////////
// adapter_detector

adapter_detector::adapter_detector(adapter_database database,
                                   const simd::instruction_set is,
                                   const double mismatch_threshold)
  : m_database(std::move(database))
  , m_adapters(collect_adapters(m_database))
  , m_aligner(adapters_to_set(m_adapters), is, mismatch_threshold)
  , m_common_prefixes(common_prefixes(m_adapters, ADAPTER_DETECT_MIN_OVERLAP))
{
  m_aligner.set_min_se_overlap(ADAPTER_DETECT_MIN_OVERLAP);
}

void
adapter_detector::detect_se(adapter_detection_stats& stats, const fastq& read)
{
  stats.m_reads_1++;
  detect_adapters(read, stats.m_mate_1);
}

void
adapter_detector::detect_pe(adapter_detection_stats& stats,
                            const fastq& read_1,
                            const fastq& read_2)
{
  // For now PE reads are processed like SE reads
  stats.m_reads_1++;
  stats.m_reads_2++;
  detect_adapters(read_1, stats.m_mate_1);
  detect_adapters(read_2, stats.m_mate_2);
}

identified_adapter_pair
adapter_detector::select_best(const adapter_detection_stats& stats) const
{
  const auto seq_1 = select_best_match(read_mate::_1,
                                       m_adapters,
                                       stats.m_mate_1,
                                       stats.m_reads_1);
  const auto seq_2 = select_best_match(read_mate::_2,
                                       m_adapters,
                                       stats.m_mate_2,
                                       stats.m_reads_2);

  return m_database.identify_exact(seq_1, seq_2);
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

    const auto adapter_len = m_adapters.at(idx).length();
    if (alignment.offset() + adapter_len >= read.length()) {
      auto length = alignment.length() - ADAPTER_DETECT_MIN_OVERLAP;

      idx_start -= m_common_prefixes.at(idx).at(length).first;
      idx_end += m_common_prefixes.at(idx).at(length).second;
    }

    for (idx = idx_start; idx <= idx_end; ++idx) {
      auto& it = stats.at(idx);

      it.hits++;
      it.aligned += alignment.length() - alignment.n_ambiguous();
      it.mismatches += alignment.n_mismatches();
    }
  }
}

} // namespace adapterremoval
