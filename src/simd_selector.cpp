// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "simd_selector.hpp" // declarations
#include "alignment.hpp"     // for sequence_aligner
#include "debug.hpp"         // for AR_REQUIRE
#include "logging.hpp"       // for log
#include "scheduler.hpp"     // for analytical_step, chunk_ptr, chunk_vec, ...
#include "sequence_sets.hpp" // for adapter_set
#include "simd.hpp"          // for supported
#include "strutils.hpp"      // for format_rough_number
#include "userconfig.hpp"    // for userconfig
#include "utilities.hpp"     // for dynamic_cast_unique
#include <chrono>            // for steady_clock
#include <limits>            // for numeric_limits

namespace adapterremoval {

namespace {

//! The minimum number of reads to benchmark before picking an instruction set
constexpr size_t MIN_READS_FOR_BENCHMARK = 10'000;
//! The maximum amount of time in ns to spend benchmarking each instruction set
constexpr std::chrono::milliseconds MAX_NS_PER_TEST{ 100 };

} // namespace

simd_selector::simd_selector(threadsafe_data<sample_set> samples,
                             threadsafe_data<simd::instruction_set> is,
                             double mismatch_threshold,
                             uint32_t shift,
                             size_t next_step)
  : analytical_step(processing_order::ordered, "simd_selector")
  , m_next_step(next_step)
  , m_samples(std::move(samples))
  , m_simd(std::move(is))
  , m_mismatch_threshold(mismatch_threshold)
  , m_shift(shift)
{
}

chunk_vec
simd_selector::process(chunk_ptr data)
{
  chunk_vec chunks;

  if (m_processed_reads < MIN_READS_FOR_BENCHMARK) {
    auto chunk = dynamic_cast_unique<fastq_chunk>(data);
    AR_REQUIRE(chunk);

    const auto adapters = m_samples.get_reader()->adapters();

    if (
      // No point in benchmarking if there's no more data
      (m_processed_reads == 0 && chunk->eof) ||
      // No point in benchmarking if we won't align anything
      (chunk->reads_2.empty() &&
       (adapters == adapter_set{} || adapters == adapter_set{ {} }))) {
      // Skip future checks, if any
      m_processed_reads = std::numeric_limits<size_t>::max();
    } else {
      const auto candidates = simd::supported();
      m_candidates.resize(candidates.size());

      bool ready = true;
      for (size_t i = 0; i < candidates.size(); ++i) {
        // A new aligner is constructed every loop to prevent automatic
        // re-ordering of sequences to affect results between runs
        sequence_aligner aligner{ adapters,
                                  candidates.at(i),
                                  m_mismatch_threshold };

        auto& candidate = m_candidates.at(i);
        if (candidate.second < MAX_NS_PER_TEST &&
            candidate.first < MIN_READS_FOR_BENCHMARK) {
          candidate.first += chunk->reads();
          candidate.second += process_chunk(aligner, *chunk);
        }

        ready &= candidate.first >= MIN_READS_FOR_BENCHMARK ||
                 candidate.second >= MAX_NS_PER_TEST;
      }

      if (ready || chunk->eof) {
        log::info() << "Selecting optimal SIMD instruction set for data:";

        auto best = std::numeric_limits<double>::min();
        simd::instruction_set selected = simd::instruction_set::none;
        for (size_t i = 0; i < candidates.size(); ++i) {
          const auto [reads, duration] = m_candidates.at(i);
          const auto seconds =
            std::chrono::duration_cast<std::chrono::duration<double>>(duration);
          const double rate = reads / std::max<double>(1e-9, seconds.count());

          log::info() << "  - " << simd::name(candidates.at(i)) << ": "
                      << format_rough_number(rate) << " reads/s";

          if (rate > best) {
            selected = candidates.at(i);
            best = rate;
          }
        }

        log::info() << "Selected SIMD instruction set " << simd::name(selected);

        *m_simd.get_writer() = selected;
        m_processed_reads = std::numeric_limits<size_t>::max();
      } else {
        m_processed_reads += chunk->reads();
      }
    }

    chunks.emplace_back(m_next_step, std::move(chunk));
  } else {
    chunks.emplace_back(m_next_step, std::move(data));
  }

  return chunks;
}

simd_selector::clock::duration
simd_selector::process_chunk(sequence_aligner& aligner,
                             const fastq_chunk& chunk) const
{
  auto start_time = clock::now();
  if (chunk.reads_2.empty()) {
    for (const auto& read : chunk.reads_1) {
      auto alignment = aligner.align_single_end(read, m_shift);
      // Insurance against the compiler eliding this entire loop
      blackbox(alignment);
    }
  } else {
    AR_REQUIRE(chunk.reads_1.size() == chunk.reads_2.size());
    auto it_1 = chunk.reads_1.begin();
    auto it_2 = chunk.reads_2.begin();

    for (; it_1 != chunk.reads_1.end(); ++it_1, ++it_2) {
      auto read_2 = *it_2;
      read_2.reverse_complement();

      auto alignment = aligner.align_paired_end(*it_1, read_2, m_shift);
      // Insurance against the compiler eliding this entire loop
      blackbox(alignment);
    }
  }

  return (clock::now() - start_time);
}

} // namespace adapterremoval
