// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "benchmarking.hpp"      // for benchmarker, benchmarker_input, ...
#include "buffer.hpp"            // for buffer
#include "debug.hpp"             // for AR_REQUIRE
#include "fastq_in.hpp"          // for INPUT_BLOCK_SIZE
#include "linereader_joined.hpp" // for joined_line_readers
#include "logging.hpp"           // for cerr, error
#include "mathutils.hpp"         // for arithmetic_mean, standard_deviation
#include "strutils.hpp"          // for format_fraction, format_thousand_sep, ..
#include <algorithm>             // for find, min, minmax_element, sort, ...
#include <array>                 // for array
#include <chrono>                // for steady_clock, time_point
#include <cmath>                 // for round
#include <cstddef>               // for size_t
#include <cstdint>               // for uint64_t
#include <iomanip>               // for setw
#include <iostream>              // for cout, endl
#include <limits>                // for numeric_limits
#include <numeric>               // for accumulate
#include <sstream>               // for ostringstream
#include <string>                // for string
#include <utility>               // for move
#include <vector>                // for vector

namespace adapterremoval {

namespace {

//! Number of loops to perform prior to benchmarking as burn-in
constexpr size_t BENCHMARK_BURN_IN = 1;
//! Benchmarks must be repeated at least this number of times
constexpr size_t BENCHMARK_MIN_LOOPS = 10;
//! Benchmarks must be repeated at most this number of times
constexpr size_t BENCHMARK_MAX_LOOPS = 1000;
//! Benchmarks must run for at this number of nano-seconds
constexpr uint64_t BENCHMARK_MIN_TIME_NS = 5'000'000'000;
//! Benchmark loops shorter than this number of nano-seconds cannot be measured
constexpr uint64_t BENCHMARK_CUTOFF_TIME_NS = 10'000;
//! Number of NS between terminal updates
constexpr uint64_t BENCHMARK_UPDATE_INTERVAL = 50'000'000;

using clock = std::chrono::steady_clock;
using time_point = std::chrono::time_point<clock>;

} // namespace

////////////////////////////////////////////////////////////////////////////////
// benchmarker_input

benchmarker_input::benchmarker_input(
  std::vector<std::string> filenames,
  const fastq_encoding& encoding,
  const uint64_t head,
  const benchmarker_input_types::type required)
  : m_head(head)
  , m_filenames(std::move(filenames))
{
  AR_REQUIRE(head != std::numeric_limits<uint64_t>::max());
  AR_REQUIRE(!(required & ~benchmarker_input_types::mask));

  if (required) {
    m_lines = read_lines();
    m_lines_loaded = true;
    m_head = std::min<uint64_t>(m_head, (m_lines.size() + 3) / 4);
  }

  if (required & benchmarker_input_types::buffers) {
    m_buffers = make_buffers();
    m_buffers_loaded = true;
  }

  if (required & (benchmarker_input_types::unsafe_fastq |
                  benchmarker_input_types::validated_fastq)) {
    m_unsafe_fastqs = make_fastqs();
    m_unsafe_fastq_loaded = true;
  }

  if ((required & benchmarker_input_types::lines) == 0) {
    // lines are no longer needed beyond this point
    m_lines_loaded = false;
    m_lines = {};
  }

  if (required & benchmarker_input_types::validated_fastq) {
    // are unsafe FASTQs needed beyond this point?
    m_unsafe_fastq_loaded = required & benchmarker_input_types::unsafe_fastq;
    m_validated_fastqs = validate_fastqs(encoding, !m_unsafe_fastq_loaded);
    m_validated_fastq_loaded = true;
  }
}

benchmarker_input::~benchmarker_input() = default;

const std::vector<std::string>&
benchmarker_input::filenames() const
{
  return m_filenames;
}

const std::vector<std::string>&
benchmarker_input::lines() const
{
  AR_REQUIRE(m_lines_loaded);
  return m_lines;
}

const std::vector<buffer>&
benchmarker_input::buffers() const
{
  AR_REQUIRE(m_buffers_loaded);
  return m_buffers;
}

const std::vector<fastq>&
benchmarker_input::unsafe_fastqs() const
{
  AR_REQUIRE(m_unsafe_fastq_loaded);
  return m_unsafe_fastqs;
}

const std::vector<fastq>&
benchmarker_input::validated_fastqs() const
{
  AR_REQUIRE(m_validated_fastq_loaded);
  return m_validated_fastqs;
}

std::vector<std::string>
benchmarker_input::read_lines()
{
  std::vector<std::string> lines;
  if (!m_filenames.empty()) {
    joined_line_readers reader(m_filenames);

    // vector and strings are reused to match actual usage via `read_fastq`
    while (lines.size() < m_head * 4) {
      lines.emplace_back();
      if (!reader.getline(lines.back())) {
        lines.pop_back();
        break;
      }
    }
  }

  return lines;
}

std::vector<buffer>
benchmarker_input::make_buffers()
{
  std::vector<buffer> buffers;
  if (m_lines.empty()) {
    return buffers;
  }

  buffer buf;
  buf.reserve(INPUT_BLOCK_SIZE);

  for (const auto& line : m_lines) {
    const char* start = line.data();
    const char* end = start + line.size();

    while (start < end) {
      const size_t to_copy =
        std::min<size_t>(INPUT_BLOCK_SIZE - buf.size(), end - start);

      buf.append(start, to_copy);
      if (buf.size() == INPUT_BLOCK_SIZE) {
        buffers.emplace_back(std::move(buf));
        buf.reserve(INPUT_BLOCK_SIZE);
      }

      start += to_copy;
    }

    buf.append_u8('\n');
    if (buf.size() == INPUT_BLOCK_SIZE) {
      buffers.emplace_back(std::move(buf));
      buf.reserve(INPUT_BLOCK_SIZE);
    }
  }

  if (buf.size()) {
    buffers.emplace_back(std::move(buf));
  }

  return buffers;
}

std::vector<fastq>
benchmarker_input::make_fastqs()
{
  std::vector<fastq> fastqs;

  if (!m_lines.empty()) {
    fastqs.reserve((m_lines.size() / 4) + 1);
    borrowed_vec_reader reader{ m_lines.begin(), m_lines.end() };
    while (fastqs.size() < (m_lines.size() / 4)) {
      fastqs.emplace_back();
      if (!fastqs.back().read_unsafe(reader)) {
        fastqs.pop_back();
        break;
      }
    }
  }

  return fastqs;
}

std::vector<fastq>
benchmarker_input::validate_fastqs(const fastq_encoding& encoding,
                                   const bool consume)
{
  std::vector<fastq> result;
  if (consume) {
    result = std::move(m_unsafe_fastqs);
  } else {
    result = m_unsafe_fastqs;
  }

  for (auto& it : result) {
    it.post_process(encoding);
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////
// benchmarker

benchmarker::benchmarker(std::string desc,
                         std::vector<std::string> toggles,
                         benchmarker_input_types::type input)
  : m_description(std::move(desc))
  , m_toggles(std::move(toggles))
  , m_input(input)
{
}

void
benchmarker::run(const benchmarker_input& input_1,
                 const benchmarker_input& input_2)
{
  static bool header = false;
  if (!header) {
    std::cout << "             Benchmark |       Min |      Mean |       Max | "
                 "SD (%) | Loops | Outliers"
              << std::endl;
    header = true;
  }

  for (size_t i = 1; i <= BENCHMARK_BURN_IN; ++i) {
    log::cerr().transient() << "\r\033[K" << std::setw(22) << m_description
                            << " burn-in: " << i << " of " << BENCHMARK_BURN_IN;

    setup(input_1, input_2);
    execute(input_1, input_2);
  }

  size_t loops = 0;
  uint64_t since_update = 0;
  do {
    uint64_t elapsed =
      std::accumulate(m_durations.begin(), m_durations.end(), uint64_t());

    do {
      setup(input_1, input_2);
      const auto start = clock::now();
      execute(input_1, input_2);
      const auto duration = (clock::now() - start).count();

      loops++;
      elapsed += duration;
      since_update += duration;
      m_durations.push_back(duration);

      if (since_update >= BENCHMARK_UPDATE_INTERVAL) {
        log::cerr() << "\r\033[K" << summarize(loops);
        since_update %= BENCHMARK_UPDATE_INTERVAL;
      }
    } while (m_durations.size() < BENCHMARK_MAX_LOOPS &&
             (m_durations.size() < BENCHMARK_MIN_LOOPS ||
              (elapsed < BENCHMARK_MIN_TIME_NS &&
               elapsed >= m_durations.size() * BENCHMARK_CUTOFF_TIME_NS)));
  } while ((loops < 2 * m_durations.size()) && grubbs_test_prune(m_durations));

  log::cerr() << "\r\033[K";
  std::cout << summarize(loops) << std::endl;
}

bool
benchmarker::enabled(const std::vector<std::string>& toggles) const
{
  // Run everything by default
  if (toggles.empty()) {
    return true;
  }

  for (const auto& toggle : m_toggles) {
    if (std::find(toggles.begin(), toggles.end(), toggle) != toggles.end()) {
      return true;
    }
  }

  return false;
}

void
benchmarker::setup(const benchmarker_input& /* input_1 */,
                   const benchmarker_input& /* input_2 */)
{
}

std::string
benchmarker::summarize(size_t loops) const
{
  const std::array<size_t, 7> COLUMN_WIDTHS{ 22, 9, 9, 9, 6, 5, 8 };
  std::array<std::string, COLUMN_WIDTHS.size()> values{
    m_description,
    "",
    "",
    "",
    "",
    stringify(m_durations.size()),
    stringify(loops - m_durations.size()),
  };

  if (!m_durations.empty()) {
    const auto min_max =
      std::minmax_element(m_durations.begin(), m_durations.end());
    const auto mean = arithmetic_mean(m_durations);

    values.at(1) = format_thousand_sep(*min_max.first / 1e3);
    values.at(2) = format_thousand_sep(std::round(mean / 1e3));
    values.at(3) = format_thousand_sep(*min_max.second / 1e3);

    if (m_durations.size() > 1) {
      const auto sd = standard_deviation(m_durations);
      values.at(4) = format_fraction(100.0 * sd, mean);
    }
  }

  std::ostringstream ss;
  for (size_t i = 0; i < values.size(); ++i) {
    if (i) {
      ss << " | ";
    }

    ss << std::setw(COLUMN_WIDTHS.at(i)) << values.at(i);
  }

  return ss.str();
}

////////////////////////////////////////////////////////////////////////////////

namespace {

std::vector<std::string>
collect_toggles(const std::vector<std::unique_ptr<benchmarker>>& benchmarks)
{
  std::vector<std::string> toggles;
  for (const auto& benchmark : benchmarks) {
    const auto& values = benchmark->toggles();

    toggles.insert(toggles.end(), values.begin(), values.end());
  }

  std::sort(toggles.begin(), toggles.end());
  auto last = std::unique(toggles.begin(), toggles.end());
  toggles.erase(last, toggles.end());

  return toggles;
}

} // namespace

bool
validate_toggles(const std::vector<std::unique_ptr<benchmarker>>& benchmarks,
                 const std::vector<std::string>& candidates)
{
  bool any_errors = false;

  const auto toggles = collect_toggles(benchmarks);
  for (const auto& toggle : candidates) {
    const auto it = std::find(toggles.begin(), toggles.end(), toggle);

    if (it == toggles.end()) {
      log::error() << "Unknown benchmarking toggle " << log_escape(toggle);
      any_errors = true;
    }
  }

  if (any_errors) {
    auto msg = log::error();
    msg << "Valid toggles are " << join_text(toggles, ", ", ", and ");
  }

  return !any_errors;
}

} // namespace adapterremoval
