// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "benchmarking.hpp" // for benchmark_toggles, benchmarker
#include "debug.hpp"        // for AR_REQUIRE
#include "logging.hpp"      // for log
#include "mathutils.hpp"    // for arithmetic_mean, standard_deviation
#include "strutils.hpp"     // for to_lower
#include <algorithm>        // for find, accumulate
#include <array>            // for array
#include <cmath>            // for round
#include <cstdint>          // for uint64_t
#include <iomanip>          // for fixed, setsetprecision
#include <iostream>         // for cout
#include <numeric>          // for accumulate
#include <sstream>          // for ostringstream
#include <string>           // for string, to_string
#include <utility>          // for move
#include <vector>           // for vector

namespace adapterremoval {

namespace {

//! Number of loops to perform prior to benchmarking as burn-in
const size_t BENCHMARK_BURN_IN = 1;
//! Benchmarks must be repeated at least this number of times
const size_t BENCHMARK_MIN_LOOPS = 10;
//! Benchmarks must be repeated at most this number of times
const size_t BENCHMARK_MAX_LOOPS = 1000;
//! Benchmarks must run for at this this number of nano-seconds
const double BENCHMARK_MIN_TIME_NS = 5'000'000'000;
//! Benchmark loops shorter than this number of nano-seconds cannot be measured
const double BENCHMARK_CUTOFF_TIME_NS = 10'000;
//! Number of NS between terminal updates
const size_t BENCHMARK_UPDATE_INTERVAL = 50'000'000;

} // namespace

benchmark_toggles::benchmark_toggles(std::vector<std::string> keys)
  : m_toggles(std::move(keys))
  , m_enabled(m_toggles.size(), false)
  , m_defaults()
{
}

bool
benchmark_toggles::update_toggles(const std::vector<std::string>& keys)
{
  bool any_errors = false;
  m_defaults = keys.empty();
  for (const auto& key : keys) {
    const auto it =
      std::find(m_toggles.begin(), m_toggles.end(), to_lower(key));

    if (it == m_toggles.end()) {
      log::error() << "Unknown benchmarking toggle '" << key << "'";
      any_errors = true;
    } else {
      m_enabled.at(it - m_toggles.begin()) = true;
    }
  }

  if (any_errors) {
    auto msg = log::error();
    msg << "Valid toggles are " << join_text(m_toggles, ", ", ", and ");
  }

  return !any_errors;
}

bool
benchmark_toggles::is_set(const std::string& key) const
{
  const auto it = std::find(m_toggles.begin(), m_toggles.end(), to_lower(key));

  AR_REQUIRE(it != m_toggles.end());
  return m_enabled.at(it - m_toggles.begin());
}

benchmarker::benchmarker(std::string desc, std::vector<std::string> toggles)
  : m_description(std::move(desc))
  , m_toggles(std::move(toggles))
{
}

void
benchmarker::run_if_toggled(const benchmark_toggles& toggles)
{
  const auto s = enabled(toggles);
  if (s != strategy::skip) {
    run(s);
  }
}

/** Called before `setup` to perform any per batch setup */
void
benchmarker::setup()
{
}

strategy
benchmarker::enabled(const benchmark_toggles& toggles) const
{
  if (toggles.defaults()) {
    return strategy::benchmark;
  }

  for (const auto& toggle : m_toggles) {
    if (toggles.is_set(toggle)) {
      return strategy::benchmark;
    }
  }

  return m_required ? strategy::passthrough : strategy::skip;
}

void
benchmarker::run(const strategy s)
{
  static bool header = false;
  if (!header) {
    std::cout << "             Benchmark |       Min |      Mean |       Max | "
                 "SD (%) | Loops | Outliers"
              << std::endl;
    header = true;
  }

  if (s != strategy::passthrough) {
    log::cerr() << "\r\033[K" << std::setw(22) << m_description << " burn-in";

    for (size_t i = 1; i <= BENCHMARK_BURN_IN; ++i) {
      setup();
      execute();

      log::cerr() << ".";
    }
  } else {
    log::cerr() << "\r\033[K" << std::setw(22) << m_description << " (setup)";
  }

  size_t loops = 0;
  uint64_t next_update = 0;
  do {
    uint64_t elapsed =
      std::accumulate(m_durations.begin(), m_durations.end(), uint64_t());

    do {
      setup();
      const auto start = clock::now();
      execute();
      const auto duration = (clock::now() - start).count();

      loops++;
      elapsed += duration;
      m_durations.push_back(duration);

      if (elapsed >= next_update) {
        log::cerr() << "\r\033[K" << summarize(loops);
        next_update += BENCHMARK_UPDATE_INTERVAL;
      }
    } while (s == strategy::benchmark &&
             m_durations.size() < BENCHMARK_MAX_LOOPS &&
             (m_durations.size() < BENCHMARK_MIN_LOOPS ||
              (elapsed < BENCHMARK_MIN_TIME_NS &&
               elapsed / m_durations.size() >= BENCHMARK_CUTOFF_TIME_NS)));
  } while (s == strategy::benchmark && (loops < 2 * m_durations.size()) &&
           grubbs_test_prune(m_durations));

  log::cerr() << "\r\033[K";

  if (s != strategy::passthrough) {
    std::cout << summarize(loops) << std::endl;
  }
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
    std::to_string(m_durations.size()),
    std::to_string(loops - m_durations.size()),
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
      values.at(4) = format_fraction(1e9 * sd, 1e7 * mean);
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

} // namespace adapterremoval
