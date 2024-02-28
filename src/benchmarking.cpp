/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
\*************************************************************************/
#include "benchmarking.hpp" // for benchmark_toggles, benchmarker
#include "debug.hpp"        // for AR_REQUIRE
#include "logging.hpp"      // for log
#include "mathutils.hpp"    // for arithmetic_mean, standard_deviation
#include "simd.hpp"         // for supported
#include "strutils.hpp"     // for to_lower
#include <algorithm>        // for find, accumulate
#include <iomanip>          // for fixed, setsetprecision
#include <numeric>          // for accumulate

namespace adapterremoval {

namespace {

//! Number of loops to perform prior to benchmarking as burn-in
const size_t BENCHMARK_BURN_IN = 1;
//! Benchmarks must be repeated at least this number of times
const size_t BENCHMARK_MIN_LOOPS = 10;
//! Benchmarks must be repeated at most this number of times
const size_t BENCHMARK_MAX_LOOPS = 1000;
//! Benchmarks must run for at this this number of nano-seconds
const double BENCHMARK_MIN_TIME_NS = 10'000'000'000;
//! Benchmark loops shorter than this number of nano-seconds cannot be measured
const double BENCHMARK_CUTOFF_TIME_NS = BENCHMARK_MIN_TIME_NS / 100'000;

} // namespace

benchmark_toggles::benchmark_toggles(const string_vec& keys)
  : m_toggles(keys)
  , m_enabled(keys.size(), false)
  , m_defaults()
{
}

bool
benchmark_toggles::update_toggles(const string_vec& keys)
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
    msg << "Valid toggles are ";
    for (size_t i = 0; i < m_toggles.size() - 1; ++i) {
      msg << m_toggles.at(i) << ", ";
    }
    msg << " and " << m_toggles.back();
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

benchmarker::benchmarker(const std::string& desc, string_vec toggles)
  : m_description(desc)
  , m_durations()
  , m_toggles(toggles)
  , m_required(false)
{
}

benchmarker::~benchmarker() {}

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
benchmarker::setup(){};

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
  log::cerr() << "\r\033[KBenchmarking: " << summarize(0);

  if (s != strategy::passthrough) {
    for (size_t i = 1; i <= BENCHMARK_BURN_IN; ++i) {
      setup();
      execute();

      log::cerr() << std::fixed << std::setprecision(5)
                  << "\rBenchmarking: " << m_description << " burn-in loop "
                  << i << " completed";
    }
  }

  const auto timer = clock();
  size_t loops = 0;
  do {
    uint64_t elapsed =
      std::accumulate(m_durations.begin(), m_durations.end(), uint64_t());

    do {
      setup();
      const auto start = timer.now();
      execute();
      const auto duration = (timer.now() - start).count();

      loops++;
      elapsed += duration;
      m_durations.push_back(duration);

      log::cerr() << "\r\033[KBenchmarking: " << summarize(loops);
    } while (s == strategy::benchmark &&
             m_durations.size() < BENCHMARK_MAX_LOOPS &&
             (m_durations.size() < BENCHMARK_MIN_LOOPS ||
              (elapsed < BENCHMARK_MIN_TIME_NS &&
               elapsed / m_durations.size() >= BENCHMARK_CUTOFF_TIME_NS)));
  } while (s == strategy::benchmark && (loops < 2 * m_durations.size()) &&
           grubbs_test_prune(m_durations));

  log::cerr() << "\r\033[K";
  log::info() << "  " << summarize(loops);
}

std::string
benchmarker::summarize(size_t loops) const
{
  std::ostringstream ss;
  if (loops) {
    ss << m_durations.size();
    if (loops > m_durations.size()) {
      ss << " + " << loops - m_durations.size();
    }

    ss << (loops != 1 ? " loops of " : " loop of ");
  }

  ss << m_description;

  if (m_durations.size()) {
    ss << " in " << std::fixed << std::setprecision(5)
       << arithmetic_mean(m_durations) / 1e9;

    if (m_durations.size() > 1) {
      ss << " +/- " << std::setprecision(6)
         << standard_deviation(m_durations) / 1e9;
    }

    ss << " seconds";
  }

  return ss.str();
}

} // namespace adapterremoval
