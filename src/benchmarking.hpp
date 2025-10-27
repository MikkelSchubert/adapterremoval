// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <chrono>      // for high_resolution_clock
#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t
#include <string>      // for string
#include <string_view> // for string_view
#include <vector>      // for vector

namespace adapterremoval {

//! Strategy for running benchmarks
enum class strategy
{
  //! Carry out full benchmarking and return the final result
  benchmark,
  //! Run the benchmark loop once; used to collect results for later benchmarks
  passthrough,
  //! Skip the benchmark entirely (non-mandatory benchmarks)
  skip,
};

/** Class for tracking and validating user toggles for individual benchmarks */
class benchmark_toggles
{
public:
  explicit benchmark_toggles(std::vector<std::string> keys);

  /** Parses user-supplied toggles and returns true if valid */
  bool update_toggles(const std::vector<std::string>& keys);

  /** Returns true if the default (most) benchmarks should be run */
  bool defaults() const { return m_defaults; }

  /** Returns true if the specified benchmark toggle is set */
  [[nodiscard]] bool is_set(std::string_view key) const;

private:
  //! Vector of supported toggles
  const std::vector<std::string> m_toggles;
  //! Vector indicating whether the user set a specific toggle
  std::vector<bool> m_enabled;
  //! Indicates if default toggles should be used (benchmark specific)
  bool m_defaults;
};

/** Base-class for benchmarks */
class benchmarker
{
public:
  benchmarker(std::string desc, std::vector<std::string> toggles);

  virtual ~benchmarker() = default;

  virtual void run_if_toggled(const benchmark_toggles& toggles);

  benchmarker(const benchmarker&) = delete;
  benchmarker(benchmarker&&) = delete;
  benchmarker& operator=(const benchmarker&) = delete;
  benchmarker& operator=(benchmarker&&) = delete;

protected:
  /** Indicate that this benchmark *must* be run at least once */
  void set_required(bool value = true) { m_required = value; }

  /** Called before `setup` to perform any per batch setup */
  virtual void setup();
  /** Called to carry out work to be benchmarked */
  virtual void execute() = 0;

  virtual strategy enabled(const benchmark_toggles& toggles) const;

private:
  size_t count() const { return m_description.size(); }

  void run(strategy s);

  std::string summarize(size_t loops) const;

  using clock = std::chrono::steady_clock;
  using time_point = std::chrono::time_point<clock>;

  const std::string m_description;
  std::vector<uint64_t> m_durations{};
  std::vector<std::string> m_toggles{};
  bool m_required = false;
};

} // namespace adapterremoval
