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
#include "alignment.hpp"         // for alignment_info, sequence_aligner
#include "debug.hpp"             // for AR_REQUIRES
#include "fastq.hpp"             // for fastq
#include "linereader_joined.hpp" // for joined_line_readers
#include "logging.hpp"           // for log
#include "mathutils.hpp"         // for arithmetic_mean, standard_deviation
#include "strutils.hpp"          // for to_lower
#include "userconfig.hpp"        // for userconfig, output_files
#include <algorithm>             // for swap
#include <chrono>                // for high_resolution_clock
#include <cstring>               // for size_t
#include <iomanip>               // for fixed, setsetprecision
#include <sstream>               // for ostringstream
#include <vector>                // for vector

namespace adapterremoval {

namespace {

const size_t BENCHMARK_BURN_IN = 1;
const size_t BENCHMARK_MIN_LOOPS = 10;
const size_t BENCHMARK_MAX_LOOPS = 1000;
const double BENCHMARK_MIN_TIME_NS = 10'000'000'000;
const double BENCHMARK_CUTOFF_TIME_NS = BENCHMARK_MIN_TIME_NS / 100'000;

//! Strategy for running benchmarks
enum class strategy
{
  //! Carry out full benchmarking and return the final result
  benchmark,
  //! Run the benchmark loop once and return the result
  passthrough,
};

/** Enum representing the available benchmarks */
enum class benchmarks : size_t
{
  //! Benchmarking of decompression and line reading
  read = 0,
  //! Benchmarking of FASTQ parsing
  parse,
  //! Benchmarking of reverse complementations
  revcompl,
  //! Benchmarking of single-end alignments
  align_se,
  //! Benchmarking of paired-end alignments
  align_pe,
  //! Benchmark SIMD alignment functions other than just the preferred one
  simd,
  //! Indicator variable for the number of items
  MAX
};

class benchmark_toggles
{
public:
  benchmark_toggles()
    : m_enabled(static_cast<size_t>(benchmarks::MAX), false)
  {
  }

  bool parse_args(const string_vec& items)
  {
    if (items.empty()) {
      m_enabled = std::vector<bool>(static_cast<size_t>(benchmarks::MAX), true);
      return true;
    }

    bool any_errors = false;
    for (const auto& it : items) {
      const auto key = to_lower(it);

      benchmarks enum_key = benchmarks::MAX;
      if (key == "read") {
        enum_key = benchmarks::read;
      } else if (key == "parse") {
        enum_key = benchmarks::parse;
      } else if (key == "revcompl") {
        enum_key = benchmarks::revcompl;
      } else if (key == "align_se") {
        enum_key = benchmarks::align_se;
      } else if (key == "align_pe") {
        enum_key = benchmarks::align_pe;
      } else if (key == "simd") {
        enum_key = benchmarks::simd;
      } else {
        log::error() << "Unknown benchmark key '" << it << "'";
        any_errors = true;
      }

      if (enum_key != benchmarks::MAX) {
        m_enabled.at(static_cast<size_t>(enum_key)) = true;
      }
    }

    return !any_errors;
  }

  strategy operator()(benchmarks b) const
  {
    return enabled(b) ? strategy::benchmark : strategy::passthrough;
  }

  strategy align(benchmarks s, simd::instruction_set is) const
  {
    if (enabled(s)) {
      // Benchmarking all instruction sets
      if (enabled(benchmarks::simd)) {
        return strategy::benchmark;
      }

      // Otherwise just benchmarking the preferred instruction set
      const auto supported = simd::supported();
      if (supported.back() == is) {
        return strategy::benchmark;
      }
    }

    return strategy::passthrough;
  }

private:
  bool enabled(benchmarks b) const
  {
    return m_enabled.at(static_cast<size_t>(b));
  }

  std::vector<bool> m_enabled;
};

template<typename T>
class benchmarker
{
public:
  benchmarker(const std::string& desc)
    : m_description(desc)
    , m_durations()
  {
  }

  virtual ~benchmarker() {}

  T run(const strategy s)
  {
    T t = T();
    const auto timer = clock();
    size_t burn_in = 0;
    size_t loops = 0;

    do {
      uint64_t elapsed =
        std::accumulate(m_durations.begin(), m_durations.end(), uint64_t());

      do {
        setup();
        const auto start = timer.now();
        auto current = execute();
        const auto duration = (timer.now() - start).count();

        if (burn_in < BENCHMARK_BURN_IN && s == strategy::benchmark) {
          burn_in++;
        } else {
          loops++;
          elapsed += duration;
          m_durations.push_back(duration);
        }

        std::swap(t, current);

        if (m_durations.size()) {
          log::cerr() << "\rRunning in " << summarize(loops);
        } else {
          log::cerr() << std::fixed << std::setprecision(5) << "\rBurning in "
                      << m_description << " (" << burn_in << " loops) ";
        }
      } while (s == strategy::benchmark &&
               m_durations.size() < BENCHMARK_MAX_LOOPS &&
               (m_durations.size() < BENCHMARK_MIN_LOOPS ||
                (elapsed < BENCHMARK_MIN_TIME_NS &&
                 elapsed / m_durations.size() >= BENCHMARK_CUTOFF_TIME_NS)));
    } while (s == strategy::benchmark && (loops < 2 * m_durations.size()) &&
             grubbs_test_prune(m_durations));

    log::cerr() << "\r\033[K";
    log::info() << "  " << summarize(loops);

    return t;
  }

protected:
  virtual void setup(){};
  virtual T execute() = 0;

private:
  size_t count() const { return m_description.size(); }

  double mean() const { return arithmetic_mean(m_durations) / 1e9; }

  double sd() const { return standard_deviation(m_durations) / 1e9; }

  std::string summarize(size_t loops) const
  {
    AR_REQUIRE(m_durations.size());
    std::ostringstream ss;

    ss << std::fixed << std::setprecision(5) << mean();
    if (m_durations.size() > 1) {
      ss << " +/- " << std::setprecision(6) << sd();
    }

    if (loops > 1) {
      ss << " seconds (" << m_durations.size();
      if (loops > m_durations.size()) {
        ss << "+" << loops - m_durations.size();
      }

      ss << " loops)";
    }

    ss << " to " << m_description;

    return ss.str();
  }

  using clock = std::chrono::high_resolution_clock;
  using time_point = std::chrono::time_point<clock>;

  const std::string m_description;
  std::vector<uint64_t> m_durations;
};

class readlines_benchmarker : public benchmarker<string_vec>
{
public:
  readlines_benchmarker(const userconfig& config, const string_vec& filenames)
    : benchmarker("read FASTQs")
    , m_filenames(filenames)
    , m_head(config.head)
  {
  }

protected:
  string_vec execute()
  {
    joined_line_readers reader(m_filenames);

    std::string line;
    string_vec lines;
    while (lines.size() / 4 < m_head && reader.getline(line)) {
      lines.push_back(line);
    }

    return lines;
  }

private:
  const string_vec m_filenames;
  const size_t m_head;
};

/** Benchmarking of FASTQ parsing excluding file IO */
class fastq_parser_benchmarker : public benchmarker<fastq_vec>
{
public:
  fastq_parser_benchmarker(const string_vec& lines)
    : benchmarker("parse FASTQs")
    , m_lines(lines)
  {
  }

protected:
  fastq_vec execute()
  {
    fastq record;
    fastq_vec records;
    vec_reader reader(m_lines);
    while (record.read(reader)) {
      records.push_back(record);
    }

    return records;
  }

private:
  const string_vec& m_lines;
};

class reverse_complement_benchmarker : public benchmarker<bool>
{
public:
  reverse_complement_benchmarker(const fastq_vec& records)
    : benchmarker("reverse complement FASTQs")
    , m_records(records)
  {
  }

protected:
  bool execute()
  {
    // Repeated to match normal procedure during PE alignments
    for (size_t i = 0; i < 2; ++i) {
      for (auto& it : m_records) {
        it.reverse_complement();
      }
    }

    return true;
  }

private:
  fastq_vec m_records;
};

/** Benchmarking of SE alignments */
class benchmarker_se_alignment : public benchmarker<alignment_info>
{
public:
  benchmarker_se_alignment(const userconfig& config,
                           const fastq_vec& reads,
                           const simd::instruction_set is)
    : benchmarker("align single-end FASTQs using " +
                  std::string(simd::name(is)))
    , m_config(config)
    , m_reads(reads)
    , m_adapters(config.adapters.get_adapter_set(0))
    , m_aligner(m_adapters, is)
  {
  }

protected:
  alignment_info execute()
  {
    alignment_info best;

    for (const auto& read : m_reads) {
      const alignment_info alignment =
        m_aligner.align_single_end(read, m_config.shift);

      if (alignment.score() > best.score()) {
        best = alignment;
      }
    }

    return best;
  }

private:
  const userconfig& m_config;
  const fastq_vec& m_reads;
  const fastq_pair_vec m_adapters;
  const sequence_aligner m_aligner;
};

/** Benchmarking of PE alignments */
class pe_alignment_benchmarker : public benchmarker<alignment_info>
{
public:
  pe_alignment_benchmarker(const userconfig& config,
                           const fastq_vec& mate_1,
                           const fastq_vec& mate_2,
                           const simd::instruction_set is)
    : benchmarker("align paired-end FASTQs using " +
                  std::string(simd::name(is)))
    , m_config(config)
    , m_mate_1(mate_1)
    , m_mate_2(mate_2)
    , m_mate_2_reversed(mate_2)
    , m_adapters(config.adapters.get_adapter_set(0))
    , m_aligner(m_adapters, is)
  {
    // Normally done as part of the alignment loop but is benchmarked separately
    for (auto& it : m_mate_2_reversed) {
      it.reverse_complement();
    }
  }

protected:
  alignment_info execute()
  {
    AR_REQUIRE(m_mate_1.size() == m_mate_2.size());

    alignment_info best;

    auto it_1 = m_mate_1.begin();
    auto it_2 = m_mate_2_reversed.begin();
    while (it_1 != m_mate_1.end()) {
      const fastq& read_1 = *it_1++;
      const fastq& read_2 = *it_2++;

      const alignment_info alignment =
        m_aligner.align_paired_end(read_1, read_2, m_config.shift);

      if (alignment.score() > best.score()) {
        best = alignment;
      }
    }

    return best;
  }

private:
  const userconfig& m_config;
  const fastq_vec& m_mate_1;
  const fastq_vec& m_mate_2;
  fastq_vec m_mate_2_reversed;
  const fastq_pair_vec m_adapters;
  const sequence_aligner m_aligner;
};

} // namespace

int
benchmark(const userconfig& config)

{
  benchmark_toggles enabled;
  if (!enabled.parse_args(config.benchmarks)) {
    return 1;
  }

  const string_vec input_1 = readlines_benchmarker(config, config.input_files_1)
                               .run(enabled(benchmarks::read));
  string_vec input_2;

  if (config.input_files_2.size()) {
    input_2 = readlines_benchmarker(config, config.input_files_2)
                .run(enabled(benchmarks::read));
  }

  const auto records_1 =
    fastq_parser_benchmarker(input_1).run(enabled(benchmarks::parse));

  fastq_vec records_2;
  if (input_2.size()) {
    records_2 =
      fastq_parser_benchmarker(input_2).run(enabled(benchmarks::parse));
  }

  reverse_complement_benchmarker(records_1).run(enabled(benchmarks::revcompl));

  for (const auto is : simd::supported()) {
    benchmarker_se_alignment(config, records_1, is)
      .run(enabled.align(benchmarks::align_se, is));
  }

  if (input_2.size()) {
    for (const auto is : simd::supported()) {
      pe_alignment_benchmarker(config, records_1, records_2, is)
        .run(enabled.align(benchmarks::align_pe, is));
    }
  }

  return 0;
}

} // namespace adapterremoval
