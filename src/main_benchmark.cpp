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
#include "benchmarking.hpp"      // for benchmarker
#include "fastq.hpp"             // for fastq
#include "linereader_joined.hpp" // for joined_line_readers
#include "logging.hpp"           // for log
#include "statistics.hpp"        // for fastq_statistics
#include "strutils.hpp"          // for to_lower
#include "userconfig.hpp"        // for userconfig
#include <cstdint>               // for size_t
#include <sstream>               // for ostringstream
#include <vector>                // for vector

namespace adapterremoval {

/*
 * Pass pointer to results to a different compilation unit, to prevent the
 * compiler from aggressivly eliding otherwise effect-free code. This will
 * probably break if LTO is enabled.
 */
void
blackbox(void* p);

namespace {

class readlines_benchmarker : public benchmarker
{
public:
  readlines_benchmarker(const string_vec& filenames_1,
                        const string_vec& filenames_2,
                        const size_t head)
    : benchmarker("FASTQ reading", { "read" })
    , m_filenames_1(filenames_1)
    , m_filenames_2(filenames_2)
    , m_head(head)
    , m_lines_1()
    , m_lines_2()
  {
    set_required();
  }

  const string_vec& lines_1() const { return m_lines_1; }

  const string_vec& lines_2() const { return m_lines_2; }

protected:
  void setup() override
  {
    m_lines_1 = string_vec();
    m_lines_2 = string_vec();
  }

  void execute() override
  {
    read_lines(m_filenames_1, m_lines_1);
    read_lines(m_filenames_2, m_lines_2);

    blackbox(&m_lines_1);
    blackbox(&m_lines_2);
  }

private:
  void read_lines(const string_vec& filenames, string_vec& lines) const
  {
    joined_line_readers reader(filenames);
    while (lines.size() / 4 < m_head) {
      lines.emplace_back(std::string());
      if (!reader.getline(lines.back())) {
        lines.pop_back();
        break;
      }
    }
  }

  const string_vec m_filenames_1;
  const string_vec m_filenames_2;
  const size_t m_head;
  // Vector containing files set of lines read
  string_vec m_lines_1;
  string_vec m_lines_2;
};

/** Benchmarking of FASTQ parsing excluding file IO */
class fastq_parser_benchmarker : public benchmarker
{
public:
  fastq_parser_benchmarker(const string_vec& lines_1, const string_vec& lines_2)
    : benchmarker("FASTQ parsing", { "parse" })
    , m_lines_1(lines_1)
    , m_lines_2(lines_2)
    , m_records_1()
    , m_records_2()
  {
    set_required();
  }

  const fastq_vec& records_1() const { return m_records_1; }

  const fastq_vec& records_2() const { return m_records_2; }

protected:
  void setup() override
  {
    m_records_1 = fastq_vec();
    m_records_2 = fastq_vec();
  }

  void execute() override
  {
    fastq record;
    {
      vec_reader reader_1(m_lines_1);
      while (record.read(reader_1)) {
        m_records_1.push_back(record);
      }
    }

    {
      vec_reader reader_2(m_lines_2);
      while (record.read(reader_2)) {
        m_records_2.push_back(record);
      }
    }

    blackbox(&m_records_1);
    blackbox(&m_records_2);
  }

private:
  const string_vec& m_lines_1;
  const string_vec& m_lines_2;
  fastq_vec m_records_1;
  fastq_vec m_records_2;
};

class reverse_complement_benchmarker : public benchmarker
{
public:
  reverse_complement_benchmarker(const fastq_vec& records_1,
                                 const fastq_vec& records_2)
    : benchmarker("reverse complement", { "revcompl" })
    , m_records_1(records_1)
    , m_records_2(records_2)
  {
  }

protected:
  void execute() override
  {
    for (auto& it : m_records_1) {
      it.reverse_complement();
    }

    for (auto& it : m_records_2) {
      it.reverse_complement();
    }
  }

private:
  fastq_vec m_records_1;
  fastq_vec m_records_2;
};

class complexity_benchmarker : public benchmarker
{
public:
  complexity_benchmarker(const fastq_vec& records_1, const fastq_vec& records_2)
    : benchmarker("read complexity", { "complexity" })
    , m_records_1(records_1)
    , m_records_2(records_2)
  {
  }

protected:
  void execute() override
  {
    double total = 0.0;
    total += complexity(m_records_1);
    total += complexity(m_records_2);

    blackbox(&total);
  }

private:
  double complexity(const fastq_vec& records) const
  {
    double total = 0.0;
    for (auto& it : records) {
      total += it.complexity();
    }

    return total;
  }

  const fastq_vec& m_records_1;
  const fastq_vec& m_records_2;
};

class trimming_benchmarker : public benchmarker
{
public:
  trimming_benchmarker(const std::string& desc,
                       const std::string& toggle,
                       const fastq_vec& records_1,
                       const fastq_vec& records_2)
    : benchmarker(desc, { "trim", "trim:" + toggle })
    , m_records_1(records_1)
    , m_records_2(records_2)
    , m_trimmed_records_1()
    , m_trimmed_records_2()
  {
  }

protected:
  void setup() override
  {
    m_trimmed_records_1 = m_records_1;
    m_trimmed_records_2 = m_records_2;
  }

  void execute() override
  {
    trim(m_trimmed_records_1);
    trim(m_trimmed_records_2);

    blackbox(&m_trimmed_records_1);
    blackbox(&m_trimmed_records_2);
  }

  virtual void trim(fastq_vec& reads) const = 0;

private:
  const fastq_vec& m_records_1;
  const fastq_vec& m_records_2;
  fastq_vec m_trimmed_records_1;
  fastq_vec m_trimmed_records_2;
};

class basic_trimming_benchmarker : public trimming_benchmarker
{
public:
  basic_trimming_benchmarker(const fastq_vec& records_1,
                             const fastq_vec& records_2)
    : trimming_benchmarker("basic trimming", "basic", records_1, records_2)
  {
  }

protected:
  void trim(fastq_vec& reads) const override
  {
    for (auto& read : reads) {
      read.trim_trailing_bases(true, 2);
    }
  }
};

class mott_trimming_benchmarker : public trimming_benchmarker
{
public:
  mott_trimming_benchmarker(const fastq_vec& records_1,
                            const fastq_vec& records_2)
    : trimming_benchmarker("mott trimming", "mott", records_1, records_2)
  {
  }

protected:
  void trim(fastq_vec& reads) const override
  {
    for (auto& read : reads) {
      read.mott_trimming(0.05);
    }
  }
};

class window_trimming_benchmarker : public trimming_benchmarker
{
public:
  window_trimming_benchmarker(const fastq_vec& records_1,
                              const fastq_vec& records_2)
    : trimming_benchmarker("window trimming", "window", records_1, records_2)
  {
  }

protected:
  void trim(fastq_vec& reads) const override
  {
    for (auto& read : reads) {
      read.trim_windowed_bases(true, 2, 0.1);
    }
  }
};

/** Class for benchmarking collection of `fastq_statistics` */
class fastq_statistics_benchmarker : public benchmarker
{
public:
  fastq_statistics_benchmarker(const fastq_vec& records_1,
                               const fastq_vec& records_2)
    : benchmarker("read statistics", { "stats" })
    , m_records_1(records_1)
    , m_records_2(records_2)
  {
  }

  void execute() override
  {
    collect_statistics(m_records_1);
    collect_statistics(m_records_2);
  }

private:
  void collect_statistics(const fastq_vec& records) const
  {
    fastq_statistics stats;
    for (const auto& it : records) {
      stats.process(it);
    }

    blackbox(&stats);
  }

  const fastq_vec& m_records_1;
  const fastq_vec& m_records_2;
};

/** Base-class for benchmarking SE/PE alignments */
class alignment_benchmarker : public benchmarker
{
public:
  alignment_benchmarker(const std::string& key, const simd::instruction_set is)
    : benchmarker(to_upper(key) + " alignment (" + simd::name(is) + ")", {})
    , m_key(key)
    , m_is(is)
  {
  }

  strategy enabled(const benchmark_toggles& toggles) const override
  {
    if (toggles.defaults() || toggles.is_set("align") ||
        toggles.is_set("align:" + m_key)) {

      if (toggles.is_set("simd") ||
          toggles.is_set("simd:" + to_lower(simd::name(m_is)))) {
        return strategy::benchmark;
      }

      // Benchmark the prefered algorithm if no algorithms were specified
      const auto supported = simd::supported();
      if (supported.size() && supported.back() == m_is) {
        for (const auto is : supported) {
          if (toggles.is_set("simd:" + to_lower(simd::name(is)))) {
            return strategy::skip;
          }
        }

        return strategy::benchmark;
      }
    }

    return strategy::skip;
  }

private:
  const std::string m_key;
  const simd::instruction_set m_is;
};

/** Benchmarking of SE alignments */
class benchmarker_se_alignment : public alignment_benchmarker

{
public:
  benchmarker_se_alignment(const userconfig& config,
                           const fastq_vec& reads,
                           const simd::instruction_set is)
    : alignment_benchmarker("se", is)
    , m_config(config)
    , m_reads(reads)
    , m_adapters(config.adapters.get_adapter_set(0))
    , m_aligner(m_adapters, is)
  {
  }

protected:
  void execute() override
  {
    alignment_info best;

    for (const auto& read : m_reads) {
      const alignment_info alignment =
        m_aligner.align_single_end(read, m_config.shift);

      if (alignment.score() > best.score()) {
        best = alignment;
      }
    }

    blackbox(&best);
  }

private:
  const userconfig& m_config;
  const fastq_vec& m_reads;
  const fastq_pair_vec m_adapters;
  const sequence_aligner m_aligner;
};

/** Benchmarking of PE alignments */
class pe_alignment_benchmarker : public alignment_benchmarker
{
public:
  pe_alignment_benchmarker(const userconfig& config,
                           const fastq_vec& mate_1,
                           const fastq_vec& mate_2,
                           const simd::instruction_set is)
    : alignment_benchmarker("pe", is)
    , m_config(config)
    , m_mate_1(mate_1)
    , m_mate_2(mate_2)
    , m_mate_2_reversed()
    , m_adapters(config.adapters.get_adapter_set(0))
    , m_aligner(m_adapters, is)
  {
  }

protected:
  void setup() override
  {
    if (m_mate_2_reversed.empty()) {
      m_mate_2_reversed = m_mate_2;

      // Done as part of the alignment loop but is benchmarked separately
      for (auto& it : m_mate_2_reversed) {
        it.reverse_complement();
      }
    }
  }

  void execute() override
  {
    AR_REQUIRE(m_mate_1.size() == m_mate_2_reversed.size());

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

    blackbox(&best);
  }

private:
  const userconfig& m_config;
  const fastq_vec& m_mate_1;
  const fastq_vec& m_mate_2;
  fastq_vec m_mate_2_reversed;
  const fastq_pair_vec m_adapters;
  const sequence_aligner m_aligner;
};

string_vec
supported_toggles()
{
  string_vec toggles = { "read",  "parse",      "revcompl",  "complexity",
                         "trim",  "trim:basic", "trim:mott", "trim:window",
                         "stats", "align",      "align:se",  "align:pe",
                         "simd" };

  for (const auto is : simd::supported()) {
    toggles.push_back(std::string("simd:") + to_lower(simd::name(is)));
  }

  return toggles;
}

} // namespace

int
benchmark(const userconfig& config)
{
  auto head = config.head;
  // Limit benchmarking to a reasonable data-set size by default
  if (config.head == std::numeric_limits<uint64_t>::max()) {
    log::warn() << "Defaulting to reading at most 1,000,000 reads/mate pairs";
    head = 1000000;
  }

  // Parse user-specified benchmarking toggles
  benchmark_toggles toggles(supported_toggles());
  if (!toggles.update_toggles(config.benchmarks)) {
    return 1;
  }

  auto lines =
    readlines_benchmarker(config.input_files_1, config.input_files_2, head);
  lines.run_if_toggled(toggles);

  auto records = fastq_parser_benchmarker(lines.lines_1(), lines.lines_2());
  records.run_if_toggled(toggles);

  reverse_complement_benchmarker(records.records_1(), records.records_2())
    .run_if_toggled(toggles);

  complexity_benchmarker(records.records_1(), records.records_2())
    .run_if_toggled(toggles);

  basic_trimming_benchmarker(records.records_1(), records.records_2())
    .run_if_toggled(toggles);

  mott_trimming_benchmarker(records.records_1(), records.records_2())
    .run_if_toggled(toggles);

  window_trimming_benchmarker(records.records_1(), records.records_2())
    .run_if_toggled(toggles);

  fastq_statistics_benchmarker(records.records_1(), records.records_2())
    .run_if_toggled(toggles);

  for (const auto is : simd::supported()) {
    benchmarker_se_alignment(config, records.records_1(), is)
      .run_if_toggled(toggles);
  }

  if (config.paired_ended_mode) {
    for (const auto is : simd::supported()) {
      pe_alignment_benchmarker(
        config, records.records_1(), records.records_2(), is)
        .run_if_toggled(toggles);
    }
  }

  return 0;
}

} // namespace adapterremoval
