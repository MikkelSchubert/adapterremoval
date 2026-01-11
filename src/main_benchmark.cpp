// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "alignment.hpp"         // for alignment_info, sequence_aligner
#include "benchmarking.hpp"      // for benchmarker
#include "debug.hpp"             // for AR_REQUIRE
#include "fastq.hpp"             // for fastq
#include "fastq_enc.hpp"         // for FASTQ_ENCODING_33
#include "linereader.hpp"        // for vec_reader
#include "linereader_joined.hpp" // for joined_line_readers
#include "logging.hpp"           // for log
#include "main.hpp"              // declarations
#include "sequence_sets.hpp"     // for adapter_set
#include "simd.hpp"              // for name, supported, instruction_set (p...
#include "statistics.hpp"        // for fastq_statistics
#include "strutils.hpp"          // for to_lower
#include "threading.hpp"         // for threadsafe_data
#include "userconfig.hpp"        // for userconfig
#include "utilities.hpp"         // for blackbox
#include <cstddef>               // for size_t
#include <cstdint>               // for uint64_t
#include <limits>                // for numeric_limits
#include <string>                // for string+, char_traits
#include <string_view>           // for string_view
#include <utility>               // for move
#include <vector>                // for vector

namespace adapterremoval {

namespace {

class readlines_benchmarker : public benchmarker
{
public:
  readlines_benchmarker(string_vec filenames_1,
                        string_vec filenames_2,
                        const size_t head)
    : benchmarker("FASTQ reading", { "read" })
    , m_filenames_1(std::move(filenames_1))
    , m_filenames_2(std::move(filenames_2))
    , m_head(head)
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

    blackbox(m_lines_1);
    blackbox(m_lines_2);
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

  string_vec m_filenames_1{};
  string_vec m_filenames_2{};
  size_t m_head = 0;
  // Vector containing files set of lines read
  string_vec m_lines_1{};
  string_vec m_lines_2{};
};

/** Benchmarking of FASTQ parsing excluding file IO */
class fastq_parser_benchmarker : public benchmarker
{
public:
  fastq_parser_benchmarker(const string_vec& lines_1, const string_vec& lines_2)
    : benchmarker("FASTQ parsing", { "parse" })
    , m_lines_1(lines_1)
    , m_lines_2(lines_2)
  {
    set_required();
  }

  const std::vector<fastq>& records_1() const { return m_records_1; }

  const std::vector<fastq>& records_2() const { return m_records_2; }

protected:
  void setup() override
  {
    m_records_1 = std::vector<fastq>();
    m_records_2 = std::vector<fastq>();
  }

  void execute() override
  {
    fastq record;
    {
      vec_reader reader_1(m_lines_1);
      while (record.read(reader_1, FASTQ_ENCODING_33)) {
        m_records_1.push_back(record);
      }
    }

    {
      vec_reader reader_2(m_lines_2);
      while (record.read(reader_2, FASTQ_ENCODING_33)) {
        m_records_2.push_back(record);
      }
    }

    blackbox(m_records_1);
    blackbox(m_records_2);
  }

private:
  const string_vec& m_lines_1;
  const string_vec& m_lines_2;
  std::vector<fastq> m_records_1{};
  std::vector<fastq> m_records_2{};
};

class reverse_complement_benchmarker : public benchmarker
{
public:
  reverse_complement_benchmarker(std::vector<fastq> records_1,
                                 std::vector<fastq> records_2)
    : benchmarker("reverse complement", { "revcompl" })
    , m_records_1(std::move(records_1))
    , m_records_2(std::move(records_2))
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
  std::vector<fastq> m_records_1{};
  std::vector<fastq> m_records_2{};
};

class complexity_benchmarker : public benchmarker
{
public:
  complexity_benchmarker(const std::vector<fastq>& records_1,
                         const std::vector<fastq>& records_2)
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

    blackbox(total);
  }

private:
  static double complexity(const std::vector<fastq>& records)
  {
    double total = 0.0;
    for (const auto& it : records) {
      total += it.complexity();
    }

    return total;
  }

  const std::vector<fastq>& m_records_1;
  const std::vector<fastq>& m_records_2;
};

class trimming_benchmarker : public benchmarker
{
public:
  trimming_benchmarker(const std::string& desc,
                       const std::string& toggle,
                       const std::vector<fastq>& records_1,
                       const std::vector<fastq>& records_2)
    : benchmarker(desc, { "trim", "trim:" + toggle })
    , m_records_1(records_1)
    , m_records_2(records_2)
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

    blackbox(m_trimmed_records_1);
    blackbox(m_trimmed_records_2);
  }

  virtual void trim(std::vector<fastq>& reads) const = 0;

private:
  const std::vector<fastq>& m_records_1;
  const std::vector<fastq>& m_records_2;
  std::vector<fastq> m_trimmed_records_1{};
  std::vector<fastq> m_trimmed_records_2{};
};

class basic_trimming_benchmarker : public trimming_benchmarker
{
public:
  basic_trimming_benchmarker(const std::vector<fastq>& records_1,
                             const std::vector<fastq>& records_2)
    : trimming_benchmarker("basic trimming", "basic", records_1, records_2)
  {
  }

protected:
  void trim(std::vector<fastq>& reads) const override
  {
    for (auto& read : reads) {
      read.trim_trailing_bases(true, 2);
    }
  }
};

class mott_trimming_benchmarker : public trimming_benchmarker
{
public:
  mott_trimming_benchmarker(const std::vector<fastq>& records_1,
                            const std::vector<fastq>& records_2)
    : trimming_benchmarker("mott trimming", "mott", records_1, records_2)
  {
  }

protected:
  void trim(std::vector<fastq>& reads) const override
  {
    for (auto& read : reads) {
      read.mott_trimming(0.05);
    }
  }
};

class window_trimming_benchmarker : public trimming_benchmarker
{
public:
  window_trimming_benchmarker(const std::vector<fastq>& records_1,
                              const std::vector<fastq>& records_2)
    : trimming_benchmarker("window trimming", "window", records_1, records_2)
  {
  }

protected:
  void trim(std::vector<fastq>& reads) const override
  {
    for (auto& read : reads) {
      read.trim_windowed_bases(true, 2, 0.1);
    }
  }
};

class poly_x_trimming_benchmarker : public trimming_benchmarker
{
public:
  poly_x_trimming_benchmarker(const std::vector<fastq>& records_1,
                              const std::vector<fastq>& records_2,
                              std::string nucleotides,
                              size_t threshold)
    : trimming_benchmarker("poly-X trimming", "polyx", records_1, records_2)
    , m_nucleotides(std::move(nucleotides))
    , m_threshold(threshold)
  {
  }

protected:
  void trim(std::vector<fastq>& reads) const override
  {
    for (auto& read : reads) {
      read.poly_x_trimming(m_nucleotides, m_threshold);
    }
  }

  std::string m_nucleotides;
  size_t m_threshold;
};

/** Class for benchmarking collection of `fastq_statistics` */
class fastq_statistics_benchmarker : public benchmarker
{
public:
  fastq_statistics_benchmarker(const std::vector<fastq>& records_1,
                               const std::vector<fastq>& records_2)
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
  static void collect_statistics(const std::vector<fastq>& records)
  {
    fastq_statistics stats;
    for (const auto& it : records) {
      stats.process(it);
    }

    blackbox(stats);
  }

  const std::vector<fastq>& m_records_1;
  const std::vector<fastq>& m_records_2;
};

/** Base-class for benchmarking SE/PE alignments */
class alignment_benchmarker : public benchmarker
{
public:
  alignment_benchmarker(std::string_view key, const simd::instruction_set is)
    : benchmarker(to_upper(key) + " alignment (" +
                    std::string{ simd::name(is) } + ")",
                  {})
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

      // Benchmark the preferred algorithm if no algorithms were specified
      const auto supported = simd::supported();
      if (!supported.empty() && supported.back() == m_is) {
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
                           const std::vector<fastq>& reads,
                           const simd::instruction_set is)
    : alignment_benchmarker("se", is)
    , m_config(config)
    , m_reads(reads)
    , m_adapters(config.samples.get_reader()->adapters())
    , m_aligner(m_adapters, is, config.mismatch_threshold)
  {
    m_aligner.set_min_se_overlap(config.min_adapter_overlap);
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

    blackbox(best);
  }

private:
  const userconfig& m_config;
  const std::vector<fastq>& m_reads;
  const adapter_set m_adapters;
  sequence_aligner m_aligner;
};

/** Benchmarking of PE alignments */
class pe_alignment_benchmarker : public alignment_benchmarker
{
public:
  pe_alignment_benchmarker(const userconfig& config,
                           const std::vector<fastq>& mate_1,
                           const std::vector<fastq>& mate_2,
                           const simd::instruction_set is)
    : alignment_benchmarker("pe", is)
    , m_config(config)
    , m_mate_1(mate_1)
    , m_mate_2(mate_2)
    , m_adapters(config.samples.get_reader()->adapters())
    , m_aligner(m_adapters, is, config.mismatch_threshold)
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

    blackbox(best);
  }

private:
  const userconfig& m_config;
  const std::vector<fastq>& m_mate_1;
  const std::vector<fastq>& m_mate_2;
  std::vector<fastq> m_mate_2_reversed{};
  adapter_set m_adapters;
  sequence_aligner m_aligner;
};

string_vec
supported_toggles()
{
  string_vec toggles = { "read",       "parse",      "revcompl",  "complexity",
                         "trim",       "trim:basic", "trim:mott", "trim:window",
                         "trim:polyx", "stats",      "align",     "align:se",
                         "align:pe",   "simd" };

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

  readlines_benchmarker lines{ config.input_files_1,
                               config.input_files_2,
                               head };
  lines.run_if_toggled(toggles);

  fastq_parser_benchmarker records{ lines.lines_1(), lines.lines_2() };
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

  {
    auto nucleotides = *config.pre_trim_poly_x.get_reader();
    if (nucleotides == "auto") {
      nucleotides = "ACGT";
    }

    poly_x_trimming_benchmarker(records.records_1(),
                                records.records_2(),
                                nucleotides,
                                config.trim_poly_x_threshold)
      .run_if_toggled(toggles);
  }

  fastq_statistics_benchmarker(records.records_1(), records.records_2())
    .run_if_toggled(toggles);

  for (const auto is : simd::supported()) {
    benchmarker_se_alignment(config, records.records_1(), is)
      .run_if_toggled(toggles);
  }

  if (config.paired_ended_mode) {
    for (const auto is : simd::supported()) {
      pe_alignment_benchmarker(config,
                               records.records_1(),
                               records.records_2(),
                               is)
        .run_if_toggled(toggles);
    }
  }

  return 0;
}

} // namespace adapterremoval
