// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "alignment.hpp"     // for alignment_info, sequence_aligner
#include "benchmarking.hpp"  // for benchmarker, benchmarker_input, ...
#include "buffer.hpp"        // for buffer
#include "buffer_reader.hpp" // for buffer_reader
#include "debug.hpp"         // for AR_REQUIRE
#include "fastq.hpp"         // for fastq
#include "fastq_enc.hpp"     // fastq_encoding
#include "fastq_in.hpp"      // for INPUT_BLOCK_SIZE, reusable_buffers, ...
#include "linereader.hpp"
#include "logging.hpp"       // for error, info, warn
#include "main.hpp"          // declarations
#include "sequence_sets.hpp" // for adapter_set
#include "simd.hpp"          // for name, supported, instruction_set, ...
#include "statistics.hpp"    // for fastq_statistics
#include "strutils.hpp"      // for to_lower, to_upper
#include "userconfig.hpp"    // for userconfig
#include "utilities.hpp"     // for blackbox
#include <algorithm>         // for find
#include <cstddef>           // for size_t
#include <cstdint>           // for uint64_t
#include <limits>            // for numeric_limits
#include <memory>            // for make_unique, unique_ptr
#include <string>            // for string
#include <string_view>       // for string_view
#include <utility>           // for move
#include <vector>            // for vector

namespace adapterremoval {

namespace {

class readlines_benchmarker : public benchmarker
{
public:
  readlines_benchmarker()
    : benchmarker("line reading", { "read" }, benchmarker_input_types::nothing)
  {
  }

protected:
  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    read_lines(input_1, m_lines_1);
    read_lines(input_2, m_lines_2);

    blackbox(m_lines_1);
    blackbox(m_lines_2);

    return ((m_lines_1.size() + 3) / 4) + ((m_lines_2.size() + 3) / 4);
  }

private:
  static void read_lines(const benchmarker_input& input,
                         std::vector<std::string>& lines)
  {
    if (!input.filenames().empty()) {
      joined_line_readers reader(input.filenames());

      const size_t head = input.head();
      if (lines.empty()) {
        // first loop (burn-in) determines the actual number of lines to read
        for (size_t i = 0; i < head * 4; i++) {
          lines.emplace_back();
          if (!reader.getline(lines.back())) {
            lines.pop_back();
            break;
          }
        }
      } else {
        // subsequent loops reuse strings to match actual usage via `read_fastq`
        for (auto& line : lines) {
          if (!reader.getline(line)) {
            break;
          }
        }
      }
    }
  }

  // Vector containing files set of lines read
  std::vector<std::string> m_lines_1{};
  std::vector<std::string> m_lines_2{};
};

/** Benchmarking line-reading from buffers */
class buffer_reader_benchmarker : public benchmarker
{
public:
  buffer_reader_benchmarker()
    : benchmarker("buffer reader",
                  { "readbuf" },
                  benchmarker_input_types::buffers)
    , m_buffers("buffers", INPUT_BLOCK_SIZE)
  {
  }

protected:
  void setup(const benchmarker_input& input_1,
             const benchmarker_input& input_2) override
  {
    copy_buffers(input_1.buffers(), m_working_buffers_1);
    copy_buffers(input_2.buffers(), m_working_buffers_2);
  }

  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    parse_buffers(m_working_buffers_1);
    parse_buffers(m_working_buffers_2);

    return input_1.head() + input_2.head();
  }

private:
  void parse_buffers(std::vector<buffer>& src)
  {
    buffer_reader reader{ { "dummy" } };
    reader.set_release_function(
      [this](buffer&& buf) { m_buffers.release(std::move(buf)); });

    std::string dst;
    for (auto&& buf : src) {
      reader.add(std::move(buf), &buf == &src.back());

      while (reader.getline(dst) == buffer_reader::result::complete) {
        blackbox(dst);
      }
    }
  }

  /** Creates disposable copy of the input buffers */
  void copy_buffers(const std::vector<buffer>& src, std::vector<buffer>& dst)
  {
    dst.clear();
    dst.reserve(src.size());
    for (const auto& buf : src) {
      buffer copy = m_buffers.acquire(0);
      copy.append(buf);

      dst.emplace_back(std::move(copy));
    }
  }

  // Reused buffers; used to eliminate overhead from malloc/free
  reusable_buffers m_buffers;

  //! Copy of buffers in m_buffers_1, to be consumed during benchmarking
  std::vector<buffer> m_working_buffers_1{};
  //! Copy of buffers in m_buffers_2, to be consumed during benchmarking
  std::vector<buffer> m_working_buffers_2{};
};

/** Benchmarking of FASTQ parsing excluding file IO */
class fastq_parser_benchmarker : public benchmarker
{
public:
  fastq_parser_benchmarker(fastq_encoding encoding)
    : benchmarker("FASTQ parsing", { "parse" }, benchmarker_input_types::lines)
    , m_encoding(std::move(encoding))
  {
  }

protected:
  void setup(const benchmarker_input& input_1,
             const benchmarker_input& input_2) override
  {
    if (m_reader_1) {
      (*m_reader_1).reset();
      (*m_reader_2).reset();
    } else {
      m_reader_1 =
        std::make_unique<borrowed_vec_reader>(input_1.lines().begin(),
                                              input_1.lines().end());
      m_reader_2 =
        std::make_unique<borrowed_vec_reader>(input_2.lines().begin(),
                                              input_2.lines().end());

      m_records_1.resize(input_1.head());
      m_records_2.resize(input_2.head());
    }
  }

  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    // FASTQ records are reused to match actual usage via `read_fastq`
    for (auto& record : m_records_1) {
      AR_REQUIRE(record.read(*m_reader_1, m_encoding));
    }

    // FASTQ records are reused to match actual usage via `read_fastq`
    for (auto& record : m_records_2) {
      AR_REQUIRE(record.read(*m_reader_2, m_encoding));
    }

    blackbox(m_records_1);
    blackbox(m_records_2);

    return input_1.head() + input_2.head();
  }

private:
  const fastq_encoding m_encoding;
  std::unique_ptr<borrowed_vec_reader> m_reader_1{};
  std::unique_ptr<borrowed_vec_reader> m_reader_2{};
  std::vector<fastq> m_records_1{};
  std::vector<fastq> m_records_2{};
};

class reverse_complement_benchmarker : public benchmarker
{
public:
  reverse_complement_benchmarker()
    : benchmarker("reverse complement",
                  { "revcompl" },
                  benchmarker_input_types::validated_fastq)
  {
  }

protected:
  void setup(const benchmarker_input& input_1,
             const benchmarker_input& input_2) override
  {
    // To reduce overhead, and since it makes no practical difference, the same
    // fastq objects are reverse complemented each loop
    if (m_records_1.empty()) {
      m_records_1 = input_1.validated_fastqs();
      m_records_2 = input_2.validated_fastqs();
    }
  }

  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    for (auto& it : m_records_1) {
      it.reverse_complement();
    }

    for (auto& it : m_records_2) {
      it.reverse_complement();
    }

    return input_1.head() + input_2.head();
  }

private:
  std::vector<fastq> m_records_1{};
  std::vector<fastq> m_records_2{};
};

class complexity_benchmarker : public benchmarker
{
public:
  complexity_benchmarker()
    : benchmarker("read complexity",
                  { "complexity" },
                  benchmarker_input_types::validated_fastq)
  {
  }

protected:
  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    double total = 0.0;
    total += complexity(input_1.validated_fastqs());
    total += complexity(input_2.validated_fastqs());

    blackbox(total);

    return input_1.head() + input_2.head();
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
};

class trimming_benchmarker : public benchmarker
{
public:
  trimming_benchmarker(const std::string& desc, const std::string& toggle)
    : benchmarker(desc,
                  { "trim", "trim:" + toggle },
                  benchmarker_input_types::validated_fastq)
  {
  }

protected:
  void setup(const benchmarker_input& input_1,
             const benchmarker_input& input_2) override
  {
    m_trimmed_records_1 = input_1.validated_fastqs();
    m_trimmed_records_2 = input_2.validated_fastqs();
  }

  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    trim(m_trimmed_records_1);
    trim(m_trimmed_records_2);

    blackbox(m_trimmed_records_1);
    blackbox(m_trimmed_records_2);

    return input_1.head() + input_2.head();
  }

  virtual void trim(std::vector<fastq>& reads) const = 0;

private:
  std::vector<fastq> m_trimmed_records_1{};
  std::vector<fastq> m_trimmed_records_2{};
};

class basic_trimming_benchmarker : public trimming_benchmarker
{
public:
  basic_trimming_benchmarker()
    : trimming_benchmarker("basic trimming", "basic")
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
  mott_trimming_benchmarker()
    : trimming_benchmarker("mott trimming", "mott")
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
  window_trimming_benchmarker()
    : trimming_benchmarker("window trimming", "window")
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
  poly_x_trimming_benchmarker(std::string nucleotides, size_t threshold)
    : trimming_benchmarker("poly-X trimming", "polyx")
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

/** Class for benchmarking collection of FASTQ statistics */
class fastq_statistics_benchmarker : public benchmarker
{
public:
  fastq_statistics_benchmarker()
    : benchmarker("read statistics",
                  { "stats" },
                  benchmarker_input_types::validated_fastq)
  {
  }

protected:
  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    collect_statistics(input_1.validated_fastqs());
    collect_statistics(input_2.validated_fastqs());

    return input_1.head() + input_2.head();
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
};

/** Base-class for benchmarking SE/PE alignments */
class alignment_benchmarker : public benchmarker
{
public:
  alignment_benchmarker(std::string key,
                        const simd::instruction_set is,
                        const bool default_is)
    : benchmarker(
        to_upper(key) + " alignment (" + std::string{ simd::name(is) } + ")",
        { "align", "align:" + key, "simd", "simd:" + to_lower(simd::name(is)) },
        benchmarker_input_types::validated_fastq)
    , m_key(std::move(key))
    , m_is(is)
    , m_default_is(default_is)
  {
  }

  [[nodiscard]] bool enabled(
    const std::vector<std::string>& toggles) const override
  {
    // Only enable the benchmarks for the preferred SIMD by default
    if (toggles.empty()) {
      return m_default_is;
    }

    const auto toggle_enabled = [&toggles](const std::string_view& toggle) {
      return std::find(toggles.begin(), toggles.end(), toggle) != toggles.end();
    };

    const auto align_key = "align:" + m_key;
    if (!toggle_enabled("align") && !toggle_enabled(align_key)) {
      return false;
    }

    // `simd` is purely used to subselect between enabled alignment benchmarks
    const auto simd_key = "simd:" + to_lower(simd::name(m_is));
    if (toggle_enabled("simd") || toggle_enabled(simd_key)) {
      return true;
    }

    // enable unless another instruction set was requested
    if (m_default_is) {
      for (const auto& toggle : toggles) {
        if (starts_with(toggle, "simd:")) {
          return false;
        }
      }

      return true;
    }

    return false;
  }

private:
  //! Key describing the kind of alignment, e.g. 'se' or 'pe'
  const std::string m_key;
  const simd::instruction_set m_is;
  const bool m_default_is;
};

/** Benchmarking of SE alignments */
class se_alignment_benchmarker : public alignment_benchmarker
{
public:
  se_alignment_benchmarker(const userconfig& config,
                           const simd::instruction_set is,
                           const bool default_is)
    : alignment_benchmarker("se", is, default_is)
    , m_config(config)
    , m_adapters(config.samples.get_reader()->adapters())
    , m_aligner(m_adapters, is, config.mismatch_threshold)
  {
    m_aligner.set_min_overlap(config.min_overlap);
  }

protected:
  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& /* input_2 */) override
  {
    alignment_info best;

    for (const auto& read : input_1.validated_fastqs()) {
      const alignment_info alignment =
        m_aligner.align_single_end(read, m_config.shift);

      if (alignment.score() > best.score()) {
        best = alignment;
      }
    }

    blackbox(best);

    return input_1.head();
  }

private:
  const userconfig& m_config;
  const adapter_set m_adapters;
  sequence_aligner m_aligner;
};

/** Benchmarking of PE alignments */
class pe_alignment_benchmarker : public alignment_benchmarker
{
public:
  pe_alignment_benchmarker(const userconfig& config,
                           const simd::instruction_set is,
                           const bool default_is)
    : alignment_benchmarker("pe", is, default_is)
    , m_config(config)
    , m_adapters(config.samples.get_reader()->adapters())
    , m_aligner(m_adapters, is, config.mismatch_threshold)
  {
  }

protected:
  void setup(const benchmarker_input& /*  input_1 */,
             const benchmarker_input& input_2) override
  {
    if (m_mate_2_reversed.empty()) {
      m_mate_2_reversed = input_2.validated_fastqs();

      // Done as part of the alignment loop but is benchmarked separately
      for (auto& it : m_mate_2_reversed) {
        it.reverse_complement();
      }
    }
  }

  uint64_t execute(const benchmarker_input& input_1,
                   const benchmarker_input& input_2) override
  {
    const auto& mate_1 = input_1.validated_fastqs();
    const auto& mate_2 = m_mate_2_reversed;
    AR_REQUIRE(mate_1.size() == mate_2.size());

    alignment_info best;

    auto it_1 = mate_1.begin();
    auto it_2 = mate_2.begin();
    while (it_1 != mate_1.end()) {
      const fastq& read_1 = *it_1++;
      const fastq& read_2 = *it_2++;

      const alignment_info alignment =
        m_aligner.align_paired_end(read_1, read_2, m_config.shift);

      if (alignment.score() > best.score()) {
        best = alignment;
      }
    }

    blackbox(best);

    return input_1.head() + input_2.head();
  }

private:
  const userconfig& m_config;
  std::vector<fastq> m_mate_2_reversed{};
  const adapter_set m_adapters;
  sequence_aligner m_aligner;
};

template<typename T, typename... Args>
void
add(std::vector<std::unique_ptr<benchmarker>>& dst, Args&&... args)
{
  static_assert(std::is_base_of<benchmarker, T>(),
                "requires benchmarker sub-class");
  dst.emplace_back(std::make_unique<T>(std::forward<Args>(args)...));
}

} // namespace

int
benchmark(const userconfig& config)
{
  auto head = config.head;
  // Limit benchmarking to a reasonable data-set size by default
  if (config.head == std::numeric_limits<uint64_t>::max()) {
    log::warn() << "Defaulting to reading at most 1,000,000 reads/mate pairs";
    head = 1'000'000;
  } else if (config.head == 0) {
    log::error() << "Cannot perform benchmarking with 0 reads/mate pairs";
    return 1;
  }

  std::vector<std::unique_ptr<benchmarker>> benchmarks;
  add<readlines_benchmarker>(benchmarks);
  add<buffer_reader_benchmarker>(benchmarks);
  add<fastq_parser_benchmarker>(benchmarks, config.io_encoding);
  add<reverse_complement_benchmarker>(benchmarks);
  add<complexity_benchmarker>(benchmarks);
  add<basic_trimming_benchmarker>(benchmarks);
  add<mott_trimming_benchmarker>(benchmarks);
  add<window_trimming_benchmarker>(benchmarks);

  {
    auto nucleotides = *config.pre_trim_poly_x.get_reader();
    if (nucleotides == "auto") {
      nucleotides = "ACGT";
    }

    add<poly_x_trimming_benchmarker>(benchmarks,
                                     nucleotides,
                                     config.trim_poly_x_threshold);
  }

  add<fastq_statistics_benchmarker>(benchmarks);

  // Benchmark the preferred algorithm if no algorithms were specified
  const auto instruction_sets = simd::supported();
  for (const auto is : instruction_sets) {
    add<se_alignment_benchmarker>(benchmarks,
                                  config,
                                  is,
                                  is == instruction_sets.back());
  }

  if (config.paired_ended_mode) {
    for (const auto is : instruction_sets) {
      add<pe_alignment_benchmarker>(benchmarks,
                                    config,
                                    is,
                                    is == instruction_sets.back());
    }
  }

  if (!validate_toggles(benchmarks, config.benchmarks)) {
    return 1;
  }

  benchmarker_input_types::type input = benchmarker_input_types::nothing;

  {
    std::vector<std::unique_ptr<benchmarker>> enabled_benchmarks;
    for (auto& benchmark : benchmarks) {
      if (benchmark->enabled(config.benchmarks)) {
        input |= benchmark->input();
        enabled_benchmarks.emplace_back(std::move(benchmark));
      }
    }

    benchmarks = std::move(enabled_benchmarks);
  }

  auto encoding = config.io_encoding;
  log::info() << "Pre-loading input data";
  benchmarker_input input_1{ config.input_files_1, encoding, head, input };
  benchmarker_input input_2{ config.input_files_2, encoding, head, input };

  if (input_1.head() != input_2.head() && !input_2.filenames().empty()) {
    log::error() << "Mismatching PE input: --in-file1 has " << input_1.head()
                 << " reads, --in-file2 has " << input_2.head() << " reads";
    return 1;
  }

  int returncode = 0;
  log::info() << "Running " << benchmarks.size() << " benchmark(s)";
  for (auto& benchmark : benchmarks) {
    if (!benchmark->run(input_1, input_2)) {
      returncode = 1;
    }

    benchmark.reset(); // free any intermediate values
  }

  return returncode;
}

} // namespace adapterremoval
