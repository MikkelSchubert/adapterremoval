// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t
#include <memory>  // for unique_ptr
#include <string>  // for string
#include <vector>  // for vector

namespace adapterremoval {

class benchmarker_input;
class buffer;
class fastq;
class fastq_encoding;

/** Specifies the input used by the benchmark */
struct benchmarker_input_types
{
  using type = size_t;

  //! Takes no input, i.e. generates data or reads from disk itself
  static constexpr type nothing = 0b0000'0000;
  //! Fixed-sized buffers lines from the input files (normalized to \n)
  static constexpr type buffers = 0b0000'0001;
  //! Unprocessed lines from the input files (normalized to \n)
  static constexpr type lines = 0b0000'0010;
  //! Parsed but not post-processed/validated FASTQ reads
  static constexpr type unsafe_fastq = 0b0000'0100;
  //! Parsed and processed/validated FASTQ reads
  static constexpr type validated_fastq = 0b0000'1000;

  //! All valid values
  static constexpr type mask = 0b0000'1111;
};

class benchmarker_input
{
public:
  /**
   * Create benchmarking input from zero or more FASTQ files, processing data to
   * provide the `required` type of input as well as all types of input required
   * to generate that input type.
   *
   * Note: Files are assumed to not contain any extraneous newlines.
   */
  benchmarker_input(std::vector<std::string> filenames,
                    const fastq_encoding& encoding,
                    uint64_t head,
                    benchmarker_input_types::type required);

  ~benchmarker_input();

  /** Returns the number of reads read from the input files */
  [[nodiscard]] uint64_t head() const { return m_head; }

  /** Returns input filenames, if any */
  [[nodiscard]] const std::vector<std::string>& filenames() const;

  /** Returns all lines from input files; asserts if not loaded */
  [[nodiscard]] const std::vector<std::string>& lines() const;

  /** Returns fixed-sized buffers from input files; asserts if not loaded */
  [[nodiscard]] const std::vector<buffer>& buffers() const;

  /** Returns un-validated FASTQs from input files; asserts if not loaded */
  [[nodiscard]] const std::vector<fastq>& unsafe_fastqs() const;

  /** Returns validated FASTQs from input files; asserts if not loaded */
  [[nodiscard]] const std::vector<fastq>& validated_fastqs() const;

  benchmarker_input(const benchmarker_input&) = delete;
  benchmarker_input(benchmarker_input&&) = default;
  benchmarker_input& operator=(const benchmarker_input&) = delete;
  benchmarker_input& operator=(benchmarker_input&&) = default;

private:
  /** Reads lines from the input files, at most `head * 4` in total */
  std::vector<std::string> read_lines();

  /** Converts lines into fixed-sized buffers */
  std::vector<buffer> make_buffers();

  /** Parses FASTQ records but does not perform post-processing */
  std::vector<fastq> make_fastqs();

  /** Post-processes and validates FASTQ records */
  std::vector<fastq> validate_fastqs(const fastq_encoding& encoding,
                                     bool consume);

  //! The maximum number of reads to load
  uint64_t m_head;

  //! Zero or more FASTQ files
  std::vector<std::string> m_filenames{};

  //! Specifies if lines were loaded from the input files
  bool m_lines_loaded = false;
  //! Lines read from the input files with trailing \n or \r\n removed
  std::vector<std::string> m_lines{};

  //! Specifies if buffers were derived from the loaded lines
  bool m_buffers_loaded = false;
  //! Fixed sizes buffers containing uncompressed FASTQ data
  std::vector<buffer> m_buffers{};

  //! Specifies if FASTQ reads were parsed from the input lines
  bool m_unsafe_fastq_loaded = false;
  //! FASTQ records loaded via `read_unsafe`, but not post-processed/validated
  std::vector<fastq> m_unsafe_fastqs{};

  //! Specifies if FASTQ reads were post-processed and validated
  bool m_validated_fastq_loaded = false;
  //! FASTQ that have been post-processed and validated
  std::vector<fastq> m_validated_fastqs{};
};

/** Base-class for benchmarks */
class benchmarker
{
public:
  benchmarker(std::string desc,
              std::vector<std::string> toggles,
              benchmarker_input_types::type input);

  virtual ~benchmarker() = default;

  /**
   * Returns true if the benchmarker is enabled by the user-specified toggles.
   * An empty vector is interpreted as the "default" selection.
   */
  [[nodiscard]] virtual bool enabled(
    const std::vector<std::string>& toggles) const;

  /** Runs the benchmark with the specified input */
  virtual void run(const benchmarker_input& input_1,
                   const benchmarker_input& input_2);

  /** Returns the toggles recognized by this benchmarker */
  [[nodiscard]] const std::vector<std::string>& toggles() const
  {
    return m_toggles;
  }

  /** Returns the input required by the benchmarker */
  [[nodiscard]] benchmarker_input_types::type input() const { return m_input; }

  benchmarker(const benchmarker&) = delete;
  benchmarker(benchmarker&&) = delete;
  benchmarker& operator=(const benchmarker&) = delete;
  benchmarker& operator=(benchmarker&&) = delete;

protected:
  /** Called before `execute` to perform any per batch setup */
  virtual void setup(const benchmarker_input& input_1,
                     const benchmarker_input& input_2);
  /**
   * Carries out the actual work to be benchmarked. Each call is expected to
   * carry out work on the same number of items to make runs comparable.
   */
  virtual void execute(const benchmarker_input& input_1,
                       const benchmarker_input& input_2) = 0;

private:
  /** Returns a summary of the current benchmarking results */
  [[nodiscard]] std::string summarize(size_t loops) const;

  //! Short description of this benchmark
  const std::string m_description;
  //! Vector of durations from executed benchmarks
  std::vector<uint64_t> m_durations{};
  //! Vector of toggles recognized by this benchmark
  std::vector<std::string> m_toggles{};
  //! Mask of types required by this benchmark
  benchmarker_input_types::type m_input;
};

bool
validate_toggles(const std::vector<std::unique_ptr<benchmarker>>& benchmarks,
                 const std::vector<std::string>& candidates);

} // namespace adapterremoval
