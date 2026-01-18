// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "argparse.hpp"      // for parser, parse_result (ptr only)
#include "commontypes.hpp"   // for string_vec, output_format, merge_strategy
#include "fastq_enc.hpp"     // for fastq_encoding, FASTQ_ENCODING_33
#include "sequence_sets.hpp" // for sample_set
#include "simd.hpp"          // for instruction_set
#include "threading.hpp"     // for threadsafe_data
#include "timer.hpp"         // for monotonic_timer
#include <cstdint>           // for uint32_t, uint64_t
#include <memory>            // for unique_ptr
#include <string>            // for basic_string, string
#include <string_view>       // for string_view
#include <utility>           // for pair
#include <vector>            // for vector

namespace adapterremoval {

namespace argparse {

class parser;
enum class parse_result;

} // namespace argparse

class adapter_database;
class alignment_info;
class output_files;
class sample_output_files;
class sample_set;
enum class progress_type;
struct output_file;

enum class ar_command
{
  trim_adapters,
  demultiplex_only,
  report_only,
  benchmark,
};

/**
 * Configuration store, containing all user-supplied options / default values,
 * as well as help-functions using these options.
 */
class userconfig
{
public:
  userconfig();
  ~userconfig();

  /** Parses a set of command-line arguments. */
  argparse::parse_result parse_args(const string_vec& argvec);

  output_files get_output_filenames() const;

  /** Returns runtime in seconds. */
  double runtime() const { return m_runtime.duration(); }

  //! Command-line arguments
  string_vec args{};

  //! Type of run to execute; see command
  ar_command run_type = ar_command::trim_adapters;

  //! Path to input file containing mate 1 reads (required)
  string_vec input_files_1{};
  //! Path to input file containing mate 2 reads (for PE reads)
  string_vec input_files_2{};
  //! Prefix used for output files for which no filename was explicitly set
  std::string out_prefix{};

  //! Name prefix for mate 1 reads
  std::string prefix_read_1{};
  //! Name prefix for mate 2 reads
  std::string prefix_read_2{};
  //! Name prefix for merged reads
  std::string prefix_merged{};

  //! Set to true if both --in-file1 and --in-file2 are set, or if either of
  //! --interleaved or --interleaved-input are set.
  bool paired_ended_mode = false;
  //! Set to true if --interleaved or --interleaved-input is set.
  bool interleaved_input = false;
  //! Set to true if --interleaved or --interleaved-output is set.
  bool interleaved_output = false;

  //! Maximum of reads/read pairs to process
  uint64_t head = 0;

  //! Character separating the mate number from the read name in FASTQ reads.
  char input_mate_separator = '\0';

  [[nodiscard]] char get_output_mate_separator(char value) const
  {
    return m_normalize_mate_separator ? m_output_mate_separator : value;
  }

  //! The minimum length of trimmed reads (ie. genomic nts) to be retained
  uint32_t min_genomic_length{};
  //! The maximum length of trimmed reads (ie. genomic nts) to be retained
  uint32_t max_genomic_length{};
  //! The minimum required overlap before trimming single-end reads.
  uint32_t min_adapter_overlap{};
  //! Rate of mismatches determining the threshold for a an acceptable
  //! alignment, depending on the length of the alignment. But see also the
  //! limits set in the function 'is_good_alignment'.
  double mismatch_threshold{};

  //! Quality format expected in input files.
  fastq_encoding io_encoding = FASTQ_ENCODING_33;

  //! Fixed number of bases to trim from 5' for mate 1 and mate 2 reads
#ifdef PRE_TRIM_5P
  std::pair<unsigned, unsigned> pre_trim_fixed_5p{};
#endif
  std::pair<unsigned, unsigned> pre_trim_fixed_3p{};
  std::pair<unsigned, unsigned> post_trim_fixed_5p{};
  std::pair<unsigned, unsigned> post_trim_fixed_3p{};

  //! Strategy used for performing quality trimming
  trimming_strategy trim = trimming_strategy::none;
  //! [mott] Error rate used for trimming using the modified Mott's algorithm.
  double trim_mott_rate{};
  //! [window] Window based trimming; a fraction / N bp size / off (negative)
  double trim_window_length{};
  //! [window/per-base] The highest quality score considered low-quality
  uint32_t trim_quality_score{};
  //! [per-base] If true, low quality bases read termini are trimmed.
  bool trim_low_quality_bases = false;
  //! [per-base] If true, ambiguous bases (N) at read termini are trimmed.
  bool trim_ambiguous_bases = false;

  //! Nucleotides to trim from poly-X tails prior to alignment/adapter trimming.
  threadsafe_data<std::string> pre_trim_poly_x{};
  //! Nucleotides to trim from poly-X tails after alignment/adapter trimming.
  threadsafe_data<std::string> post_trim_poly_x{};
  //! Minimum number of bases in poly-X tails.
  uint32_t trim_poly_x_threshold{};

  //! The maximum number of ambiguous bases (N) in an read; reads exceeding
  //! this number following trimming (optionally) are discarded.
  uint32_t max_ambiguous_bases{};
  //! The maximum fraction of ambiguous bases (N) in an read.
  double max_ambiguous_base_fraction{};

  //! The minimum average phred score of non-empty reads
  double min_mean_quality = 0.0;
  //! The minimum complexity score for FASTQ reads (see FASTQ::complexity()).
  double min_complexity{};

  //! If true, only the 3p is trimmed for low quality bases (if enabled)
  bool preserve5p = false;

  //! If true, PE reads overlapping at least 'merge_threshold' are
  //! merged to generate a higher quality consensus sequence.
  merge_strategy merge = merge_strategy::none;
  //! The minimum required genomic overlap before merging reads into one.
  uint32_t merge_threshold{};
  //! The maximum quality allowed when recalculating quality scores
  uint32_t merge_quality_max{};

  // Allow for slipping base-pairs by allowing missing bases in adapter
  uint32_t shift{};

  //! The maximum number of threads used by the program
  uint32_t max_threads{};

  //! Automatically determine most suitable SIMD instruction set
  bool simd_auto_select = false;
  //! SIMD instruction set used for alignments
  threadsafe_data<simd::instruction_set> simd{ simd::instruction_set::none };

  //! The format in which in which output reads are written to files
  output_format out_file_format = output_format::fastq_gzip;
  //! The format in which in which output reads are written to STDOUT
  output_format out_stdout_format = output_format::fastq;
  //! Compression level used for output reads where appropriate
  uint32_t compression_level{};

  //! Maximum number of mismatches (considering both barcodes for PE)
  uint32_t barcode_mm{};
  //! Maximum number of mismatches (considering both barcodes for PE)
  uint32_t barcode_mm_r1{};
  //! Maximum number of mismatches (considering both barcodes for PE)
  uint32_t barcode_mm_r2{};
  //! Normalize the orientation of merged reads to forward
  bool normalize_orientation{};

  //! Strategy for selecting what adapters to trim
  adapter_selection adapter_selection_strategy = adapter_selection::manual;
  //! Fallback strategy if the adapters could not be detected automatically
  adapter_fallback adapter_fallback_strategy = adapter_fallback::abort;
  //! Sample specific barcodes and adapters. In non-demultiplexing mode this
  //! set contains a single unnamed sample with empty barcodes
  threadsafe_data<sample_set> samples{};

  /** For auto-adapter selection; returns known and user adapters sequences */
  [[nodiscard]] adapter_database known_adapters() const;

  //! Title used for HTML report
  std::string report_title{};
  //! Fraction of reads used for quality/content curves, etc.
  double report_sample_rate{};
  //! Number of reads used to estimate duplication in input files
  uint64_t report_duplication{};

  //! The kind of progress indicator to use
  progress_type log_progress{};

  //! Sink for benchmark
  string_vec benchmarks{};

  /* Helper functions for logging / reporting */
  bool is_adapter_trimming_enabled() const;
  bool is_demultiplexing_enabled() const;
  bool is_read_merging_enabled() const;

  bool is_any_quality_trimming_enabled() const;
  bool is_low_quality_trimming_enabled() const;
  bool is_terminal_base_pre_trimming_enabled() const;
  bool is_terminal_base_post_trimming_enabled() const;
  bool is_poly_x_tail_pre_trimming_enabled() const;
  bool is_poly_x_tail_post_trimming_enabled() const;

  bool is_any_filtering_enabled() const;
  bool is_short_read_filtering_enabled() const;
  bool is_long_read_filtering_enabled() const;
  bool is_ambiguous_base_filtering_enabled() const;
  bool is_mean_quality_filtering_enabled() const;
  bool is_low_complexity_filtering_enabled() const;

  userconfig(const userconfig&) = delete;
  userconfig(userconfig&&) = delete;
  userconfig& operator=(const userconfig&) = delete;
  userconfig& operator=(userconfig&&) = delete;

private:
  /** Sets up adapter sequences based on user settings. */
  bool setup_adapter_sequences();
  /** Sets up demultiplexing settings: Mismatch, samples, and barcodes */
  bool setup_demultiplexing();

  /** Returns the file format in which to write the give file */
  [[nodiscard]] output_format infer_output_format(
    std::string_view filename) const;

  /**
   * Generates filename + format for a output (FQ/SAM/BAM/HTML/JSON) file.
   * sample_name may optionally contain a sample name, otherwise the file is
   * expected to be used by all samples (if any); keys contain filename parts
   * such as "r1", from which the final filename is constructed. The last key is
   * ignored by formats that aggregate output by default; ext is the mandatory
   * file extension.
   */
  [[nodiscard]] output_file new_output_file(std::string_view key,
                                            std::string_view sample,
                                            std::vector<std::string_view> keys,
                                            std::string_view ext) const;

  //! Argument parser setup to parse the arguments expected by AR
  std::unique_ptr<argparse::parser> m_argparser{};

  //! Sink for --adapter1, adapter sequence expected at 3' of mate 1 reads
  std::string adapter_1{};
  //! Sink for --adapter2, adapter sequence expected at 3' of mate 2 reads
  std::string adapter_2{};
  //! Sink for --adapter-table; list of adapter #1 and #2 sequences
  std::string adapter_table{};

  //! Sink for --barcode-table; list of barcode #1 (and #2 sequences)
  std::string barcode_table{};

  //! Sink for user-supplied quality score formats; use quality_input_fmt.
  std::string quality_input_base{};
  //! Sink for --interleaved
  bool interleaved = false;

  //! Sink for --read-group
  string_vec read_group{};

  //! Sinks for --pre-trim5p/--pre-trimp3p
#ifdef PRE_TRIM_5P
  string_vec pre_trim5p{};
#endif
  string_vec pre_trim3p{};
  //! Sinks for --post-trim5p/--post-trim3p
  string_vec post_trim5p{};
  string_vec post_trim3p{};
  //! Sink for --pre-trim-polyx
  string_vec pre_trim_poly_x_sink{};
  //! Sink for --post-trim-polyx
  string_vec post_trim_poly_x_sink{};

  //! Sink for log color (on/off/auto)
  std::string log_color{};
  //! Sink for log levels
  std::string log_level{};

  //! Measures runtime since the program was started
  monotonic_timer m_runtime{};

  //! Should the output mate separator be normalized
  bool m_normalize_mate_separator = false;
  //! Mate separator to use when writing FASTQ files
  char m_output_mate_separator = '\0';
};

} // namespace adapterremoval
