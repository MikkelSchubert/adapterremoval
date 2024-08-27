/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#pragma once

#include "argparse.hpp"      // for parse_result, parser
#include "commontypes.hpp"   // for string_vec, read_type, merge_strategy
#include "fastq_enc.hpp"     // for fastq_encoding
#include "sequence_sets.hpp" // for adapter_set
#include "serializer.hpp"    // for read_group
#include "simd.hpp"          // for size_t, instruction_set
#include "timer.hpp"         // for monotonic_timer
#include <cstddef>           // for size_t
#include <cstdint>           // for uint64_t
#include <string>            // for string
#include <utility>           // for pair

namespace adapterremoval {

enum class progress_type;
class output_files;
class sample_output_files;
struct alignment_info;
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
  ~userconfig() = default;

  /** Parses a set of command-line arguments. */
  argparse::parse_result parse_args(const string_vec& argvec);

  output_files get_output_filenames() const;

  /** Characterize an alignment based on user settings. */
  bool is_good_alignment(const alignment_info& alignment) const;

  /** Returns true if the alignment is sufficient for merge. */
  bool can_merge_alignment(const alignment_info& alignment) const;

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
  char mate_separator{};

  //! The minimum length of trimmed reads (ie. genomic nts) to be retained
  unsigned min_genomic_length{};
  //! The maximum length of trimmed reads (ie. genomic nts) to be retained
  unsigned max_genomic_length{};
  //! The minimum required overlap before trimming single-end reads.
  unsigned min_adapter_overlap{};
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
  unsigned trim_quality_score{};
  //! [per-base] If true, low quality bases read termini are trimmed.
  bool trim_low_quality_bases = false;
  //! [per-base] If true, ambiguous bases (N) at read termini are trimmed.
  bool trim_ambiguous_bases = false;

  //! Nucleotides to trim from poly-X tails prior to alignment/adapter trimming.
  std::string pre_trim_poly_x{};
  //! Nucleotides to trim from poly-X tails after alignment/adapter trimming.
  std::string post_trim_poly_x{};
  //! Minimum number of bases in poly-X tails.
  unsigned trim_poly_x_threshold{};

  //! The maximum number of ambiguous bases (N) in an read; reads exceeding
  //! this number following trimming (optionally) are discarded.
  unsigned max_ambiguous_bases{};
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
  unsigned merge_threshold{};
  //! The maximum quality allowed when recalculating quality scores
  unsigned merge_quality_max{};

  // Allow for slipping base-pairs by allowing missing bases in adapter
  unsigned shift{};

  //! The maximum number of threads used by the program
  unsigned max_threads{};
  //! SIMD instruction set used for alignments
  simd::instruction_set simd = simd::instruction_set::none;

  //! The format in which in which output reads are written to files
  output_format out_file_format = output_format::fastq_gzip;
  //! The format in which in which output reads are written to STDOUT
  output_format out_stdout_format = output_format::fastq;
  //! Compression level used for output reads where appropriate
  unsigned int compression_level{};
  //! Read group for SAM/BAM output
  read_group output_read_group{};

  //! Maximum number of mismatches (considering both barcodes for PE)
  unsigned barcode_mm{};
  //! Maximum number of mismatches (considering both barcodes for PE)
  unsigned barcode_mm_r1{};
  //! Maximum number of mismatches (considering both barcodes for PE)
  unsigned barcode_mm_r2{};

  //! Adapter sequences expected to be found in the input (without barcodes)
  adapter_set adapters{};
  //! Sample names and barcodes for demultiplexing
  barcode_set samples{};

  //! Fraction of reads used for quality/content curves, etc.
  double report_sample_rate{};
  //! Number of reads used to estimate duplication in input files
  unsigned report_duplication{};

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
    const std::string& filename) const;

  [[nodiscard]] output_file new_output_file(const std::string& key,
                                            const string_vec& values) const;

  //! Argument parser setup to parse the arguments expected by AR
  argparse::parser argparser{};

  //! Sink for --adapter1, adapter sequence expected at 3' of mate 1 reads
  std::string adapter_1{};
  //! Sink for --adapter2, adapter sequence expected at 3' of mate 2 reads
  std::string adapter_2{};
  //! Sink for --adapter-list; list of adapter #1 and #2 sequences
  std::string adapter_list{};

  //! Sink for --barcode-list; list of barcode #1 (and #2 sequences)
  std::string barcode_list{};

  //! Sink for user-supplied quality score formats; use quality_input_fmt.
  std::string quality_input_base{};
  //! Sink for the mate separator character; use mate separator
  std::string mate_separator_str{};
  //! Sink for --interleaved
  bool interleaved = false;

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
};

} // namespace adapterremoval
