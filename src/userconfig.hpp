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

#include <array>    // for array
#include <stddef.h> // for size_t
#include <string>   // for string
#include <utility>  // for pair
#include <vector>   // for vector

#include "adapterset.hpp"  // for adapter_set
#include "argparse.hpp"    // for parse_result, parser
#include "commontypes.hpp" // for string_vec, read_type, read_type::max
#include "fastq_enc.hpp"   // for fastq_encoding
#include "timer.hpp"       // for monotonic_timer

namespace adapterremoval {

struct alignment_info;
class trimming_statistics;
enum class progress_type;

//! Path used to indicate that a file is not needed
const std::string DEV_NULL = "/dev/null";

enum class ar_command
{
  trim_adapters,
  identify_adapters,
  demultiplex_sequences,
  report_only,
};

/** Per sample output filenames / steps  */
class output_sample_files
{
public:
  output_sample_files();

  /** Sets the output filename for a given read type. */
  void set_filename(read_type rtype, const std::string& filename);
  /** Pushes a pipeline step for a given output filename. There can  */
  void push_pipeline_step(size_t step);

  /** Unique output filenames, indexed using `offset` */
  const string_vec& filenames() const;
  /** Unique pipeline steps, indexed using `offset` */
  const std::vector<size_t>& pipeline_steps() const;

  /** Returns the offset to the pipeline step/filename for a given read type. */
  size_t offset(read_type value) const;

  //! Constant used to represent disabled output files/steps.
  static const size_t disabled;

private:
  //! Unique output filenames. Multiple read types may be mapped to a filename
  string_vec m_filenames;
  //! Unique pipeline steps IDs. Multiple read types may be mapped to a step
  std::vector<size_t> m_pipeline_steps;
  //! Mapping of read types to filenames/steps.
  std::array<size_t, static_cast<size_t>(read_type::max)> m_offsets;
};

/** Class used to organize filenames of output files. */
class output_files
{
public:
  output_files();

  //! JSON file containing settings / statistics
  std::string settings_json;
  //! HTML file containing settings / statistics / plots
  std::string settings_html;

  //! Filename for unidentified mate 1 reads (demultiplexing)
  std::string unidentified_1;
  //! Filename for unidentified mate 1 reads (demultiplexing)
  std::string unidentified_2;

  std::vector<output_sample_files> samples;
};

/**
 * Configuration store, containing all user-supplied options / default values,
 * as well as help-functions using these options.
 */
class userconfig
{
public:
  /**
   * @param name Name of program.
   * @param version Version string excluding program name.
   * @param help Help text describing program.
   */
  userconfig(const std::string& name,
             const std::string& version,
             const std::string& help);

  /** Parses a set of commandline arguments. */
  argparse::parse_result parse_args(int argc, char* argv[]);

  output_files get_output_filenames() const;

  /** Characterize an alignment based on user settings. */
  bool is_good_alignment(const alignment_info& alignment) const;

  /** Returns true if the alignment is sufficient for merge. */
  bool can_merge_alignment(const alignment_info& alignment) const;

  /** Returns runtime in seconds. */
  double runtime() const;

  //! Command-line arguments
  string_vec args;

  //! Type of run to execute; see command
  ar_command run_type;

  //! Path to input file containing mate 1 reads (required)
  string_vec input_files_1;
  //! Path to input file containing mate 2 reads (for PE reads)
  string_vec input_files_2;
  //! Prefix used for output files for which no filename was explicitly set
  std::string out_basename;
  //! Template filename used for writing JSON report
  std::string out_json;
  //! Template filename used for writing HTML report
  std::string out_html;
  //! Template filename used for writing mate 1 reads
  std::string out_pair_1;
  //! Template filename used for writing mate 2 reads
  std::string out_pair_2;
  //! Template filename used for writing merged mate 1/2 reads
  std::string out_merged;
  //! Template filename used for writing discarded reads
  std::string out_discarded;
  //! Template filename used for writing singleton reads
  std::string out_singleton;

  //! Name prefix for mate 1 reads
  std::string prefix_read_1;
  //! Name prefix for mate 2 reads
  std::string prefix_read_2;
  //! Name prefix for merged reads
  std::string prefix_merged;

  //! Set to true if both --input1 and --input2 are set, or if either of
  //! --interleaved or --interleaved-input are set.
  bool paired_ended_mode;
  //! Set to true if --interleaved or --interleaved-input is set.
  bool interleaved_input;
  //! Set to true if --interleaved or --interleaved-output is set.
  bool interleaved_output;

  //! Maximum of reads/read pairs to process
  uint64_t head;

  //! Character separating the mate number from the read name in FASTQ reads.
  char mate_separator;

  //! The minimum length of trimmed reads (ie. genomic nts) to be retained
  unsigned min_genomic_length;
  //! The maximum length of trimmed reads (ie. genomic nts) to be retained
  unsigned max_genomic_length;
  //! The minimum required overlap before trimming single-end reads.
  unsigned min_adapter_overlap;
  //! Rate of mismatches determining the threshold for a an acceptable
  //! alignment, depending on the length of the alignment. But see also the
  //! limits set in the function 'is_good_alignment'.
  double mismatch_threshold;

  //! Quality format expected in input files.
  fastq_encoding io_encoding;

  //! Fixed number of bases to trim from 5' for mate 1 and mate 2 reads
  std::pair<unsigned, unsigned> pre_trim_fixed_5p;
  std::pair<unsigned, unsigned> pre_trim_fixed_3p;
  std::pair<unsigned, unsigned> post_trim_fixed_5p;
  std::pair<unsigned, unsigned> post_trim_fixed_3p;

  //! Error rate used for quality trimming using the modified Mott's algorithm.
  double trim_error_rate;
  //! Nucleotides to trim from poly-X tails prior to alignment/adapter trimming.
  std::string pre_trim_poly_x;
  //! Nucleotides to trim from poly-X tails after alignment/adapter trimming.
  std::string post_trim_poly_x;
  //! Minimum number of bases in poly-X tails.
  unsigned trim_poly_x_threshold;

  //! DEPRECATED: If true, read termini are trimmed for low-quality bases.
  bool trim_by_quality;
  //! DEPRECATED: Window based trimming; a fraction / N bp size / off (negative)
  double trim_window_length;
  //! DEPRECATED: The highest quality score which is considered low-quality
  unsigned low_quality_score;
  //! DEPRECATED: If true, ambiguous bases (N) at read termini are trimmed.
  bool trim_ambiguous_bases;

  //! The maximum number of ambiguous bases (N) in an read; reads exceeding
  //! this number following trimming (optionally) are discarded.
  unsigned max_ambiguous_bases;
  //! The minimum complexity score for FASTQ reads (see FASTQ::complexity()).
  double min_complexity;

  //! If true, only the 3p is trimmed for low quality bases (if enabled)
  bool preserve5p;

  //! If true, PE reads overlapping at least 'merge_threshold' are
  //! merged to generate a higher quality consensus sequence.
  bool merge;
  //! The minimum required genomic overlap before merging reads into one.
  unsigned merge_threshold;

  // Allow for slipping basepairs by allowing missing bases in adapter
  unsigned shift;

  //! The maximum number of threads used by the program
  unsigned max_threads;

  //! GZip compression enabled / disabled
  bool gzip;
  //! GZip compression level used for output reads
  unsigned int gzip_level;

  //! Maximum number of mismatches (considering both barcodes for PE)
  unsigned barcode_mm;
  //! Maximum number of mismatches (considering both barcodes for PE)
  unsigned barcode_mm_r1;
  //! Maximum number of mismatches (considering both barcodes for PE)
  unsigned barcode_mm_r2;

  adapter_set adapters;

  //! Fraction of reads used for quality/content curves, etc.
  double report_sample_rate;
  //! Number of reads used to estimate duplication in input files
  unsigned report_duplication;

  //! The kind of progress indicator to use
  progress_type log_progress;

  /* Helper functions for logging / reporting */
  bool is_adapter_trimming_enabled() const;

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
  bool is_low_complexity_filtering_enabled() const;

  //! Copy construction not supported
  userconfig(const userconfig&) = delete;
  //! Assignment not supported
  userconfig& operator=(const userconfig&) = delete;

private:
  /** Sets up adapter sequences based on user settings. */
  bool setup_adapter_sequences();

  std::string new_filename(const std::string& key,
                           const std::string& first,
                           const std::string& second = std::string()) const;

  //! Argument parser setup to parse the arguments expected by AR
  argparse::parser argparser;

  //! Sink for --adapter1, adapter sequence expected at 3' of mate 1 reads
  std::string adapter_1;
  //! Sink for --adapter2, adapter sequence expected at 3' of mate 2 reads
  std::string adapter_2;
  //! Sink for --adapter-list; list of adapter #1 and #2 sequences
  std::string adapter_list;

  //! Sink for --barcode-list; list of barcode #1 (and #2 sequences)
  std::string barcode_list;

  //! Sink for user-supplied quality score formats; use quality_input_fmt.
  std::string quality_input_base;
  //! Sink for the mate separator character; use mate separator
  std::string mate_separator_str;
  //! Sink for --interleaved
  bool interleaved;

  //! Sinks for --pre-trim5p/--pre-trimp3p
  string_vec pre_trim5p;
  string_vec pre_trim3p;
  //! Sinks for --post-trim5p/--post-trimp3p
  string_vec post_trim5p;
  string_vec post_trim3p;
  //! Sink for --pre-trim-polyx
  string_vec pre_trim_poly_x_sink;
  //! Sink for --post-trim-polyx
  string_vec post_trim_poly_x_sink;

  //! Sink for log color (on/off/auto)
  std::string log_color;
  //! Sink for log levels
  std::string log_level;
  //! Sink for progress indicators
  std::string log_progress_sink;

  //! Sink for head
  std::string head_sink;

  //! Measures runtime since the program was started
  monotonic_timer m_runtime;

  //! Sink for deprecated knobs
  unsigned m_deprecated_knobs;
};

} // namespace adapterremoval
