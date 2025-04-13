// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "userconfig.hpp"  // declarations
#include "alignment.hpp"   // for alignment_info
#include "commontypes.hpp" // for string_vec, DEV_STDOUT, DEV_STDERR, ...
#include "debug.hpp"       // for AR_REQUIRE, AR_FAIL
#include "errors.hpp"      // for fastq_error
#include "fastq.hpp"       // for ACGT, ACGT::indices, ACGT::values
#include "fastq_enc.hpp"   // for PHRED_SCORE_MAX
#include "licenses.hpp"    // for LICENSES
#include "logging.hpp"     // for log_stream, error, set_level, set_colors, info
#include "main.hpp"        // for HELPTEXT, NAME, VERSION
#include "output.hpp"      // for DEV_NULL, output_files, output_file
#include "progress.hpp"    // for progress_type, progress_type::simple, progr...
#include "sequence.hpp"    // for dna_sequence
#include "simd.hpp"        // for size_t, name, supported, instruction_set
#include "strutils.hpp"    // for shell_escape, str_to_u32
#include <algorithm>       // for find, max, min
#include <cerrno>          // for errno
#include <cmath>           // for pow
#include <cstdlib>         // for getenv
#include <cstring>         // for size_t, strerror, strcmp
#include <filesystem>      // for weakly_canonical
#include <limits>          // for numeric_limits
#include <stdexcept>       // for invalid_argument
#include <string>          // for string, basic_string, operator==, operator+
#include <string_view>     // for string_view
#include <tuple>           // for get, tuple
#include <unistd.h>        // for access, isatty, R_OK, STDERR_FILENO

namespace adapterremoval {

namespace {

const char* HELPTEXT =
  "This program searches for and removes remnant adapter sequences, poly-X "
  "tails and low-quality base from FASTQ reads. For detailed explanation of "
  "the parameters, please refer to the man page. For comments, suggestions "
  "and feedback please use\n"
  "\n"
  "  https://github.com/MikkelSchubert/adapterremoval/issues/new\n"
  "\n"
  "If you use the program, please cite the paper\n"
  "\n"
  "  Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid\n"
  "  adapter trimming, identification, and read merging. BMC Research\n"
  "  Notes, 12;9(1):88. https://doi.org/10.1186/s13104-016-1900-2\n"
  "\n"
  "Use the filename '-' to read from STDIN or to write to STDOUT. If the same "
  "filenames are used for two or more of the --out-* options (excluding "
  "--out-json and --out-html), then the combined output is written to that "
  "file in interleaved mode.\n";

////////////////////////////////////////////////////////////////////////////////
// Helper functions

std::pair<unsigned, unsigned>
parse_trim_argument(const string_vec& values)
{
  unsigned mate_1 = 0;
  unsigned mate_2 = 0;

  switch (values.size()) {
    case 1:
      mate_1 = str_to_u32(values.front());
      mate_2 = mate_1;
      break;

    case 2:
      mate_1 = str_to_u32(values.front());
      mate_2 = str_to_u32(values.back());
      break;

    default:
      throw std::invalid_argument("please specify exactly one or two values");
  }

  return { mate_1, mate_2 };
}

bool
parse_poly_x_option(const std::string& key,
                    const string_vec& values,
                    std::string& out)
{
  out.clear();
  if (values.empty()) {
    out = "ACGT";
    return true;
  }

  std::array<bool, ACGT::indices> enabled = {};
  for (const auto& value : values) {
    for (const auto nuc : to_upper(value)) {
      switch (nuc) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
          enabled.at(ACGT::to_index(nuc)) = true;
          break;

        default:
          log::error() << "Option " << key << " called with invalid value "
                       << shell_escape(value) << ". Only A, C, G, and T are "
                       << "permitted!";

          return false;
      }
    }
  }

  for (const auto nuc : ACGT::values) {
    if (enabled.at(ACGT::to_index(nuc))) {
      out.push_back(nuc);
    }
  }

  return true;
}

bool
parse_counts(const argparse::parser& args,
             const std::string& key,
             uint64_t& out)
{
  const auto sink = args.value(std::string{ key });
  if (sink.empty()) {
    return true;
  }

  uint64_t unit = 1;
  std::string sink_without_unit = sink;
  if (sink.back() < '0' || sink.back() > '9') {
    switch (sink.back()) {
      case 'k':
      case 'K':
        unit = 1000;
        break;

      case 'm':
      case 'M':
        unit = 1000'000;
        break;

      case 'g':
      case 'G':
        unit = 1000'000'000;
        break;

      default:
        log::error() << "Invalid unit in command-line option " << key
                     << shell_escape(sink);
        return false;
    }

    sink_without_unit.pop_back();
  }

  try {
    // This should not be able to overflow as log2(2^32 * 1e9) ~= 62,
    // but will need to be changed if we want to allow large raw numbers
    out = static_cast<uint64_t>(str_to_u32(sink_without_unit)) * unit;
  } catch (const std::invalid_argument&) {
    log::error() << "Invalid value in command-line option --sink "
                 << shell_escape(sink);
    return false;
  }

  return true;
}

bool
check_no_clobber(const std::string& label,
                 const string_vec& in_files,
                 const output_file& out_file)
{
  for (const auto& in_file : in_files) {
    if (in_file == out_file.name && in_file != DEV_NULL) {
      log::error() << "Input file would be overwritten: " << label << " "
                   << in_file;
      return false;
    }
  }

  return true;
}

/** Replace the STDIN pseudo-filename with the device path */
void
normalize_input_file(std::string& filename)
{
  if (filename == DEV_PIPE) {
    filename = DEV_STDIN;
  }
}

/** Replace the STDIN pseudo-filename with the device path */
void
normalize_output_file(std::string& filename)
{
  if (filename == DEV_PIPE) {
    filename = DEV_STDOUT;
  }
}

void
append_normalized_input_files(string_pair_vec& out, const string_vec& filenames)
{
  for (const auto& filename : filenames) {
    try {
      out.emplace_back(std::filesystem::weakly_canonical(filename), filename);
    } catch (const std::filesystem::filesystem_error&) {
      // Permission errors are handled by the explicit access checks below
      out.emplace_back(filename, filename);
    }
  }
}

bool
check_input_files(const string_vec& filenames_1, const string_vec& filenames_2)
{
  string_pair_vec filenames;
  append_normalized_input_files(filenames, filenames_1);
  append_normalized_input_files(filenames, filenames_2);
  std::sort(filenames.begin(), filenames.end());

  bool any_errors = false;
  for (size_t i = 1; i < filenames.size(); ++i) {
    const auto& it_0 = filenames.at(i - 1);
    const auto& it_1 = filenames.at(i);

    if (it_0.second == it_1.second) {
      log::error() << "Input file " << log_escape(it_0.second)
                   << " has been specified multiple times using --in-file1 "
                      "and/or --in-file2";
      any_errors = true;
    } else if (it_0.first == it_1.first) {
      log::error() << "The path of input file " << log_escape(it_0.second)
                   << " and the path of input " << "file "
                   << log_escape(it_1.second) << " both point to the file "
                   << log_escape(it_0.first);
      any_errors = true;
    }
  }

  for (const auto& it : filenames) {
    if (access(it.second.c_str(), R_OK)) {
      log::error() << "Cannot access input file " << log_escape(it.second)
                   << ": " << std::strerror(errno);
      any_errors = true;
    }
  }

  return !any_errors;
}

bool
check_output_files(const std::string& label,
                   const string_vec& filenames,
                   const output_files& output_files)
{

  if (!check_no_clobber(label, filenames, output_files.unidentified_1)) {
    return false;
  }

  if (!check_no_clobber(label, filenames, output_files.unidentified_2)) {
    return false;
  }

  for (const auto& sample : output_files.samples()) {
    for (size_t i = 0; i < sample.size(); ++i) {
      if (!check_no_clobber(label, filenames, sample.file(i))) {
        return false;
      }
    }
  }

  return true;
}

/**
 * Tries to parse a simple command-line argument while ignoring the validity
 * of the overall command-line. This is only intended to make pre-configured
 * logging output consistent with post-configured output if possible.
 */
std::string
try_parse_argument(const string_vec& args,
                   const std::string& key,
                   const std::string& fallback)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it != args.end() && (it + 1) != args.end()) {
    return *(it + 1);
  }

  return fallback;
}

/** Returns vector of keys for output files that have been set by the user. */
string_vec
user_supplied_keys(const argparse::parser& argparser, const string_vec& keys)
{
  string_vec result;
  for (const auto& key : keys) {
    if (argparser.is_set(key)) {
      result.push_back(key);
    }
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////

bool
fancy_output_allowed()
{
  if (::isatty(STDERR_FILENO)) {
    // NO_COLOR is checked as suggested by https://no-color.org/
    const char* no_color = std::getenv("NO_COLOR");
    const char* term = std::getenv("TERM");

    return !(no_color && no_color[0] != '\0') &&
           !(term && strcmp(term, "dumb") == 0);
  }

  return false;
}

void
configure_log_levels(const std::string& value, bool fallible = false)
{
  const auto log_level = to_lower(value);

  if (log_level == "debug") {
    log::set_level(log::level::debug);
  } else if (log_level == "info") {
    log::set_level(log::level::info);
  } else if (log_level == "warning") {
    log::set_level(log::level::warning);
  } else if (log_level == "error") {
    log::set_level(log::level::error);
  } else {
    AR_REQUIRE(fallible, "unhandled log_level value");
  }
}

void
configure_log_colors(const std::string& colors, bool fallible = false)
{
  if (colors == "always") {
    log::set_colors(true);
  } else if (colors == "never") {
    log::set_colors(false);
  } else if (colors == "auto") {
    log::set_colors(fancy_output_allowed());
  } else {
    AR_REQUIRE(fallible, "unhandled log_colors value");
  }
}

progress_type
configure_log_progress(const std::string& progress)
{
  if (progress == "never") {
    return progress_type::none;
  } else if (progress == "spin") {
    return progress_type::spinner;
  } else if (progress == "log") {
    return progress_type::simple;
  } else if (progress == "auto") {
    if (fancy_output_allowed()) {
      return progress_type::spinner;
    } else {
      return progress_type::simple;
    }
  }

  AR_FAIL("unhandled log_progress value");
}

fastq_encoding
configure_encoding(const std::string& value,
                   degenerate_encoding degenerate,
                   uracil_encoding uracils)
{
  if (value == "33") {
    return fastq_encoding{ quality_encoding::phred_33, degenerate, uracils };
  } else if (value == "64") {
    return fastq_encoding{ quality_encoding::phred_64, degenerate, uracils };
  } else if (value == "solexa") {
    return fastq_encoding{ quality_encoding::solexa, degenerate, uracils };
  } else if (value == "sam") {
    return fastq_encoding{ quality_encoding::sam, degenerate, uracils };
  }

  AR_FAIL("unhandled qualitybase value");
}

bool
parse_output_formats(const argparse::parser& argparser,
                     output_format& file_format,
                     output_format& stdout_format)
{
  if (argparser.is_set("--gzip")) {
    file_format = stdout_format = output_format::fastq_gzip;
    return true;
  }

  auto format_s = argparser.value("--out-format");
  if (!output_files::parse_format(format_s, file_format)) {
    log::error() << "Invalid output format " + log_escape(format_s);
    return false;
  }

  // Default to writing uncompressed output to STDOUT
  if (!argparser.is_set("--stdout-format")) {
    switch (file_format) {
      case output_format::fastq:
      case output_format::fastq_gzip:
        stdout_format = output_format::fastq;
        return true;
      case output_format::sam:
      case output_format::sam_gzip:
        stdout_format = output_format::sam;
        return true;
      case output_format::bam:
      case output_format::ubam:
        stdout_format = output_format::ubam;
        return true;
      default:
        AR_FAIL("invalid output format");
    }
  }

  format_s = argparser.value("--stdout-format");
  if (!output_files::parse_format(format_s, stdout_format)) {
    log::error() << "Invalid output format " + log_escape(format_s);
    return false;
  }

  return true;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// Implementations for `userconfig`

std::string userconfig::start_time = timestamp("%FT%T%z");

userconfig::userconfig()
{
  argparser.set_name(NAME);
  argparser.set_version(VERSION);
  argparser.set_preamble(HELPTEXT);
  argparser.set_licenses(LICENSES);
  argparser.set_terminal_width(log::get_terminal_width());

  //////////////////////////////////////////////////////////////////////////////
  argparser.add("--threads", "N")
    .help("Maximum number of threads")
    .bind_u32(&max_threads)
    .with_default(2)
    .with_minimum(1);

  {
    std::vector<std::string> choices;
    for (const auto is : simd::supported()) {
      choices.emplace_back(simd::name(is));
    }

    AR_REQUIRE(!choices.empty());
    argparser.add("--simd", "NAME")
      .help("SIMD instruction set to use; defaults to the most advanced "
            "instruction set supported by this computer")
      .bind_str(nullptr)
      .with_choices(choices)
      .with_default(choices.back());
  }

  argparser.add("--benchmark")
    .help("Carry out benchmarking of AdapterRemoval sub-systems")
    .conflicts_with("--demultiplex-only")
    .conflicts_with("--report-only")
    .conflicts_with("--interleaved")
    .conflicts_with("--interleaved-input")
#if !defined(DEBUG)
    .hidden()
#endif
    .bind_vec(&benchmarks)
    .with_min_values(0);

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("INPUT FILES:");

  argparser.add("--in-file1", "FILE")
    .help("One or more input files containing mate 1 reads [REQUIRED]")
    .deprecated_alias("--file1")
    .bind_vec(&input_files_1)
    .with_preprocessor(normalize_input_file);
  argparser.add("--in-file2", "FILE")
    .help("Input files containing mate 2 reads; if used, then the same number "
          "of files as --in-file1 must be listed [OPTIONAL]")
    .deprecated_alias("--file2")
    .bind_vec(&input_files_2)
    .with_preprocessor(normalize_input_file);
  argparser.add("--head", "N")
    .help("Process only the first N reads in single-end mode or the first N "
          "read-pairs in paired-end mode. Accepts suffixes K (thousands), M "
          "(millions), and G (billions) [default: all reads]")
    .bind_str(nullptr);

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("OUTPUT FILES:");

  argparser.add("--out-prefix", "PREFIX")
    .help("Prefix for output files for which the corresponding --out option "
          "was not set [default: not set]")
    .deprecated_alias("--basename")
    .bind_str(&out_prefix)
    .with_default(DEV_NULL);

  argparser.add_separator();
  argparser.add("--out-file1", "FILE")
    .help("Output file containing trimmed mate 1 reads. Setting this value in "
          "in demultiplexing mode overrides --out-prefix for this file")
    .deprecated_alias("--output1")
    .bind_str(nullptr)
    .with_default("{prefix}[.sample].r1.fastq")
    .with_preprocessor(normalize_output_file);
  argparser.add("--out-file2", "FILE")
    .help("Output file containing trimmed mate 2 reads. Setting this value in "
          "in demultiplexing mode overrides --out-prefix for this file")
    .deprecated_alias("--output2")
    .bind_str(nullptr)
    .with_default("{prefix}[.sample].r2.fastq")
    .with_preprocessor(normalize_output_file);
  argparser.add("--out-merged", "FILE")
    .help("Output file that, if --merge is set, contains overlapping "
          "read-pairs that have been merged into a single read (PE mode only). "
          "Setting this value in demultiplexing mode overrides --out-prefix "
          "for this file")
    .deprecated_alias("--outputcollapsed")
    .bind_str(nullptr)
    .with_default("{prefix}[.sample].merged.fastq")
    .with_preprocessor(normalize_output_file);
  argparser.add("--out-singleton", "FILE")
    .help("Output file containing paired reads for which the mate "
          "has been discarded. This file is only created if filtering is "
          "enabled. Setting this value in demultiplexing mode overrides "
          "--out-prefix for this file")
    .deprecated_alias("--singleton")
    .bind_str(nullptr)
    .with_default("{prefix}[.sample].singleton.fastq")
    .with_preprocessor(normalize_output_file);

  argparser.add_separator();
  argparser.add("--out-unidentified1", "FILE")
    .help("In demultiplexing mode, contains mate 1 reads that could not be "
          "assigned to a single sample")
    .bind_str(nullptr)
    .with_default("{prefix}.unidentified.r1.fastq")
    .with_preprocessor(normalize_output_file);
  argparser.add("--out-unidentified2", "FILE")
    .help("In demultiplexing mode, contains mate 2 reads that could not be "
          "assigned to a single sample")
    .bind_str(nullptr)
    .with_default("{prefix}.unidentified.r2.fastq")
    .with_preprocessor(normalize_output_file);
  argparser.add("--out-discarded", "FILE")
    .help("Output file containing filtered reads. Setting this value in "
          "demultiplexing mode overrides --out-prefix for this file [default: "
          "not saved]")
    .deprecated_alias("--discarded")
    .bind_str(nullptr)
    .with_preprocessor(normalize_output_file);

  argparser.add_separator();
  argparser.add("--out-json", "FILE")
    .help("Output file containing statistics about input files, trimming, "
          "merging, and more in JSON format")
    .bind_str(nullptr)
    .with_default("{prefix}.json")
    .with_preprocessor(normalize_output_file);
  argparser.add("--out-html", "FILE")
    .help("Output file containing statistics about input files, trimming, "
          "merging, and more in HTML format")
    .bind_str(nullptr)
    .with_default("{prefix}.html")
    .with_preprocessor(normalize_output_file);

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("FASTQ OPTIONS:");

  argparser.add("--quality-format", "N")
    .help("Format used to encode Phred scores in input")
    .deprecated_alias("--qualitybase")
    .bind_str(&quality_input_base)
    .with_choices({ "33", "64", "solexa", "sam" })
    .with_default("33");
  argparser.add("--mate-separator", "CHAR")
    .help("Character separating the mate number (1 or 2) from the read name in "
          "FASTQ records. Will be determined automatically if not specified")
    .bind_str(&mate_separator_str);

  argparser.add("--interleaved-input")
    .help("The (single) input file provided contains both the mate 1 and mate "
          "2 reads, one pair after the other, with one mate 1 reads followed "
          "by one mate 2 read. This option is implied by the --interleaved "
          "option")
    .conflicts_with("--in-file2")
    .bind_bool(&interleaved_input);
  argparser.add("--interleaved-output")
    .help("If set, trimmed paired-end reads are written to a single file "
          "containing mate 1 and mate 2 reads, one pair after the other. This "
          "option is implied by the --interleaved option")
    .conflicts_with("--out-file2")
    .bind_bool(&interleaved_output);
  argparser.add("--interleaved")
    .help("This option enables both the --interleaved-input option and the "
          "--interleaved-output option")
    .conflicts_with("--in-file2")
    .conflicts_with("--out-file2")
    .bind_bool(&interleaved);

  argparser.add("--mask-degenerate-bases")
    .help("Mask degenerate/ambiguous bases (B/D/H/K/M/N/R/S/V/W/Y) in the "
          "input by replacing them with an 'N'; if this option is not used, "
          "AdapterRemoval will abort upon encountering degenerate bases");
  argparser.add("--convert-uracils")
    .help("Convert uracils (U) to thymine (T) in input reads; if this option "
          "is not used, AdapterRemoval will abort upon encountering uracils");

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("OUTPUT FORMAT:");

  argparser.add("--gzip")
    .hidden()
    .deprecated()
    .conflicts_with("--out-format")
    .conflicts_with("--stdout-format");
  argparser.add("--out-format", "X")
    .help("Selects the default output format; either 'fastq' for uncompressed "
          "FASTQ reads, 'fastq.gz' for gzip compressed FASTQ reads, 'sam' for "
          "uncompressed SAM records, 'sam.gz' for gzip compressed SAM records, "
          "'bam' for BGZF compressed BAM records, and 'ubam' for uncompressed "
          "BAM records. Setting an `--out-*` option overrides this option "
          "based on the filename used (except .ubam)")
    .bind_str(nullptr)
    .with_choices({ "fastq", "fastq.gz", "sam", "sam.gz", "bam", "ubam" })
    .with_default("fastq.gz");
  argparser.add("--stdout-format", "X")
    .help("Selects the output format for data written to STDOUT; choices are "
          "the same as for --out-format [default: the same format as "
          "--out-format, but uncompressed]")
    .bind_str(nullptr)
    .with_choices({ "fastq", "fastq.gz", "sam", "sam.gz", "bam", "ubam" });
  argparser.add("--read-group", "RG")
    .help("Add read-group to SAM/BAM output. Takes zero or more arguments in "
          "the form 'tag:value' where tag consists of two alphanumerical "
          "characters and where value is one or more characters. An argument "
          "may also contain multiple, tab-separated tag/value pairs. An ID tag "
          "is automatically generated if no ID tag is specified")
    .bind_vec(&read_group)
    .with_min_values(0);
  argparser.add("--compression-level", "N")
    .help(
      "Sets the compression level for compressed output. Valid values are 0 to "
      "13: Level 0 is uncompressed but includes gzip headers/checksums, level "
      "1 is streamed for SAM/FASTQ output (this may be required in rare cases "
      "for compatibility), and levels 2 to 13 are block compressed using the "
      "BGZF format")
    .deprecated_alias("--gzip-level")
    .bind_u32(&compression_level)
    .with_maximum(13)
    .with_default(5);

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("PROCESSING:");

  argparser.add("--adapter1", "SEQ")
    .help("Adapter sequence expected to be found in mate 1 reads. Any 'N' in "
          "this sequence is treated as a wildcard")
    .bind_str(&adapter_1)
    .with_default("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
  argparser.add("--adapter2", "SEQ")
    .help("Adapter sequence expected to be found in mate 2 reads. Any 'N' in "
          "this sequence is treated as a wildcard")
    .bind_str(&adapter_2)
    .with_default("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");
  argparser.add("--adapter-list", "FILE")
    .help("Read adapter pairs from the first two columns of a white-space "
          "separated table. AdapterRemoval will then select the best matching "
          "adapter pair for each pair of input reads when trimming. Only the "
          "first column is required for single-end trimming")
    .conflicts_with("--adapter1")
    .conflicts_with("--adapter2")
    .bind_str(&adapter_list);

  argparser.add_separator();
  argparser.add("--min-adapter-overlap", "N")
    .help("In single-end mode, reads are only trimmed if the overlap between "
          "read and the adapter is at least X bases long, not counting "
          "ambiguous nucleotides (Ns)")
    .deprecated_alias("--minadapteroverlap")
    .bind_u32(&min_adapter_overlap)
    .with_default(1)
    .with_minimum(1);
  argparser.add("--mismatch-rate", "X")
    .help("Max error-rate when aligning reads and/or adapters. If > 1, the max "
          "error-rate is set to 1 / X; if < 0, the defaults are used, "
          "otherwise the user-supplied value is used directly [default: 1/6 "
          "for trimming; 1/10 when identifying adapters]")
    .deprecated_alias("--mm")
    .bind_double(&mismatch_threshold)
    .with_default(-1.0);
  argparser.add("--shift", "N")
    .help("Consider alignments where up to N nucleotides are missing from the "
          "5' termini")
    .bind_u32(&shift)
    .with_default(2);

  argparser.add_separator();
  argparser.add("--merge")
    .help("When set, paired ended read alignments of --merge-threshold or "
          "more bases are merged into a single consensus sequence. Merged "
          "reads are written to prefix.merged by default. Has no effect "
          "in single-end mode")
    .deprecated_alias("--collapse");
  argparser.add("--merge-threshold", "N")
    .help("Paired reads must overlap at least this many bases to be considered "
          "overlapping for the purpose of read merging. Overlapping bases "
          "where one or both bases are ambiguous (N) are not counted")
    .deprecated_alias("--minalignmentlength")
    .bind_u32(&merge_threshold)
    .with_default(11);
  argparser.add("--merge-strategy", "X")
    .help(
      "The 'maximum' strategy uses Q=max(Q1,Q2) for matches while the "
      "'additive' strategy uses Q=Q1+Q2. Both strategies use Q=abs(Q1-Q2) for "
      "mismatches and picks the highest quality base, unless the qualities are "
      "the same in which case 'N' is used. Setting this option implies --merge")
    .bind_str(nullptr)
    .with_choices({ "maximum", "additive" })
    .with_default("maximum");
  argparser.add("--merge-quality-max", "N")
    .help("Sets the maximum Phred score for re-calculated quality scores when "
          "read merging is enabled with the 'additive' merging strategy. The "
          "value must be in the range 0 to 93, corresponding to Phred+33 "
          "encoded values of '!' to '~'")
    .deprecated_alias("--qualitymax")
    .bind_u32(&merge_quality_max)
    .with_maximum(PHRED_SCORE_MAX)
    .with_default(41);
  argparser.add("--collapse-deterministic")
    .conflicts_with("--collapse-conservatively")
    .conflicts_with("--merge-strategy")
    .deprecated();
  argparser.add("--collapse-conservatively")
    .conflicts_with("--collapse-deterministic")
    .conflicts_with("--merge-strategy")
    .deprecated();

  argparser.add_separator();
  argparser.add("--prefix-read1", "X")
    .help("Adds the specified prefix to read 1 names [default: no prefix]")
    .bind_str(&prefix_read_1);
  argparser.add("--prefix-read2", "X")
    .help("Adds the specified prefix to read 2 names [default: no prefix]")
    .bind_str(&prefix_read_2);
  argparser.add("--prefix-merged", "X")
    .help("Adds the specified prefix to merged read names [default: no prefix]")
    .bind_str(&prefix_merged);

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("QUALITY TRIMMING:");

#ifdef PRE_TRIM_5P
  argparser.add("--pre-trim5p", "N")
    .help("Trim the 5' of reads by a fixed amount after demultiplexing (if "
          "enabled) but before trimming adapters and low quality bases. "
          "Specify one value to trim mate 1 and mate 2 reads the same amount, "
          "or two values separated by a space to trim each mate a different "
          "amount [default: no trimming]")
    .bind_vec(&pre_trim5p)
    .with_max_values(2);
#endif
  argparser.add("--pre-trim3p", "N")
    .help("Trim the 3' of reads by a fixed amount after demultiplexing (if "
          "enabled) but before trimming adapters and low quality bases. "
          "Specify one value to trim mate 1 and mate 2 reads the same amount, "
          "or two values separated by a space to trim each mate a different "
          "amount [default: no trimming]")
    .bind_vec(&pre_trim3p)
    .with_max_values(2);

  argparser.add("--post-trim5p", "N")
    .help("Trim the 5' by a fixed amount after removing adapters, but before "
          "carrying out quality based trimming [default: no trimming]")
    .deprecated_alias("--trim5p")
    .bind_vec(&post_trim5p)
    .with_max_values(2);
  argparser.add("--post-trim3p", "N")
    .deprecated_alias("--trim3p")
    .help("Trim the 3' by a fixed amount after removing adapters, but before "
          "carrying out quality based trimming [default: no trimming]")
    .bind_vec(&post_trim3p)
    .with_max_values(2);

  argparser.add_separator();
  argparser.add("--quality-trimming", "method")
    .help("Strategy for trimming low quality bases: 'mott' for the modified "
          "Mott's algorithm; 'window' for window based trimming; 'per-base' "
          "for a per-base trimming of low quality base; and 'none' for no "
          "trimming of low quality bases")
    .deprecated_alias("--trim-strategy") // name used during v3 alpha 1
    .bind_str(nullptr)
    .with_choices({ "mott", "window", "per-base", "none" })
    .with_default("mott");

  argparser.add("--trim-mott-quality", "N")
    .help("The inclusive threshold value used when performing quality based "
          "trimming using the modified Mott's algorithm. The value must be in "
          "the range 0 to 93, corresponding to Phred+33 encoded values of '!' "
          "to '~'")
    .deprecated_alias("--trim-mott-rate")
    .conflicts_with("--trim-windows")
    .conflicts_with("--trim-ns")
    .conflicts_with("--trim-qualities")
    .conflicts_with("--trim-min-quality")
    .bind_double(&trim_mott_rate)
    .with_minimum(0)
    .with_maximum(93)
    .with_default(13);
  argparser.add("--trim-windows", "X")
    .help("Specifies the size of the window used for '--quality-trimming "
          "window': If >= 1, this value will be used as the window size; if "
          "the value is < 1, window size is the read length times this value. "
          "If the resulting window size is 0 or larger than the read length, "
          "the read length is used as the window size")
    .deprecated_alias("--trimwindows")
    .conflicts_with("--trim-mott-quality")
    .conflicts_with("--trim-qualities")
    .bind_double(&trim_window_length)
    .with_minimum(0.0)
    .with_default(0.1);
  argparser.add("--trim-min-quality", "N")
    .help("Inclusive minimum quality used when trimming low-quality bases with "
          "--quality-trimming options 'window' and 'per-base'. The value must "
          "be in the range 0 to 93, corresponding to Phred+33 encoded values "
          "of '!' to '~'")
    .deprecated_alias("--minquality")
    .conflicts_with("--trim-mott-quality")
    .bind_u32(&trim_quality_score)
    .with_maximum(PHRED_SCORE_MAX)
    .with_default(2);
  argparser.add("--trim-ns")
    .help("If set, trim ambiguous bases (N) at 5'/3' termini when using the "
          "'window' or the 'per-base' trimming strategy")
    .conflicts_with("--trim-mott-quality")
    .deprecated_alias("--trimns")
    .bind_bool(&trim_ambiguous_bases);
  argparser.add("--trim-qualities")
    .help("If set, trim low-quality bases (< --trim-min-quality) when using "
          "the 'per-base' trimming strategy")
    .deprecated_alias("--trimqualities")
    .conflicts_with("--trim-mott-quality")
    .conflicts_with("--trim-windows")
    .bind_bool(&trim_low_quality_bases);

  argparser.add_separator();
  argparser.add("--pre-trim-polyx", "X")
    .help("Enable trimming of poly-X tails prior to read alignment and adapter "
          "trimming. Zero or more nucleotides (A, C, G, T) may be specified. "
          "Zero or more nucleotides may be specified after the option "
          "separated by spaces, with zero nucleotides corresponding to all of "
          "A, C, G, and T")
    .bind_vec(&pre_trim_poly_x_sink)
    .with_min_values(0);
  argparser.add("--post-trim-polyx", "X")
    .help("Enable trimming of poly-X tails after read alignment and adapter "
          "trimming/merging, but before trimming of low-quality bases. Merged "
          "reads are not trimmed by this option (both ends are 5'). Zero or "
          "more nucleotides (A, C, G, T) may be specified. Zero or more "
          "nucleotides may be specified after the option separated by spaces, "
          "with zero nucleotides corresponding to all of A, C, G, and T")
    .bind_vec(&post_trim_poly_x_sink)
    .with_min_values(0);
  argparser.add("--trim-polyx-threshold", "N")
    .help("The minimum number of bases in a poly-X tail")
    .bind_u32(&trim_poly_x_threshold)
    .with_default(10);

  argparser.add_separator();
  argparser.add("--preserve5p")
    .help("If set, bases at the 5p will not be trimmed by when performing "
          "quality based trimming of reads. Merged reads will not be quality "
          "trimmed when this option is enabled [default: 5p bases are trimmed]")
    .bind_bool(&preserve5p);

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("FILTERING:");

  argparser.add("--max-ns", "N")
    .help("Reads containing more ambiguous bases (N) than this number after "
          "trimming are discarded [default: no maximum]")
    .deprecated_alias("--maxns")
    .bind_u32(&max_ambiguous_bases)
    .with_default(std::numeric_limits<uint32_t>::max());

  argparser.add("--min-length", "N")
    .help("Reads shorter than this length following trimming are discarded")
    .deprecated_alias("--minlength")
    .bind_u32(&min_genomic_length)
    .with_default(15);
  argparser.add("--max-length", "N")
    .help("Reads longer than this length following trimming are discarded "
          "[default: no maximum]")
    .deprecated_alias("--maxlength")
    .bind_u32(&max_genomic_length)
    .with_default(std::numeric_limits<uint32_t>::max());

  argparser.add("--min-mean-quality", "N")
    .help("Reads with a mean Phred quality score less than this value "
          "following trimming are discarded. The value must be in the range 0 "
          "to 93, corresponding to Phred+33 encoded values of '!' to '~' "
          "[default: no minimum]")
    .bind_double(&min_mean_quality)
    .with_minimum(0.0)
    .with_maximum(PHRED_SCORE_MAX)
    .with_default(0.0);

  argparser.add("--min-complexity", "X")
    .help(
      "Filter reads with a complexity score less than this value. Complexity "
      "is measured as the fraction of positions that differ from the previous "
      "position. A suggested value is 0.3 [default: no minimum]")
    .bind_double(&min_complexity)
    .with_minimum(0.0)
    .with_maximum(1.0)
    .with_default(0);

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("DEMULTIPLEXING:");

  argparser.add("--barcode-list", "FILE")
    .help("List of barcodes or barcode pairs for single or double-indexed "
          "demultiplexing. Note that both indexes should be specified for "
          "both single-end and paired-end trimming, if double-indexed "
          "multiplexing was used, in order to ensure that the demultiplexed "
          "reads can be trimmed correctly")
    .bind_str(&barcode_list);
  argparser.add("--multiple-barcodes")
    .help("Allow for more than one barcode (pair) for each sample. If this "
          "option is not specified, AdapterRemoval will abort if multiple "
          "barcodes/barcode pairs identify the same sample");
  argparser.add("--mixed-orientation", "X")
    .help("Process barcodes be sequences in both the barcode1-insert-barcode2 "
          "(forward) orientation and barcode2-insert-barcode1 (reverse) "
          "orientation. Takes an optional argument specifying the orientation "
          "of the barcodes in the `--barcode-list`, defaulting to `forward`")
    .deprecated_alias("--reversible-barcodes")
    .depends_on("--barcode-list")
    .bind_str(nullptr)
    .with_default("unspecified")
    .with_implicit_argument("forward")
    .with_choices({ "unspecified", "forward", "reverse", "explicit" });

  argparser.add_separator();
  argparser.add("--barcode-mm", "N")
    .help("Maximum number of mismatches allowed when counting mismatches in "
          "both the mate 1 and the mate 2 barcode for paired reads")
    .bind_u32(&barcode_mm)
    .with_default(0);
  argparser.add("--barcode-mm-r1", "N")
    .help("Maximum number of mismatches allowed for the mate 1 barcode. "
          "Cannot be higher than the --barcode-mm value [default: same value "
          "as --barcode-mm]")
    .bind_u32(&barcode_mm_r1)
    .with_default(0);
  argparser.add("--barcode-mm-r2", "N")
    .help("Maximum number of mismatches allowed for the mate 2 barcode. "
          "Cannot be higher than the --barcode-mm value [default: same value "
          "as --barcode-mm]")
    .bind_u32(&barcode_mm_r2)
    .with_default(0);
  argparser.add("--demultiplex-only")
    .help("Only carry out demultiplexing using the list of barcodes "
          "supplied with --barcode-list. No other processing is done")
    .depends_on("--barcode-list")
    .conflicts_with("--report-only");

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("REPORTS:");

  argparser.add("--report-only")
    .help("Write a report of the input data without performing any processing "
          "of the FASTQ reads. Adapter sequence inference is performed for PE "
          "data based on overlapping mate reads. A report including read "
          "processing, but without output, can be generated by setting "
          "--output options to /dev/null")
    .deprecated_alias("--identify-adapters")
    .conflicts_with("--barcode-list")
    .conflicts_with("--benchmark")
    .conflicts_with("--demultiplex-only");

  argparser.add("--report-title", "X")
    .help("Title used for HTML report")
    .bind_str(&report_title)
    .with_default(NAME + " " + VERSION);
  argparser.add("--report-sample-rate", "X")
    .help("Fraction of reads to use when generating base quality/composition "
          "curves for trimming reports. Using all data (--report-sample-nth "
          "1.0) results in an 10-30% decrease in throughput")
    .bind_double(&report_sample_rate)
    .with_minimum(0.0)
    .with_maximum(1.0)
    .with_default(0.1);
  argparser.add("--report-duplication", "N")
    .help("FastQC based duplicate detection, based on the frequency of the "
          "first N unique sequences observed. If no value is given, an N of "
          "100k is used, corresponding to FastQC defaults; a value of 0 "
          "disables the analysis. Accepts suffixes K, M, and G")
    .bind_str(nullptr)
    .with_implicit_argument("100k");

  //////////////////////////////////////////////////////////////////////////////
  argparser.add_header("LOGGING:");

  argparser.add("--log-level", "X")
    .help("The minimum severity of messages to be written to STDERR")
    .bind_str(&log_level)
    .with_choices({ "debug", "info", "warning", "error" })
    .with_default("info");

  argparser.add("--log-colors", "X")
    .help("Enable/disable the use of colors when writing log messages. If set "
          "to auto, colors will only be enabled if STDERR is a terminal and "
          "the NO_COLORS is environmental variable is not set")
    .bind_str(&log_color)
    .with_choices({ "auto", "always", "never" })
    .with_default("auto");
  argparser.add("--log-progress", "X")
    .help("Specify the type of progress reports used. If set to auto, then a "
          "spinner will be used if STDERR is a terminal and the NO_COLORS "
          "environmental variable is not set, otherwise logging will be used")
    .bind_str(nullptr)
    .with_choices({ "auto", "log", "spin", "never" })
    .with_default("auto");
}

argparse::parse_result
userconfig::parse_args(const string_vec& argvec)
{
  args = argvec;
  if (args.size() <= 1) {
    argparser.print_help();
    return argparse::parse_result::error;
  }

  // ad-hoc arg parsing to make argparse output consistent with rest of run
  configure_log_colors(try_parse_argument(args, "--log-color", "auto"), true);
  configure_log_levels(try_parse_argument(args, "--log-level", "info"), true);

  const argparse::parse_result result = argparser.parse_args(args);
  if (result != argparse::parse_result::ok) {
    return result;
  }

  configure_log_colors(log_color);
  configure_log_levels(log_level);
  log_progress = configure_log_progress(argparser.value("--log-progress"));

  {
    const auto degenerate = argparser.is_set("--mask-degenerate-bases")
                              ? degenerate_encoding::mask
                              : degenerate_encoding::reject;
    const auto uracils = argparser.is_set("--convert-uracils")
                           ? uracil_encoding::convert
                           : uracil_encoding::reject;

    io_encoding = configure_encoding(quality_input_base, degenerate, uracils);
  }

  if (argparser.is_set("--mate-separator")) {
    if (mate_separator_str.size() != 1) {
      log::error() << "The argument for --mate-separator must be "
                      "exactly one character long, not "
                   << mate_separator_str.size() << " characters!";
      return argparse::parse_result::error;
    } else {
      mate_separator = mate_separator_str.at(0);
    }
  }

  if (argparser.is_set("--demultiplex-only")) {
    run_type = ar_command::demultiplex_only;
  } else if (argparser.is_set("--report-only")) {
    run_type = ar_command::report_only;
  } else if (argparser.is_set("--benchmark")) {
    run_type = ar_command::benchmark;
  }

  {
    const auto strategy = argparser.value("--quality-trimming");
    if (strategy == "mott") {
      trim = trimming_strategy::mott;
      trim_mott_rate = std::pow(10.0, trim_mott_rate / -10.0);
    } else if (strategy == "window") {
      trim = trimming_strategy::window;
    } else if (strategy == "per-base") {
      trim = trimming_strategy::per_base;

      if (!trim_low_quality_bases && !trim_ambiguous_bases) {
        log::error() << "The per-base quality trimming strategy is enabled, "
                     << "but neither trimming of low-quality bases (via "
                     << "--trim-qualities) nor trimming of Ns (via --trim-ns) "
                     << "is enabled.";
        return argparse::parse_result::error;
      }
    } else if (strategy == "none") {
      trim = trimming_strategy::none;
    } else {
      AR_FAIL(shell_escape(strategy));
    }
  }

  // Check for invalid combinations of settings
  if (input_files_1.empty() && input_files_2.empty()) {
    log::error()
      << "No input files (--in-file1 / --in-file2) specified.\n"
      << "Please specify at least one input file using --in-file1 FILENAME.";

    return argparse::parse_result::error;
  } else if (!input_files_2.empty() &&
             (input_files_1.size() != input_files_2.size())) {
    log::error()
      << "Different number of files specified for --in-file1 and --in-file2.";

    return argparse::parse_result::error;
  } else if (!input_files_2.empty()) {
    paired_ended_mode = true;
  }

  interleaved_input |= interleaved;
  interleaved_output |= interleaved;

  if (interleaved_input) {
    // Enable paired end mode .. other than the FASTQ reader, all other
    // parts of the pipeline simply run in paired-end mode.
    paired_ended_mode = true;
  }

  if (paired_ended_mode) {
    min_adapter_overlap = 0;

    // merge related options implies --merge
    if (argparser.is_set("--collapse-deterministic")) {
      merge = merge_strategy::additive;
    } else if (argparser.is_set("--collapse-conservatively")) {
      merge = merge_strategy::maximum;
    } else if (argparser.is_set("--merge") ||
               argparser.is_set("--merge-strategy")) {
      const auto strategy = argparser.value("--merge-strategy");
      if (strategy == "maximum") {
        merge = merge_strategy::maximum;
      } else if (strategy == "additive") {
        merge = merge_strategy::additive;
      } else {
        AR_FAIL(strategy);
      }
    }
  }

  // (Optionally) read adapters from file and validate
  if (!setup_adapter_sequences()) {
    return argparse::parse_result::error;
  }

  // (Optionally) read barcodes from file and validate
  if (!setup_demultiplexing()) {
    return argparse::parse_result::error;
  }

  if (argparser.is_set("--read-group")) {
    auto merged_tags = join_text(read_group, "\t");

    try {
      samples.set_read_group(merged_tags);
    } catch (const std::invalid_argument& error) {
      log::error() << "Invalid argument --read-group "
                   << log_escape(merged_tags) << ": " << error.what();

      return argparse::parse_result::error;
    }
  }

  // Set mismatch threshold
  if (mismatch_threshold > 1) {
    mismatch_threshold = 1.0 / mismatch_threshold;
  } else if (mismatch_threshold < 0) {
    if (run_type == ar_command::report_only) {
      mismatch_threshold = 1.0 / 10.0;
    } else {
      // Defaults for PE / SE trimming (changed in v3)
      mismatch_threshold = 1.0 / 6.0;
    }
  }

  {
    bool found = false;
    const auto simd_choice = argparser.value("--simd");
    for (const auto is : simd::supported()) {
      if (simd_choice == simd::name(is)) {
        simd = is;
        found = true;
        break;
      }
    }

    AR_REQUIRE(found);
  }

  using fixed_trimming =
    std::tuple<const char*, const string_vec&, std::pair<unsigned, unsigned>&>;

  const std::vector<fixed_trimming> fixed_trimming_options = {
#ifdef PRE_TRIM_5P
    { "--pre-trim5p", pre_trim5p, pre_trim_fixed_5p },
#endif
    { "--pre-trim3p", pre_trim3p, pre_trim_fixed_3p },
    { "--post-trim5p", post_trim5p, post_trim_fixed_5p },
    { "--post-trim3p", post_trim3p, post_trim_fixed_3p },
  };

  for (const auto& it : fixed_trimming_options) {
    try {
      if (argparser.is_set(std::get<0>(it))) {
        std::get<2>(it) = parse_trim_argument(std::get<1>(it));
      }
    } catch (const std::invalid_argument& error) {
      log::error() << "Could not parse " << std::get<0>(it)
                   << " argument(s): " << error.what();

      return argparse::parse_result::error;
    }
  }

  if (!parse_output_formats(argparser, out_file_format, out_stdout_format)) {
    return argparse::parse_result::error;
  }

  // An empty prefix or directory would results in the creation of dot-files
  if (out_prefix.empty()) {
    log::error() << "--out-prefix must be a non-empty value.";

    return argparse::parse_result::error;
  } else if (out_prefix.back() == '/') {
    log::error() << "--out-prefix must not be a directory: "
                 << shell_escape(out_prefix);

    return argparse::parse_result::error;
  } else if (out_prefix == DEV_NULL && run_type != ar_command::benchmark) {
    // Relevant output options depend on input files and other settings
    const std::vector<std::pair<std::string, bool>> output_keys = {
      { "--out-file1",
        is_adapter_trimming_enabled() || is_demultiplexing_enabled() },
      { "--out-file2",
        is_adapter_trimming_enabled() || is_demultiplexing_enabled() },
      { "--out-singleton", is_any_filtering_enabled() },
      { "--out-merged", is_read_merging_enabled() },
      { "--out-discarded", is_any_filtering_enabled() },
      { "--out-unidentified1", is_demultiplexing_enabled() },
      { "--out-unidentified2", is_demultiplexing_enabled() },
      { "--out-json", true },
      { "--out-html", true },
      { "--out-prefix", true },
    };

    string_vec required_keys;
    for (const auto& it : output_keys) {
      if (it.second) {
        required_keys.push_back(it.first);
      }
    }

    const auto user_keys = user_supplied_keys(argparser, required_keys);
    if (user_keys.empty()) {
      auto error = log::error();
      error << "No output would be generated; at least one of the options "
            << join_text(required_keys, ", ", ", or ")
            << " must be used. The --out-prefix option automatically enables "
               "all relevant --out options.";

      return argparse::parse_result::error;
    }
  }

  {
    const std::string key = "--pre-trim-polyx";
    if (argparser.is_set(key) &&
        !parse_poly_x_option(key, pre_trim_poly_x_sink, pre_trim_poly_x)) {
      return argparse::parse_result::error;
    }
  }

  {
    const std::string key = "--post-trim-polyx";
    if (argparser.is_set(key) &&
        !parse_poly_x_option(key, post_trim_poly_x_sink, post_trim_poly_x)) {
      return argparse::parse_result::error;
    }
  }

  if (!min_genomic_length) {
    log::warn() << "--min-length is set to 0. This may produce FASTQ files "
                   "that are incompatible with some tools!";
  }

  // Default to all reads, but don't print value with --help
  head = std::numeric_limits<uint64_t>::max();
  if (!parse_counts(argparser, "--head", head)) {
    return argparse::parse_result::error;
  }

  if (!parse_counts(argparser, "--report-duplication", report_duplication)) {
    return argparse::parse_result::error;
  }

  return argparse::parse_result::ok;
}

bool
userconfig::is_good_alignment(const alignment_info& alignment) const
{
  if (!alignment.length || alignment.score() <= 0) {
    return false;
  }

  // Only pairs of called bases are considered part of the alignment
  const size_t n_aligned = alignment.length - alignment.n_ambiguous;
  if (n_aligned < min_adapter_overlap && !paired_ended_mode) {
    return false;
  }

  auto mm_threshold = static_cast<size_t>(mismatch_threshold * n_aligned);
  if (n_aligned < 6) {
    mm_threshold = 0;
  } else if (n_aligned < 10) {
    // Allow at most 1 mismatch, possibly set to 0 by the user
    mm_threshold = std::min<size_t>(1, mm_threshold);
  }

  return alignment.n_mismatches <= mm_threshold;
}

bool
userconfig::can_merge_alignment(const alignment_info& alignment) const
{
  if (alignment.length < alignment.n_ambiguous) {
    throw std::invalid_argument("#ambiguous bases > read length");
  }

  return alignment.length - alignment.n_ambiguous >= merge_threshold;
}

output_format
userconfig::infer_output_format(const std::string& filename) const
{
  if (filename == DEV_STDOUT) {
    return out_stdout_format;
  }

  output_format result = out_file_format;
  // Parse failures are ignored here; default to --out-format
  output_files::parse_extension(filename, result);

  return result;
}

output_files
userconfig::get_output_filenames() const
{
  output_files files;

  files.settings_json = new_output_file("--out-json", {}, {}, ".json").name;
  files.settings_html = new_output_file("--out-html", {}, {}, ".html").name;

  auto ext = output_files::file_extension(out_file_format);
  std::string_view out1 = interleaved_output ? "" : ".r1";
  std::string_view out2 = interleaved_output ? "" : ".r2";

  if (is_demultiplexing_enabled()) {
    files.unidentified_1 = new_output_file("--out-unidentified1",
                                           {},
                                           { ".unidentified", out1 },
                                           ext);

    if (paired_ended_mode) {
      if (interleaved_output) {
        files.unidentified_2 = files.unidentified_1;
      } else {
        files.unidentified_2 = new_output_file("--out-unidentified2",
                                               {},
                                               { ".unidentified", out2 },
                                               ext);
      }
    }
  }

  for (const auto& sample : samples) {
    const auto& name = sample.name();
    sample_output_files map;

    const auto mate_1 = new_output_file("--out-file1", name, { out1 }, ext);
    map.set_file(read_file::mate_1, mate_1);

    if (paired_ended_mode) {
      if (interleaved_output) {
        map.set_file(read_file::mate_2, mate_1);
      } else {
        map.set_file(read_file::mate_2,
                     new_output_file("--out-file2", name, { out2 }, ext));
      }
    }

    if (run_type == ar_command::trim_adapters) {
      if (is_any_filtering_enabled()) {
        map.set_file(
          read_file::discarded,
          new_output_file("--out-discarded", name, { ".discarded" }, ext));
      }

      if (paired_ended_mode) {
        if (is_any_filtering_enabled()) {
          map.set_file(
            read_file::singleton,
            new_output_file("--out-singleton", name, { ".singleton" }, ext));
        }

        if (is_read_merging_enabled()) {
          map.set_file(
            read_file::merged,
            new_output_file("--out-merged", name, { ".merged" }, ext));
        }
      }
    }

    files.add_sample(std::move(map));
  }

  return files;
}

output_file
userconfig::new_output_file(const std::string& key,
                            std::string_view sample,
                            std::vector<std::string_view> keys,
                            std::string_view ext) const
{
  AR_REQUIRE(!ext.empty());
  const auto default_is_fastq = out_file_format == output_format::fastq ||
                                out_file_format == output_format::fastq_gzip;

  std::string out;
  if (argparser.is_set(key)) {
    out = argparser.value(key);

    // global files, e.g. reports and unidentified reads
    if (sample.empty()) {
      return { out, infer_output_format(out) };
    }
  } else if (default_is_fastq && key == "--out-discarded") {
    // Discarded reads are dropped by default for non-archival formats
    out = DEV_NULL;
  } else {
    out = out_prefix;
  }

  if (out == DEV_NULL) {
    return { out, output_format::fastq };
  }

  if (!(default_is_fastq || keys.empty())) {
    // SAM/BAM files are combined by default
    keys.pop_back();
  }

  if (!sample.empty()) {
    keys.insert(keys.begin(), sample);
  }

  keys.emplace_back(ext);

  for (const auto& value : keys) {
    if (!value.empty() && value.front() != '.') {
      out.push_back('.');
    }

    out.append(value);
  }

  return output_file{ out, infer_output_format(out) };
}

bool
check_and_set_barcode_mm(const argparse::parser& argparser,
                         const std::string& key,
                         uint32_t barcode_mm,
                         uint32_t& dst)
{
  if (!argparser.is_set(key)) {
    dst = barcode_mm;
  } else if (dst > barcode_mm) {
    log::error()
      << "The maximum number of errors for " << key
      << " is set \n"
         "to a higher value than the total number of mismatches allowed\n"
         "for barcodes (--barcode-mm). Please correct these settings.";
    return false;
  }

  return true;
}

bool
userconfig::is_adapter_trimming_enabled() const
{
  return run_type == ar_command::trim_adapters;
}

bool
userconfig::is_demultiplexing_enabled() const
{
  return !barcode_list.empty();
}

bool
userconfig::is_read_merging_enabled() const
{
  return is_adapter_trimming_enabled() && merge != merge_strategy::none;
}

bool
userconfig::is_any_quality_trimming_enabled() const
{
  return is_adapter_trimming_enabled() &&
         (is_low_quality_trimming_enabled() ||
          is_terminal_base_pre_trimming_enabled() ||
          is_terminal_base_post_trimming_enabled() ||
          is_poly_x_tail_pre_trimming_enabled() ||
          is_poly_x_tail_post_trimming_enabled());
}

bool
userconfig::is_low_quality_trimming_enabled() const
{
  return trim != trimming_strategy::none;
}

bool
userconfig::is_terminal_base_pre_trimming_enabled() const
{
  return
#ifdef PRE_TRIM_5P
    pre_trim_fixed_5p.first || pre_trim_fixed_5p.second ||
#endif
    pre_trim_fixed_3p.first || pre_trim_fixed_3p.second;
}

bool
userconfig::is_terminal_base_post_trimming_enabled() const
{
  return post_trim_fixed_5p.first || post_trim_fixed_5p.second ||
         post_trim_fixed_3p.first || post_trim_fixed_3p.second;
}

bool
userconfig::is_poly_x_tail_pre_trimming_enabled() const
{
  return !pre_trim_poly_x.empty();
}

bool
userconfig::is_poly_x_tail_post_trimming_enabled() const
{
  return !post_trim_poly_x.empty();
}

bool
userconfig::is_any_filtering_enabled() const
{
  return is_adapter_trimming_enabled() &&
         (is_short_read_filtering_enabled() ||
          is_long_read_filtering_enabled() ||
          is_ambiguous_base_filtering_enabled() ||
          is_mean_quality_filtering_enabled() ||
          is_low_complexity_filtering_enabled());
}

bool
userconfig::is_short_read_filtering_enabled() const
{
  return min_genomic_length > 0;
}

bool
userconfig::is_long_read_filtering_enabled() const
{
  return max_genomic_length !=
         std::numeric_limits<decltype(max_genomic_length)>::max();
}

bool
userconfig::is_ambiguous_base_filtering_enabled() const
{
  return max_ambiguous_bases !=
         std::numeric_limits<decltype(max_ambiguous_bases)>::max();
}

bool
userconfig::is_mean_quality_filtering_enabled() const
{
  return min_mean_quality > 0;
}

bool
userconfig::is_low_complexity_filtering_enabled() const
{
  return min_complexity > 0;
}

bool
userconfig::setup_adapter_sequences()
{
  adapter_set adapters;
  if (argparser.is_set("--adapter-list")) {
    try {
      adapters.load(adapter_list, paired_ended_mode);
    } catch (const std::exception& error) {
      log::error() << "Error reading adapters from " << log_escape(adapter_list)
                   << ": " << error.what();
      return false;
    }

    log::info() << "Read " << adapters.size()
                << " adapters / adapter pairs from '" << adapter_list << "'";
  } else {
    try {
      adapters.add(dna_sequence{ adapter_1 }, dna_sequence{ adapter_2 });
    } catch (const fastq_error& error) {
      log::error() << "Error parsing adapter sequence(s):\n"
                   << "   " << error.what();

      return false;
    }
  }

  samples.set_adapters(std::move(adapters));

  return true;
}

bool
userconfig::setup_demultiplexing()
{
  if (!argparser.is_set("--barcode-mm")) {
    barcode_mm = barcode_mm_r1 + barcode_mm_r2;
  }

  if (!check_and_set_barcode_mm(argparser,
                                "--barcode-mm-r1",
                                barcode_mm,
                                barcode_mm_r1)) {
    return false;
  }

  if (!check_and_set_barcode_mm(argparser,
                                "--barcode-mm-r2",
                                barcode_mm,
                                barcode_mm_r2)) {
    return false;
  }

  if (argparser.is_set("--barcode-list")) {
    const auto orientation =
      parse_table_orientation(argparser.value("--mixed-orientation"));

    barcode_config config;
    config.paired_end_mode(paired_ended_mode)
      .allow_multiple_barcodes(argparser.is_set("--multiple-barcodes"))
      .orientation(orientation);

    try {
      samples.load(barcode_list, config);
    } catch (const std::exception& error) {
      log::error() << "Error reading barcodes from " << log_escape(barcode_list)
                   << ": " << error.what();
      return false;
    }

    log::info() << "Read " << samples.size() << " sets of barcodes from "
                << shell_escape(barcode_list);
  }

  const auto& output_files = get_output_filenames();

  return check_input_files(input_files_1, input_files_2) &&
         check_output_files("--in-file1", input_files_1, output_files) &&
         check_output_files("--in-file2", input_files_2, output_files);
}

} // namespace adapterremoval
