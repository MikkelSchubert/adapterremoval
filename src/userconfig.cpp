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
#include <algorithm> // for find, min, max
#include <cmath>     // for pow
#include <cstring>   // for size_t, strerror
#include <errno.h>   // for errno
#include <limits>    // for numeric_limits
#include <stdexcept> // for invalid_argument
#include <string>    // for string, operator<<, char_traits, operator==
#include <tuple>     // for tuple
#include <unistd.h>  // for access, isatty, R_OK, STDERR_FILENO

#include "alignment.hpp"  // for alignment_info
#include "debug.hpp"      // for AR_FAIL
#include "logging.hpp"    // for log
#include "progress.hpp"   // for progress_type
#include "strutils.hpp"   // for str_to_unsigned, toupper
#include "userconfig.hpp" // declarations

namespace adapterremoval {

const size_t output_sample_files::disabled = std::numeric_limits<size_t>::max();

////////////////////////////////////////////////////////////////////////////////
// Helper functions

std::pair<unsigned, unsigned>
parse_trim_argument(const string_vec& values)
{
  unsigned mate_1 = 0;
  unsigned mate_2 = 0;

  switch (values.size()) {
    case 1:
      mate_1 = str_to_unsigned(values.front());
      mate_2 = mate_1;
      break;

    case 2:
      mate_1 = str_to_unsigned(values.front());
      mate_2 = str_to_unsigned(values.back());
      break;

    default:
      throw std::invalid_argument("please specify exactly one or two values");
  }

  return std::pair<unsigned, unsigned>(mate_1, mate_2);
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

  std::array<bool, ACGT::size> enabled = {};
  for (const auto& value : values) {
    for (const auto nuc : toupper(value)) {
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
check_no_clobber(const std::string& label,
                 const string_vec& in_files,
                 const std::string& out_file)
{
  for (const auto& in_file : in_files) {
    if (in_file == out_file && in_file != DEV_NULL) {
      log::error() << "Input file would be overwritten: " << label << " "
                   << in_file;
      return false;
    }
  }

  return true;
}

bool
check_input_and_output(const std::string& label,
                       const string_vec& filenames,
                       const output_files& output_files)
{
  for (const auto& filename : filenames) {
    if (access(filename.c_str(), R_OK)) {
      log::error() << "Cannot read file: " << label << " " << filename
                   << "': " << std::strerror(errno);

      return false;
    }
  }

  if (!check_no_clobber(label, filenames, output_files.unidentified_1)) {
    return false;
  }

  if (!check_no_clobber(label, filenames, output_files.unidentified_2)) {
    return false;
  }

  for (const auto& sample : output_files.samples) {
    for (const auto& out_file : sample.filenames()) {
      if (!check_no_clobber(label, filenames, out_file)) {
        return false;
      }
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `output_files`

output_files::output_files()
  : settings_json()
  , settings_html()
  , unidentified_1()
  , unidentified_2()
  , samples()
{
}

output_sample_files::output_sample_files()
  : m_filenames()
  , m_pipeline_steps()
  , m_offsets()
{
  m_offsets.fill(output_sample_files::disabled);
}

void
output_sample_files::set_filename(const read_type rtype,
                                  const std::string& filename)
{
  const auto index = static_cast<size_t>(rtype);
  AR_REQUIRE(m_offsets.at(index) == output_sample_files::disabled);

  // If the file type isn't being saved, then there is no need to process the
  // reads. This saves time especially when output compression is enabled.
  if (filename != DEV_NULL) {
    // FIXME: This assumes that filesystem is case sensitive
    auto it = std::find(m_filenames.begin(), m_filenames.end(), filename);
    if (it == m_filenames.end()) {
      m_filenames.push_back(filename);

      m_offsets.at(index) = m_filenames.size() - 1;
    } else {
      m_offsets.at(index) = it - m_filenames.begin();
    }
  }
}

const string_vec&
output_sample_files::filenames() const
{
  return m_filenames;
}

const std::vector<size_t>&
output_sample_files::pipeline_steps() const
{
  return m_pipeline_steps;
}

void
output_sample_files::push_pipeline_step(size_t step)
{
  AR_REQUIRE(m_pipeline_steps.size() < m_filenames.size());
  m_pipeline_steps.push_back(step);
}

size_t
output_sample_files::offset(read_type value) const
{
  return m_offsets.at(static_cast<size_t>(value));
}

////////////////////////////////////////////////////////////////////////////////

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
output_files_set(const argparse::parser& argparser)
{
  const string_vec keys = {
    "--out-json",      "--out-html",   "--out-file1",     "--out-file2",
    "--out-singleton", "--out-merged", "--out-discarded",
  };

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
    const char* no_color = ::getenv("NO_COLOR");

    return !(no_color && no_color[0] != '\0');
  }

  return false;
}

void
configure_log_levels(const std::string& value, bool fallible = false)
{
  const auto log_level = tolower(value);

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
configure_encoding(const std::string& value)
{
  if (value == "33") {
    return FASTQ_ENCODING_33;
  } else if (value == "64") {
    return FASTQ_ENCODING_64;
  } else if (value == "solexa") {
    return FASTQ_ENCODING_SOLEXA;
  } else if (value == "sam") {
    return FASTQ_ENCODING_SAM;
  }

  AR_FAIL("unhandled qualitybase value");
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `userconfig`

userconfig::userconfig(const std::string& name,
                       const std::string& version,
                       const std::string& help)
  : args()
  , run_type(ar_command::trim_adapters)
  , input_files_1()
  , input_files_2()
  , out_basename()
  , out_json()
  , out_html()
  , out_pair_1()
  , out_pair_2()
  , out_merged()
  , out_discarded()
  // FIXME: Support both .1 and .2
  , out_singleton()
  , prefix_read_1()
  , prefix_read_2()
  , prefix_merged()
  , paired_ended_mode()
  , interleaved_input()
  , interleaved_output()
  , head()
  , mate_separator()
  , min_genomic_length()
  , max_genomic_length()
  , min_adapter_overlap()
  , min_alignment_length()
  , mismatch_threshold()
  , io_encoding(FASTQ_ENCODING_33)
  , pre_trim_fixed_5p()
  , pre_trim_fixed_3p()
  , post_trim_fixed_5p()
  , post_trim_fixed_3p()
  , trim_error_rate()
  , pre_trim_poly_x()
  , post_trim_poly_x()
  , trim_poly_x_threshold()
  , trim_by_quality()
  , trim_window_length()
  , low_quality_score()
  , trim_ambiguous_bases()
  , max_ambiguous_bases()
  , min_complexity()
  , preserve5p()
  , merge()
  , shift()
  , max_threads()
  , gzip()
  , gzip_level()
  , barcode_mm()
  , barcode_mm_r1()
  , barcode_mm_r2()
  , adapters()
  , report_sample_rate()
  , report_duplication()
  , log_progress()
  , argparser(name, version, help)
  , adapter_1()
  , adapter_2()
  , adapter_list()
  , barcode_list()
  , quality_input_base()
  , mate_separator_str()
  , interleaved()
  , pre_trim5p()
  , pre_trim3p()
  , post_trim5p()
  , post_trim3p()
  , pre_trim_poly_x_sink()
  , post_trim_poly_x_sink()
  , log_color()
  , log_level()
  , log_progress_sink()
  , m_runtime()
  , m_deprecated_knobs()
{
  argparser.add("--file1", "FILE")
    .help("Input files containing mate 1 reads or single-ended reads; one or "
          "more files may be listed [REQUIRED]")
    .bind_vec(&input_files_1);
  argparser.add("--file2", "FILE")
    .help("Input files containing mate 2 reads; if used, then the same number "
          "of files as --file1 must be listed [OPTIONAL]")
    .bind_vec(&input_files_2);
  argparser.add("--head", "N")
    .help("Process only the first N reads in single-end mode or the first N "
          "read-pairs in paired-end mode [default: all reads]")
    .bind_uint(&head)
    .with_default(std::numeric_limits<unsigned>::max());

  argparser.add("--identify-adapters")
    .help("Attempt to identify the adapter pair of PE reads, by searching for "
          "overlapping mate reads")
    .conflicts_with("--demultiplex-only")
    .conflicts_with("--report-only");

  argparser.add("--threads", "N")
    .help("Maximum number of threads")
    .bind_uint(&max_threads)
    .with_default(1);

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

  argparser.add("--interleaved")
    .help("This option enables both the --interleaved-input option and the "
          "--interleaved-output option")
    .conflicts_with("--file2")
    .conflicts_with("--out-file2")
    .bind_bool(&interleaved);
  argparser.add("--interleaved-input")
    .help("The (single) input file provided contains both the mate 1 and mate "
          "2 reads, one pair after the other, with one mate 1 reads followed "
          "by one mate 2 read. This option is implied by the --interleaved "
          "option")
    .conflicts_with("--file2")
    .bind_bool(&interleaved_input);
  argparser.add("--interleaved-output")
    .help("If set, trimmed paired-end reads are written to a single file "
          "containing mate 1 and mate 2 reads, one pair after the other. This "
          "option is implied by the --interleaved option")
    .conflicts_with("--out-file2")
    .bind_bool(&interleaved_output);

  argparser.add_header("OUTPUT READ NAMES:");

  argparser.add("--prefix-read1")
    .help("Adds the specified prefix to read 1 names")
    .bind_str(&prefix_read_1);
  argparser.add("--prefix-read2")
    .help("Adds the specified prefix to read 2 names")
    .bind_str(&prefix_read_2);
  argparser.add("--prefix-merged")
    .help("Adds the specified prefix to merged read names")
    .bind_str(&prefix_merged);

  argparser.add_header("OUTPUT FILES:");
  argparser.add("--basename", "PREFIX")
    .help("Default prefix for all output files for which no filename was "
          "explicitly set")
    .bind_str(&out_basename)
    .with_default("your_output");
  argparser.add("--out-json", "FILE")
    .help("Output file containing statistics about trimming, merging, and more "
          "in JSON format")
    .bind_str(&out_json)
    .with_default("{basename}[.sample].json");
  argparser.add("--out-html", "FILE")
    .help("Output report containing statistics about trimming, merging, and "
          "more in HTML format")
    .bind_str(&out_html)
    .with_default("{basename}[.sample].html");

  argparser.add("--out-file1", "FILE")
    .help("Output file containing trimmed mate 1 reads. When performing "
          "demultiplexing, this path is treated as a prefix for output files ")
    .deprecated_alias("--output1")
    .bind_str(&out_pair_1)
    .with_default("{basename}[.sample].r1.fastq");
  argparser.add("--out-file2", "FILE")
    .help("Output file containing trimmed mate 2 reads")
    .deprecated_alias("--output2")
    .bind_str(&out_pair_2)
    .with_default("{basename}[.sample].r2.fastq");
  argparser.add("--out-singleton", "FILE")
    .help("Output file to which containing paired reads for which the mate "
          "has been discarded")
    .deprecated_alias("--singleton")
    .bind_str(&out_singleton)
    .with_default("{basename}[.sample].singleton.fastq");
  argparser.add("--out-merged", "FILE")
    .help("If --merge is set, contains overlapping mate-pairs which "
          "have been merged into a single read (PE mode) or reads for which "
          "the adapter was identified by a minimum overlap, indicating that "
          "the entire template molecule is present. This does not include "
          "which have subsequently been trimmed due to low-quality or "
          "ambiguous nucleotides")
    .deprecated_alias("--outputcollapsed")
    .bind_str(&out_merged)
    .with_default("{basename}[.sample].merged.fastq");
  argparser.add("--out-discarded", "FILE")
    .help("Contains reads discarded due to the --minlength, --maxlength or "
          "--maxns options")
    .deprecated_alias("--discarded")
    .bind_str(&out_discarded)
    .with_default("{basename}[.sample].discarded.fastq");

  argparser.add_header("OUTPUT COMPRESSION:");
  argparser.add("--gzip").help("Enable gzip-compression.").bind_bool(&gzip);
  argparser
    .add("--gzip-level", "N")
#ifdef USE_LIBDEFLATE
    .help("GZip compression level, 0 - 9. For compression levels 4 - 9, output "
          "consist of concatenated GZip blocks, which may cause compatibility "
          "problems in some rare cases. See the documentaiton for more "
          "information")
    .bind_uint(&gzip_level)
    .with_default(6);
#else
    .help("GZip compression level, 0 - 3")
    .bind_uint(&gzip_level)
    .with_default(3);
#endif

  argparser.add_header("TRIMMING SETTINGS:");
  argparser.add("--adapter1", "SEQ")
    .help("Adapter sequence expected to be found in mate 1 reads")
    .bind_str(&adapter_1)
    .with_default("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
  argparser.add("--adapter2", "SEQ")
    .help("Adapter sequence expected to be found in mate 2 reads")
    .bind_str(&adapter_2)
    .with_default("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");
  argparser.add("--adapter-list", "FILE")
    .help("Read table of white-space separated adapters pairs, used as if the "
          "first column was supplied to --adapter1, and the second column was "
          "supplied to --adapter2; only the first adapter in each pair is "
          "required SE trimming mode")
    .bind_str(&adapter_list);

  argparser.add_separator();
  argparser.add("--minadapteroverlap", "N")
    .help("In single-end mode, reads are only trimmed if the overlap between "
          "read and the adapter is at least X bases long, not counting "
          "ambiguous nucleotides (N); this is independent of the "
          "--minalignmentlength when using --merge, allowing a conservative "
          "selection of putative complete inserts while ensuring that all "
          "possible adapter contamination is trimmed")
    .bind_uint(&min_adapter_overlap)
    .with_default(0);
  argparser.add("--mm", "X")
    .help("Max error-rate when aligning reads and/or adapters. If > 1, the max "
          "error-rate is set to 1 / X; if < 0, the defaults are used, "
          "otherwise the user-supplied value is used directly [default: 1/3 "
          "for trimming; 1/10 when identifying adapters]")
    .bind_double(&mismatch_threshold)
    .with_default(-1.0);
  argparser.add("--shift", "N")
    .help("Consider alignments where up to N nucleotides are missing from the "
          "5' termini")
    .bind_uint(&shift)
    .with_default(2);

  argparser.add_separator();
  argparser.add("--pre-trim5p", "N")
    .help("Trim the 5' of reads by a fixed amount after demultiplexing (if "
          "enabled) but before trimming adapters and low quality bases. "
          "Specify one value to trim mate 1 and mate 2 reads the same amount, "
          "or two values separated by a space to trim each mate different "
          "amounts [default: no trimming]")
    .bind_vec(&pre_trim5p)
    .with_max_values(2);
  argparser.add("--pre-trim3p", "N")
    .help("Trim the 3' of reads by a fixed amount. See --pre-trim5p [default: "
          "no trimming]")
    .bind_vec(&pre_trim3p)
    .with_max_values(2);

  argparser.add("--post-trim5p", "N")
    .help("Trim the 5' of reads by a fixed amount after removing adapters, "
          "but before carrying out quality based trimming. See --pre-trim5p "
          "for more information [default: no trimming]")
    .deprecated_alias("--trim5p")
    .bind_vec(&post_trim5p)
    .with_max_values(2);
  argparser.add("--post-trim3p", "N")
    .deprecated_alias("--trim3p")
    .help("Trim the 3' of reads by a fixed amount. See --post-trim5p [default: "
          "no trimming]")
    .bind_vec(&post_trim3p)
    .with_max_values(2);

  argparser.add("--trim-error-rate", "X")
    .help("The threshold value used when performing trimming quality based "
          "trimming using the modified Mott's algorithm. A value of zero or "
          "less disables trimming; a value greater than one is assumed to be "
          "a Phred encoded error rate (e.g. 13 ~= 0.05)")
    .conflicts_with("--trimwindows")
    .conflicts_with("--trimqualities")
    .conflicts_with("--trimns")
    .bind_double(&trim_error_rate)
    .with_default(0);

  argparser.add("--pre-trim-polyx", "X")
    .help("Enable trimming of poly-X tails prior to read alignment and adapter "
          "trimming. Zero or more nucleotides (A, C, G, T) may be specified. "
          "Zero or more nucleotides may be specified after the option "
          "seperated by spaces, with zero nucleotides corresponding to all of "
          "A, C, G, and T")
    .bind_vec(&pre_trim_poly_x_sink)
    .with_min_values(0);
  argparser.add("--post-trim-polyx", "X")
    .help("Enable trimming of poly-X tails after read alignment and adapter "
          "trimming/merging, but before trimming of low-quality bases. Merged "
          "reads are not trimmed by this option (both ends are 5'). Zero or "
          "more nucleotides (A, C, G, T) may be specified. Zero or more "
          "nucleotides may be specified after the option seperated by spaces, "
          "with zero nucleotides corresponding to all of A, C, G, and T")
    .bind_vec(&post_trim_poly_x_sink)
    .with_min_values(0);
  argparser.add("--trim-polyx-threshold", "N")
    .help("Enable trimming of poly-G tails")
    .bind_uint(&trim_poly_x_threshold)
    .with_default(10);

  argparser.add_header("FILTERING SETTINGS:");
  argparser.add("--maxns", "N")
    .help("Reads containing more ambiguous bases (N) than this number after "
          "trimming are discarded")
    .bind_uint(&max_ambiguous_bases)
    .with_default(1000);
  argparser.add("--preserve5p")
    .help("If set, bases at the 5p will not be trimmed by when performing "
          "quality based trimming of reads. Merged reads will not be quality "
          "trimmed when this option is enabled [default: 5p bases are trimmed]")
    .bind_bool(&preserve5p);

  argparser.add_separator();
  argparser.add("--minlength", "N")
    .help("Reads shorter than this length are discarded following trimming")
    .bind_uint(&min_genomic_length)
    .with_default(15);
  argparser.add("--maxlength", "N")
    .help("Reads longer than this length are discarded following trimming")
    .bind_uint(&max_genomic_length)
    .with_default(std::numeric_limits<unsigned>::max());

  argparser.add("--min-complexity", "X")
    .help("Filter low-complexity reads after carrying out adapter trimming and "
          "trimming of low-quality bases. Complexity is a value in the range 0 "
          "to 1, with 0 being the least complex reads")
    .bind_double(&min_complexity)
    .with_default(0);

  argparser.add_header("READ MERGING:");
  argparser.add("--merge")
    .help("When set, paired ended read alignments of --minalignmentlength or "
          "more bases are merged into a single consensus sequence. Merged "
          "reads are written to basename.merged by default. Has no effect "
          "in single-end mode")
    .deprecated_alias("--collapse")
    .bind_bool(&merge);
  argparser.add("--minalignmentlength", "N")
    .help("If --merge is set, paired reads must overlap at least this "
          "number of bases to be merged, and single-ended reads must "
          "overlap at least this number of bases with the adapter to be "
          "considered complete template molecules")
    .bind_uint(&min_alignment_length)
    .with_default(11);

  argparser.add_header("DEMULTIPLEXING:");
  argparser.add("--barcode-list", "FILE")
    .help("List of barcodes or barcode pairs for single or double-indexed "
          "demultiplexing. Note that both indexes should be specified for "
          "both single-end and paired-end trimming, if double-indexed "
          "multiplexing was used, in order to ensure that the demultiplexed "
          "reads can be trimmed correctly")
    .depends_on("--basename")
    .bind_str(&barcode_list);
  argparser.add("--barcode-mm", "N")
    .help("Maximum number of mismatches allowed when counting mismatches in "
          "both the mate 1 and the mate 2 barcode for paired reads")
    .bind_uint(&barcode_mm)
    .with_default(0);
  argparser.add("--barcode-mm-r1", "N")
    .help("Maximum number of mismatches allowed for the mate 1 barcode. "
          "Cannot be higher than the --barcode-mm value [default: same value "
          "as --barcode-mm]")
    .bind_uint(&barcode_mm_r1)
    .with_default(0);
  argparser.add("--barcode-mm-r2", "N")
    .help("Maximum number of mismatches allowed for the mate 2 barcode. "
          "Cannot be higher than the --barcode-mm value [default: same value "
          "as --barcode-mm]")
    .bind_uint(&barcode_mm_r2)
    .with_default(0);
  argparser.add("--demultiplex-only")
    .help("Only carry out demultiplexing using the list of barcodes "
          "supplied with --barcode-list. No other processing is done")
    .depends_on("--barcode-list")
    .conflicts_with("--identify-adapters")
    .conflicts_with("--report-only");

  argparser.add_header("REPORTS:");
  argparser.add("--report-only")
    .help("Write a report of the input data without performing any processing "
          "of the FASTQ reads. Report-only post-trimming/demultiplexing runs "
          "can be accomplished by setting --output options to /dev/null")
    .conflicts_with("--demultiplex-only")
    .conflicts_with("--identify-adapters");

  argparser.add("--report-sample-rate", "X")
    .help("Fraction of reads to use when generating base quality/composition "
          "curves for trimming reports. Using all data (--report-sample-nth "
          "1.0) results in an 10-30% decrease in throughput")
    .bind_double(&report_sample_rate)
    .with_default(0.1);
  argparser.add("--report-duplication", "N")
    .help("FastQC based duplicate detection, based on the frequency of the "
          "first N unique sequences observed. A value of 100,000 corresponds "
          "to FastQC defaults; a value of 0 disables the analysis")
    .bind_uint(&report_duplication)
    .with_default(0);

  argparser.add_header("LOGGING:");
  argparser.add("--log-level", "X")
    .help("The minimum severity of messagest to be written to stderr")
    .bind_str(&log_level)
    .with_choices({ "debug", "info", "warning", "error" })
    .with_default("info");

  argparser.add("--log-colors", "X")
    .help("Enable/disable the use of colors when writing log messages. If set "
          "to auto, colors will only be enabled if STDOUT is a terminal and "
          "the NO_COLORS is environmetal variable is not set")
    .bind_str(&log_color)
    .with_choices({ "auto", "always", "never" })
    .with_default("auto");
  argparser.add("--log-progress", "X")
    .help("Specify the type of progress reports used. If set to auto, then a "
          "spinner will be used if STDERR is a terminal and the NO_COLORS "
          "environmetal variable is not set, otherwise logging will be used")
    .bind_str(&log_progress_sink)
    .with_choices({ "auto", "log", "spin", "never" })
    .with_default("auto");

  argparser.add("--trimwindows", "X")
    .deprecated()
    .bind_double(&trim_window_length)
    .with_default(-1.0);
  argparser.add("--trimns").deprecated().bind_bool(&trim_ambiguous_bases);
  argparser.add("--trimqualities").deprecated().bind_bool(&trim_by_quality);
  argparser.add("--minquality", "N")
    .deprecated()
    .bind_uint(&low_quality_score)
    .with_default(2);

  argparser.add("--seed").deprecated().bind_uint(&m_deprecated_knobs);
  argparser.add("--collapse-deterministic").deprecated();
  argparser.add("--collapse-conservatively").deprecated();
  argparser.add("--qualitymax", "N")
    .deprecated()
    .bind_uint(&m_deprecated_knobs);
}

argparse::parse_result
userconfig::parse_args(int argc, char* argv[])
{
  if (argc <= 1) {
    argparser.print_help();
    return argparse::parse_result::error;
  }

  args = string_vec(argv, argv + argc);

  // ad-hoc arg parsing to make argparse output consistent with rest of run
  configure_log_colors(try_parse_argument(args, "--log-color", "auto"), true);
  configure_log_levels(try_parse_argument(args, "--log-level", "info"), true);

  const argparse::parse_result result = argparser.parse_args(argc, argv);
  if (result != argparse::parse_result::ok) {
    return result;
  }

  configure_log_colors(log_color);
  configure_log_levels(log_level);
  log_progress = configure_log_progress(log_progress_sink);
  io_encoding = configure_encoding(quality_input_base);

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

  if (argparser.is_set("--identify-adapters")) {
    // By default quality scores are ignored when inferring adapter
    // sequences. However, arguments are still checked above.
    io_encoding = FASTQ_ENCODING_SAM;
    run_type = ar_command::identify_adapters;
  } else if (argparser.is_set("--demultiplex-only")) {
    run_type = ar_command::demultiplex_sequences;
  } else if (argparser.is_set("--report-only")) {
    run_type = ar_command::report_only;
  }

  if (low_quality_score > static_cast<unsigned>(PHRED_SCORE_MAX)) {
    log::error() << "Invalid value for --minquality: " << low_quality_score
                 << "\n"
                 << "   must be in the range 0 .. " << PHRED_SCORE_MAX;
    return argparse::parse_result::error;
  } else if (trim_window_length >= 0) {
    trim_by_quality = true;
  } else if (argparser.is_set("--trimwindows")) {
    log::error() << "Invalid value for --trimwindows (" << trim_window_length
                 << "); value must be >= 0.";
    return argparse::parse_result::error;
  }

  // Check for invalid combinations of settings
  if (input_files_1.empty() && input_files_2.empty()) {
    log::error()
      << "No input files (--file1 / --file2) specified.\n"
      << "Please specify at least one input file using --file1 FILENAME.";

    return argparse::parse_result::error;
  } else if (!input_files_2.empty() &&
             (input_files_1.size() != input_files_2.size())) {
    log::error()
      << "Different number of files specified for --file1 and --file2.";

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
    merge |= argparser.is_set("--collapse-deterministic");
    merge |= argparser.is_set("--collapse-conservatively");
  } else {
    merge = false;
  }

  if (run_type == ar_command::identify_adapters && !paired_ended_mode) {
    log::error() << "Both input files (--file1 / --file2) must be "
                 << "specified when using --identify-adapters, or input must "
                 << "be interleaved FASTQ reads (requires --interleaved).";

    return argparse::parse_result::error;
  }

  // (Optionally) read adapters from file and validate
  if (!setup_adapter_sequences()) {
    return argparse::parse_result::error;
  }

  // Set mismatch threshold
  if (mismatch_threshold > 1) {
    mismatch_threshold = 1.0 / mismatch_threshold;
  } else if (mismatch_threshold < 0) {
    if (run_type == ar_command::identify_adapters) {
      mismatch_threshold = 1.0 / 10.0;
    } else {
      // Defaults for PE / SE trimming (changed in v2)
      mismatch_threshold = 1.0 / 3.0;
    }
  }

  // The actual max depends on the presence of libdeflate, but it does not seem
  // justified to blow up even if the user specifies a level above that
  // supported by isa-l.
  if (gzip_level > 9) {
    log::error() << "--gzip-level must be in the range 0 to 9, not "
                 << gzip_level;
    return argparse::parse_result::error;
  }

  if (!max_threads) {
    log::error() << "--threads must be at least 1!";
    return argparse::parse_result::error;
  }

  using fixed_trimming =
    std::tuple<const char*, const string_vec&, std::pair<unsigned, unsigned>&>;

  const std::vector<fixed_trimming> fixed_trimming_options = {
    { "--pre-trim5p", pre_trim5p, pre_trim_fixed_5p },
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

  if (trim_error_rate > 1) {
    trim_error_rate = std::pow(10.0, trim_error_rate / -10.0);
  }

  // Required since a missing filename part results in the creation of dot-files
  if (out_basename.empty()) {
    log::error() << "--basename must be a non-empty value.";

    return argparse::parse_result::error;
  } else if (out_basename.back() == '/') {
    log::error() << "--basename must not be a directory: "
                 << shell_escape(out_basename);

    return argparse::parse_result::error;
  }

  if (adapters.barcode_count()) {
    bool any_illegal_keys = false;
    for (const auto& key : output_files_set(argparser)) {
      if (argparser.to_str(key) != DEV_NULL) {
        log::error() << "Command-line option " << key << " can only be set to "
                     << "/dev/null when demultiplexing!";
        any_illegal_keys = true;
      }
    }

    if (any_illegal_keys) {
      return argparse::parse_result::error;
    }
  }

  if (run_type != ar_command::identify_adapters &&
      !argparser.is_set("--basename") && output_files_set(argparser).empty()) {
    log::error() << "No output files were specified. Use either --basename or "
                    "one or more --out-* arguments. If you actually don't want "
                    "any output, then use '--basename /dev/null'";

    return argparse::parse_result::error;
  }

  {
    const std::string key = "--pre-trim-polyx";
    if (argparser.is_set(key)) {
      if (!parse_poly_x_option(key, pre_trim_poly_x_sink, pre_trim_poly_x)) {
        return argparse::parse_result::error;
      }
    }
  }

  {
    const std::string key = "--post-trim-polyx";
    if (argparser.is_set(key)) {
      if (!parse_poly_x_option(key, post_trim_poly_x_sink, post_trim_poly_x)) {
        return argparse::parse_result::error;
      }
    }
  }

  return argparse::parse_result::ok;
}

bool
userconfig::is_good_alignment(const alignment_info& alignment) const
{
  if (!alignment.length || alignment.score <= 0) {
    return false;
  }

  // Only pairs of called bases are considered part of the alignment
  const size_t n_aligned =
    static_cast<size_t>(alignment.length - alignment.n_ambiguous);
  size_t mm_threshold = static_cast<size_t>(mismatch_threshold * n_aligned);

  if (n_aligned < min_adapter_overlap && !paired_ended_mode) {
    return false;
  }

  if (n_aligned < 6) {
    mm_threshold = 0;
  } else if (n_aligned < 10) {
    // Allow at most 1 mismatch, possibly set to 0 by the user
    mm_threshold = std::min<size_t>(1, mm_threshold);
  }

  if (alignment.n_mismatches > mm_threshold) {
    return false;
  }

  return true;
}

bool
userconfig::can_merge_alignment(const alignment_info& alignment) const
{
  if (alignment.length < alignment.n_ambiguous) {
    throw std::invalid_argument("#ambiguous bases > read length");
  }

  return alignment.length - alignment.n_ambiguous >= min_alignment_length;
}

output_files
userconfig::get_output_filenames() const
{
  output_files files;

  files.settings_json = new_filename("--out-json", ".json");
  files.settings_html = new_filename("--out-html", ".html");

  const std::string ext = gzip ? ".fastq.gz" : ".fastq";
  const std::string out1 = (interleaved_output ? "" : ".r1") + ext;
  const std::string out2 = (interleaved_output ? "" : ".r2") + ext;

  files.unidentified_1 = new_filename("--out-file1", ".unidentified" + out1);
  files.unidentified_2 = new_filename("--out-file2", ".unidentified" + out2);

  const bool demultiplexing = adapters.barcode_count();
  files.samples.resize(adapters.adapter_set_count());

  for (size_t i = 0; i < files.samples.size(); ++i) {
    const auto sample = demultiplexing ? adapters.get_sample_name(i) : "";
    auto& map = files.samples.at(i);

    const auto mate_1_filename = new_filename("--out-file1", sample, out1);
    map.set_filename(read_type::mate_1, mate_1_filename);

    if (paired_ended_mode) {
      if (interleaved_output) {
        map.set_filename(read_type::mate_2, mate_1_filename);
      } else {
        map.set_filename(read_type::mate_2,
                         new_filename("--out-file2", sample, out2));
      }
    }

    if (run_type == ar_command::trim_adapters) {
      map.set_filename(
        read_type::discarded,
        new_filename("--out-discarded", sample, ".discarded" + ext));

      if (paired_ended_mode) {
        map.set_filename(
          read_type::singleton,
          new_filename("--out-singleton", sample, ".singleton" + ext));

        if (merge) {
          map.set_filename(
            read_type::merged,
            new_filename("--out-merged", sample, ".merged" + ext));
        }
      }
    }
  }

  return files;
}

std::string
userconfig::new_filename(const std::string& key,
                         const std::string& first,
                         const std::string& second) const
{
  if (argparser.is_set(key)) {
    return argparser.to_str(key);
  } else if (!argparser.is_set("--basename")) {
    return DEV_NULL;
  } else if (out_basename == DEV_NULL) {
    // Special case to allow dry runs with no output
    return DEV_NULL;
  }

  std::string out = out_basename;
  if (first.size() && first.front() != '.') {
    out.push_back('.');
  }

  out.append(first);
  out.append(second);

  return out;
}

bool
check_and_set_barcode_mm(const argparse::parser& argparser,
                         const std::string& key,
                         unsigned barcode_mm,
                         unsigned& dst)
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
userconfig::is_low_quality_trimming_enabled() const
{
  return trim_error_rate > 0 || trim_window_length >= 0 ||
         trim_ambiguous_bases || trim_by_quality || trim_error_rate > 0;
}

bool
userconfig::is_terminal_base_pre_trimming_enabled() const
{
  return pre_trim_fixed_5p.first || pre_trim_fixed_5p.second ||
         pre_trim_fixed_3p.first || pre_trim_fixed_3p.second;
}

bool
userconfig::is_terminal_base_post_trimming_enabled() const
{
  return post_trim_fixed_5p.first || post_trim_fixed_5p.second ||
         post_trim_fixed_3p.first || post_trim_fixed_3p.second;
}

bool
userconfig::is_short_read_filtering_enabled() const
{
  return min_genomic_length > 0;
}

bool
userconfig::is_long_read_filtering_enabled() const
{
  return max_genomic_length != std::numeric_limits<unsigned>::max();
}

bool
userconfig::is_ambiguous_base_filtering_enabled() const
{
  return max_ambiguous_bases != std::numeric_limits<unsigned>::max();
}

bool
userconfig::is_low_complexity_filtering_enabled() const
{
  return min_complexity > 0;
}

bool
userconfig::setup_adapter_sequences()
{
  const bool adapters_is_set =
    argparser.is_set("--adapter1") || argparser.is_set("--adapter2");
  const bool adapter_list_is_set = argparser.is_set("--adapter-list");

  if (adapters_is_set && adapter_list_is_set) {
    log::error() << "Use either --adapter1 and --adapter2, or "
                 << "--adapter-list, not both!";

    return false;
  }

  if (adapter_list_is_set) {
    if (!adapters.load_adapters(adapter_list, paired_ended_mode)) {
      return false;
    } else if (adapters.adapter_count()) {
      log::info() << "Read " << adapters.adapter_count()
                  << " adapters / adapter pairs from '" << adapter_list << "'";
    } else {
      log::error() << "No adapter sequences found in table!";
      return false;
    }
  } else {
    try {
      adapters.add_adapters(adapter_1, adapter_2);
    } catch (const fastq_error& error) {
      log::error() << "Error parsing adapter sequence(s):\n"
                   << "   " << error.what();

      return false;
    }
  }

  if (!argparser.is_set("--barcode-mm")) {
    barcode_mm = barcode_mm_r1 + barcode_mm_r2;
  }

  if (!check_and_set_barcode_mm(
        argparser, "--barcode-mm-r1", barcode_mm, barcode_mm_r1)) {
    return false;
  }

  if (!check_and_set_barcode_mm(
        argparser, "--barcode-mm-r2", barcode_mm, barcode_mm_r2)) {
    return false;
  }

  if (argparser.is_set("--barcode-list")) {
    if (!adapters.load_barcodes(barcode_list, paired_ended_mode)) {
      return false;
    } else if (adapters.adapter_count()) {
      log::info() << "Read " << adapters.barcode_count()
                  << " sets of barcodes from " << shell_escape(barcode_list);
    } else {
      log::error() << "No barcodes sequences found in table!";
      return false;
    }
  }

  const auto& output_files = get_output_filenames();
  if (!check_input_and_output("--file1", input_files_1, output_files)) {
    return false;
  } else if (!check_input_and_output("--file2", input_files_2, output_files)) {
    return false;
  }

  return true;
}

double
userconfig::runtime() const
{
  return m_runtime.duration();
}

} // namespace adapterremoval
