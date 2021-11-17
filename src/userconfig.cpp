/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <cstring>   // for size_t, strerror
#include <errno.h>   // for errno
#include <iostream>  // for operator<<, basic_ostream, endl, cerr, ostream
#include <limits>    // for numeric_limits
#include <stdexcept> // for invalid_argument
#include <string>    // for string, operator<<, char_traits, operator==
#include <unistd.h>  // for access, R_OK

#include "alignment.hpp" // for alignment_info
#include "strutils.hpp"  // for template_replace, str_to_unsigned, toupper
#include "userconfig.hpp"

const size_t output_sample_files::disabled = std::numeric_limits<size_t>::max();

////////////////////////////////////////////////////////////////////////////////
// Helper functions

bool
select_encoding(const std::string& value,
                size_t quality_max,
                fastq_encoding& out)
{
  quality_encoding encoding;

  const std::string uppercase_value = toupper(value);
  if (uppercase_value == "33") {
    encoding = quality_encoding::phred_33;
  } else if (uppercase_value == "64") {
    encoding = quality_encoding::phred_64;
  } else if (uppercase_value == "SOLEXA") {
    encoding = quality_encoding::solexa;
  } else {
    std::cerr << "Error: Invalid value for --qualitybase: '" << value << "'\n"
              << "   expected values 33, 64, or solexa." << std::endl;

    return false;
  }

  out = fastq_encoding(encoding, quality_max);

  return true;
}

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
check_no_clobber(const std::string& label,
                 const string_vec& in_files,
                 const std::string& out_file)
{
  for (const auto& in_file : in_files) {
    if (in_file == out_file) {
      std::cerr << "ERROR: Input file would be overwritten: " << label << " "
                << in_file << std::endl;
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
      std::cerr << "ERROR: Cannot read file: " << label << " " << filename
                << "': " << std::strerror(errno) << std::endl;

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
    for (const auto& out_file : sample.filenames) {
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
  : settings()
  , unidentified_1()
  , unidentified_2()
  , samples()
{}

output_sample_files::output_sample_files()
  : filenames()
  , steps()
  , offsets()
{
  offsets.fill(output_sample_files::disabled);
}

size_t
output_sample_files::add(const std::string& filename)
{
  // If discarded then no post processing is needed; this saves time especially
  // when output compression is enabled.
  if (filename == "/dev/null") {
    return output_sample_files::disabled;
  }

  auto it = std::find(filenames.begin(), filenames.end(), filename);
  if (it == filenames.end()) {
    filenames.push_back(filename);

    return filenames.size() - 1;
  }

  return it - filenames.begin();
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
  , out_basename("your_output")
  , out_settings("{basename}{.sample}.json")
  , out_interleaved("{basename}{.sample}.fastq")
  , out_pair_1("{basename}{.sample}.r1.fastq")
  , out_pair_2("{basename}{.sample}.r2.fastq")
  , out_merged("{basename}{.sample}.merged.fastq")
  , out_discarded("{basename}{.sample}.discarded.fastq")
  // FIXME: Support both .1 and .2
  , out_singleton("{basename}{.sample}.singleton.fastq")
  , paired_ended_mode(false)
  , interleaved_input(false)
  , interleaved_output(false)
  , mate_separator(MATE_SEPARATOR)
  , min_genomic_length(15)
  , max_genomic_length(std::numeric_limits<unsigned>::max())
  , min_adapter_overlap(0)
  , min_alignment_length(11)
  , mismatch_threshold(-1.0)
  , io_encoding(FASTQ_ENCODING_33)
  , quality_max(MAX_PHRED_SCORE_DEFAULT)
  , trim_fixed_5p(0, 0)
  , trim_fixed_3p(0, 0)
  , trim_by_quality(false)
  , trim_window_length(std::numeric_limits<double>::quiet_NaN())
  , low_quality_score(2)
  , trim_ambiguous_bases(false)
  , max_ambiguous_bases(1000)
  , preserve5p(false)
  , merge(false)
  , merge_conservatively(false)
  , shift(2)
  , max_threads(1)
  , gzip(false)
  , gzip_stream(false)
  , gzip_level(6)
  , barcode_mm(0)
  , barcode_mm_r1(0)
  , barcode_mm_r2(0)
  , adapters()
  , report_sample_rate(0.1)
  , argparser(name, version, help)
  , adapter_1("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
  , adapter_2("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
  , adapter_list()
  , barcode_list()
  , quality_input_base("33")
  , mate_separator_str(1, MATE_SEPARATOR)
  , interleaved(false)
  , trim5p()
  , trim3p()
  , m_runtime()
  , m_deprecated_knobs()
  , m_deprecated_flags()
{
  argparser["--file1"] = new argparse::many(
    &input_files_1,
    "FILE [FILE ...]",
    "Input files containing mate 1 reads or single-ended reads; "
    "one or more files may be listed [REQUIRED].");
  argparser["--file2"] = new argparse::many(
    &input_files_2,
    "[FILE ...]",
    "Input files containing mate 2 reads; if used, then the same "
    "number of files as --file1 must be listed [OPTIONAL].");

  argparser["--identify-adapters"] = new argparse::flag(
    nullptr,
    "Attempt to identify the adapter pair of PE reads, by searching "
    "for overlapping mate reads [default: %default].");
  argparser["--threads"] = new argparse::knob(
    &max_threads, "THREADS", "Maximum number of threads [default: %default]");

  argparser.add_header("FASTQ OPTIONS:");
  argparser["--qualitybase"] = new argparse::any(
    &quality_input_base,
    "BASE",
    "Quality base used to encode Phred scores in input; either 33 or 64 "
    "[default: %default].");
  argparser["--qualitymax"] = new argparse::knob(
    &quality_max,
    "BASE",
    "Specifies the maximum Phred score expected in input files, and "
    "used when writing output. ASCII encoded values are limited to "
    "the characters '!' (ASCII = 33) to '~' (ASCII = 126), meaning "
    "that possible scores are 0 - 93 with offset 33, and 0 - 62 "
    "for offset 64 scores [default: %default].");
  argparser["--mate-separator"] = new argparse::any(
    &mate_separator_str,
    "CHAR",
    "Character separating the mate number (1 or 2) from the read name "
    "in FASTQ records [default: '%default'].");

  argparser["--interleaved"] = new argparse::flag(
    &interleaved,
    "This option enables both the --interleaved-input option and the "
    "--interleaved-output option [default: %default].");
  argparser["--interleaved-input"] = new argparse::flag(
    &interleaved_input,
    "The (single) input file provided contains both the mate 1 and "
    "mate 2 reads, one pair after the other, with one mate 1 reads "
    "followed by one mate 2 read. This option is implied by the "
    "--interleaved option [default: %default].");
  argparser["--interleaved-output"] = new argparse::flag(
    &interleaved_output,
    "If set, trimmed paired-end reads are written to a single file "
    "containing mate 1 and mate 2 reads, one pair after the other. "
    "This option is implied by the --interleaved option [default: "
    "%default].");

  argparser.add_header("OUTPUT FILES:");
  argparser["--basename"] = new argparse::any(
    &out_basename,
    "BASENAME",
    "Default prefix for all output files for which no filename was "
    "explicitly set [default: %default].");
  argparser["--settings"] = new argparse::any(
    &out_settings,
    "FILE",
    "Output file containing information on the parameters used in the "
    "run as well as overall statistics on the reads after trimming "
    "[default: %default]");
  argparser["--output1"] = new argparse::any(
    &out_pair_1,
    "FILE",
    "Output file containing trimmed mate1 reads [default: %default]");
  argparser["--output2"] = new argparse::any(
    &out_pair_2,
    "FILE",
    "Output file containing trimmed mate 2 reads [default: %default]");
  argparser["--singleton"] = new argparse::any(
    &out_singleton,
    "FILE",
    "Output file to which containing paired reads for which the mate "
    "has been discarded [default: %default]");
  argparser["--outputmerged"] = new argparse::any(
    &out_merged,
    "FILE",
    "If --merge is set, contains overlapping mate-pairs which "
    "have been merged into a single read (PE mode) or reads for which "
    "the adapter was identified by a minimum overlap, indicating that "
    "the entire template molecule is present. This does not include "
    "which have subsequently been trimmed due to low-quality or "
    "ambiguous nucleotides [default: %default]");
  argparser["--discarded"] = new argparse::any(
    &out_discarded,
    "FILE",
    "Contains reads discarded due to the --minlength, --maxlength or "
    "--maxns options [default: %default]");

  argparser.add_header("OUTPUT COMPRESSION:");
  argparser["--gzip"] =
    new argparse::flag(&gzip, "Enable gzip compression [default: %default]");
  argparser["--gzip-stream"] = new argparse::flag(
    &gzip_stream,
    "Compress output using GZip streams instead of compressing independent "
    "blocks of 64kb data. Block based compression "
#ifdef USE_LIBDEFLATE
    "uses libdeflate for greater throughput and "
#endif
    "allows even greater throughput in threaded mode at the cost of about 3% "
    "greater file size and possible incompatibility with a few programs. "
    "Implies --gzip [default: %default]");
  argparser["--gzip-level"] = new argparse::knob(
    &gzip_level, "LEVEL", "Compression level, 0 - 9 [default: %default]");

  argparser.add_header("TRIMMING SETTINGS:");
  argparser["--adapter1"] =
    new argparse::any(&adapter_1,
                      "SEQUENCE",
                      "Adapter sequence expected to be found in mate 1 reads "
                      "[default: %default].");
  argparser["--adapter2"] =
    new argparse::any(&adapter_2,
                      "SEQUENCE",
                      "Adapter sequence expected to be found in mate 2 reads "
                      "[default: %default].");
  argparser["--adapter-list"] = new argparse::any(
    &adapter_list,
    "FILENAME",
    "Read table of white-space separated adapters pairs, used as if "
    "the first column was supplied to --adapter1, and the second "
    "column was supplied to --adapter2; only the first adapter in "
    "each pair is required SE trimming mode [default: %default].");

  argparser.add_seperator();
  argparser["--minadapteroverlap"] = new argparse::knob(
    &min_adapter_overlap,
    "LENGTH",
    "In single-end mode, reads are only trimmed if the overlap "
    "between read and the adapter is at least X bases long, not "
    "counting ambiguous nucleotides (N); this is independent of the "
    "--minalignmentlength when using --merge, allowing a "
    "conservative selection of putative complete inserts while "
    "ensuring that all possible adapter contamination is trimmed "
    "[default: %default].");
  argparser["--mm"] = new argparse::floaty_knob(
    &mismatch_threshold,
    "MISMATCH_RATE",
    "Max error-rate when aligning reads and/or adapters. If > 1, the "
    "max error-rate is set to 1 / MISMATCH_RATE; if < 0, the defaults "
    "are used, otherwise the user-supplied value is used directly "
    "[default: 1/3 for trimming; 1/10 when identifying adapters].");
  argparser["--shift"] = new argparse::knob(
    &shift,
    "N",
    "Consider alignments where up to N nucleotides are missing from "
    "the 5' termini [default: %default].");

  argparser.add_seperator();
  argparser["--trim5p"] = new argparse::many(
    &trim5p,
    "N [N]",
    "Trim the 5' of reads by a fixed amount after removing adapters, "
    "but before carrying out quality based trimming. Specify one "
    "value to trim mate 1 and mate 2 reads the same amount, or two "
    "values separated by a space to trim each mate different amounts "
    "[default: no trimming].");
  argparser["--trim3p"] = new argparse::many(
    &trim3p, "N [N]", "Trim the 3' of reads by a fixed amount. See --trim5p.");

  argparser["--trimns"] =
    new argparse::flag(&trim_ambiguous_bases,
                       "If set, trim ambiguous bases (N) at 5'/3' termini "
                       "[default: %default]");
  argparser["--maxns"] = new argparse::knob(
    &max_ambiguous_bases,
    "MAX",
    "Reads containing more ambiguous bases (N) than this number after "
    "trimming are discarded [default: %default].");
  argparser["--trimqualities"] = new argparse::flag(
    &trim_by_quality,
    "If set, trim bases at 5'/3' termini with quality scores <= to "
    "--minquality value [default: %default]");
  argparser["--trimwindows"] = new argparse::floaty_knob(
    &trim_window_length,
    "INT",
    "If set, quality trimming will be carried out using window based "
    "approach, where windows with an average quality less than "
    "--minquality will be trimmed. If >= 1, this value will be used "
    "as the window size. If the value is < 1, the value will be "
    "multiplied with the read length to determine a window size per "
    "read. If the resulting window size is 0 or larger than the read "
    "length, the read length is used as the window size. This option "
    "implies --trimqualities [default: %default].");
  argparser["--minquality"] =
    new argparse::knob(&low_quality_score,
                       "PHRED",
                       "Inclusive minimum; see --trimqualities for details "
                       "[default: %default]");
  argparser["--preserve5p"] = new argparse::flag(
    &preserve5p,
    "If set, bases at the 5p will not be trimmed by --trimns, "
    "--trimqualities, and --trimwindows. Merged reads will "
    "not be quality trimmed when this option is enabled "
    "[default: 5p bases are trimmed]");

  argparser.add_seperator();
  argparser["--minlength"] =
    new argparse::knob(&min_genomic_length,
                       "LENGTH",
                       "Reads shorter than this length are discarded "
                       "following trimming [default: %default].");
  argparser["--maxlength"] =
    new argparse::knob(&max_genomic_length,
                       "LENGTH",
                       "Reads longer than this length are discarded "
                       "following trimming [default: %default].");

  argparser.add_header("READ MERGING:");
  argparser["--merge"] = new argparse::flag(
    &merge,
    "When set, paired ended read alignments of --minalignmentlength or more "
    "bases are merged into a single consensus sequence. Merged reads are "
    "written to basename.merged by default. Has no effect in single-end "
    "mode [default: %default].");
  argparser["--merge-conservatively"] = new argparse::flag(
    &merge_conservatively,
    "Enables a more conservative merging algorithm inspired by fastq-join, "
    "in which the higher quality score is picked for matching bases and the "
    "max score minus the min score is picked for mismatching bases. For more "
    "details, see the documentation. Setting this option also sets --merge "
    "[default: %default].");
  argparser["--minalignmentlength"] = new argparse::knob(
    &min_alignment_length,
    "LENGTH",
    "If --merge is set, paired reads must overlap at least this "
    "number of bases to be merged, and single-ended reads must "
    "overlap at least this number of bases with the adapter to be "
    "considered complete template molecules [default: %default].");
  argparser["--seed"] = new argparse::knob(&m_deprecated_knobs, "", "HIDDEN");

  argparser.add_header("DEMULTIPLEXING:");
  argparser["--barcode-list"] = new argparse::any(
    &barcode_list,
    "FILENAME",
    "List of barcodes or barcode pairs for single or double-indexed "
    "demultiplexing. Note that both indexes should be specified for "
    "both single-end and paired-end trimming, if double-indexed "
    "multiplexing was used, in order to ensure that the demultiplexed "
    "reads can be trimmed correctly [default: %default].");
  argparser["--barcode-mm"] = new argparse::knob(
    &barcode_mm,
    "N",
    "Maximum number of mismatches allowed when counting mismatches in "
    "both the mate 1 and the mate 2 barcode for paired reads.");
  argparser["--barcode-mm-r1"] = new argparse::knob(
    &barcode_mm_r1,
    "N",
    "Maximum number of mismatches allowed for the mate 1 barcode; "
    "if not set, this value is equal to the '--barcode-mm' value; "
    "cannot be higher than the '--barcode-mm value'.");
  argparser["--barcode-mm-r2"] = new argparse::knob(
    &barcode_mm_r2,
    "N",
    "Maximum number of mismatches allowed for the mate 2 barcode; "
    "if not set, this value is equal to the '--barcode-mm' value; "
    "cannot be higher than the '--barcode-mm value'.");
  argparser["--demultiplex-only"] = new argparse::flag(
    nullptr,
    "Only carry out demultiplexing using the list of barcodes "
    "supplied with --barcode-list. No other processing is done.");

  argparser.add_header("REPORTS:");
  argparser["--report-only"] = new argparse::flag(
    nullptr,
    "Write a report of the input data without performing any processing of the "
    "FASTQ reads. To generate a post-trimming/demultiplexing report without "
    "writing FASTQ files, set --output options to /dev/null.");

  argparser["--report-sample-rate"] = new argparse::floaty_knob(
    &report_sample_rate,
    "X",
    "Fraction of reads to use when generating base quality/composition curves "
    "for trimming reports. Using all data (--report-sample-nth 1.0) results in "
    "an about 10-30% decrease in throughput depending on settings [%default]");

  // Deprecated command-line options
  argparser["--collapse-deterministic"] =
    new argparse::flag(&m_deprecated_flags, "HIDDEN");

  // Aliases for backwards compatibility
  argparser.create_alias("--merge", "--collapse");
  argparser.create_alias("--merge-conservatively", "--collapse-conservatively");
  argparser.create_alias("--outputmerged", "--outputcollapsed");

  // Required options
  argparser.option_requires("--demultiplex-only", "--barcode-list");

  // Probibited combinations
  argparser.option_prohibits("--demultiplex-only", "--identify-adapters");
  argparser.option_prohibits("--demultiplex-only", "--report-only");
  argparser.option_prohibits("--identify-adapters", "--report-only");
  argparser.option_prohibits("--interleaved", "--file2");
  argparser.option_prohibits("--interleaved-input", "--file2");
}

argparse::parse_result
userconfig::parse_args(int argc, char* argv[])
{
  if (argc <= 1) {
    argparser.print_help();
    return argparse::parse_result::error;
  }

  args = string_vec(argv, argv + argc);
  const argparse::parse_result result = argparser.parse_args(argc, argv);
  if (result != argparse::parse_result::ok) {
    return result;
  }

  // --qualitybase is not always used (e.g. when identifying adapters), but is
  // always checked in order to catch invalid argument.
  if (!select_encoding(quality_input_base, quality_max, io_encoding)) {
    return argparse::parse_result::error;
  }

  if (mate_separator_str.size() != 1) {
    std::cerr << "Error: The argument for --mate-separator must be "
                 "exactly one character long, not "
              << mate_separator_str.size() << " characters!" << std::endl;
    return argparse::parse_result::error;
  } else {
    mate_separator = mate_separator_str.at(0);
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

  if (low_quality_score > static_cast<unsigned>(MAX_PHRED_SCORE)) {
    std::cerr << "Error: Invalid value for --minquality: " << low_quality_score
              << "\n"
              << "   must be in the range 0 .. " << MAX_PHRED_SCORE
              << std::endl;
    return argparse::parse_result::error;
  } else if (trim_window_length >= 0) {
    trim_by_quality = true;
  } else if (trim_window_length < 0.0) {
    std::cerr << "Error: Invalid value for --trimwindows ("
              << trim_window_length << "); value must be >= 0." << std::endl;
    return argparse::parse_result::error;
  }

  // Check for invalid combinations of settings
  if (input_files_1.empty() && input_files_2.empty()) {
    std::cerr
      << "Error: No input files (--file1 / --file2) specified.\n"
      << "Please specify at least one input file using --file1 FILENAME."
      << std::endl;

    return argparse::parse_result::error;
  } else if (!input_files_2.empty() &&
             (input_files_1.size() != input_files_2.size())) {
    std::cerr
      << "Error: Different number of files specified for --file1 and --file2."
      << std::endl;

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
    // --collapse-deterministic implies --merge
    merge |= argparser.is_set("--collapse-deterministic");
    // --merge-conservatively implies --merge
    merge |= merge_conservatively;
  } else {
    merge = false;
    merge_conservatively = false;
  }

  if (run_type == ar_command::identify_adapters && !paired_ended_mode) {
    std::cerr << "Error: Both input files (--file1 / --file2) must be "
              << "specified when using --identify-adapters, or input must "
              << "be interleaved FASTQ reads (requires --interleaved)."
              << std::endl;

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

  gzip |= gzip_stream;

  if (gzip_level > 9) {
    std::cerr << "Error: --gzip-level must be in the range 0 to 9, not "
              << gzip_level << std::endl;
    return argparse::parse_result::error;
  }

  if (!max_threads) {
    std::cerr << "Error: --threads must be at least 1!" << std::endl;
    return argparse::parse_result::error;
  }

  try {
    if (argparser.is_set("--trim5p")) {
      trim_fixed_5p = parse_trim_argument(trim5p);
    }
  } catch (const std::invalid_argument& error) {
    std::cerr << "Error: Could not parse --trim5p argument(s): " << error.what()
              << std::endl;

    return argparse::parse_result::error;
  }

  try {
    if (argparser.is_set("--trim3p")) {
      trim_fixed_3p = parse_trim_argument(trim3p);
    }
  } catch (const std::invalid_argument& error) {
    std::cerr << "Error: Could not parse --trim3p argument(s): " << error.what()
              << std::endl;

    return argparse::parse_result::error;
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

  const size_t n_aligned = alignment.length - alignment.n_ambiguous;
  if (n_aligned < min_alignment_length) {
    return false;
  }

  return merge || run_type == ar_command::identify_adapters;
}

output_files
userconfig::get_output_filenames() const
{
  output_files files;

  files.settings = template_replace(out_settings, "basename", out_basename);
  files.settings = template_replace(files.settings, "sample", "");

  const std::string out1 = interleaved_output ? out_interleaved : out_pair_1;
  const std::string out2 = interleaved_output ? out_interleaved : out_pair_2;

  files.unidentified_1 = get_output_filename("--output1", out1, "unidentified");
  files.unidentified_2 = get_output_filename("--output2", out2, "unidentified");

  const bool demultiplexing = adapters.barcode_count();

  for (size_t i = 0; i < adapters.adapter_set_count(); ++i) {
    const std::string name = demultiplexing ? adapters.get_sample_name(i) : "";

    files.samples.emplace_back();
    auto& map = files.samples.back();

    map.offset(read_type::mate_1) =
      map.add(get_output_filename("--output1", out1, name));

    if (paired_ended_mode) {
      if (interleaved_output) {
        map.offset(read_type::mate_2) = map.offset(read_type::mate_1);
      } else {
        map.offset(read_type::mate_2) =
          map.add(get_output_filename("--output2", out2, name));
      }
    }

    if (run_type == ar_command::trim_adapters) {
      map.offset(read_type::discarded_1) = map.offset(read_type::discarded_2) =
        map.add(get_output_filename("--discarded", out_discarded, name));

      if (paired_ended_mode) {
        map.offset(read_type::singleton_1) =
          map.offset(read_type::singleton_2) =
            map.add(get_output_filename("--singleton", out_singleton, name));

        if (merge) {
          map.offset(read_type::merged) =
            map.add(get_output_filename("--outputmerged", out_merged, name));
        }
      }
    }
  }

  return files;
}

std::string
userconfig::get_output_filename(const std::string& key,
                                const std::string& filename,
                                const std::string& sample) const
{
  auto tmp = template_replace(filename, "basename", out_basename);
  tmp = template_replace(tmp, "sample", sample);

  if (gzip && !argparser.is_set(key)) {
    tmp.append(".gz");
  }

  return tmp;
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
    std::cerr
      << "The maximum number of errors for " << key
      << " is set \n"
         "to a higher value than the total number of mismatches allowed\n"
         "for barcodes (--barcode-mm). Please correct these settings."
      << std::endl;
    return false;
  }

  return true;
}

bool
userconfig::setup_adapter_sequences()
{
  const bool adapters_is_set =
    argparser.is_set("--adapter1") || argparser.is_set("--adapter2");
  const bool adapter_list_is_set = argparser.is_set("--adapter-list");

  if (adapters_is_set && adapter_list_is_set) {
    std::cerr << "ERROR: "
              << "Use either --adapter1 and --adapter2, or "
              << "--adapter-list, not both!" << std::endl;

    return false;
  }

  if (adapter_list_is_set) {
    if (!adapters.load_adapters(adapter_list, paired_ended_mode)) {
      return false;
    } else if (adapters.adapter_count()) {
      std::cerr << "Read " << adapters.adapter_count()
                << " adapters / adapter pairs from '" << adapter_list << "'..."
                << std::endl;
    } else {
      std::cerr << "Error: No adapter sequences found in table!" << std::endl;
      return false;
    }
  } else {
    try {
      adapters.add_adapters(adapter_1, adapter_2);
    } catch (const fastq_error& error) {
      std::cerr << "Error parsing adapter sequence(s):\n"
                << "   " << error.what() << std::endl;

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
      std::cerr << "Read " << adapters.barcode_count()
                << " barcodes / barcode pairs from '" << barcode_list << "' ..."
                << std::endl;
    } else {
      std::cerr << "Error: No barcodes sequences found in table!" << std::endl;
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
