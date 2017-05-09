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
#include <iostream>
#include <cstring>
#include <cerrno>
#include <string>
#include <stdexcept>
#include <sys/time.h>
#include <limits>

#include "alignment.h"
#include "debug.h"
#include "fastq.h"
#include "strutils.h"
#include "userconfig.h"


namespace ar
{


size_t get_seed()
{
    struct timeval timestamp;
    gettimeofday(&timestamp, NULL);

    return (timestamp.tv_sec << 20) | timestamp.tv_usec;
}


fastq_encoding_ptr select_encoding(const std::string& name,
                                   const std::string& value,
                                   size_t quality_max = MAX_PHRED_SCORE_DEFAULT)
{
    fastq_encoding_ptr ptr;

    const std::string uppercase_value = toupper(value);
    if (uppercase_value == "33") {
        ptr.reset(new fastq_encoding(PHRED_OFFSET_33, quality_max));
    } else if (uppercase_value == "64") {
        ptr.reset(new fastq_encoding(PHRED_OFFSET_64, quality_max));
    } else if (uppercase_value == "SOLEXA") {
        ptr.reset(new fastq_encoding_solexa(quality_max));
    } else {
        std::cerr << "Error: Invalid value for " << name << ": '"
                  << value << "'\n" << "   expected values 33, 64, or solexa."
                  << std::endl;
    }

    return ptr;
}


userconfig::userconfig(const std::string& name,
                       const std::string& version,
                       const std::string& help)
    : run_type(ar_trim_adapters)
    , basename("your_output")
    , input_files_1()
    , input_files_2()
    , paired_ended_mode(false)
    , interleaved_input(false)
    , interleaved_output(false)
    , combined_output(false)
    , mate_separator(MATE_SEPARATOR)
    , min_genomic_length(15)
    , max_genomic_length(std::numeric_limits<unsigned>::max())
    , min_adapter_overlap(0)
    , min_alignment_length(11)
    , mismatch_threshold(-1.0)
    , quality_input_fmt()
    , quality_output_fmt()
    , trim_by_quality(false)
    , trim_window_length(std::numeric_limits<double>::quiet_NaN())
    , low_quality_score(2)
    , trim_ambiguous_bases(false)
    , max_ambiguous_bases(1000)
    , collapse(false)
    , shift(2)
    , seed(get_seed())
    , max_threads(1)
    , gzip(false)
    , gzip_level(6)
    , bzip2(false)
    , bzip2_level(9)
    , barcode_mm(0)
    , barcode_mm_r1(0)
    , barcode_mm_r2(0)
    , adapters()
    , argparser(name, version, help)
    , adapter_1("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
    , adapter_2("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")
    , adapter_list()
    , barcode_list()
    , quality_input_base("33")
    , quality_output_base("33")
    , quality_max(MAX_PHRED_SCORE_DEFAULT)
    , mate_separator_str(1, MATE_SEPARATOR)
    , interleaved(false)
    , identify_adapters(false)
    , demultiplex_sequences(false)
{
    argparser["--file1"] =
        new argparse::many(&input_files_1, "FILE [FILE ...]",
            "Input files containing mate 1 reads or single-ended reads; "
            "one or more files may be listed [REQUIRED].");
    argparser["--file2"] =
        new argparse::many(&input_files_2, "[FILE ...]",
            "Input files containing mate 2 reads; if used, then the same "
            "number of files as --file1 must be listed [OPTIONAL].");

    argparser.add_header("FASTQ OPTIONS:");
    argparser["--qualitybase"] =
        new argparse::any(&quality_input_base, "BASE",
            "Quality base used to encode Phred scores in input; either 33, "
            "64, or solexa [current: %default].");
    argparser["--qualitybase-output"] =
        new argparse::any(&quality_output_base, "BASE",
            "Quality base used to encode Phred scores in output; either 33, "
            "64, or solexa. By default, reads will be written in the same "
            "format as the that specified using --qualitybase.");
    argparser["--qualitymax"] =
        new argparse::knob(&quality_max, "BASE",
            "Specifies the maximum Phred score expected in input files, and "
            "used when writing output. ASCII encoded values are limited to "
            "the characters '!' (ASCII = 33) to '~' (ASCII = 126), meaning "
            "that possible scores are 0 - 93 with offset 33, and 0 - 62 "
            "for offset 64 and Solexa scores [default: %default].");
    argparser["--mate-separator"] =
        new argparse::any(&mate_separator_str, "CHAR",
            "Character separating the mate number (1 or 2) from the read name "
            "in FASTQ records [default: '%default'].");

    argparser["--interleaved"] =
        new argparse::flag(&interleaved,
            "This option enables both the --interleaved-input option and the "
            "--interleaved-output option [current: %default].");
    argparser["--interleaved-input"] =
        new argparse::flag(&interleaved_input,
            "The (single) input file provided contains both the mate 1 and "
            "mate 2 reads, one pair after the other, with one mate 1 reads "
            "followed by one mate 2 read. This option is implied by the "
            "--interleaved option [current: %default].");
    argparser["--interleaved-output"] =
        new argparse::flag(&interleaved_output,
            "If set, trimmed paired-end reads are written to a single file "
            "containing mate 1 and mate 2 reads, one pair after the other. "
            "This option is implied by the --interleaved option [current: "
            "%default].");
    argparser["--combined-output"] =
        new argparse::flag(&combined_output,
            "If set, all reads are written to the same file(s), specified by "
            "--output1 and --output2 (--output1 only if --interleaved-output "
            "is not set). Each read is further marked by either a \"PASSED\" "
            "or a \"FAILED\" flag, and any read that has been FAILED "
            "(including the mate for collapsed reads) are replaced with a "
            "single 'N' with Phred score 0 [current: %default].");

    argparser.add_header("OUTPUT FILES:");
    argparser["--basename"] =
        new argparse::any(&basename, "BASENAME",
            "Default prefix for all output files for which no filename was "
            "explicitly set [current: %default].");
    argparser["--settings"] =
        new argparse::any(NULL, "FILE",
            "Output file containing information on the parameters used in the "
            "run as well as overall statistics on the reads after trimming "
            "[default: BASENAME.settings]");
    argparser["--output1"] =
        new argparse::any(NULL, "FILE",
            "Output file containing trimmed mate1 reads [default: "
            "BASENAME.pair1.truncated (PE), BASENAME.truncated (SE), or "
            "BASENAME.paired.truncated (interleaved PE)]");
    argparser["--output2"] =
        new argparse::any(NULL, "FILE",
            "Output file containing trimmed mate 2 reads [default: "
            "BASENAME.pair2.truncated (only used in PE mode, but not if "
            "--interleaved-output is enabled)]");
    argparser["--singleton"] =
        new argparse::any(NULL, "FILE",
            "Output file to which containing paired reads for which the mate "
            "has been discarded [default: BASENAME.singleton.truncated]");
    argparser["--outputcollapsed"] =
        new argparse::any(NULL, "FILE",
            "If --collapsed is set, contains overlapping mate-pairs which "
            "have been merged into a single read (PE mode) or reads for which "
            "the adapter was identified by a minimum overlap, indicating that "
            "the entire template molecule is present. This does not include "
            "which have subsequently been trimmed due to low-quality or "
            "ambiguous nucleotides [default: BASENAME.collapsed]");
    argparser["--outputcollapsedtruncated"] =
        new argparse::any(NULL, "FILE",
            "Collapsed reads (see --outputcollapsed) which were trimmed due "
            "the presence of low-quality or ambiguous nucleotides "
            "[default: BASENAME.collapsed.truncated]");
    argparser["--discarded"] =
        new argparse::any(NULL, "FILE",
            "Contains reads discarded due to the --minlength, --maxlength or "
            "--maxns options [default: BASENAME.discarded]");

#if defined(AR_GZIP_SUPPORT) || defined(AR_BZIP2_SUPPORT)
   argparser.add_header("OUTPUT COMPRESSION:");
#endif

#ifdef AR_GZIP_SUPPORT
    argparser["--gzip"] =
        new argparse::flag(&gzip,
            "Enable gzip compression [current: %default]");
    argparser["--gzip-level"] =
        new argparse::knob(&gzip_level, "LEVEL",
            "Compression level, 0 - 9 [current: %default]");
#endif
#ifdef AR_BZIP2_SUPPORT
    argparser["--bzip2"] =
        new argparse::flag(&bzip2,
            "Enable bzip2 compression [current: %default]");
    argparser["--bzip2-level"] =
        new argparse::knob(&bzip2_level, "LEVEL",
            "Compression level, 0 - 9 [current: %default]");
#endif

    argparser.add_header("TRIMMING SETTINGS:");
    // Backwards compatibility with AdapterRemoval v1; not recommended due to
    // schematicts that differ from most other adapter trimming programs,
    // namely requiring that the --pcr2 sequence is that which is observed in
    // the reverse complement of mate 2, rather than in the raw reads.
    argparser["--pcr1"] = new argparse::any(&adapter_1, "SEQUENCE", "HIDDEN");
    argparser["--pcr2"] = new argparse::any(&adapter_2, "SEQUENCE", "HIDDEN");

    argparser["--adapter1"] =
        new argparse::any(&adapter_1, "SEQUENCE",
            "Adapter sequence expected to be found in mate 1 reads "
            "[current: %default].");
    argparser["--adapter2"] =
        new argparse::any(&adapter_2, "SEQUENCE",
            "Adapter sequence expected to be found in mate 2 reads "
            "[current: %default].");
    argparser["--adapter-list"] =
        new argparse::any(&adapter_list, "FILENAME",
            "Read table of white-space separated adapters pairs, used as if "
            "the first column was supplied to --adapter1, and the second "
            "column was supplied to --adapter2; only the first adapter in "
            "each pair is required SE trimming mode [current: %default].");

    argparser.add_seperator();
    argparser["--mm"]
        = new argparse::floaty_knob(&mismatch_threshold, "MISMATCH_RATE",
            "Max error-rate when aligning reads and/or adapters. If > 1, the "
            "max error-rate is set to 1 / MISMATCH_RATE; if < 0, the defaults "
            "are used, otherwise the user-supplied value is used directly "
            "[defaults: 1/3 for trimming; 1/10 when identifying adapters].");
    argparser["--maxns"] =
        new argparse::knob(&max_ambiguous_bases, "MAX",
            "Reads containing more ambiguous bases (N) than this number after "
            "trimming are discarded [current: %default].");
    argparser["--shift"] =
        new argparse::knob(&shift, "N",
            "Consider alignments where up to N nucleotides are missing from "
            "the 5' termini [current: %default].");

    argparser.add_seperator();
    argparser["--trimns"] =
        new argparse::flag(&trim_ambiguous_bases,
            "If set, trim ambiguous bases (N) at 5'/3' termini "
            "[current: %default]");
    argparser["--trimqualities"] =
        new argparse::flag(&trim_by_quality,
            "If set, trim bases at 5'/3' termini with quality scores <= to "
            "--minquality value [current: %default]");
    argparser["--trimwindows"] =
        new argparse::floaty_knob(&trim_window_length, "INT",
            "If set, quality trimming will be carried out using window based "
            "approach, where windows with an average quality less than "
            "--minquality will be trimmed. If >= 1, this value will be used "
            "as the window size. If the value is < 1, the value will be "
            "multiplied with the read length to determine a window size per "
            "read. If the resulting window size is 0 or larger than the read "
            "length, the read length is used as the window size. This option "
            "implies --trimqualities [current: %default].");
    argparser["--minquality"] =
        new argparse::knob(&low_quality_score, "PHRED",
            "Inclusive minimum; see --trimqualities for details "
            "[current: %default]");

    argparser["--minlength"] =
        new argparse::knob(&min_genomic_length, "LENGTH",
            "Reads shorter than this length are discarded "
            "following trimming [current: %default].");
    argparser["--maxlength"] =
        new argparse::knob(&max_genomic_length, "LENGTH",
            "Reads longer than this length are discarded "
            "following trimming [current: %default].");
    argparser["--collapse"] =
        new argparse::flag(&collapse,
            "When set, paired ended read alignments of --minalignmentlength "
            "or more bases are combined into a single consensus sequence, "
            "representing the complete insert, and written to either "
            "basename.collapsed or basename.collapsed.truncated (if trimmed "
            "due to low-quality bases following collapse); for single-ended "
            "reads, putative complete inserts are identified as having at "
            "least --minalignmentlength bases overlap with the adapter "
            "sequence, and are written to the the same files "
            "[current: %default].");
    argparser["--minalignmentlength"] =
        new argparse::knob(&min_alignment_length, "LENGTH",
            "If --collapse is set, paired reads must overlap at least this "
            "number of bases to be collapsed, and single-ended reads must "
            "overlap at least this number of bases with the adapter to be "
            "considered complete template molecules [current: %default].");
    argparser["--minadapteroverlap"] =
        new argparse::knob(&min_adapter_overlap, "LENGTH",
            "In single-end mode, reads are only trimmed if the overlap "
            "between read and the adapter is at least X bases long, not "
            "counting ambiguous nucleotides (N); this is independent of the "
            "--minalignmentlength when using --collapse, allowing a "
            "conservative selection of putative complete inserts while "
            "ensuring that all possible adapter contamination is trimmed "
            "[current: %default].");

    argparser.add_header("DEMULTIPLEXING:");
    argparser["--barcode-list"] =
        new argparse::any(&barcode_list, "FILENAME",
            "List of barcodes or barcode pairs for single or double-indexed "
            "demultiplexing. Note that both indexes should be specified for "
            "both single-end and paired-end trimming, if double-indexed "
            "multiplexing was used, in order to ensure that the demultiplexed "
            "reads can be trimmed correctly [current: %default].");
    argparser["--barcode-mm"] =
        new argparse::knob(&barcode_mm, "N",
            "Maximum number of mismatches allowed when counting mismatches in "
            "both the mate 1 and the mate 2 barcode for paired reads.");
    argparser["--barcode-mm-r1"] =
        new argparse::knob(&barcode_mm_r1, "N",
            "Maximum number of mismatches allowed for the mate 1 barcode; "
            "if not set, this value is equal to the '--barcode-mm' value; "
            "cannot be higher than the '--barcode-mm value'.");
    argparser["--barcode-mm-r2"] =
        new argparse::knob(&barcode_mm_r2, "N",
            "Maximum number of mismatches allowed for the mate 2 barcode; "
            "if not set, this value is equal to the '--barcode-mm' value; "
            "cannot be higher than the '--barcode-mm value'.");
    argparser["--demultiplex-only"] =
        new argparse::flag(&demultiplex_sequences,
            "Only carry out demultiplexing using the list of barcodes "
            "supplied with --barcode-list; do not attempt to trim adapters "
            "to carry out other processing.");

    argparser.add_header("MISC:");
    argparser["--identify-adapters"] =
        new argparse::flag(&identify_adapters,
            "Attempt to identify the adapter pair of PE reads, by searching "
            "for overlapping reads [current: %default].");
    argparser["--seed"] =
        new argparse::knob(&seed, "SEED",
            "Sets the RNG seed used when choosing between bases with equal "
            "Phred scores when collapsing. Note that runs are not "
            "deterministic if more than one thread is used. If not specified, "
            "a seed is generated using the current time.");
    argparser["--threads"] =
        new argparse::knob(&max_threads, "THREADS",
            "Maximum number of threads [current: %default]");
}


argparse::parse_result userconfig::parse_args(int argc, char *argv[])
{
    if (argc <= 1) {
        argparser.print_help();
        return argparse::pr_error;
    }

    const argparse::parse_result result = argparser.parse_args(argc, argv);
    if (result != argparse::pr_ok) {
        return result;
    }

    quality_input_fmt = select_encoding("--qualitybase", quality_input_base, quality_max);
    if (!quality_input_fmt.get()) {
        return argparse::pr_error;
    }

    if (argparser.is_set("--qualitybase-output")) {
        quality_output_fmt = select_encoding("--qualitybase-out", quality_output_base, quality_max);
        if (!quality_output_fmt.get()) {
            return argparse::pr_error;
        }
    } else {
        // Default to using the same output encoding as the input
        quality_output_fmt = select_encoding("--qualitybase", quality_input_base, quality_max);
        if (!quality_output_fmt.get()) {
            return argparse::pr_error;
        }
    }

    if (mate_separator_str.size() != 1) {
        std::cerr << "Error: The argument for --mate-separator must be "
                     "exactly one character long, not "
                     << mate_separator_str.size()
                     << " characters!" << std::endl;
        return argparse::pr_error;
    } else {
        mate_separator = mate_separator_str.at(0);
    }

    // Previous settings are overwritten; this ensures that bad arguments are still caught
    if (identify_adapters) {
        // By default quality scores are ignored when inferring adapter
        // sequences. However, arguments are still checked above.
        quality_input_fmt.reset(new fastq_encoding(PHRED_OFFSET_33, MAX_PHRED_SCORE));
        quality_output_fmt.reset(new fastq_encoding(PHRED_OFFSET_33, MAX_PHRED_SCORE));
        run_type = ar_identify_adapters;
    }

    if (demultiplex_sequences) {
        if (identify_adapters) {
            std::cerr << "Error: Cannot use --identify-adapters and "
                      << "--demultiplex-only at the same time!"
                      << std::endl;

            return argparse::pr_error;
        } else if (!argparser.is_set("--barcode-list")) {
            std::cerr << "Error: Cannot use --demultiplex-only without specifying "
                      << "a list of barcodes using --barcode-list!"
                      << std::endl;

            return argparse::pr_error;
        }

        run_type = ar_demultiplex_sequences;
    }

    if (low_quality_score > static_cast<unsigned>(MAX_PHRED_SCORE)) {
        std::cerr << "Error: Invalid value for --minquality: "
                  << low_quality_score << "\n"
                  << "   must be in the range 0 .. " << MAX_PHRED_SCORE
                  << std::endl;
        return argparse::pr_error;
    } else if (trim_window_length >= 0) {
        trim_by_quality = true;
    } else if (trim_window_length < 0.0) {
        std::cerr << "Error: Invalid value for --trimwindows ("
                  << trim_window_length << "); value must be >= 0."
                  << std::endl;
        return argparse::pr_error;
    }

    // Check for invalid combinations of settings
    if (input_files_1.empty() && input_files_2.empty()) {
        std::cerr << "Error: No input files (--file1 / --file2) specified.\n"
                  << "Please specify at least one input file using --file1 FILENAME."
                  << std::endl;

        return argparse::pr_error;
    } else if (!input_files_2.empty() && (input_files_1.size() != input_files_2.size())) {
        std::cerr << "Error: Different number of files specified for --file1 and --file2." << std::endl;

        return argparse::pr_error;
    } else if (!input_files_2.empty()) {
        paired_ended_mode = true;
        min_adapter_overlap = 0;
    }

    interleaved_input |= interleaved;
    interleaved_output |= interleaved;

    if (interleaved_input) {
        if (!input_files_2.empty()) {
            std::cerr << "Error: The options --interleaved and "
                      << "--interleaved-input cannot be used "
                      << "together with the --file2 option; only --file1 must "
                      << "be specified!"
                      << std::endl;
            return argparse::pr_error;
        }

        // Enable paired end mode .. other than the FASTQ reader, all other
        // parts of the pipeline simply run in paired-end mode.
        paired_ended_mode = true;
    }

    if (identify_adapters && !paired_ended_mode) {
        std::cerr << "Error: Both input files (--file1 / --file2) must be "
                  << "specified when using --identify-adapters, or input must "
                  << "be interleaved FASTQ reads (requires --interleaved)."
                  << std::endl;

        return argparse::pr_error;
    }

    // (Optionally) read adapters from file and validate
    if (!setup_adapter_sequences()) {
        return argparse::pr_error;
    }

    // Set mismatch threshold
    if (mismatch_threshold > 1) {
        mismatch_threshold = 1.0 / mismatch_threshold;
    } else if (mismatch_threshold < 0) {
        if (identify_adapters) {
            mismatch_threshold = 1.0 / 10.0;
        } else {
            // Defaults for PE / SE trimming (changed in v2)
            mismatch_threshold = 1.0 / 3.0;
        }
    }

    if (gzip_level > 9) {
        std::cerr << "Error: --gzip-level must be in the range 0 to 9, not "
                  << gzip_level << std::endl;
        return argparse::pr_error;

    }

#ifdef AR_BZIP2_SUPPORT
    if (bzip2_level < 1 || bzip2_level > 9) {
        std::cerr << "Error: --bzip2-level must be in the range 1 to 9, not "
                  << bzip2_level << std::endl;
        return argparse::pr_error;
    } else if (bzip2 && gzip) {
        std::cerr << "Error: Cannot enable --gzip and --bzip2 at the same time!"
                  << std::endl;
        return argparse::pr_error;
    }
#endif

    if (!max_threads) {
        std::cerr << "Error: --threads must be at least 1!" << std::endl;
        return argparse::pr_error;
    } else if (max_threads > 1 && argparser.is_set("--seed")) {
        std::cerr << "Warning: The option --seed should not be used when "
                  << "using multiple threads; multi-threaded behavior is not "
                  << "deterministic even with a fixed seed!"
                  << std::endl;
    }

    return argparse::pr_ok;
}


statistics_ptr userconfig::create_stats() const
{
    statistics_ptr stats(new statistics());
    stats->number_of_reads_with_adapter.resize(adapters.adapter_count());
    return stats;
}


bool userconfig::is_good_alignment(const alignment_info& alignment) const
{
    if (!alignment.length || alignment.score <= 0) {
        return false;
    }

    // Only pairs of called bases are considered part of the alignment
    const size_t n_aligned = static_cast<size_t>(alignment.length - alignment.n_ambiguous);
    size_t mm_threshold = static_cast<size_t>(mismatch_threshold * n_aligned);

    if (n_aligned < min_adapter_overlap) {
        return false;
    }

    if (n_aligned < 6) {
        mm_threshold = 0;
    } else if (n_aligned < 10) {
        // --mm may imply fewer allowed mismatches than 1, so always compare
        mm_threshold = std::min<size_t>(1, mm_threshold);
    }

    if (alignment.n_mismatches > mm_threshold) {
        return false;
    }

    return true;
}


bool userconfig::is_alignment_collapsible(const alignment_info& alignment) const
{
    if (alignment.length < alignment.n_ambiguous) {
        throw std::invalid_argument("#ambiguous bases > read length");
    }

    const size_t n_aligned = alignment.length - alignment.n_ambiguous;
    if (n_aligned < min_alignment_length) {
        return false;
    }

    return collapse || identify_adapters;
}


bool userconfig::is_acceptable_read(const fastq& seq) const
{
    const size_t seq_len = seq.length();

    return seq_len >= min_genomic_length
        && seq_len <= max_genomic_length
        && (max_ambiguous_bases >= seq_len
            || seq.count_ns() <= max_ambiguous_bases);
}


std::string userconfig::get_output_filename(const std::string& key,
                                            size_t nth) const
{
    std::string filename = basename;

    if (key == "demux_stats") {
        return filename += ".settings";
    } else if (key == "demux_unknown") {
        filename += ".unidentified";

        AR_DEBUG_ASSERT(nth <= 9);
        if (!interleaved_output) {
            if (nth) {
                filename.push_back('_');
                filename.push_back('0' + nth);
            }
        } else {
            filename += ".paired";
        }

        // Currently only when demultiplexing; for backwards compatibility
        if (run_type == ar_demultiplex_sequences) {
            filename += ".fastq";
        }

        if (gzip) {
            filename += ".gz";
        } else if (bzip2) {
            filename += ".bz2";
        }

        return filename;
    } else if (adapters.barcode_count()) {
        filename.push_back('.');
        filename.append(adapters.get_sample_name(nth));
    } else if (argparser.is_set(key)) {
        return argparser.at(key)->to_str();
    }

    if (key == "--settings") {
        return filename + ".settings";
    } else if (key == "--outputcollapsed") {
        filename += ".collapsed";
    } else if (key == "--outputcollapsedtruncated") {
        filename += ".collapsed.truncated";
    } else if (key == "--discarded") {
        filename += ".discarded";
    } else if (paired_ended_mode) {
        if (key == "--output1") {
            if (interleaved_output) {
                filename += ".paired";
            } else {
                filename += ".pair1";
            }
        } else if (key == "--output2") {
            filename += ".pair2";
        } else if (key == "--singleton") {
            filename += ".singleton";
        } else {
            throw std::invalid_argument("invalid read-type in userconfig::get_output_filename constructor: " + key);
        }

        if (run_type != ar_demultiplex_sequences) {
            filename += ".truncated";
        }
    } else if (key == "--output1") {
        if (run_type != ar_demultiplex_sequences) {
            filename += ".truncated";
        }
    } else if (key != "demux_unknown") {
        throw std::invalid_argument("invalid read-type in userconfig::get_output_filename constructor: " + key);
    }

    // Currently only when demultiplexing; for backwards compatibility
    if (run_type == ar_demultiplex_sequences) {
        filename += ".fastq";
    }

    if (gzip) {
        filename += ".gz";
    } else if (bzip2) {
        filename += ".bz2";
    }

    return filename;
}


fastq::ntrimmed userconfig::trim_sequence_by_quality_if_enabled(fastq& read) const
{
    if (trim_window_length >= 0) {
        return read.trim_windowed_bases(trim_ambiguous_bases,
                                        low_quality_score,
                                        trim_window_length);
    } else if (trim_ambiguous_bases || trim_by_quality) {
        const char quality_score = trim_by_quality ? low_quality_score : -1;

        return read.trim_trailing_bases(trim_ambiguous_bases,
                                        quality_score);
    }

    return fastq::ntrimmed();
}


bool check_and_set_barcode_mm(const argparse::parser& argparser,
                              const std::string& key,
                              unsigned barcode_mm,
                              unsigned& dst)
{
    if (!argparser.is_set(key)) {
        dst = barcode_mm;
    } else if (dst > barcode_mm) {
        std::cerr << "The maximum number of errors for " << key << " is set \n"
                     "to a higher value than the total number of mismatches allowed\n"
                     "for barcodes (--barcode-mm). Please correct these settings."
                  << std::endl;
        return false;
    }

    return true;
}


bool userconfig::setup_adapter_sequences()
{
    const bool pcr_is_set
        = argparser.is_set("--pcr1") || argparser.is_set("--pcr2");
    const bool adapters_is_set
        = argparser.is_set("--adapter1") || argparser.is_set("--adapter2");
    const bool adapter_list_is_set = argparser.is_set("--adapter-list");

    if (pcr_is_set) {
        std::cerr << "WARNING: Command-line options --pcr1 and --pcr2 are deprecated.\n"
                  << "         Using --adapter1 and --adapter2 is recommended.\n"
                  << "         Please see documentation for more information.\n"
                  << std::endl;
    }

    if (pcr_is_set && (adapters_is_set || adapter_list_is_set)) {
        std::cerr << "ERROR: "
                  << "Either use --pcr1 and --pcr2, or use --adapter1 and "
                  << "--adapter2 / --adapter-list, not both!\n\n"
                  << std::endl;

        return false;
    } else if (adapters_is_set && adapter_list_is_set) {
        std::cerr << "ERROR: "
                  << "Use either --adapter1 and --adapter2, or "
                  << "--adapter-list, not both!"
                  << std::endl;

        return false;
    }

    if (adapter_list_is_set) {
        if (!adapters.load_adapters(adapter_list, paired_ended_mode)) {
            return false;
        } else if (adapters.adapter_count()) {
            std::cerr << "Read " << adapters.adapter_count()
                      << " adapters / adapter pairs from '" << adapter_list
                      << "'..." << std::endl;
        } else {
            std::cerr << "Error: No adapter sequences found in table!" << std::endl;
            return false;
        }
    } else {
        try {
            adapters.add_adapters(adapter_1, adapter_2, !pcr_is_set);
        } catch (const fastq_error& error) {
            std::cerr << "Error parsing adapter sequence(s):\n"
                      << "   " << error.what() << std::endl;

            return false;
        }
    }

    if (!argparser.is_set("--barcode-mm")) {
        barcode_mm = barcode_mm_r1 + barcode_mm_r2;
    }

    if (!check_and_set_barcode_mm(argparser, "--barcode-mm-r1", barcode_mm, barcode_mm_r1)) {
        return false;
    }

    if (!check_and_set_barcode_mm(argparser, "--barcode-mm-r2", barcode_mm, barcode_mm_r2)) {
        return false;
    }

    if (argparser.is_set("--barcode-list")) {
        if (!adapters.load_barcodes(barcode_list, paired_ended_mode)) {
            return false;
        } else if (adapters.adapter_count()) {
            std::cerr << "Read " << adapters.barcode_count()
                      << " barcodes / barcode pairs from '" << barcode_list
                      << "' ..." << std::endl;
        } else {
            std::cerr << "Error: No barcodes sequences found in table!" << std::endl;
            return false;
        }
    }

    return true;
}

} // namespace ar
