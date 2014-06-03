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
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <string>

#include "userconfig.h"
#include "fastq.h"
#include "alignment.h"


bool cleanup_and_validate_sequence(std::string& sequence,
                                   const std::string& desc)
{
    try {
        fastq::clean_sequence(sequence);
    } catch (const fastq_error&) {
        std::cerr << "Error: Invalid nucleotide sequence supplied to " << desc
                  << ": '" << sequence << "'" << std::endl;
        return false;
    }

    return true;
}


fastq read_adapter_sequence(std::stringstream& instream)
{
    std::string adapter_seq;
    instream >> adapter_seq;

    if (instream.eof() && adapter_seq.empty()) {
        throw fastq_error("Adapter table contains partial record");
    } else if (instream.fail()) {
        throw fastq_error("IO error reading adapter-table");
    } else if (adapter_seq.empty()) {
        throw fastq_error("Adapter table contains empty entries");
    }

    fastq::clean_sequence(adapter_seq);
    return fastq("adapter", adapter_seq,
                 std::string(adapter_seq.length(), 'I'));
}


userconfig::userconfig(const std::string& name,
                             const std::string& version,
                             const std::string& help)
    : argparser(name, version, help)
    , basename("output")
    , input_file_1()
    , input_file_2()
    , paired_ended_mode(false)
    , adapters()
    , trim_barcodes_mode(false)
    , barcodes()
    , min_genomic_length(15)
    , min_alignment_length(11) // The minimum required genomic overlap before collapsing reads into one
    , mismatch_threshold(-1.0)
    , quality_input_fmt(phred_33)
    , quality_output_fmt(phred_33)
    , trim_by_quality(false)
    , low_quality_score(2)
    , trim_ambiguous_bases(false)
    , max_ambiguous_bases(1000)
    , collapse(false)
    , shift(2)
    , seed(static_cast<size_t>(time(NULL)))
    , identify_adapters(false)
    , PCR1("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
    , PCR2("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
    , barcode()
    , quality_input_base("33")
    , quality_output_base("33")
{
    argparser["--file1"] = new argparse::any(&input_file_1, "FILE",
        "Input file containing mate 1 reads or single-ended reads [REQUIRED].");
    argparser["--file2"] = new argparse::any(&input_file_2, "FILE",
        "Input file containing mate 2 reads [OPTIONAL].");

    // Output files
    argparser.add_seperator();
    argparser["--basename"] = new argparse::any(&basename, "BASENAME",
        "Default prefix for all output files for which no filename was explicitly set [current: %default].");
    argparser["--settings"] = new argparse::any(NULL, "FILE", "BASENAME.settings");
    argparser["--output1"] = new argparse::any(NULL, "FILE",
        "BASENAME.pair1.truncated (PE) or BASENAME.truncated (SE)");
    argparser["--output2"] = new argparse::any(NULL, "FILE",
        "BASENAME.pair2.truncated (only used in PE mode).");
    argparser["--singleton"] = new argparse::any(NULL, "FILE", "BASENAME.singleton.truncated");
    argparser["--outputcollapsed"] = new argparse::any(NULL, "FILE", "BASENAME.collapsed");
    argparser["--outputcollapsedtruncated"] = new argparse::any(NULL, "FILE", "BASENAME.collapsed.truncated");
    argparser["--discarded"] = new argparse::any(NULL, "FILE", "BASENAME.discarded");

    argparser.add_seperator();
    argparser["--pcr1"] = new argparse::any(&PCR1, "SEQUENCE",
        "Adapter sequence expected to be found in mate 1 reads [current: %default].");
    argparser["--pcr2"] = new argparse::any(&PCR2, "SEQUENCE",
        "Adapter sequence expected to be found in reverse complemented mate 2 reads [current: %default].");
    argparser["--pcr-list"] = new argparse::any(&PCR_list, "FILENAME",
        "List of adapters pairs, used as if supplied to --pcr1 / --pcr2; only the first adapter "
        "in each pair is required / used in SE mode [current: none].");
    argparser["--mm"] = new argparse::floaty_knob(&mismatch_threshold, "MISMATCH_RATE",
        "Max error-rate when aligning reads and/or adapters "
        "[default: 1/3 for single-ended; 1/10 for paired-ended; 0 when identifing adapters].");
    argparser["--maxns"] = new argparse::knob(&max_ambiguous_bases, "MAX",
        "Reads containing more ambiguous bases (N) than this number after trimming are discarded [current: %default].");
    argparser["--shift"] = new argparse::knob(&shift, "N",
        "Consider alignments where up to N nucleotides are missing from the 5' termini [current: %default].");

    argparser.add_seperator();
    argparser["--qualitybase"] = new argparse::any(&quality_input_base, "BASE",
        "Quality base used to encode Phred scores in input; either 33, 64, or solexa [current: %default].");
    argparser["--qualitybase-output"] = new argparse::any(&quality_output_base, "BASE",
        "Quality base used to encode Phred scores in output; either 33, 64 [current: %default].");
    argparser["--5prime"] = new argparse::any(&barcode, "BARCODE",
        "If set, the NT barcode is detected (max 1 mismatch) in and trimmed from the 5' of mate 1 reads [current: %default].");
    argparser["--5prime-list"] = new argparse::any(&barcode_list, "FILENAME",
        "List of barcode sequences, used as if supplied to --5prime; is "
        "parsed as SE adapters supplied to --pcr-list. [current: none].");
    argparser["--trimns"] = new argparse::flag(&trim_ambiguous_bases,
        "If set, trim ambiguous bases (N) at 5'/3' termini [current: %default]");
    argparser["--trimqualities"] = new argparse::flag(&trim_by_quality,
        "If set, trim bases at 5'/3' termini with quality scores <= to --minquality value [current: %default]");
    argparser["--minquality"] = new argparse::knob(&low_quality_score, "PHRED",
        "Inclusive minimum; see --trimqualities for details [current: %default]");
    argparser["--minlength"] = new argparse::knob(&min_genomic_length, "LENGTH",
        "Reads shorter than this length are written to BASENAME.discarded following trimming [current: %default].");

    argparser["--minalignmentlength"] = new argparse::knob(&min_alignment_length, "LENGTH",
        "If --collapse is set, reads must overlap at least this number of bases to be collapsed [current: %default].");
    argparser["--collapse"] = new argparse::flag(&collapse,
        "If set, paired ended reads which overlapp at least --minalignmentlength bases are combined into a single consensus read [current: %default].");

    argparser.add_seperator();
    argparser["--identify-adapters"] = new argparse::flag(&identify_adapters,
        "Attempt to identify the adapter pair of PE reads, by searching for overlapping reads [current: %default].");
    argparser["--seed"] = new argparse::knob(&seed, "SEED",
        "Sets the RNG seed used when choosing between bases with equal Phred scores when collapsing [current: %default].");
#warning FIXME: argparser["--debug"] = new argparse::flag(&DEBUG);
}


bool userconfig::parse_args(int argc, char *argv[])
{
    if (argc <= 1) {
        argparser.print_help();
        return false;
    } else if (!argparser.parse_args(argc, argv)) {
        return false;
    }

    if (quality_input_base == "33") {
        quality_input_fmt = phred_33;
    } else if (quality_input_base == "64") {
        quality_input_fmt = phred_64;
    } else if (quality_input_base == "solexa") {
        quality_input_fmt = solexa;
    } else {
        std::cerr << "Error: Invalid value for --qualitybase: '"
                  << quality_input_base << "'\n"
                  << "   expected values 33, 64, or solexa." << std::endl;
        return false;
    }

    if (quality_output_base == "33") {
        quality_output_fmt = phred_33;
    } else if (quality_output_base == "64") {
        quality_output_fmt = phred_64;
    } else {
        std::cerr << "Error: Invalid value for --qualitybase-out: '"
                  << quality_output_base << "'\n"
                  << "   expected values 33 or 64." << std::endl;
        return false;
    }

    if (low_quality_score > MAX_PHRED_SCORE) {
        std::cerr << "Error: Invalid value for --minquality: "
                  << low_quality_score << "\n"
                  << "   must be in the range 0 .. " << MAX_PHRED_SCORE
                  << std::endl;
        return false;
    }

    if (!cleanup_and_validate_sequence(PCR1, "--pcr1")) {
        return false;
    } else if (!cleanup_and_validate_sequence(PCR2, "--pcr2")) {
        return false;
    } else if (!cleanup_and_validate_sequence(barcode, "--5prime")) {
        return false;
    }

    // Check for invalid combinations of settings
    const bool file_1_set = argparser.is_set("--file1");
    const bool file_2_set = argparser.is_set("--file2");

    if (collapse && !file_2_set) {
        std::cerr << "You tried to collapse a single file. Collapse has been reset to false." << std::endl;
        collapse = false;
    } else if (!(file_1_set || file_2_set)) {
        std::cerr << "Error: No input files (--file1 / --file2) specified.\n"
                  << "Please specify at least one input file using --file1 FILENAME.\n"
                  << std::endl;
        argparser.print_help();
        return false;
    } else if (file_2_set && !file_2_set) {
        std::cerr << "Error: --file2 specified, but --file1 is not specified\n" << std::endl;
        argparser.print_help();
        return false;
    } else if ((argparser.is_set("--pcr1") || argparser.is_set("--pcr2"))
               && argparser.is_set("--pcr-list")) {
        std::cerr << "Error: Use either --pcr1 / --pcr2, or --pcr-list, not both!\n" << std::endl;
        return false;
    } else if (argparser.is_set("--5prime") && argparser.is_set("--5prime-list")) {
        std::cerr << "Error: Use either --5prime or --5prime-list, not both!\n" << std::endl;
        return false;
    }

    paired_ended_mode = file_2_set;

    if (argparser.is_set("--5prime")) {
        trim_barcodes_mode = true;
        barcodes.push_back(fastq_pair(fastq("Barcode", barcode, std::string(barcode.length(), 'J')), fastq()));
    } else if (argparser.is_set("--5prime-list")) {
        trim_barcodes_mode = true;
        if (!read_adapters_sequences(barcode_list, barcodes, paired_ended_mode)) {
            return false;
        } else if (adapters.empty()) {
            std::cerr << "Error: No barcode sequences found in table!" << std::endl;
            return false;
        }
    }

    if (argparser.is_set("--pcr-list")) {
        if (!read_adapters_sequences(PCR_list, adapters, paired_ended_mode)) {
            return false;
        } else if (adapters.empty()) {
            std::cerr << "Error: No adapter sequences found in table!" << std::endl;
            return false;
        }
    } else {
        adapters.push_back(fastq_pair(fastq("PCR1", PCR1, std::string(PCR1.length(), 'J')),
                                      fastq("PCR2", PCR2, std::string(PCR2.length(), 'J'))));
    }

    // Set mismatch threshold
    if (mismatch_threshold > 1) {
        mismatch_threshold = 1.0 / mismatch_threshold;
    } else if (mismatch_threshold < 0) {
        // Default values
        if (identify_adapters) {
            mismatch_threshold = 0.0;
        } else if (paired_ended_mode) {
            mismatch_threshold = 1.0 / 10.0;
        } else {
            mismatch_threshold = 1.0 / 3.0;
        }
    }

    // Set seed for RNG; rand is used in collapse_paired_ended_sequences()
    srandom(seed);

    return true;
}


statistics userconfig::create_stats() const
{
    statistics stats;
    stats.number_of_barcodes_trimmed.resize(barcodes.size());
    stats.number_of_reads_with_adapter.resize(adapters.size());
    return stats;
}



userconfig::alignment_type userconfig::evaluate_alignment(const alignment_info& alignment) const
{
    if (!alignment.length) {
        return not_aligned;
    }

    // Only pairs of called bases are considered part of the alignment
    const size_t n_aligned = static_cast<size_t>(alignment.length - alignment.n_ambiguous);

    size_t mm_threshold = static_cast<size_t>(mismatch_threshold * n_aligned);
    if (n_aligned < 6) {
        mm_threshold = 0;
    } else if (n_aligned < 10) {
        // --mm may imply fewer allowed mismatches than 1, so always compare
        mm_threshold = std::min<size_t>(1, mm_threshold);
    }

    // If too many mismatches
    if (alignment.n_mismatches > mm_threshold) {
        return not_aligned;
    } else if (collapse || identify_adapters) {
        if (n_aligned < static_cast<size_t>(min_alignment_length)) {
            // If the aligned part is too short to collapse the reads,
            // treat them as unaligned. This is also done when attempting to
            // identify adapter sequences, to avoid very short overlaps
            // expected due between the ends of the sequences
            return not_aligned;
        }
    } else if (alignment.score <= 0) {
        // Very poor alignment, will not be considered
        return poor_alignment;
    }

    return valid_alignment;
}


bool userconfig::is_acceptable_read(const fastq& seq) const
{
    return seq.length() >= min_genomic_length
        && seq.count_ns() <= max_ambiguous_bases;
}


void userconfig::open_with_default_filename(std::ofstream& stream,
                                            const std::string& key,
                                            const std::string& postfix) const
{
    std::string filename = basename + postfix;
    if (argparser.is_set(key)) {
        filename = argparser.at(key)->to_str();
    }

    try {
        stream.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        stream.open(filename.c_str(), std::ofstream::out);
    } catch (const std::ofstream::failure&) {
        std::string message = std::string("Failed to open file '") + filename + "': ";
        throw std::ofstream::failure(message + std::strerror(errno));
    }
}


void userconfig::open_ifstream(std::ifstream& stream, const std::string& filename) const
{
    try {
        // failbit is not set, since we expect to try reading past EOF
        stream.exceptions(std::ifstream::badbit);
        stream.open(filename.c_str(), std::ifstream::in);
        if (!stream.is_open()) {
            throw std::ifstream::failure("failed to open file");
        }
    } catch (const std::ifstream::failure&) {
        std::string message = std::string("Failed to open file '") + filename + "': ";
        throw std::ifstream::failure(message + std::strerror(errno));
    }
}


void userconfig::trim_barcodes_if_enabled(fastq& read, statistics& stats) const
{
    if (trim_barcodes_mode) {
        const alignment_info alignment = trim_barcodes(read, barcodes, shift);
        if (alignment.length) {
            stats.number_of_barcodes_trimmed.at(alignment.adapter_id)++;
        }
    }
}


fastq::ntrimmed userconfig::trim_sequence_by_quality_if_enabled(fastq& read) const
{
    fastq::ntrimmed trimmed;
    if (trim_ambiguous_bases || trim_by_quality) {
        char quality_score = trim_by_quality ? low_quality_score : -1;
        trimmed = read.trim_low_quality_bases(trim_ambiguous_bases,
                                              quality_score);
    }

    return trimmed;
}


bool userconfig::read_adapters_sequences(const std::string& filename,
                                         fastq_pair_vec& adapters,
                                         bool paired_ended)
{
    size_t line_num = 1;
    std::ifstream adapter_file;
    try {
        open_ifstream(adapter_file, filename);

        std::string line;
        while (std::getline(adapter_file, line)) {
            const size_t index = line.find_first_not_of(" \t");
            if (index == std::string::npos || line.at(index) == '#') {
                line_num++;
                continue;
            }

            std::stringstream instream(line);

            fastq adapter_5p = read_adapter_sequence(instream);
            fastq adapter_3p;
            if (paired_ended) {
                adapter_3p = read_adapter_sequence(instream);
            }

            adapters.push_back(fastq_pair(adapter_5p, adapter_3p));
            line_num++;
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error reading adapter sequences (line " << line_num
                  << "); aborting:\n    " << error.what() << std::endl;
        return false;
    } catch (const fastq_error& error) {
        std::cerr << "Error parsing adapter sequences (line " << line_num
                  << "); aborting:\n    " << error.what() << std::endl;
        return false;
    }

    return true;
}
