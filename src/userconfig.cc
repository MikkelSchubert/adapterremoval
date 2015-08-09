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
#include <cstring>
#include <cerrno>
#include <string>
#include <stdexcept>
#include <sys/time.h>

#include "userconfig.h"
#include "fastq.h"
#include "alignment.h"
#include "linereader.h"


size_t get_seed()
{
    struct timeval timestamp;
    gettimeofday(&timestamp, NULL);

    return (timestamp.tv_sec << 20) | timestamp.tv_usec;
}


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
    , basename("your_output")
    , input_file_1()
    , input_file_2()
    , paired_ended_mode(false)
    , adapters()
    , barcodes()
    , min_genomic_length(15)
    , max_genomic_length(std::numeric_limits<unsigned>::max())
    , min_alignment_length(11)
    , mismatch_threshold(-1.0)
    , quality_input_fmt(fastq::phred_33)
    , quality_output_fmt(fastq::phred_33)
    , trim_by_quality(false)
    , low_quality_score(2)
    , trim_ambiguous_bases(false)
    , max_ambiguous_bases(1000)
    , collapse(false)
    , shift(2)
    , seed(get_seed())
    , identify_adapters(false)
    , quiet(false)
    , max_threads(1)
    , gzip(false)
    , gzip_level(6)
    , bzip2(false)
    , bzip2_level(9)
    , adapter_1("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
    , adapter_2("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")
    , adapter_list()
    , barcode()
    , barcode_list()
    , quality_input_base("33")
    , quality_output_base("NA")
{
    argparser["--file1"] =
        new argparse::any(&input_file_1, "FILE",
            "Input file containing mate 1 reads or single-ended reads "
            "[REQUIRED].");
    argparser["--file2"] =
        new argparse::any(&input_file_2, "FILE",
            "Input file containing mate 2 reads [OPTIONAL].");

    argparser.add_seperator();
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
            "BASENAME.pair1.truncated (PE) or BASENAME.truncated (SE)]");
    argparser["--output2"] =
        new argparse::any(NULL, "FILE",
            "Output file containing trimmed mate 2 reads [default: "
            "BASENAME.pair2.truncated (only used in PE mode)]");
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
            "the presence of low-quality or ambiguous nucleotides"
            "[default: BASENAME.collapsed.truncated]");
    argparser["--discarded"] =
        new argparse::any(NULL, "FILE",
            "Contains reads discarded due to the --minlength, --maxlength or "
            "--maxns options [default: BASENAME.discarded]");

    argparser.add_seperator();
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
            "List of adapters pairs, used as if supplied to --pcr1 / --pcr2; "
            "only the first adapter in each pair is required / used in SE "
            "mode [current: %default].");

    argparser.add_seperator();
    argparser["--mm"]
        = new argparse::floaty_knob(&mismatch_threshold, "MISMATCH_RATE",
            "Max error-rate when aligning reads and/or adapters. If > 1, the "
            "max error-rate is set to 1 / MISMATCH_RATE; if < 0, the defaults "
            "are used, otherwise the user-supplied value is used directly. "
            "[defaults: 1/3 for trimming; 1/10 when identifing adapters].");
    argparser["--maxns"] =
        new argparse::knob(&max_ambiguous_bases, "MAX",
            "Reads containing more ambiguous bases (N) than this number after "
            "trimming are discarded [current: %default].");
    argparser["--shift"] =
        new argparse::knob(&shift, "N",
            "Consider alignments where up to N nucleotides are missing from "
            "the 5' termini [current: %default].");

    argparser.add_seperator();
    argparser["--qualitybase"] =
        new argparse::any(&quality_input_base, "BASE",
            "Quality base used to encode Phred scores in input; either 33, "
            "64, or solexa [current: %default].");
    argparser["--qualitybase-output"] =
        new argparse::any(&quality_output_base, "BASE",
            "Quality base used to encode Phred scores in output; either 33, "
            "64. By default, reads in Phred+33 format will be written as "
            "Phred+33, while Phred+64 / Solexa reads will be written as "
            "Phred+64.");
    argparser["--5prime"] =
        new argparse::any(&barcode, "BARCODE",
            "If set, the NT barcode is detected (max 1 mismatch) in and "
            "trimmed from the 5' of mate 1 reads [current: %default].");
    argparser["--5prime-list"] =
        new argparse::any(&barcode_list, "FILENAME",
            "List of barcode sequences, with one barcode per line. The best "
            "barcode is selected using the criteria of --5prime, and trimmed "
            "from the 5' end of mate 1 reads [current: %default].");
    argparser["--trimns"] =
        new argparse::flag(&trim_ambiguous_bases,
            "If set, trim ambiguous bases (N) at 5'/3' termini "
            "[current: %default]");
    argparser["--trimqualities"] =
        new argparse::flag(&trim_by_quality,
            "If set, trim bases at 5'/3' termini with quality scores <= to "
            "--minquality value [current: %default]");
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
            "If set, paired ended reads which overlap at least "
            "--minalignmentlength bases are combined into a single consensus "
            "read; for single-ended reads, full inserts are identified by "
            "requiring --minalignmentlength overlap with the adapter sequence "
            "[current: %default].");
    argparser["--minalignmentlength"] =
        new argparse::knob(&min_alignment_length, "LENGTH",
            "If --collapse is set, paired reads must overlap at least this "
            "number of bases to be collapsed, and single-ended reads must "
            "overlap at least this number of bases with the adapter to be "
            "considered complete template molecules [current: %default].");

    argparser.add_seperator();
    argparser["--identify-adapters"] =
        new argparse::flag(&identify_adapters,
            "Attempt to identify the adapter pair of PE reads, by searching "
            "for overlapping reads [current: %default].");
    argparser["--seed"] =
        new argparse::knob(&seed, "SEED",
            "Sets the RNG seed used when choosing between bases with equal "
            "Phred scores when collapsing [current: %default].");

    argparser.add_seperator();
    argparser["--gzip"] =
        new argparse::flag(&gzip,
            "Enable gzip compression [current: %default]");
    argparser["--gzip-level"] =
        new argparse::knob(&gzip_level, "LEVEL",
            "Compression level, 0 - 9 [current: %default]");
#ifdef AR_BZIP2_SUPPORT
    argparser["--bzip2"] =
        new argparse::flag(&bzip2,
            "Enable bzip2 compression [current: %default]");
    argparser["--bzip2-level"] =
        new argparse::knob(&bzip2_level, "LEVEL",
            "Compression level, 0 - 9 [current: %default]");
#endif
    argparser["--quiet"] =
        new argparse::flag(&quiet,
            "Only print warnings / errors to STDERR [current: %default]");

#ifdef AR_PTHREAD_SUPPORT
    argparser["--threads"] =
        new argparse::knob(&max_threads, "THREADS",
            "Maximum number of threads [current: %default]");
#endif
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

    if (argparser.is_set("--qualitybase")) {
        if (quality_input_base == "33") {
            quality_input_fmt = fastq::phred_33;
        } else if (quality_input_base == "64") {
            quality_input_fmt = fastq::phred_64;
        } else if (quality_input_base == "solexa") {
            quality_input_fmt = fastq::solexa;
        } else {
            std::cerr << "Error: Invalid value for --qualitybase: '"
                      << quality_input_base << "'\n"
                      << "   expected values 33, 64, or solexa." << std::endl;
            return argparse::pr_error;
        }
    } else if (identify_adapters) {
        // By default quality scores are ignored when inferring adapter sequences
        quality_input_fmt = fastq::ignored;
    }

    if (quality_output_base == "33") {
        quality_output_fmt = fastq::phred_33;
    } else if (quality_output_base == "64") {
        quality_output_fmt = fastq::phred_64;
    } else if (quality_output_base == "NA") {
        if (quality_input_fmt == fastq::phred_33) {
            quality_output_fmt = fastq::phred_33;
        } else {
            quality_output_fmt = fastq::phred_64;
        }
    } else {
        std::cerr << "Error: Invalid value for --qualitybase-out: '"
                  << quality_output_base << "'\n"
                  << "   expected values 33 or 64." << std::endl;
        return argparse::pr_error;
    }

    if (low_quality_score > MAX_PHRED_SCORE) {
        std::cerr << "Error: Invalid value for --minquality: "
                  << low_quality_score << "\n"
                  << "   must be in the range 0 .. " << MAX_PHRED_SCORE
                  << std::endl;
        return argparse::pr_error;
    }

    if (!setup_barcode_sequences()) {
        return argparse::pr_error;
    }

    // Check for invalid combinations of settings
    const bool file_1_set = argparser.is_set("--file1");
    const bool file_2_set = argparser.is_set("--file2");

    if (!(file_1_set || file_2_set)) {
        std::cerr << "Error: No input files (--file1 / --file2) specified.\n"
                  << "Please specify at least one input file using --file1 FILENAME."
                  << std::endl;
        return argparse::pr_error;
    } else if (file_2_set && !file_1_set) {
        std::cerr << "Error: --file2 specified, but --file1 is not specified." << std::endl;
        return argparse::pr_error;
    } else if (identify_adapters && !(file_1_set && file_2_set)) {
        std::cerr << "Error: Both input files (--file1 / --file2) must be "
                  << "specified when using --identify-adapters."
                  << std::endl;
        return argparse::pr_error;
    }

    paired_ended_mode = file_2_set;

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
    }

    return argparse::pr_ok;
}


std::auto_ptr<statistics> userconfig::create_stats() const
{
    std::auto_ptr<statistics> stats(new statistics());
    stats->number_of_barcodes_trimmed.resize(barcodes.size());
    stats->number_of_reads_with_adapter.resize(adapters.size());
    return stats;
}



userconfig::alignment_type userconfig::evaluate_alignment(const alignment_info& alignment) const
{
    if (!alignment.length) {
        return not_aligned;
    } else if (alignment.score <= 0) {
        // Very poor alignment, will not be considered
        return poor_alignment;
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

    if (alignment.n_mismatches > mm_threshold) {
        return not_aligned;
    }

    return valid_alignment;
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
    return seq.length() >= min_genomic_length
        && seq.length() <= max_genomic_length
        && seq.count_ns() <= max_ambiguous_bases;
}


std::auto_ptr<std::ostream> userconfig::open_with_default_filename(
                                            const std::string& key,
                                            const std::string& postfix,
                                            bool compressed) const
{
    std::string filename = basename + postfix;
    if (argparser.is_set(key)) {
        filename = argparser.at(key)->to_str();
    } else if (compressed && gzip && gzip_level) {
        filename += ".gz";
    } else if (compressed && bzip2) {
        filename += ".bz2";
    }

    std::auto_ptr<std::ofstream> stream(new std::ofstream(filename.c_str(),
                                                           std::ofstream::out));

    if (!stream->is_open()) {
        std::string message = std::string("Failed to open file '") + filename + "': ";
        throw std::ofstream::failure(message + std::strerror(errno));
    }

    stream->exceptions(std::ofstream::failbit | std::ofstream::badbit);
    return std::auto_ptr<std::ostream>(stream);
}


void userconfig::trim_barcodes_if_enabled(fastq& read, statistics& stats) const
{
    if (!barcodes.empty()) {
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


bool userconfig::setup_barcode_sequences()
{
    const bool barcode_is_set = argparser.is_set("--5prime");
    const bool barcode_list_is_set = argparser.is_set("--5prime-list");

    if (barcode_is_set && barcode_list_is_set) {
        std::cerr << "Error: Use either --5prime or --5prime-list, not both!" << std::endl;
        return argparse::pr_error;
    } else if (barcode_list_is_set) {
        if (!read_adapter_sequences(barcode_list, barcodes, "barcode", false)) {
            return false;
        } else if (barcodes.empty()) {
            std::cerr << "Error: No barcodes sequences found in table!" << std::endl;
            return false;
        }
    } else if (barcode_is_set) {
        if (!cleanup_and_validate_sequence(barcode, "--5prime")) {
            return false;
        }

        fastq barcode_fq = fastq("PCR1", barcode, std::string(barcode.length(), 'J'));

        barcodes.push_back(fastq_pair(barcode_fq, fastq("DUMMMY", "", "")));
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
        if (!read_adapter_sequences(adapter_list, adapters, "adapter", paired_ended_mode)) {
            return false;
        } else if (adapters.empty()) {
            std::cerr << "Error: No adapter sequences found in table!" << std::endl;
            return false;
        }
    } else {
        const char* label_1 = pcr_is_set ? "--pcr1" : "--adapter1";
        const char* label_2 = pcr_is_set ? "--pcr2" : "--adapter2";

        if (!cleanup_and_validate_sequence(adapter_1, label_1)) {
            return false;
        } else if (!cleanup_and_validate_sequence(adapter_2, label_2)) {
            return false;
        }

        fastq adapter1 = fastq("PCR1", adapter_1, std::string(adapter_1.length(), 'J'));
        fastq adapter2 = fastq("PCR2", adapter_2, std::string(adapter_2.length(), 'J'));
        if (!pcr_is_set) {
            // --pcr2 is expected to already be reverse completed; whereas
            // --adapter2 should correspond to the sequences observed directly
            // in the .fastq files.
            adapter2.reverse_complement();
        }

        adapters.push_back(fastq_pair(adapter1, adapter2));
    }

    return true;

}


bool userconfig::read_adapter_sequences(const std::string& filename,
                                         fastq_pair_vec& adapters,
                                         const std::string& name,
                                         bool paired_ended) const
{
    size_t line_num = 1;
    try {
        line_reader adapter_file(filename);

        std::string line;
        while (adapter_file.getline(line)) {
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
                adapter_3p.reverse_complement();
            }

            adapters.push_back(fastq_pair(adapter_5p, adapter_3p));
            line_num++;
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error reading " << name << " sequences (line " << line_num
                  << "); aborting:\n    " << error.what() << std::endl;
        return false;
    } catch (const fastq_error& error) {
        std::cerr << "Error parsing " << name << " sequences (line " << line_num
                  << "); aborting:\n    " << error.what() << std::endl;
        return false;
    }

    return true;
}
