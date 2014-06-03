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
#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include "argparse.h"
#include "fastq.h"
#include "alignment.h"
#include "statistics.h"


struct alignment_info;


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

    /**
     * Parses a set of commandline arguments.
     *
     * The function returns false on failure, or if --help / --version or no
     * arguments were supplied by the user.
     */
    bool parse_args(int argc, char *argv[]);


    statistics create_stats() const;


    enum alignment_type
    {
        valid_alignment,
        poor_alignment,
        not_aligned
    };

    alignment_type evaluate_alignment(const alignment_info& alignment) const;

    /** Returns true if the read matches the quality criteria set by the user. **/
    bool is_acceptable_read(const fastq& seq) const;


    void open_with_default_filename(std::ofstream& stream,
                                    const std::string& key,
                                    const std::string& postfix) const;


    void open_ifstream(std::ifstream& stream, const std::string& filename) const;


    void trim_barcodes_if_enabled(fastq& read, statistics& stats) const;


    fastq::ntrimmed trim_sequence_by_quality_if_enabled(fastq& read) const;


    //! Argument parser setup to parse the arguments expected by AR
    argparse::parser argparser;

    //! Prefix used for output files for which no filename was explicitly set
    std::string basename;
    //! Path to input file containing mate 1 reads (required)
    std::string input_file_1;
    //! Path to input file containing mate 2 reads (for PE reads)
    std::string input_file_2;

    //! Set to true if both --input1 and --input2 are set.
    bool paired_ended_mode;

    //! Pairs of adapters; may only contain the first value in SE enabled
    fastq_pair_vec adapters;

    //! Set to true if a nucleotide barcode has been supplied by the user.
    bool trim_barcodes_mode;
    //! Nucleotide barcodes to be trimmed from the 5' termini of mate 1 reads
    //! Only the first value in the pair is defined.
    fastq_pair_vec barcodes;

    //! The minimum length of trimmed reads (ie. genomic nts) to be retained
    unsigned min_genomic_length;
    //! The minimum required genomic overlap before collapsing reads into one.
    unsigned min_alignment_length;
    //! Rate of mismatches determining the threshold for a an acceptable
    //! alignment, depending on the length of the alignment. But see also the
    //! limits set in the function 'evaluate_alignment'.
    double mismatch_threshold;

    //! Quality format expected in input files.
    quality_format quality_input_fmt;
    //! Quality format to write FASTQ records in; either phred_33 and phred_64.
    quality_format quality_output_fmt;

    //! If true, read termini are trimmed for low-quality bases.
    bool trim_by_quality;
    //! The highest quality score which is considered low-quality
    unsigned low_quality_score;

    //! If true, ambiguous bases (N) at read termini are trimmed.
    bool trim_ambiguous_bases;
    //! The maximum number of ambiguous bases (N) in an read; reads exceeding
    //! this number following trimming (optionally) are discarded.
    unsigned max_ambiguous_bases;

    //! If true, PE reads overlapping at least 'min_alignment_length' are
    //! collapsed to generate a higher quality consensus sequence.
    bool collapse;
    // Allow for slipping basepairs by allowing missing bases in adapter
    unsigned shift;

    //! RNG seed for randomly selecting between to bases with the same quality
    //! when collapsing overllapping PE reads.
    unsigned seed;

    //! If true, the program attempts to identify the adapter pair of PE reads
    bool identify_adapters;

private:
    //! Sink for --pcr1, adapter sequence expected at 3' of mate 1 reads
    std::string PCR1;
    //! Sink for --pcr2, adapter sequence expected at 3' of mate 2 reads
    std::string PCR2;
    //! Sink for --5prime, barcode to be trimmed at 5' of mate 1 reads
    std::string barcode;

    //! Sink for user-supplied quality score formats; use quality_input_fmt.
    std::string quality_input_base;
    //! Sink for user-supplied quality score formats; use quality_output_fmt.
    std::string quality_output_base;
};


#endif
