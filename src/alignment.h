/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
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
#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>

class fastq;



char p_to_phred_33(double p);


struct alignment_info
{
    alignment_info();

    size_t n_ambiguous;
    size_t n_mismatches;
    int score;
    int offset;
    size_t length;
};



bool truncate_barcode(fastq& read, const std::string& barcode, int shift);



alignment_info align_single_ended_sequence(const fastq& read,
                                           const fastq& adapter,
                                           int shift);


/**
 *
 *
 *
 * Note the returned offset is relative read1, not to adapter2 + read1,
 * and can be used to directly infer the alignment between read1 and read2.
 */
alignment_info align_paired_ended_sequences(const fastq& read1,
                                            const fastq& read2,
                                            const fastq& adapter1,
                                            const fastq& adapter2,
                                            int shift);


void truncate_single_ended_sequence(const alignment_info& alignment,
                                    fastq& read);


/**
 *
 * @return The number of sequences (0 .. 2) which were truncated.
 */
size_t truncate_paired_ended_sequences(const alignment_info& alignment,
                                       fastq& read1,
                                       fastq& read2);


fastq collapse_paired_ended_sequences(const alignment_info& alignment,
                                      const fastq& read1,
                                      const fastq& read2);


bool extract_adapter_sequences(const alignment_info& alignment,
                               fastq& pcr1,
                               fastq& pcr2);

#endif
