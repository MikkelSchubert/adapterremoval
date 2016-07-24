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
#include <random>

#include "fastq.h"

namespace ar
{

/**
 * Summarizes an alignment.
 *
 * A single offset value is used to represent the alignment between two
 * sequences, with values ranging from -inf to seq1.len() - 1. The alignment
 * of the first base in each sequence against each other is defined as having
 * the offset 0, with each other offset defined as the relative position of
 * seq 2 base 0 to seq 1 base 0:
 *
 * Seq 1:           aaaaaaaaa
 * seq 2:           bbbbbbbbbbb
 * Offset: 0
 *
 * Seq 1:            aaaaaaaaa
 * Seq 2: bbbbbbbbbbb
 * Offset: -11
 *
 * Seq 1:           aaaaaaaaa
 * Seq 2:                   bbbbbbbbbbb
 * Offset: 8
 *
 * The meaning of the offset is slightly different in SE and PE mode; in SE
 * mode seq2 is the adapter sequence, and the offset therefore unambigiously
 * shows the starting position of the adapter, regardless of the size of the
 * adapter sequence.
 *
 * In PE mode, while the alignment is between seq1+PCR2 and PCR1+seq2, the
 * offset returned is relative to to seq1 and seq2 only. Thus, if '1'
 * represents PCR1 and '2' represents PCR2, the following alignment results
 * from PE mode (ie. offset = -9 rather than 3):
 *
 * Seq 1:    22222222222aaaaaaaaa
 * Seq 2:      bbbbbbbbbbb1111111111
 * Offset: -9
 *
 * Note that the offset can never be greater than len(read 1) - 1, but can be
 * less than -len(seq 2) + 1, if a positive shift is set in PE mode. In PE
 * mode, an offset <= -len(seq 2) indicates that no non-adapter sequence was
 * found, while an offset <= 0 indicates the same for SE mode. Offsets less
 * than -len(seq 2) for PE or less than 0 for SE indicates that bases have been
 * skipped during sequencing, and are discoverable if a shift is set:
 *
 * Read 1:           ...22222222222aaaaaaaaa
 * Read 2:             bbbbbbbbbbb1111111111...
 * Offset: -12
 *
 */
struct alignment_info
{
    /** Defaults to unaligned (len = 0), for adapter_id -1. **/
    alignment_info();

    /**
     * Returns true if this is a better alignment than other.
     *
     * When selecting among multiple alignments, the follow criteria are used:
     * 1. The alignment with the highest score is preferred.
     * 2. If score is equal, the longest alignment is preferred.
     * 3. If score and length is equal, the alignment with fewest Ns is preferred.
     */
    bool is_better_than(const alignment_info& other) const;


    //! Alignment score; equal to length - n_ambiguous - 2 * n_mismatches;
    int score;
    //! Zero based id of the adapter which offered the best alignment. Is less
    //! than zero if no alignment was found.
    int offset;
    //! The number of base-pairs included in the alignment. This number
    //! includes both bases aligned between the two mates (in PE mode) and the
    //! number of bases aligned between mates and adapter sequences.
    size_t length;
    //! Number of positions in the alignment in which the two sequences were
    //! both called (not N) but differed
    size_t n_mismatches;
    //! Number of positions in the alignment where one or both bases were N.
    size_t n_ambiguous;
    //! Offset describing the alignment between the two sequences (see above).
    int adapter_id;
};


/**
 * Attempts to align adapters sequences against a SE read.
 *
 * @param read A read potentially containing adapter sequences
 * @param adapters A set of adapter pairs; only the first adapters are used.
 * @param max_shift Allow up to this number of missing bases at the 5' end of
 *                  the read, when aligning the adapter.
 * @return The best alignment, or a length 0 alignment if not aligned.
 *
 * The best alignment is selected using alignment_info::is_better_than.
 */
alignment_info align_single_ended_sequence(const fastq& read,
                                           const fastq_pair_vec& adapters,
                                           int max_shift);


/**
 * Attempts to align PE mates, along with any adapter pairs.
 *
 * @param read1 A mate 1 read potentially containing adapter sequences
 * @param read2 A mate 2 read potentially containing adapter sequences
 * @param adapters A set of adapter pairs; both in each pair adapters are used.
 * @param max_shift Allow up to this number of missing bases at the 5' end of
 *                  both mate reads.
 * @return The best alignment, or a length 0 alignment if not aligned.
 *
 * The alignment is carried out following the concatenation of pcr2 and read1,
 * and the concatenation of read2 and pcr1, resulting in this alignment:
 *
 *                pcr2-read1
 *                read2-pcr1
 *
 * Note the returned offset is relative read1, not to adapter2 + read1,
 * and can be used to directly infer the alignment between read1 and read2.
 */
alignment_info align_paired_ended_sequences(const fastq& read1,
                                            const fastq& read2,
                                            const fastq_pair_vec& adapters,
                                            int max_shift);


/**
 * Truncates a SE read according to the alignment, such that the second read
 * used in the alignment (assumed to represent adapter sequence) is excluded
 * from the read passed to this function.
 */
void truncate_single_ended_sequence(const alignment_info& alignment,
                                    fastq& read);


/**
 * Truncate a pair of PE reads, such that any adapter sequence inferred from
 * the alignment is excluded from both mates.
 *
 * @return The number of sequences (0 .. 2) which contained adapter sequence.
 */
size_t truncate_paired_ended_sequences(const alignment_info& alignment,
                                       fastq& read1,
                                       fastq& read2);


/**
 * Collapses two overlapping PE mates into a single sequence, recalculating the
 * quality scores to reflect the added support offered by two reads at the same
 * nucleotides. In the case of different bases at the same position, the
 * highest quality base is selected; if each base have the same quality score,
 * a random base is selected. In both cases, the quality score is updated to
 * reflect the lower quality implied by these observations.
 *
 * @return A single FASTQ record representing the collapsed sequence.
 *
 * Note that the sequences are assumed to have been trimmed using the
 * truncate_paired_ended_sequences function, and will produce undefined
 * results if this is not the case!
 */
fastq collapse_paired_ended_sequences(const alignment_info& alignment,
                                      const fastq& read1,
                                      const fastq& read2,
                                      std::mt19937& rng,
                                      const char mate_sep=MATE_SEPARATOR);


/**
 * Truncates reads such that only adapter sequence remains.
 *
 * @return True if either or both reads containted adapter sequence.
 *
 * Reads that do not contain any adapter sequence are completely truncated,
 * such no bases remain of the original sequence.
 */
bool extract_adapter_sequences(const alignment_info& alignment,
                               fastq& pcr1,
                               fastq& pcr2);

} // namespace ar

#endif
