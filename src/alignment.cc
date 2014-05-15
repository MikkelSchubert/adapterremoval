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

#include <cstdlib>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include <cstring>

#include "alignment.h"
#include "fastq.h"


#ifndef NO_SSE
#include <xmmintrin.h>

//! Mask representing those (sparse) bits used when comparing multiple
//! nucleotides. These are simply the lowest significane bit in each byte.
const __m128i BIT_MASK_128 = _mm_set1_epi8(1);
//! Zero'd 128b integer.
const __m128i ZERO_128 = _mm_set1_epi8(0);

/** Counts the number of bits set in a __m128i. **/
inline size_t COUNT_BITS_128(__m128i value)
{
    // Calculates the abs. difference between each pair of bytes in the upper
    // and lower 64bit integers, and places the sum of these differences in
    // the 0th and 4th shorts (16b).
    value = _mm_sad_epu8(ZERO_128, value);
    // Return the 0th and 4th shorts containing the sums calculated above
    return _mm_extract_epi16(value, 0) + _mm_extract_epi16(value, 4);
}
#else
#warning SEE optimizations disabled!
#endif


/**
 * Compares two subsequences in an alignment to a previous (best) alignment.
 *
 * @param best The currently best alignment, used for evaluating this alignment
 * @param current The current alignment to be evaluated (counts are assumed to be zero'd!)
 * @param seq_1_ptr Pointer to the first base in the first sequence in the alignment.
 * @param seq_2_ptr Pointer to the first base in the second sequence in the alignment.
 * @return True if the current alignment is at least as good as the best alignment, false otherwise.
 *
 * If the function returns false, the current alignment cannot be assumed to
 * have been completely evaluated (due to early termination), and hence counts
 * and scores are not reliable. The function furthermore assumes that the
 * sequences contain only the characters "ACGTN", and will produce undefined
 * results if this is not the case!
 */
bool compare_subsequences(const alignment_info& best, alignment_info& current,
                          const char* seq_1_ptr, const char* seq_2_ptr)
{
    for (int remaining_bases = current.length; remaining_bases > 0;) {
#ifndef NO_SSE
        if (remaining_bases >= 8) {
            size_t nbases = 0;
            __m128i s1, s2;
            if (remaining_bases >= 16) {
                s1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_1_ptr));
                s2 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_2_ptr));
                nbases = 16;
            } else {
                s1 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(seq_1_ptr));
                s2 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(seq_2_ptr));
                nbases = 8;
            }

            seq_1_ptr += nbases;
            seq_2_ptr += nbases;

            remaining_bases -= nbases;

            // Equivalent to (s1 | s2) >> 3
            // The third bit is uniquely set for N, but not for A/C/G/T
            const __m128i ns_mask = _mm_srli_epi16(_mm_or_si128(s1, s2), 3);
            // Calculate the 
            const __m128i mm_mask = ~_mm_or_si128(_mm_cmpeq_epi8(s1, s2), ns_mask);

            current.n_ambiguous += COUNT_BITS_128(_mm_and_si128(ns_mask, BIT_MASK_128));
            current.n_mismatches += COUNT_BITS_128(_mm_and_si128(mm_mask, BIT_MASK_128));
        } else
#endif
        {
            const char nt_1 = *seq_1_ptr++;
            const char nt_2 = *seq_2_ptr++;

            if (nt_1 == 'N' || nt_2 == 'N') {
                current.n_ambiguous++;
            } else if (nt_1 != nt_2) {
                current.n_mismatches++;
            }

            remaining_bases -= 1;
        }

        // Matches count for 1, Ns for 0, and mismatches for -1
        current.score = (current.length - remaining_bases) - current.n_ambiguous - (current.n_mismatches * 2);
        // Terminate early if the remaining alignments cannot possibly
        // offer an alignment that would be accepted (see below).
        if (current.score + remaining_bases < best.score) {
            return false;
        }
    }


    return (current.score >= best.score);
}


alignment_info pairwise_align_sequences(const std::string& seq1,
                                        const std::string& seq2,
                                        int min_offset = std::numeric_limits<int>::min(),
                                        int max_offset = std::numeric_limits<int>::max())
{
    const int start_offset = std::min<int>(max_offset, static_cast<int>(seq1.length()) - 1);
    const int end_offset = std::max<int>(min_offset, -static_cast<int>(seq2.length()) + 1);

    alignment_info best;
    for (int offset = start_offset; offset >= end_offset; --offset) {
        const size_t initial_seq1_offset = std::max<int>(0,  offset);
        const size_t initial_seq2_offset = std::max<int>(0, -offset);

        alignment_info current;
        current.offset = offset;
        current.length = std::min(seq1.length() - initial_seq1_offset,
                                  seq2.length() - initial_seq2_offset);

        const char* seq_1_ptr = seq1.data() + initial_seq1_offset;
        const char* seq_2_ptr = seq2.data() + initial_seq2_offset;

        if (compare_subsequences(best, current, seq_1_ptr, seq_2_ptr)) {
            best = current;
        }
    }

    return best;
}


struct phred_scores
{
    phred_scores()
      : identical_nts(std::numeric_limits<char>::min())
      , different_nts(std::numeric_limits<char>::min())
    {
    }

    //! Phred score to assign if the two nucleotides are identical
    char identical_nts;
    //! Phred score to assign if the two nucleotides differ
    char different_nts;
};


/**
 * Calculates the phred score to be assigned to a consensus base, depending on
 * the 
 *
 * The function assumes that the 
 */
std::vector<phred_scores> calculate_phred_score()
{
    std::vector<double> Perror(MAX_PHRED_SCORE + 1, 0.0);
    std::vector<double> Ptrue(MAX_PHRED_SCORE + 1, 0.0);
    for (size_t i = 0; i <= MAX_PHRED_SCORE; ++i) {
        const double p_err = std::min<double>(0.75, std::pow(10, static_cast<double>(i) / -10.0));

        Perror.at(i) = std::log(p_err / 3.0);
        Ptrue.at(i) = std::log(1.0 - p_err);
    }

    std::vector<phred_scores> new_scores((MAX_PHRED_SCORE + 1) * (MAX_PHRED_SCORE + 1));
    for (size_t i = 0; i <= MAX_PHRED_SCORE; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            const size_t index = (i * MAX_PHRED_SCORE) + j;
            phred_scores& scores = new_scores.at(index);

            {   // When two nucleotides are identical
                const double ptrue = Ptrue.at(i) + Ptrue.at(j);
                const double perror = Perror.at(i) + Perror.at(j);
                const double normconstant = 1.0 + 3.0 * std::exp(perror - ptrue);
                scores.identical_nts = p_to_phred_33(1.0 - 1.0 / normconstant);
            }

            {   // When two nucleotides differ
                const double ptrue = Ptrue.at(i) + Perror.at(j);
                const double perror_one = Perror.at(i) + Ptrue.at(j);
                const double perror_both = Perror.at(i) + Perror.at(j);
                const double normconstant = 1.0 + 2.0 * std::exp(perror_both - ptrue) + std::exp(perror_one - ptrue);
                scores.different_nts = p_to_phred_33(1.0 - 1.0 / normconstant);
            }
        }
    }

    return new_scores;
}


const phred_scores& get_updated_phred_scores(char qual_1, char qual_2)
{
    if (qual_1 < qual_2) {
        throw std::invalid_argument("qual_1 must be superior or equal to qual_2");
    }

    // Cache of pre-calculated Phred scores for consensus bases; see above.
    static const std::vector<phred_scores> updated_phred_scores = calculate_phred_score();

    const size_t index = (static_cast<size_t>(qual_1 - PHRED_OFFSET_33) * MAX_PHRED_SCORE) + static_cast<size_t>(qual_2 - PHRED_OFFSET_33);
    return updated_phred_scores.at(index);
}


typedef std::pair<std::string, std::string> string_pair;

string_pair collapse_sequence(const std::string& sequence1,
                              const std::string& sequence2,
                              const std::string& qualities1,
                              const std::string& qualities2)
{
    if (sequence1.length() != sequence2.length() ||
        sequence1.length() != qualities1.length() ||
        sequence2.length() != qualities2.length()) {
        throw std::invalid_argument("length mismatch between sequences and/or qualities");
    }

    std::string collapsed_seq(sequence1.length(), 'X');
    std::string collapsed_qual(sequence1.length(), '\0');

    // Find new consensus sequence
    for (size_t i = 0; i < collapsed_seq.size(); ++i) {
        char nt_1 = sequence1.at(i);
        char nt_2 = sequence2.at(i);
        char qual_1 = qualities1.at(i);
        char qual_2 = qualities2.at(i);

        // Handle the case of different nucleotides with identical quality scores
        if (nt_1 == 'N' || nt_2 == 'N') {
            // If one of the bases are N, then we suppose that we just have (at
            // most) a single read at that site.
            if (nt_1 != 'N') {
                collapsed_seq.at(i) = nt_1;
                collapsed_qual.at(i) = qual_1;
            } else if (nt_2 != 'N') {
                collapsed_seq.at(i) = nt_2;
                collapsed_qual.at(i) = qual_2;
            } else {
                collapsed_seq.at(i) = 'N';
                collapsed_qual.at(i) = PHRED_OFFSET_33;
            }
        } else if (nt_1 != nt_2 && qual_1 == qual_2) {
            // Randomly select one of the reads
            const int shuffle = random() % 2;
            collapsed_seq.at(i) = shuffle ? nt_1 : nt_2;

            const phred_scores& new_scores = get_updated_phred_scores(qual_1, qual_2);
            collapsed_qual.at(i) = new_scores.different_nts;
        } else {
            // Ensure that nt_1 / qual_1 always contains the preferred nt / score
            // This is an assumption of the g_updated_phred_scores cache.
            if (qual_1 < qual_2) {
                std::swap(nt_1, nt_2);
                std::swap(qual_1, qual_2);
            }

            const phred_scores& new_scores = get_updated_phred_scores(qual_1, qual_2);

            collapsed_seq.at(i) = nt_1;
            collapsed_qual.at(i) = (nt_1 == nt_2) ? new_scores.identical_nts : new_scores.different_nts;
        }
    }

    return string_pair(collapsed_seq, collapsed_qual);
}


///////////////////////////////////////////////////////////////////////////////
// Public functions


char p_to_phred_33(double p)
{
    const int raw_score = static_cast<int>(-10.0 * std::log10(p));
    const char phred_score = static_cast<char>(std::min<int>(MAX_PHRED_SCORE, raw_score));
    return phred_score + PHRED_OFFSET_33;
}



alignment_info::alignment_info()
    : n_ambiguous(0)
    , n_mismatches(0)
    , score(0)
    , offset(0)
    , length(0)
{
}


bool truncate_barcode(fastq& read, const std::string& barcode, int shift)
{
    // Nothing to truncate
    if (!read.length() || barcode.empty()) {
        return false;
    }

    const alignment_info alignment = pairwise_align_sequences(read.sequence(),
                                                              barcode,
                                                              -shift,
                                                              0);

    if (alignment.length - alignment.n_ambiguous && alignment.n_mismatches <= 1) {
        read.truncate(static_cast<size_t>(barcode.length() + alignment.offset));
        return true;
    }

    return false;
}


alignment_info align_single_ended_sequence(const fastq& read,
                                           const fastq& adapter,
                                           int shift)
{
    return pairwise_align_sequences(read.sequence(), adapter.sequence(), -shift);
}


alignment_info align_paired_ended_sequences(const fastq& read1,
                                            const fastq& read2,
                                            const fastq& adapter1,
                                            const fastq& adapter2,
                                            int shift)
{
    const std::string sequence1 = adapter2.sequence() + read1.sequence();
    const std::string sequence2 = read2.sequence() + adapter1.sequence();
    // Only consider alignments where at least one nucleotide from each read
    // is aligned against the other, included shifted alignments to account
    // for missing bases at the 5' ends of the reads.
    #warning validate!
    const int min_offset = adapter2.length() - read2.length() - shift;
    alignment_info alignment = pairwise_align_sequences(sequence1, sequence2, min_offset);
    if (alignment.length) {
        // Convert the alignment into an alignment involving only read 1/2
        alignment.offset -= adapter2.length();
    }

    return alignment;
}


void truncate_single_ended_sequence(const alignment_info& alignment,
                                    fastq& read)
{
    // Given a shift, the alignment of the adapter may start one or more
    // bases before the start of the sequence, leading to a negative offset
    const size_t len = std::max<int>(0, alignment.offset);

    return read.truncate(0, len);
}



size_t truncate_paired_ended_sequences(const alignment_info& alignment,
                                     fastq& read1,
                                     fastq& read2)
{
    size_t had_adapter = 0;
    const int template_length = std::max<int>(0, static_cast<int>(read2.length()) + alignment.offset);
    if (alignment.offset > static_cast<int>(read1.length())) {
        throw std::invalid_argument("invalid offset");
    } else if (alignment.offset >= 0) {
        // Read1 can potentially extend past read2, but by definition read2
        // cannot extend past read1 when the offset is not negative, so there
        // is no need to edit read2.
        had_adapter += static_cast<size_t>(template_length) < read1.length();
        read1.truncate(0, static_cast<size_t>(template_length));
    } else {
        had_adapter += static_cast<size_t>(template_length) < read1.length();
        had_adapter += static_cast<size_t>(template_length) < read2.length();

        read1.truncate(0, static_cast<size_t>(template_length));
        read2.truncate(static_cast<size_t>(static_cast<int>(read2.length()) - template_length));
    }

    return had_adapter;
}


fastq collapse_paired_ended_sequences(const alignment_info& alignment,
                                      const fastq& read1,
                                      const fastq& read2)
{
    if (alignment.offset > 0) {
        // Offset to the first base overlapping read 2
        const size_t read_1_offset = static_cast<size_t>(alignment.offset);
        const std::string read_1_seq = read1.sequence().substr(0, read_1_offset);
        const std::string read_1_qual = read1.qualities().substr(0, read_1_offset);

        // Offset to the last base overlapping read 1
        const size_t read_2_offset = static_cast<size_t>(static_cast<int>(read1.length()) - alignment.offset);
        const std::string read_2_qual = read2.qualities().substr(read_2_offset);
        const std::string read_2_seq = read2.sequence().substr(read_2_offset);

        string_pair collapsed = collapse_sequence(read1.sequence().substr(read_1_offset),
                                                  read2.sequence().substr(0, read_2_offset),
                                                  read1.qualities().substr(read_1_offset),
                                                  read2.qualities().substr(0, read_2_offset));

        return fastq(read1.header(),
                     read_1_seq + collapsed.first + read_2_seq,
                     read_1_qual + collapsed.second + read_2_qual);
    } else {
        string_pair collapsed = collapse_sequence(read1.sequence(),
                                                  read2.sequence(),
                                                  read1.qualities(),
                                                  read2.qualities());

        return fastq(read1.header(),
                     collapsed.first,
                     collapsed.second);
    }
}

bool extract_adapter_sequences(const alignment_info& alignment,
                               fastq& read1,
                               fastq& read2)
{
    const int template_length = std::max(0, static_cast<int>(read2.length()) + alignment.offset);
    if (alignment.offset > static_cast<int>(read1.length())) {
        throw std::invalid_argument("invalid offset");
    } else if (alignment.offset >= 0) {
        read1.truncate(std::min<size_t>(read1.length(), template_length));
        read2.truncate(0, 0);
    } else {
        read1.truncate(template_length);
        read2.truncate(0, static_cast<size_t>(static_cast<int>(read2.length()) - template_length));
    }

    return read1.sequence().length() || read2.sequence().length();
}
