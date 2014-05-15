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
#include <limits>
#include <sstream>
#include <gtest/gtest.h>

#include "alignment.h"
#include "fastq.h"


alignment_info new_aln(int offset, int score = 0, size_t length = 0, size_t nmm = 0, size_t nn = 0)
{
    alignment_info aln;
    aln.offset = offset;
    aln.score = score;
    aln.length = length;
    aln.n_mismatches = nmm;
    aln.n_ambiguous = nn;

    return aln;
}

bool operator==(const alignment_info& first, const alignment_info& second)
{
    return (first.offset == second.offset)
        && (first.score == second.score)
        && (first.length == second.length)
        && (first.n_mismatches == second.n_mismatches)
        && (first.n_ambiguous == second.n_ambiguous);
}


std::ostream& operator<<(std::ostream& stream, const alignment_info& aln)
{
    stream << "alignment_info(" << aln.offset << ", "
                                << aln.score << ", "
                                << aln.length << ", "
                                << aln.n_mismatches << ", "
                                << aln.n_ambiguous << ")";
    return stream;
}


void ASSERT_TRUNCATED_PE_IS_UNCHANGED(const alignment_info& alignment,
                                      const fastq& record1,
                                      const fastq& record2)
{
    fastq tmp_record1 = record1;
    fastq tmp_record2 = record2;
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, tmp_record1, tmp_record2));
    ASSERT_EQ(record1, tmp_record1);
    ASSERT_EQ(record2, tmp_record2);
}

void compare_subsequences(const alignment_info& best, alignment_info& current,
                          const char* seq_1_ptr, const char* seq_2_ptr);





///////////////////////////////////////////////////////////////////////////////
// Cases for SE alignments (a = read 1, b = adapter, o = overlap):
//  1. No overlap = aaaaaa bbbbbb
//  2. Partial overlap = aaaaaooobbbbbb
//  3. Complete overlap = ooooooooo
//  4. a contains b = aaaaoooooo
//  5. b contains a = ooooobbbbb
//  6. a extends past b = aaaaoooooaaaa
//  7. b extends past a = bbbbooooobbbb
//  8. both extends past other = bbbbooooaaaa


///////////////////////////////////////////////////////////////////////////////
// Case 1: No overlap between sequences:
//         AAAAAAAAAAA
//                       BBBBBBBBBBBBBBB
// Expected result = suboptimal alignment or no overlap

TEST(alignment_se, unalignable_sequence)
{
    const fastq record("Rec",  "AAAA", "!!!!");
    const fastq adapter("Rec", "TTTT", "!!!!");
    const alignment_info expected;
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);
}


TEST(alignment_se, no_expected_overlap)
{
    const fastq record("Rec",  "ACGTAGTA",  "!!!!!!!!");
    const fastq adapter("Rec", "TGAGACGGT", "!!!!!!!!!");
    const alignment_info expected = new_aln(6, 0, 2, 1, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);
}


///////////////////////////////////////////////////////////////////////////////
// Case 2: Partial overlap between sequences:
//         AAAAAAAAAAA
//                BBBBBBBBBBBBBBB
// Expected result = optimal alignment between ends
TEST(alignment_se, partial_overlap)
{
    const fastq record("Rec",  "ACGTAGTAA", "123457890");
    const fastq adapter("Rec", "AGTAAGGT",  "!!!!!!!!");
    const alignment_info expected = new_aln(4, 5, 5, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACGT", "1234"), tmp_record);
}


TEST(alignment_se, partial_overlap_with_mismatch)
{
    const fastq record("Rec",  "ACGTAGTAA", "123457890");
    const fastq adapter("Rec", "AGGAAGGT",  "!!!!!!!!");
    const alignment_info expected = new_aln(4, 3, 5, 1, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACGT", "1234"), tmp_record);
}


TEST(alignment_se, partial_overlap_with_n)
{
    const fastq record("Rec",  "ACGTAGTAA", "123457890");
    const fastq adapter("Rec", "AGNAAGGT",  "!!!!!!!!");
    const alignment_info expected = new_aln(4, 4, 5, 0, 1);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACGT", "1234"), tmp_record);
}


///////////////////////////////////////////////////////////////////////////////
// Case 3: Complete overlap between sequences:
//         AAAAAAAAAAA
//         BBBBBBBBBBB
// Expected result = Optimal alignment involving all bases

TEST(alignment_se, completely_overlapping_sequences)
{
    const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
    const fastq adapter = record;
    const alignment_info expected = new_aln(0, 8, 8, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, completely_overlapping_sequences_with_1_mismatch)
{
    const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
    const fastq adapter("Rec", "GCGTAGTA", "!!!!!!!!");
    const alignment_info expected = new_aln(0, 6, 8, 1, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, completely_overlapping_sequences_with_1_mismatch_and_1_n)
{
    const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
    const fastq adapter("Rec", "GCGTAGTN", "!!!!!!!!");
    const alignment_info expected = new_aln(0, 5, 8, 1, 1);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


///////////////////////////////////////////////////////////////////////////////
// Cases 4 and 5: Sequence A completely contains sequence B (and vice versa)
//         AAAAAAAAAAA  AAAA
//               BBBBB  BBBBBBBBBB
// Expected result = Optimal alignment involving all bases of the shortest read

TEST(alignment_se, sequence_a_contains_b)
{
    const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
    const fastq adapter("Adp",   "TAGTA", "!!!!!");
    const alignment_info expected = new_aln(3, 5, 5, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACG", "ABC"), tmp_record);
}


TEST(alignment_se, sequence_a_contains_b__with_1_mismatch)
{
    const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
    const fastq adapter("Adp", "TATTA", "!!!!!");
    const alignment_info expected = new_aln(3, 3, 5, 1, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACG", "ABC"), tmp_record);
}


TEST(alignment_se, sequence_a_contains_b__with_1_n)
{
    const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
    const fastq adapter("Adp", "TAGNA", "!!!!!");
    const alignment_info expected = new_aln(3, 4, 5, 0, 1);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACG", "ABC"), tmp_record);
}


TEST(alignment_se, sequence_b_contains_a)
{
    const fastq record("Rec",  "ACGT", "!!!!");
    const fastq adapter("Adp", "ACGTAGTA", "!!!!!!!!");
    const alignment_info expected = new_aln(0, 4, 4, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, sequence_b_contains_a__with_1_mismatch)
{
    const fastq record("Rec",  "ACGT",     "!!!!");
    const fastq adapter("Adp", "GCGTAGTA", "!!!!!!!!");
    const alignment_info expected = new_aln(0, 2, 4, 1, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, sequence_b_contains_a__with_1_n)
{
    const fastq record("Rec",  "ACGT", "!!!!");
    const fastq adapter("Adp", "ACGNAGTA", "!!!!!!!!");
    const alignment_info expected = new_aln(0, 3, 4, 0, 1);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


///////////////////////////////////////////////////////////////////////////////
// Cases 6 and 7: Sequence A extends past sequence B (and vice versa)
//         AAAAAAAAAAAAAA    AAAA
//               BBBBB     BBBBBBBBBBBBB

TEST(alignment_se, sequence_a_extends_past_b)
{
    const fastq record("Rec",  "ACGTAGTATA", "0123456789");
    const fastq adapter("Adp",     "AGTA",   "!!!!");
    const alignment_info expected = new_aln(4, 4, 4, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACGT", "0123"), tmp_record);
}


TEST(alignment_se, sequence_b_extends_past_a__no_shift)
{
    const fastq record("Rec",   "CGTA",      "#!%%");
    const fastq adapter("Adp", "ACGTAGTATA", "!!!!!!!!!!");
    const alignment_info expected = new_aln(3, 1, 1, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "CGT", "#!%"), tmp_record);
}


TEST(alignment_se, sequence_b_extends_past_a__shift_of_1)
{
    const fastq record("Rec",   "CGTA",      "#!%%");
    const fastq adapter("Adp", "ACGTAGTATA", "!!!!!!!!!!");
    const alignment_info expected = new_aln(-1, 4, 4, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 1);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


///////////////////////////////////////////////////////////////////////////////
// Cases 8: Both sequences extends past the other
//                      AAAAAAAAAAAAAAAAAA
//               BBBBBBBBBBBBBBBBBB
// Expected result = Optimal alignment involving all overlapping bases

TEST(alignment_se, sequences_extend_past_mate)
{
    const fastq record("Rec",      "ACGTAGTATATAGT", "!!!!!!!!!!!!!!");
    const fastq adapter("Adp", "CCGAACGTAGTATA",     "!!!!!!!!!!!!!!");
    const alignment_info expected = new_aln(-4, 10, 10, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapter, 4);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


///////////////////////////////////////////////////////////////////////////////
// Misc

TEST(alignment_se, shift_is_lower_than_possible)
{
    const fastq record("Rec",  "AAAA", "!!!!");
    const fastq adapter("Adp", "TTTT", "!!!!");
    const alignment_info expected;
    const alignment_info result = align_single_ended_sequence(record, adapter, -10);
    ASSERT_EQ(expected, result);
}


///////////////////////////////////////////////////////////////////////////////
// Case 1 PE: No overlap between sequences:
//         AAAAAAAAAAA
//                       BBBBBBBBBBBBBBB
// Expected result = suboptimal alignment or no overlap

TEST(alignment_pe, unalignable_sequence)
{
    const fastq record1("Rec",  "AAAA", "!!!!");
    const fastq record2("Rec", "TTTT", "!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected;
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}


TEST(alignment_pe, no_expected_overlap)
{
    const fastq record1("Rec",  "ACGTAGTA",  "!!!!!!!!");
    const fastq record2("Rec", "TGAGACGGT", "!!!!!!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(6, 0, 2, 1, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}


///////////////////////////////////////////////////////////////////////////////
// Case 2 PE: Partial overlap between sequences:
//         AAAAAAAAAAA
//                BBBBBBBBBBBBBBB
// Expected result = optimal alignment between ends

TEST(alignment_pe, partial_overlap)
{
    const fastq record1("Rec", "ACGTAGTAA", "!!!!!!!!!");
    const fastq record2("Rec", "AGTAAGGT",  "!!!!!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(4, 5, 5, 0, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}


///////////////////////////////////////////////////////////////////////////////
// Case 3 PE: Complete overlap between sequences:
//         AAAAAAAAAAA
//         BBBBBBBBBBB
// Expected result = Optimal alignment involving all bases

TEST(alignment_pe, completely_overlapping_sequences)
{
    const fastq record1("Rec", "ACGTAGTA", "!!!!!!!!");
    const fastq record2 = record1;
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(0, 8, 8, 0, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}


///////////////////////////////////////////////////////////////////////////////
// Cases 4 and 5 PE: Sequence A completely contains sequence B (and vice versa)
//         AAAAAAAAAAA  AAAA
//               BBBBB  BBBBBBBBBB
// Expected result = Optimal alignment involving all bases of the shortest read

TEST(alignment_pe, sequence_a_contains_b)
{
    const fastq record1("Rec1", "ACGTAGTA", "!!!!!!!!");
    const fastq record2("Rec2", "TAGTA", "!!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(3, 5, 5, 0, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

TEST(alignment_pe, sequence_b_contains_a)
{
    const fastq record1("Rec1", "ACGT", "!!!!");
    const fastq record2("Rec2", "ACGTAGTA", "!!!!!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(0, 4, 4, 0, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}


///////////////////////////////////////////////////////////////////////////////
// Cases 6 and 7 PE: Sequence A extends past sequence B (and vice versa)
//         AAAAAAAAAAAAAA    AAAA
//               BBBBB     BBBBBBBBBBBBB

TEST(alignment_pe, sequence_a_extends_past_b)
{
    const fastq record1("Rec1",  "ACGTAGTACG", "!!!!!!!!!!");
    const fastq record2("Rec2",      "AGTA",   "!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(4, 6, 6, 0, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record1 = record1;
    fastq tmp_record2 = record2;
    ASSERT_EQ(1, truncate_paired_ended_sequences(result, tmp_record1, tmp_record2));
    ASSERT_EQ(fastq("Rec1", "ACGTAGTA", "!!!!!!!!"), tmp_record1);
    ASSERT_EQ(fastq("Rec2", "AGTA", "!!!!"), tmp_record2);
}

TEST(alignment_pe, sequence_b_extends_past_a)
{
    const fastq record1("Rec1",   "CGTA",      "!!!!");
    const fastq record2("Rec2", "ACCGTAGTAT", "!!!!!!!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(-2, 6, 6, 0, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record1 = record1;
    fastq tmp_record2 = record2;
    ASSERT_EQ(1, truncate_paired_ended_sequences(result, tmp_record1, tmp_record2));
    ASSERT_EQ(fastq("Rec1", "CGTA", "!!!!"), tmp_record1);
    ASSERT_EQ(fastq("Rec2", "CGTAGTAT", "!!!!!!!!"), tmp_record2);
}


///////////////////////////////////////////////////////////////////////////////
// Cases 8 PE: Both sequences extends past the other
//                      AAAAAAAAAAAAAAAAAA
//               BBBBBBBBBBBBBBBBBB
// Expected result = Optimal alignment involving all overlapping bases

TEST(alignment_pe, sequences_extend_past_mate)
{
    const fastq record1("Rec1",     "ACGTAGTATACGCT", "!!!!!!!!!!!!!!");
    const fastq record2("Rec2", "GTACACGTAGTATA",     "!!!!!!!!!!!!!!");
    const fastq pcr1("PCR1", "CGCTGA", "!!!!!!");
    const fastq pcr2("PCR2", "TGTAC",  "!!!!!");
    const alignment_info expected = new_aln(-4, 18, 18, 0, 0);
    const alignment_info result = align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record1 = record1;
    fastq tmp_record2 = record2;
    ASSERT_EQ(2, truncate_paired_ended_sequences(result, tmp_record1, tmp_record2));
    ASSERT_EQ(fastq("Rec1", "ACGTAGTATA", "!!!!!!!!!!"), tmp_record1);
    ASSERT_EQ(fastq("Rec2", "ACGTAGTATA", "!!!!!!!!!!"), tmp_record2);
}


TEST(alignment_pe, sequences_extend_past_mate__missing_base__no_shift)
{
    // Test the case where both reads are adapters, but are missing a single base
    // Normally, alignments that do not invovle read1 vs read2 are skipped, but
    // missing bases may cause some alignments to be missed.
    const fastq record1("Rec1", "CCGACC", "!!!!!!");
    const fastq record2("Rec2", "ATGCCT", "!!!!!!");
    const fastq pcr1("PCR1", "CCCGACCCGT", "!!!!!!!!!!");
    const fastq pcr2("PCR2", "AAGATGCCTT", "!!!!!!!!!!");

    // Sub-optimal alignment:
    //   aagatgccttCCGACC
    //          ATGCCTcccgacccgt
    ASSERT_EQ(new_aln(-3, 1, 9, 4, 0), align_paired_ended_sequences(record1, record2, pcr1, pcr2, 0));
    // Optimal alignment, only possible with shift
    //   aagatgccttCCGACC
    //      ATGCCTcccgacccgt
    ASSERT_EQ(new_aln(-7, 11, 13, 1, 0), align_paired_ended_sequences(record1, record2, pcr1, pcr2, 1));
}



///////////////////////////////////////////////////////////////////////////////
// Collapsing of reads


TEST(collapsing, partial_overlap)
{
    fastq record1("Rec1", "ATATTATA", "01234567");
    fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(4);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATAACGT", "01234567EFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, complete_overlap_both_directions)
{
    fastq record1("Rec1", "ATATTATAA", "JJJJJJJJJ");
    fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
    const alignment_info alignment = new_aln(-1);
    ASSERT_EQ(2, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATA", "JJJJJJJJ");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, complete_overlap_mate_1)
{
    fastq record1("Rec1", "ATATTATAG", "JJJJJJJJJ");
    fastq record2("Rec2", "ATATTATA",  "JJJJJJJJ");
    const alignment_info alignment = new_aln(0);
    ASSERT_EQ(1, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATA", "JJJJJJJJ");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, complete_overlap_mate_2)
{
    fastq record1("Rec1", "ATATTATA", "JJJJJJJJ");
    fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
    const alignment_info alignment = new_aln(-1);
    ASSERT_EQ(1, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATA", "JJJJJJJJ");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, unequal_sequence_length__mate_1_shorter)
{
    fastq record1("Rec1", "ATA", "012");
    fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(3);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATANNNNACGT", "012ABCDEFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, unequal_sequence_length__mate_2_shorter)
{
    fastq record1("Rec1", "ATATTATA", "01234567");
    fastq record2("Rec2", "ACG", "EFG");
    const alignment_info alignment = new_aln(8);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATAACG", "01234567EFG");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, ambiguous_sites_are_filled_from_the_mate)
{
    fastq record1("Rec1", "NNNNNNTATA", "0123456789");
    fastq record2("Rec2", "ACGTNNNNNN", "ABCDEFGHIJ");
    const alignment_info alignment;
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ACGTNNTATA", "ABCD!!6789");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__identical_nucleotides)
{
    fastq record1("Rec1", "GCATGATATA", "012345!0:A");
    fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
    const alignment_info alignment = new_aln(6);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "GCATGATATATACAAC", "012345(FBJEFGHIJ");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__identical_nucleotides__scores_are_capped_at_41)
{
    fastq record1("Rec1", "GCATGATATA", "0123456789");
    fastq record2("Rec2", "TATATACAAC", "ABCDEFGHIJ");
    const alignment_info alignment = new_aln(6);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "GCATGATATATACAAC", "012345JJJJEFGHIJ");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__different_nucleotides)
{
    fastq record1("Rec1", "GCATGAGCAT", "012345!0:A");
    fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
    const alignment_info alignment = new_aln(6);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "GCATGATAATTACAAC", "012345(%5%EFGHIJ");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__different_nucleotides__same_quality_1)
{
    const fastq record1("Rec1", "G", "1");
    const fastq record2("Rec2", "T", "1");
    const alignment_info alignment;
    srandom(1);
    const fastq collapsed_expected = fastq("Rec1", "G", "#");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__different_nucleotides__same_quality_2)
{
    const fastq record1("Rec1", "G", "1");
    const fastq record2("Rec2", "T", "1");
    const alignment_info alignment;
    srandom(2);
    const fastq collapsed_expected = fastq("Rec1", "T", "#");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


///////////////////////////////////////////////////////////////////////////////
// Barcode trimming

TEST(trim_barcode, empty_barcode)
{
    const fastq expected("Read", "ACGTAG", "103459");
    fastq record = expected;
    ASSERT_FALSE(truncate_barcode(record, "", 0));
    ASSERT_EQ(expected, record);
}


TEST(trim_barcode, empty_sequence)
{
    const fastq expected("Read", "", "");
    fastq record = expected;
    ASSERT_FALSE(truncate_barcode(record, "ACGTA", 0));
    ASSERT_EQ(expected, record);
}


TEST(trim_barcode, all_empty_all_the_time)
{
    const fastq expected("Read", "", "");
    fastq record = expected;
    ASSERT_FALSE(truncate_barcode(record, "", 0));
    ASSERT_EQ(expected, record);
}


TEST(trim_barcode, barcode_not_found)
{
    const fastq expected("Read", "CATCATACGTAG", "!!!!!!103459");
    fastq record = expected;
    ASSERT_FALSE(truncate_barcode(record, "TGCTGC", 0));
    ASSERT_EQ(expected, record);
}


TEST(trim_barcode, barcode_found)
{
    const fastq expected("Read", "ACGTAG", "103459");
    fastq record("Read", "TGCTGCACGTAG", "!!!!!!103459");
    ASSERT_TRUE(truncate_barcode(record, "TGCTGC", 0));
    ASSERT_EQ(expected, record);
}

TEST(trim_barcode, barcode_found_one_mismatch)
{
    const fastq expected("Read", "ACGTAG", "103459");
    fastq record("Read", "TGCAGCACGTAG", "!!!!!!103459");
    ASSERT_TRUE(truncate_barcode(record, "TGCTGC", 0));
    ASSERT_EQ(expected, record);
}

TEST(trim_barcode, barcode_not_found_two_mismatches)
{
    const fastq expected("Read", "TCCAGCACGTAG", "!!!!!!103459");
    fastq record = expected;
    ASSERT_FALSE(truncate_barcode(record, "TGCTGC", 0));
    ASSERT_EQ(expected, record);
}

TEST(trim_barcode, barcode_shift)
{
    const fastq expected_1("Read", "GCTGCACGTAG", "!!!!!103459");
    const fastq expected_2("Read", "ACGTAG", "103459");
    fastq record = expected_1;
    ASSERT_FALSE(truncate_barcode(record, "TGCTGC", 0));
    ASSERT_EQ(expected_1, record);
    ASSERT_TRUE(truncate_barcode(record, "TGCTGC", 1));
    ASSERT_EQ(expected_2, record);
}
