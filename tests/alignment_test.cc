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
#include <vector>
#include <gtest/gtest.h>

#include "testing.h"
#include "alignment.h"
#include "fastq.h"

namespace ar
{

alignment_info new_aln(int score = 0, int offset = 0, size_t length = 0,
                       size_t nmm = 0, size_t nn = 0, int adapter = 0)
{
    alignment_info aln;
    aln.offset = offset;
    aln.score = score;
    aln.length = length;
    aln.n_mismatches = nmm;
    aln.n_ambiguous = nn;
    aln.adapter_id = adapter;

    return aln;
}

bool operator==(const alignment_info& first, const alignment_info& second)
{
    return (first.offset == second.offset)
        && (first.score == second.score)
        && (first.length == second.length)
        && (first.n_mismatches == second.n_mismatches)
        && (first.n_ambiguous == second.n_ambiguous)
        && (first.adapter_id == second.adapter_id);
}


std::ostream& operator<<(std::ostream& stream, const alignment_info& aln)
{
    stream << "alignment_info(" << aln.score << ", "
                                << aln.offset << ", "
                                << aln.length << ", "
                                << aln.n_mismatches << ", "
                                << aln.n_ambiguous << ", "
                                << aln.adapter_id << ")";
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


fastq_pair_vec create_adapter_vec(const fastq& pcr1, const fastq& pcr2 = fastq())
{
    fastq_pair_vec adapters;
    adapters.push_back(fastq_pair(pcr1, pcr2));
    return adapters;
}


std::random_device g_seed;
std::mt19937 g_rng(g_seed());




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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Rec", "TTTT", "!!!!"));
    const alignment_info expected;
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);
}


TEST(alignment_se, no_expected_overlap)
{
    const fastq record("Rec",  "ACGTAGTA",  "!!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Rec", "TGAGACGGT", "!!!!!!!!!"));
    const alignment_info expected = new_aln(0, 6, 2, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Rec", "AGTAAGGT",  "!!!!!!!!"));
    const alignment_info expected = new_aln(5, 4, 5);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACGT", "1234"), tmp_record);
}


TEST(alignment_se, partial_overlap_with_mismatch)
{
    const fastq record("Rec",  "ACGTAGTAA", "123457890");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Rec", "AGGAAGGT",  "!!!!!!!!"));
    const alignment_info expected = new_aln(3, 4, 5, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACGT", "1234"), tmp_record);
}


TEST(alignment_se, partial_overlap_with_n)
{
    const fastq record("Rec",  "ACGTAGTAA", "123457890");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Rec", "AGNAAGGT",  "!!!!!!!!"));
    const alignment_info expected = new_aln(4, 4, 5, 0, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(record);
    const alignment_info expected = new_aln(8, 0, 8);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, completely_overlapping_sequences_with_1_mismatch)
{
    const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Rec", "GCGTAGTA", "!!!!!!!!"));
    const alignment_info expected = new_aln(6, 0, 8, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, completely_overlapping_sequences_with_1_mismatch_and_1_n)
{
    const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
    const fastq_pair_vec adapter = create_adapter_vec(fastq("Rec", "GCGTAGTN", "!!!!!!!!"));
    const alignment_info expected = new_aln(5, 0, 8, 1, 1);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp",   "TAGTA", "!!!!!"));
    const alignment_info expected = new_aln(5, 3, 5);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACG", "ABC"), tmp_record);
}


TEST(alignment_se, sequence_a_contains_b__with_1_mismatch)
{
    const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "TATTA", "!!!!!"));
    const alignment_info expected = new_aln(3, 3, 5, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACG", "ABC"), tmp_record);
}


TEST(alignment_se, sequence_a_contains_b__with_1_n)
{
    const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "TAGNA", "!!!!!"));
    const alignment_info expected = new_aln(4, 3, 5, 0, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACG", "ABC"), tmp_record);
}


TEST(alignment_se, sequence_b_contains_a)
{
    const fastq record("Rec",  "ACGT", "!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "ACGTAGTA", "!!!!!!!!"));
    const alignment_info expected = new_aln(4, 0, 4);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, sequence_b_contains_a__with_1_mismatch)
{
    const fastq record("Rec",  "ACGT",     "!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "GCGTAGTA", "!!!!!!!!"));
    const alignment_info expected = new_aln(2, 0, 4, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


TEST(alignment_se, sequence_b_contains_a__with_1_n)
{
    const fastq record("Rec",  "ACGT", "!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "ACGNAGTA", "!!!!!!!!"));
    const alignment_info expected = new_aln(3, 0, 4, 0, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
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
    const fastq record("Rec", "ACGTAGTATA", "0123456789");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "AGTA", "!!!!"));
    const alignment_info expected = new_aln(4, 4, 4);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "ACGT", "0123"), tmp_record);
}


TEST(alignment_se, sequence_b_extends_past_a__no_shift)
{
    const fastq record("Rec", "CGTA", "#!%%");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "ACGTAGTATA", "!!!!!!!!!!"));
    const alignment_info expected = new_aln(1, 3, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "CGT", "#!%"), tmp_record);
}


TEST(alignment_se, sequence_b_extends_past_a__shift_of_1)
{
    const fastq record("Rec", "CGTA", "#!%%");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "ACGTAGTATA", "!!!!!!!!!!"));
    const alignment_info expected = new_aln(4, -1, 4, 0, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapters, 1);
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
    const fastq record("Rec", "ACGTAGTATATAGT", "!!!!!!!!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "CCGAACGTAGTATA", "!!!!!!!!!!!!!!"));
    const alignment_info expected = new_aln(10, -4, 10, 0, 0, 0);
    const alignment_info result = align_single_ended_sequence(record, adapters, 4);
    ASSERT_EQ(expected, result);

    fastq tmp_record = record;
    truncate_single_ended_sequence(result, tmp_record);
    ASSERT_EQ(fastq("Rec", "", ""), tmp_record);
}


///////////////////////////////////////////////////////////////////////////////
// Empty sequence or adapter

TEST(alignment_se, empty_sequence)
{
    const fastq record;
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "CCGAACGTAGTATA", "!!!!!!!!!!!!!!"));
    const alignment_info expected;
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);
}


TEST(alignment_se, empty_adapter)
{
    const fastq record("Rec", "ACGTAGTATATAGT", "!!!!!!!!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq());
    const alignment_info expected;
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);
}


///////////////////////////////////////////////////////////////////////////////
// Misc

TEST(alignment_se, shift_is_lower_than_possible)
{
    const fastq record("Rec",  "AAAA", "!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("Adp", "TTTT", "!!!!"));
    const alignment_info expected;
    const alignment_info result = align_single_ended_sequence(record, adapters, -10);
    ASSERT_EQ(expected, result);
}


TEST(alignment_se, only_adapter_1_is_used)
{
    const fastq_pair_vec adapters = create_adapter_vec(fastq("barcode", "AAA", "JJJ"),
                                                       fastq("barcode", "TTTAAA", "JJJJJJ"));
    const fastq record("Rec",  "CCCCTTTAAA", "0987654321");
    const alignment_info expected = new_aln(3, 7, 3);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);
}


TEST(alignment_se, prefer_best_alignement__first)
{
    fastq_pair_vec adapters;
    adapters.push_back(fastq_pair(fastq("adapter", "TGCTGC", "JJJJJJ"), fastq()));
    adapters.push_back(fastq_pair(fastq("adapter", "TGCTGA", "JJJJJJ"), fastq()));

    const fastq record("Read", "TAGTCGCTATGCTGC", "!!!!!!!!!103459");
    const alignment_info expected = new_aln(6, 9, 6);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);
}


TEST(alignment_se, prefer_best_alignement__second)
{
    fastq_pair_vec adapters;
    adapters.push_back(fastq_pair(fastq("adapter", "TGCTGA", "JJJJJJ"), fastq()));
    adapters.push_back(fastq_pair(fastq("adapter", "TGCTGC", "JJJJJJ"), fastq()));

    const fastq record("Read", "TAGTCGCTATGCTGC", "!!!!!!!!!103459");
    const alignment_info expected = new_aln(6, 9, 6, 0, 0, 1);
    const alignment_info result = align_single_ended_sequence(record, adapters, 0);
    ASSERT_EQ(expected, result);
}


TEST(alignment_se, prefer_best_alignement__neither)
{
    fastq_pair_vec barcodes;
    barcodes.push_back(fastq_pair(fastq("barcode", "AAAAAA", "JJJJJJ"), fastq()));
    barcodes.push_back(fastq_pair(fastq("barcode", "CCCCCC", "JJJJJJ"), fastq()));

    const fastq record = fastq("Read", "AACTGTACGTAGTT", "!!!!!!10345923");
    const alignment_info result = align_single_ended_sequence(record, barcodes, 0);
    ASSERT_EQ(alignment_info(), result);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected;
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}


TEST(alignment_pe, no_expected_overlap)
{
    const fastq record1("Rec",  "ACGTAGTA",  "!!!!!!!!");
    const fastq record2("Rec", "TGAGACGGT", "!!!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(0, 6, 2, 1);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(5, 4, 5);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(8, 0, 8);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(5, 3, 5);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}


TEST(alignment_pe, sequence_b_contains_a)
{
    const fastq record1("Rec1", "ACGT", "!!!!");
    const fastq record2("Rec2", "ACGTAGTA", "!!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(4, 0, 4);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(6, 4, 6);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(6, -2, 6);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
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
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(18, -4, 18);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
    ASSERT_EQ(expected, result);

    fastq tmp_record1 = record1;
    fastq tmp_record2 = record2;
    ASSERT_EQ(2, truncate_paired_ended_sequences(result, tmp_record1, tmp_record2));
    ASSERT_EQ(fastq("Rec1", "ACGTAGTATA", "!!!!!!!!!!"), tmp_record1);
    ASSERT_EQ(fastq("Rec2", "ACGTAGTATA", "!!!!!!!!!!"), tmp_record2);
}


TEST(alignment_pe, only_adadapter_sequence)
{
    const fastq record1("Rec1", "CCCGAC", "!!!!!!");
    const fastq record2("Rec2", "ATGCCTT", "!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CCCGACCCGT", "!!!!!!!!!!"),
                                                       fastq("PCR2", "AAGATGCCTT", "!!!!!!!!!!"));

    ASSERT_EQ(new_aln(13, -7, 13), align_paired_ended_sequences(record1, record2, adapters, 0));
}


TEST(alignment_pe, only_adadapter_sequence__missing_base__shift)
{
    // Test the case where both reads are adapters, but are missing a single base
    // Normally, alignments that do not invovle read1 vs read2 are skipped, but
    // missing bases may cause some alignments to be missed.
    const fastq record1("Rec1", "CCGACC", "!!!!!!");
    const fastq record2("Rec2", "ATGCCT", "!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CCCGACCCGT", "!!!!!!!!!!"),
                                                       fastq("PCR2", "AAGATGCCTT", "!!!!!!!!!!"));

    // Sub-optimal alignment:
    //   aagatgccttCCGACC
    //          ATGCCTcccgacccgt
    ASSERT_EQ(new_aln( 1, -3, 9, 4), align_paired_ended_sequences(record1, record2, adapters, 0));
    // Optimal alignment, only possible with shift
    //   aagatgccttCCGACC
    //      ATGCCTcccgacccgt
    ASSERT_EQ(new_aln(11, -7, 13, 1), align_paired_ended_sequences(record1, record2, adapters, 1));
}


///////////////////////////////////////////////////////////////////////////////

TEST(alignment_pe, empty_mate_1)
{
    const fastq record1("Rec", "ACGTAGTAA", "!!!!!!!!!");
    const fastq record2("Rec", "AGTAAGGT",  "!!!!!!!!");
    const fastq_pair_vec adapters = create_adapter_vec(fastq("PCR1", "CGCTGA", "!!!!!!"),
                                                       fastq("PCR2", "TGTAC",  "!!!!!"));
    const alignment_info expected = new_aln(5, 4, 5);
    const alignment_info result = align_paired_ended_sequences(record1, record2, adapters, 0);
    ASSERT_EQ(expected, result);
    ASSERT_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);

}


///////////////////////////////////////////////////////////////////////////////
// Collapsing of reads

TEST(collapsing, partial_overlap)
{
    fastq record1("Rec1", "ATATTATA", "01234567");
    fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(0, 4);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATAACGT", "01234567EFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, complete_overlap_both_directions)
{
    fastq record1("Rec1",  "ATATTATAA", "JJJJJJJJJ");
    fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
    const alignment_info alignment = new_aln(0, -1);
    ASSERT_EQ(2, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATA", "wwwwwwww", FASTQ_ENCODING_SAM);
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, complete_overlap_mate_1)
{
    fastq record1("Rec1", "ATATTATAG", "JJJJJJJJJ");
    fastq record2("Rec2", "ATATTATA",  "JJJJJJJJ");
    const alignment_info alignment = new_aln();
    ASSERT_EQ(1, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATA", "wwwwwwww", FASTQ_ENCODING_SAM);
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, complete_overlap_mate_2)
{
    fastq record1("Rec1",  "ATATTATA", "JJJJJJJJ");
    fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
    const alignment_info alignment = new_aln(0, -1);
    ASSERT_EQ(1, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATA", "wwwwwwww", FASTQ_ENCODING_SAM);
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, unequal_sequence_length__mate_1_shorter)
{
    fastq record1("Rec1", "ATA", "012");
    fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(0, 3);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATANNNNACGT", "012ABCDEFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, unequal_sequence_length__mate_1_shorter__mate_2_extends_past)
{
    fastq record1("Rec1", "ATA", "012");
    fastq record2("Rec2", "AANNNNACGT", "90ABCDEFGH");
    const alignment_info alignment = new_aln(0, -2);
    ASSERT_EQ(1, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATANACGT", "012DEFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, unequal_sequence_length__mate_2_shorter)
{
    fastq record1("Rec1", "ATATTATA", "01234567");
    fastq record2("Rec2", "ACG", "EFG");
    const alignment_info alignment = new_aln(0, 8);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ATATTATAACG", "01234567EFG");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, ambiguous_sites_are_filled_from_the_mate)
{
    fastq record1("Rec1", "NNNNNNTATA", "0123456789");
    fastq record2("Rec2", "ACGTNNNNNN", "ABCDEFGHIJ");
    const alignment_info alignment;
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "ACGTNNTATA", "ABCD!!6789");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__identical_nucleotides)
{
    fastq record1("Rec1", "GCATGATATA", "012345!0:A");
    fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
    const alignment_info alignment = new_aln(0, 6);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "GCATGATATATACAAC", "012345(FBcEFGHIJ", FASTQ_ENCODING_SAM);
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__identical_nucleotides__scores_are_capped_at_41)
{
    fastq record1("Rec1", "GCATGATATA", "0123456789");
    fastq record2("Rec2", "TATATACAAC", "ABCDEFGHIJ");
    const alignment_info alignment = new_aln(0, 6);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "GCATGATATATACAAC", "012345Z\\^`EFGHIJ", FASTQ_ENCODING_SAM);
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__different_nucleotides)
{
    fastq record1("Rec1", "GCATGAGCAT", "012345!0:A");
    fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
    const alignment_info alignment = new_aln(0, 6);
    ASSERT_EQ(0, truncate_paired_ended_sequences(alignment, record1, record2));
    const fastq collapsed_expected = fastq("Rec1", "GCATGATAATTACAAC", "012345(%5%EFGHIJ");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__different_nucleotides__same_quality_1)
{
    const fastq record1("Rec1", "G", "1");
    const fastq record2("Rec2", "T", "1");
    const alignment_info alignment;
    std::seed_seq seed{1};
    std::mt19937 rng(seed);

    const fastq collapsed_expected = fastq("Rec1", "G", "#");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, consensus_bases__different_nucleotides__same_quality_2)
{
    const fastq record1("Rec1", "G", "1");
    const fastq record2("Rec2", "T", "1");
    const alignment_info alignment;
    std::seed_seq seed{2};
    std::mt19937 rng(seed);

    const fastq collapsed_expected = fastq("Rec1", "T", "#");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, rng);
    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, offset_past_the_end)
{
    const fastq record1("Rec1", "G", "1");
    const fastq record2("Rec2", "T", "1");
    const alignment_info alignment = new_aln(0, 2);
    ASSERT_THROW(collapse_paired_ended_sequences(alignment, record1, record2, g_rng), std::invalid_argument);
}


TEST(collapsing, partial_overlap__mate_numbers_removed)
{
    const fastq record1("Read/1", "ATATTATA", "01234567");
    const fastq record2("Read/2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(0, 4);
    const fastq collapsed_expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);

    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, partial_overlap__mate_numbers_removed__mate_1_meta_kept)
{
    const fastq record1("Read/1 Meta1", "ATATTATA", "01234567");
    const fastq record2("Read/2 Meta2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(0, 4);
    const fastq collapsed_expected = fastq("Read Meta1", "ATATTATAACGT", "01234567EFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);

    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, partial_overlap__mate_numbers_removed__non_std_mate_sep)
{
    const fastq record1("Read:1", "ATATTATA", "01234567");
    const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(0, 4);
    const fastq collapsed_expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng, ':');

    ASSERT_EQ(collapsed_expected, collapsed_result);
}


TEST(collapsing, partial_overlap__mate_numbers_removed__non_std_mate_sep__not_set)
{
    const fastq record1("Read:1", "ATATTATA", "01234567");
    const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
    const alignment_info alignment = new_aln(0, 4);
    const fastq collapsed_expected = fastq("Read:1", "ATATTATAACGT", "01234567EFGH");
    const fastq collapsed_result = collapse_paired_ended_sequences(alignment, record1, record2, g_rng);

    ASSERT_EQ(collapsed_expected, collapsed_result);
}


///////////////////////////////////////////////////////////////////////////////
// Barcode extraction

TEST(extract_adapter_sequences, empty_sequences)
{
    const fastq expected_1 = fastq("read1", "", "");
    const fastq expected_2 = fastq("read2", "", "");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(alignment_info(), read1, read2);
    ASSERT_EQ(expected_1, read1);
    ASSERT_EQ(expected_2, read2);
}


TEST(extract_adapter_sequences, empty_sequence_1)
{
    const fastq expected_1 = fastq("read1", "", "");
    const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
    fastq read1 = expected_1;
    fastq read2 = fastq("read2", "GGGGCC", "!!!!!!");
    extract_adapter_sequences(alignment_info(), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, empty_sequence_2)
{
    const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
    const fastq expected_2 = fastq("read2", "", "");
    fastq read1 = fastq("read1", "", "");
    fastq read2 = expected_2;
    extract_adapter_sequences(alignment_info(), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, empty_sequence_3)
{
    const fastq expected_1 = fastq("read1", "", "");
    const fastq expected_2 = fastq("read2", "", "");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(alignment_info(), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, case_1_no_alignment)
{
    const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
    const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(alignment_info(), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, case_2_partial_overlap)
{
    const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
    const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(new_aln(0, 2), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);

}


TEST(extract_adapter_sequences, case_3_complete_overlap)
{
    const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
    const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(new_aln(), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, case_4_read1_contains_read2)
{
    const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
    const fastq expected_2 = fastq("read2", "GGCC", "!!!!");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(new_aln(0, 2), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, case_5_read2_contains_read1)
{
    const fastq expected_1 = fastq("read1", "AATT", "!!!!");
    const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(new_aln(), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, case_6_read1_extends_past_read2)
{
    const fastq expected_1 = fastq("read1", "AATTTTCC", "12345678");
    const fastq expected_2 = fastq("read2", "GGGGGG", "!!!!!!");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(new_aln(), read1, read2);
    ASSERT_EQ(fastq("read1", "CC", "78"), read1);
    ASSERT_EQ(fastq("read2", "", ""), read2);
}


TEST(extract_adapter_sequences, case_7_read2_extends_past_read1)
{
    const fastq expected_1 = fastq("read1", "TTTTTT", "!!!!!!");
    const fastq expected_2 = fastq("read2", "AAGGGGGG", "12345678");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(new_aln(0, -2), read1, read2);
    ASSERT_EQ(fastq("read1", "", ""), read1);
    ASSERT_EQ(fastq("read2", "AA", "12"), read2);
}


TEST(extract_adapter_sequences, case_8_reads_extends_past_each_pther)
{
    const fastq expected_1 = fastq("read1", "TTTTTTCCC", "ABCDEFGHI");
    const fastq expected_2 = fastq("read2", "AAGGGGGG",  "12345678");
    fastq read1 = expected_1;
    fastq read2 = expected_2;
    extract_adapter_sequences(new_aln(0, -2), read1, read2);
    ASSERT_EQ(fastq("read1", "CCC", "GHI"), read1);
    ASSERT_EQ(fastq("read2", "AA", "12"), read2);
}


///////////////////////////////////////////////////////////////////////////////
// Brute-force checking of alignment calculations
// Simply check all combinations involving 3 bases varying, for a range of
// sequence lengths to help catch corner cases with the optimizations

// The function is not exposed, so a declaration is required
bool compare_subsequences(const alignment_info& best, alignment_info& current,
                          const char* seq_1_ptr, const char* seq_2_ptr);


/** Naive reimplementation of alignment calculation. **/
void update_alignment(alignment_info& aln,
                      const std::string& a,
                      const std::string& b,
                      size_t nbases)
{
    if (a.length() != b.length()) {
        throw std::invalid_argument("length does not match");
    }

    for (size_t i = 0; i < nbases; ++i) {
        const char nt1 = a.at(i);
        const char nt2 = b.at(i);

        if (nt1 == 'N' || nt2 == 'N') {
            aln.n_ambiguous++;
        } else if (nt1 == nt2) {
            aln.score++;
        } else {
            aln.n_mismatches++;
            aln.score--;
        }
    }
}


/** Returns all 3 nt combinations of the bases ACGTN. **/
std::vector<std::string> get_combinations()
{
    std::vector<std::string> result;
    const std::string nts = "ACGTN";
    for (size_t i = 0; i < nts.length(); ++i) {
        for (size_t j = 0; j < nts.length(); ++j) {
            for (size_t k = 0; k < nts.length(); ++k) {
                std::string combination(3, 'A');
                combination.at(0) = nts.at(i);
                combination.at(1) = nts.at(j);
                combination.at(2) = nts.at(k);
                result.push_back(combination);
            }
        }
    }

    return result;
}


TEST(compare_subsequences, brute_force_validation)
{
    const alignment_info best;
    const std::vector<std::string> combinations = get_combinations();
    for (size_t seqlen = 10; seqlen <= 20; ++seqlen) {
        for (size_t pos = 0; pos < seqlen; ++pos) {
            const size_t nbases = std::min<int>(3, seqlen - pos);

            for (size_t i = 0; i < combinations.size(); ++i) {
                for (size_t j = 0; j < combinations.size(); ++j) {
                    alignment_info expected;
                    expected.length = seqlen;
                    expected.score = seqlen - nbases;
                    update_alignment(expected, combinations.at(i), combinations.at(j), nbases);

                    std::string mate1 = std::string(seqlen, 'A');
                    mate1.replace(pos, nbases, combinations.at(i).substr(0, nbases));
                    std::string mate2 = std::string(seqlen, 'A');
                    mate2.replace(pos, nbases, combinations.at(j).substr(0, nbases));

                    alignment_info current;
                    current.length = seqlen;
                    compare_subsequences(best, current, mate1.c_str(), mate2.c_str());

                    if (!(expected == current)) {
                        std::cerr << "seqlen = " << seqlen << "\n"
                                  << "pos    = " << pos << "\n"
                                  << "nbases = " << nbases << "\n"
                                  << "mate1  = " << mate1 << "\n"
                                  << "mate2  = " << mate2 << std::endl;
                        ASSERT_EQ(expected, current);
                    }
                }
            }
        }
    }
}

} // namespace ar
