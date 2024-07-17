/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "alignment.hpp"   // for alignment_info, sequence_merger, extract_...
#include "commontypes.hpp" // for merge_strategy, merge_strategy::determini...
#include "debug.hpp"       // for assert_failed
#include "errors.hpp"      // for assert_failed
#include "fastq.hpp"       // for fastq, fastq_pair_vec, fastq_pair
#include "fastq_enc.hpp"   // for FASTQ_ENCODING_SAM
#include "simd.hpp"        // for size_t, instruction_set, supported, get_c...
#include "testing.hpp"     // for catch.hpp, StringMaker
#include <algorithm>       // for min
#include <cstddef>         // for size_t
#include <cstdint>         // for int64_t
#include <random>          // for mt19937, seed_seq, random_device
#include <sstream>         // for operator<<, ostream, basic_ostream, char_...
#include <string>          // for string, basic_string, operator<<
#include <vector>          // for vector

namespace adapterremoval {

struct ALN;

#define TEST_ALIGNMENT_SETTER(TYPE, NAME)                                      \
  ALN& NAME(TYPE value)                                                        \
  {                                                                            \
    info.adapter_id = 0;                                                       \
    info.NAME = value;                                                         \
    return *this;                                                              \
  }

// Parameterize tests over supported SIMD instruction sets
#define PARAMETERIZE_IS GENERATE(from_range(simd::supported()))

bool
operator==(const alignment_info& first, const alignment_info& second)
{
  return (first.offset == second.offset) && (first.length == second.length) &&
         (first.n_mismatches == second.n_mismatches) &&
         (first.n_ambiguous == second.n_ambiguous) &&
         (first.adapter_id == second.adapter_id);
}

struct ALN
{
  ALN()
    : info()
  {
  }

  TEST_ALIGNMENT_SETTER(int, offset);
  TEST_ALIGNMENT_SETTER(size_t, length);
  TEST_ALIGNMENT_SETTER(size_t, n_mismatches);
  TEST_ALIGNMENT_SETTER(size_t, n_ambiguous);
  TEST_ALIGNMENT_SETTER(int, adapter_id);

  operator alignment_info() { return info; }

  alignment_info info;
};

bool
operator==(const alignment_info& first, const ALN& second)
{
  return first == second.info;
}

fastq_pair_vec
create_adapter_vec(const fastq& adapter1, const fastq& adapter2 = fastq())
{
  fastq_pair_vec adapters;
  adapters.push_back(fastq_pair(adapter1, adapter2));
  return adapters;
}

alignment_info
align_single_ended_sequence(const fastq& read,
                            const fastq_pair_vec& adapters,
                            int max_shift)
{
  return sequence_aligner(adapters, PARAMETERIZE_IS)
    .align_single_end(read, max_shift);
}

alignment_info
align_paired_ended_sequences(const fastq& read1,
                             const fastq& read2,
                             const fastq_pair_vec& adapters,
                             int max_shift)
{
  return sequence_aligner(adapters, PARAMETERIZE_IS)
    .align_paired_end(read1, read2, max_shift);
}

void
truncate_single_ended_sequence(const alignment_info& alignment, fastq& read)
{
  alignment.truncate_single_end(read);
}

size_t
truncate_paired_ended_sequences(const alignment_info& alignment,
                                fastq& read1,
                                fastq& read2)
{
  return alignment.truncate_paired_end(read1, read2);
}

void
REQUIRE_TRUNCATED_PE_IS_UNCHANGED(const alignment_info& alignment,
                                  const fastq& record1,
                                  const fastq& record2)
{
  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(
    truncate_paired_ended_sequences(alignment, tmp_record1, tmp_record2) == 0);
  REQUIRE(tmp_record1 == record1);
  REQUIRE(tmp_record2 == record2);
}

std::random_device g_seed;
std::mt19937 g_rng_instance(g_seed());
std::mt19937* g_rng(&g_rng_instance);

///////////////////////////////////////////////////////////////////////////////

TEST_CASE("alignment_info::score")
{
  REQUIRE(ALN().info.score() == 0);
  REQUIRE(ALN().length(13).info.score() == 13);
  REQUIRE(ALN().length(13).n_ambiguous(5).info.score() == 8);
  REQUIRE(ALN().length(13).n_mismatches(5).info.score() == 3);
  REQUIRE(ALN().length(13).n_ambiguous(5).info.score() == 8);
  REQUIRE(ALN().length(13).n_ambiguous(6).n_mismatches(5).info.score() == -3);
}

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

TEST_CASE("SE: Unalignable sequence yields default alignment",
          "[alignment::single_end]")
{
  const fastq record("Rec", "AAAA", "!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Rec", "TTTT", "!!!!"));

  REQUIRE(align_single_ended_sequence(record, adapters, 0) == ALN());
}

TEST_CASE("SE: Random sequences yields suboptimal alignment",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Rec", "TGAGACGGT", "!!!!!!!!!"));

  REQUIRE(align_single_ended_sequence(record, adapters, 0) ==
          ALN().offset(6).length(2).n_mismatches(1));
}

///////////////////////////////////////////////////////////////////////////////
// Case 2: Partial overlap between sequences:
//         AAAAAAAAAAA
//                BBBBBBBBBBBBBBB
// Expected result = optimal alignment between ends
TEST_CASE("SE: Partial alignment between ends", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTAA", "123457890");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Rec", "AGTAAGGT", "!!!!!!!!"));
  const alignment_info expected = ALN().offset(4).length(5);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACGT", "1234"));
}

TEST_CASE("SE: Partial alignment with mismatches between ends",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTAA", "123457890");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Rec", "AGGAAGGT", "!!!!!!!!"));
  const alignment_info expected = ALN().offset(4).length(5).n_mismatches(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACGT", "1234"));
}

TEST_CASE("SE: Partial alignment with ambiguous between ends",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTAA", "123457890");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Rec", "AGNAAGGT", "!!!!!!!!"));
  const alignment_info expected = ALN().offset(4).length(5).n_ambiguous(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACGT", "1234"));
}

///////////////////////////////////////////////////////////////////////////////
// Case 3: Complete overlap between sequences:
//         AAAAAAAAAAA
//         BBBBBBBBBBB
// Expected result = Optimal alignment involving all bases

TEST_CASE("SE: Completely overlapping sequences", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(record);
  const alignment_info expected = ALN().length(8);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("SE: Completely overlapping sequences with mismatches",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Rec", "GCGTAGTA", "!!!!!!!!"));
  const alignment_info expected = ALN().length(8).n_mismatches(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("SE: Completely overlapping sequences with mismatches and ambiguous",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
  const fastq_pair_vec adapter =
    create_adapter_vec(fastq("Rec", "GCGTAGTN", "!!!!!!!!"));
  const alignment_info expected =
    ALN().length(8).n_mismatches(1).n_ambiguous(1);
  const alignment_info result = align_single_ended_sequence(record, adapter, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

///////////////////////////////////////////////////////////////////////////////
// Cases 4 and 5: Sequence A completely contains sequence B (and vice versa)
//         AAAAAAAAAAA  AAAA
//               BBBBB  BBBBBBBBBB
// Expected result = Optimal alignment involving all bases of the shortest read

TEST_CASE("Complete adapter inside sequence", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "TAGTA", "!!!!!"));
  const alignment_info expected = ALN().offset(3).length(5);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACG", "ABC"));
}

TEST_CASE("Complete adapter inside sequence with mismatch",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "TATTA", "!!!!!"));
  const alignment_info expected = ALN().offset(3).length(5).n_mismatches(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACG", "ABC"));
}

TEST_CASE("Complete adapter inside sequence with ambiguous",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "TAGNA", "!!!!!"));
  const alignment_info expected = ALN().offset(3).length(5).n_ambiguous(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACG", "ABC"));
}

TEST_CASE("Complete sequence inside adapter", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGT", "!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "ACGTAGTA", "!!!!!!!!"));
  const alignment_info expected = ALN().length(4);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("Complete sequence inside adapter with mismatches",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGT", "!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "GCGTAGTA", "!!!!!!!!"));
  const alignment_info expected = ALN().length(4).n_mismatches(1);
  const auto result = align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("Complete sequence inside adapter with ambiguous",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGT", "!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "ACGNAGTA", "!!!!!!!!"));
  const alignment_info expected = ALN().length(4).n_ambiguous(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

///////////////////////////////////////////////////////////////////////////////
// Cases 6 and 7: Sequence A extends past sequence B (and vice versa)
//         AAAAAAAAAAAAAA    AAAA
//               BBBBB     BBBBBBBBBBBBB

TEST_CASE("Sequence extends past adapter", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTATA", "0123456789");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "AGTA", "!!!!"));
  const alignment_info expected = ALN().offset(4).length(4);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACGT", "0123"));
}

TEST_CASE("Sequence extends past adapter, no shift", "[alignment::single_end]")
{
  const fastq record("Rec", "CGTA", "#!%%");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "ACGTAGTATA", "!!!!!!!!!!"));
  const alignment_info expected = ALN().offset(3).length(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "CGT", "#!%"));
}

TEST_CASE("Sequence extends past adapter, shift 1", "[alignment::single_end]")
{
  const fastq record("Rec", "CGTA", "#!%%");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "ACGTAGTATA", "!!!!!!!!!!"));
  const alignment_info expected = ALN().offset(-1).length(4);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 1);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

///////////////////////////////////////////////////////////////////////////////
// Cases 8: Both sequences extends past the other
//                      AAAAAAAAAAAAAAAAAA
//               BBBBBBBBBBBBBBBBBB
// Expected result = Optimal alignment involving all overlapping bases

TEST_CASE("Sequence and adapter extends past each other",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTATATAGT", "!!!!!!!!!!!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "CCGAACGTAGTATA", "!!!!!!!!!!!!!!"));
  const alignment_info expected = ALN().offset(-4).length(10);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 4);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  truncate_single_ended_sequence(result, tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

///////////////////////////////////////////////////////////////////////////////
// Empty sequence or adapter

TEST_CASE("Empty sequence alignment", "[alignment::single_end]")
{
  const fastq record;
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "CCGAACGTAGTATA", "!!!!!!!!!!!!!!"));
  const alignment_info expected;
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Empty adapter alignment", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTATATAGT", "!!!!!!!!!!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(fastq());
  const alignment_info expected;
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Max error rates

TEST_CASE("Longest valid alignment is returned", "[alignment::single_end]")
{
  const fastq record("Rec", "AAATAAAAA");
  const auto adapters = create_adapter_vec(fastq("Rec", "AAAAAAAAA"));

  sequence_aligner aligner(adapters, PARAMETERIZE_IS);
  const auto result = aligner.align_single_end(record, 0);
  const alignment_info expected = ALN().length(9).n_mismatches(1);

  REQUIRE(result == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Misc

TEST_CASE("Lower than possible shift is allowed", "[alignment::single_end]")
{
  const fastq record("Rec", "AAAA", "!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("Adp", "TTTT", "!!!!"));
  const alignment_info expected;
  const alignment_info result =
    align_single_ended_sequence(record, adapters, -10);
  REQUIRE(result == expected);
}

TEST_CASE("Only adapter 1 is used", "[alignment::single_end]")
{
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("barcode", "AAA", "JJJ"), fastq("barcode", "TTTAAA", "JJJJJJ"));
  const fastq record("Rec", "CCCCTTTAAA", "0987654321");
  const alignment_info expected = ALN().offset(7).length(3);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Best matching adapter is returned: First", "[alignment::single_end]")
{
  fastq_pair_vec adapters;
  adapters.push_back(fastq_pair(fastq("adapter", "TGCTGC", "JJJJJJ"), fastq()));
  adapters.push_back(fastq_pair(fastq("adapter", "TGCTGA", "JJJJJJ"), fastq()));

  const fastq record("Read", "TAGTCGCTATGCTGC", "!!!!!!!!!103459");
  const alignment_info expected = ALN().offset(9).length(6);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Best matching adapter returned: Second", "[alignment::single_end]")
{
  fastq_pair_vec adapters;
  adapters.push_back(fastq_pair(fastq("adapter", "TGCTGA", "JJJJJJ"), fastq()));
  adapters.push_back(fastq_pair(fastq("adapter", "TGCTGC", "JJJJJJ"), fastq()));

  const fastq record("Read", "TAGTCGCTATGCTGC", "!!!!!!!!!103459");
  const alignment_info expected = ALN().offset(9).length(6).adapter_id(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Best matching adapter returned: Neither", "[alignment::single_end]")
{
  fastq_pair_vec barcodes;
  barcodes.push_back(fastq_pair(fastq("barcode", "AAAAAA", "JJJJJJ"), fastq()));
  barcodes.push_back(fastq_pair(fastq("barcode", "CCCCCC", "JJJJJJ"), fastq()));

  const fastq record = fastq("Read", "AACTGTACGTAGTT", "!!!!!!10345923");
  const alignment_info result =
    align_single_ended_sequence(record, barcodes, 0);
  REQUIRE(result == ALN());
}

///////////////////////////////////////////////////////////////////////////////
// Case 1 PE: No overlap between sequences:
//         AAAAAAAAAAA
//                       BBBBBBBBBBBBBBB
// Expected result = suboptimal alignment or no overlap

TEST_CASE("Unalignable sequence pair", "[alignment::paired_end]")
{
  const fastq record1("Rec", "AAAA", "!!!!");
  const fastq record2("Rec", "TTTT", "!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected;
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

TEST_CASE("No overlap in sequence pair", "[alignment::paired_end]")
{
  const fastq record1("Rec", "ACGTAGTA", "!!!!!!!!");
  const fastq record2("Rec", "TGAGACGGT", "!!!!!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().offset(6).length(2).n_mismatches(1);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

///////////////////////////////////////////////////////////////////////////////
// Case 2 PE: Partial overlap between sequences:
//         AAAAAAAAAAA
//                BBBBBBBBBBBBBBB
// Expected result = optimal alignment between ends

TEST_CASE("Partial overlap in sequence pair", "[alignment::paired_end]")
{
  const fastq record1("Rec", "ACGTAGTAA", "!!!!!!!!!");
  const fastq record2("Rec", "AGTAAGGT", "!!!!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().offset(4).length(5);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

///////////////////////////////////////////////////////////////////////////////
// Case 3 PE: Complete overlap between sequences:
//         AAAAAAAAAAA
//         BBBBBBBBBBB
// Expected result = Optimal alignment involving all bases

TEST_CASE("Completely overlapping sequence pair", "[alignment::paired_end]")
{
  const fastq record1("Rec", "ACGTAGTA", "!!!!!!!!");
  const fastq record2 = record1;
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().length(8);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

///////////////////////////////////////////////////////////////////////////////
// Cases 4 and 5 PE: Sequence A completely contains sequence B (and vice versa)
//         AAAAAAAAAAA  AAAA
//               BBBBB  BBBBBBBBBB
// Expected result = Optimal alignment involving all bases of the shortest read

TEST_CASE("Sequence A contains sequence B", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "ACGTAGTA", "!!!!!!!!");
  const fastq record2("Rec2", "TAGTA", "!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().offset(3).length(5);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

TEST_CASE("Sequence B contains sequence A", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "ACGT", "!!!!");
  const fastq record2("Rec2", "ACGTAGTA", "!!!!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().length(4);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

///////////////////////////////////////////////////////////////////////////////
// Cases 6 and 7 PE: Sequence A extends past sequence B (and vice versa)
//         AAAAAAAAAAAAAA    AAAA
//               BBBBB     BBBBBBBBBBBBB
TEST_CASE("Sequence A extends past sequence B", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "ACGTAGTACG", "!!!!!!!!!!");
  const fastq record2("Rec2", "AGTA", "!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().offset(4).length(6);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(truncate_paired_ended_sequences(result, tmp_record1, tmp_record2) ==
          1);
  REQUIRE(tmp_record1 == fastq("Rec1", "ACGTAGTA", "!!!!!!!!"));
  REQUIRE(tmp_record2 == fastq("Rec2", "AGTA", "!!!!"));
}

TEST_CASE("Sequence B extends past sequence A", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "CGTA", "!!!!");
  const fastq record2("Rec2", "ACCGTAGTAT", "!!!!!!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().offset(-2).length(6);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(truncate_paired_ended_sequences(result, tmp_record1, tmp_record2) ==
          1);
  REQUIRE(tmp_record1 == fastq("Rec1", "CGTA", "!!!!"));
  REQUIRE(tmp_record2 == fastq("Rec2", "CGTAGTAT", "!!!!!!!!"));
}

///////////////////////////////////////////////////////////////////////////////
// Cases 8 PE: Both sequences extends past the other
//                      AAAAAAAAAAAAAAAAAA
//               BBBBBBBBBBBBBBBBBB
// Expected result = Optimal alignment involving all overlapping bases

TEST_CASE("Sequences extends past each other", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "ACGTAGTATACGCT", "!!!!!!!!!!!!!!");
  const fastq record2("Rec2", "GTACACGTAGTATA", "!!!!!!!!!!!!!!");
  const fastq_pair_vec adapters = create_adapter_vec(
    fastq("PCR1", "CGCTGA", "!!!!!!"), fastq("PCR2", "TGTAC", "!!!!!"));
  const alignment_info expected = ALN().offset(-4).length(18);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(truncate_paired_ended_sequences(result, tmp_record1, tmp_record2) ==
          2);
  REQUIRE(tmp_record1 == fastq("Rec1", "ACGTAGTATA", "!!!!!!!!!!"));
  REQUIRE(tmp_record2 == fastq("Rec2", "ACGTAGTATA", "!!!!!!!!!!"));
}

TEST_CASE("Adapter only sequences", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "CCCGAC", "!!!!!!");
  const fastq record2("Rec2", "ATGCCTT", "!!!!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("PCR1", "CCCGACCCGT", "!!!!!!!!!!"),
                       fastq("PCR2", "AAGATGCCTT", "!!!!!!!!!!"));
  const alignment_info expected = ALN().offset(-7).length(13);
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);

  REQUIRE(result == expected);
}

TEST_CASE("Adapter only sequences, with missing base",
          "[alignment::paired_end]")
{
  // Test the case where both reads are adapters, but are missing a single base
  // Normally, alignments that do not involve read1 vs read2 are skipped, but
  // missing bases may cause some alignments to be missed.
  const fastq record1("Rec1", "CCGACC", "!!!!!!");
  const fastq record2("Rec2", "ATGCCT", "!!!!!!");
  const fastq_pair_vec adapters =
    create_adapter_vec(fastq("PCR1", "CCCGACCCGT", "!!!!!!!!!!"),
                       fastq("PCR2", "AAGATGCCTT", "!!!!!!!!!!"));

  // Sub-optimal alignment:
  //   aagatgccttCCGACC
  //          ATGCCTcccgacccgt
  const alignment_info expected_1 = ALN().offset(-3).length(9).n_mismatches(4);
  const alignment_info result_1 =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result_1 == expected_1);
  // Optimal alignment, only possible with shift
  //   aagatgccttCCGACC
  //      ATGCCTcccgacccgt
  const alignment_info expected_2 = ALN().offset(-7).length(13).n_mismatches(1);
  const alignment_info result_2 =
    align_paired_ended_sequences(record1, record2, adapters, 1);
  REQUIRE(result_2 == expected_2);
}

///////////////////////////////////////////////////////////////////////////////
//

TEST_CASE("Invalid alignment", "[alignment::paired_end]")
{
  fastq record1("Rec", "", "");
  fastq record2("Rec", "", "");
  const alignment_info alignment = ALN().offset(1);

  REQUIRE_THROWS_AS(
    truncate_paired_ended_sequences(alignment, record1, record2),
    assert_failed);
}

///////////////////////////////////////////////////////////////////////////////
// Merging of reads using the conservative method implemented in AdapterRemoval

TEST_CASE("Merge partial overlap [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "ATATTATA", "01234567");
  fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "ATATTATAA", "JJJJJJJJJ");
  fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
  const alignment_info alignment = ALN().offset(-1);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 2);
  const fastq expected =
    fastq("Rec1", "ATATTATA", "wwwwwwww", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap for mate 1 [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "ATATTATAG", "JJJJJJJJJ");
  fastq record2("Rec2", "ATATTATA", "JJJJJJJJ");
  const alignment_info alignment = ALN();
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 1);
  const fastq expected =
    fastq("Rec1", "ATATTATA", "wwwwwwww", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap for mate 2 [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "ATATTATA", "JJJJJJJJ");
  fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
  const alignment_info alignment = ALN().offset(-1);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 1);
  const fastq expected =
    fastq("Rec1", "ATATTATA", "wwwwwwww", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 1 shorter [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "ATA", "012");
  fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(3);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATANNNNACGT", "012ABCDEFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 1 shorter, mate 2 extends past "
          "[deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "ATA", "012");
  fastq record2("Rec2", "AANNNNACGT", "90ABCDEFGH");
  const alignment_info alignment = ALN().offset(-2);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 1);
  const fastq expected = fastq("Rec1", "ATANACGT", "012DEFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 2 shorter [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "ATATTATA", "01234567");
  fastq record2("Rec2", "ACG", "EFG");
  const alignment_info alignment = ALN().offset(8);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATATTATAACG", "01234567EFG");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Ambiguous sites are filled from mate [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "NNNNNNTATA", "0123456789");
  fastq record2("Rec2", "ACGTNNNNNN", "ABCDEFGHIJ");
  const alignment_info alignment;
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ACGTNNTATA", "ABCD!!6789");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Identical nucleotides gets higher qualities [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "GCATGATATA", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected =
    fastq("Rec1", "GCATGATATATACAAC", "012345(FBcEFGHIJ", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Identical nucleotides gets higher qualities, no more than 41 "
          "[deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "GCATGATATA", "0123456789");
  fastq record2("Rec2", "TATATACAAC", "ABCDEFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected =
    fastq("Rec1", "GCATGATATATACAAC", "012345Z\\^`EFGHIJ", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Higher quality nucleotide is selected [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "GCATGAGCAT", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "GCATGATAATTACAAC", "012345(%5%EFGHIJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Randomly select between different nucleotides with same quality #1",
          "[alignment::collapse]")
{
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment;
  std::seed_seq seed{ 1 };
  std::mt19937 rng(seed);
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::original);
  merger.set_rng(&rng);

  const fastq collapsed_expected = fastq("Rec1", "G", "#");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == collapsed_expected);
}

TEST_CASE("Randomly select between different nucleotides with same quality #2",
          "[alignment::collapse]")
{
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment;
  std::seed_seq seed{ 2 };
  std::mt19937 rng(seed);
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::original);
  merger.set_rng(&rng);

  const fastq collapsed_expected = fastq("Rec1", "T", "#");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == collapsed_expected);
}

TEST_CASE("Set conflicts to N/! [deterministic]")
{
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment;
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);

  const fastq expected = fastq("Rec1", "N", "!");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Offsets past the end throws [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment = ALN().offset(2);
  REQUIRE_THROWS_AS(merger.merge(alignment, record1, record2), assert_failed);
}

TEST_CASE("Mate numbering is removed [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Read/1", "ATATTATA", "01234567");
  const fastq record2("Read/2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate numbering removed, meta from 1 kept [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Read/1 Meta1", "ATATTATA", "01234567");
  const fastq record2("Read/2 Meta2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read Meta1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate numbering removed, custom separator [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  merger.set_mate_separator(':');

  fastq record1("Read:1", "ATATTATA", "01234567");
  const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate number removed, custom separator, not set [deterministic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Read:1", "ATATTATA", "01234567");
  const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read:1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Merging of reads using the deterministic merging method

TEST_CASE("Merge partial overlap")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATATTATA", "01234567");
  fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATATTATAA", "JJJJJJJJJ");
  fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
  const alignment_info alignment = ALN().offset(-1);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 2);
  const fastq expected = fastq("Rec1", "ATATTATA", "JJJJJJJJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap for mate 1")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATATTATAG", "JJJJJJJJJ");
  fastq record2("Rec2", "ATATTATA", "JJJJJJJJ");
  const alignment_info alignment = ALN();
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 1);
  const fastq expected = fastq("Rec1", "ATATTATA", "JJJJJJJJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap for mate 2")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATATTATA", "JJJJJJJJ");
  fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
  const alignment_info alignment = ALN().offset(-1);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 1);
  const fastq expected = fastq("Rec1", "ATATTATA", "JJJJJJJJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 1 shorter")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATA", "012");
  fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(3);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATANNNNACGT", "012ABCDEFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 1 shorter, mate 2 extends past")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATA", "012");
  fastq record2("Rec2", "AANNNNACGT", "90ABCDEFGH");
  const alignment_info alignment = ALN().offset(-2);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 1);
  const fastq expected = fastq("Rec1", "ATANACGT", "012DEFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 2 shorter")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATATTATA", "01234567");
  fastq record2("Rec2", "ACG", "EFG");
  const alignment_info alignment = ALN().offset(8);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATATTATAACG", "01234567EFG");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Ambiguous sites are filled from mate")
{
  sequence_merger merger;
  fastq record1("Rec1", "NNNNNNTATA", "0123456789");
  fastq record2("Rec2", "ACGTNNNNNN", "ABCDEFGHIJ");
  const alignment_info alignment;
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ACGTNNTATA", "ABCD!!6789");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Identical nucleotides gets higher qualities")
{
  sequence_merger merger;
  fastq record1("Rec1", "GCATGATATA", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "GCATGATATATACAAC", "012345(3:AEFGHIJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Identical nucleotides gets higher qualities [determinsitic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "GCATGATATA", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected =
    fastq("Rec1", "GCATGATATATACAAC", "012345(FBcEFGHIJ", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Identical nucleotides gets higher qualities, max 41 [determinsitic]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::deterministic);
  fastq record1("Rec1", "GCATGATATA", "0123456789");
  fastq record2("Rec2", "TATATACAAC", "ABCDEFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected =
    fastq("Rec1", "GCATGATATATACAAC", "012345Z\\^`EFGHIJ", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Higher quality nucleotide is selected")
{
  sequence_merger merger;
  fastq record1("Rec1", "GCATGAGCAT", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(truncate_paired_ended_sequences(alignment, record1, record2) == 0);
  const fastq expected = fastq("Rec1", "GCATGATAATTACAAC", "012345($5#EFGHIJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Set conflicts to N/!")
{
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment;
  sequence_merger merger;

  const fastq expected = fastq("Rec1", "N", "!");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Offsets past the end throws")
{
  sequence_merger merger;
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment = ALN().offset(2);
  REQUIRE_THROWS_AS(merger.merge(alignment, record1, record2), assert_failed);
}

TEST_CASE("Mate numbering is removed")
{
  sequence_merger merger;
  fastq record1("Read/1", "ATATTATA", "01234567");
  const fastq record2("Read/2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate numbering removed, meta from 1 kept")
{
  sequence_merger merger;
  fastq record1("Read/1 Meta1", "ATATTATA", "01234567");
  const fastq record2("Read/2 Meta2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read Meta1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate numbering removed, non-standard separator")
{
  sequence_merger merger;
  merger.set_mate_separator(':');

  fastq record1("Read:1", "ATATTATA", "01234567");
  const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate numbering removed, non-standard separator, not set")
{
  sequence_merger merger;
  fastq record1("Read:1", "ATATTATA", "01234567");
  const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read:1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Adapter extraction

TEST_CASE("Extracting empty sequences yields empty sequences #1",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "", "");
  const fastq expected_2 = fastq("read2", "", "");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == expected_1);
  REQUIRE(read2 == expected_2);
}

TEST_CASE("Extracting empty sequences yields empty sequences #2",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "", "");
  const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
  fastq read1 = expected_1;
  fastq read2 = fastq("read2", "GGGGCC", "!!!!!!");
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting empty sequences yields empty sequences #3",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
  const fastq expected_2 = fastq("read2", "", "");
  fastq read1 = fastq("read1", "", "");
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting empty sequences yields empty sequences #4",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "", "");
  const fastq expected_2 = fastq("read2", "", "");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting with no alignment yields empty sequences",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
  const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting with partial overlap yields empty sequences",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
  const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN().offset(2), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting with complete overlap yields empty sequences",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
  const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting sequence 2 inside sequence 1 yields empty sequences",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "AATTTT", "!!!!!!");
  const fastq expected_2 = fastq("read2", "GGCC", "!!!!");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN().offset(2), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting sequence 1 inside sequence 2 yields empty sequences",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "AATT", "!!!!");
  const fastq expected_2 = fastq("read2", "GGGGCC", "!!!!!!");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting sequence 1 extending past sequence 2",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "AATTTTCC", "12345678");
  const fastq expected_2 = fastq("read2", "GGGGGG", "!!!!!!");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN(), read1, read2);
  REQUIRE(read1 == fastq("read1", "CC", "78"));
  REQUIRE(read2 == fastq("read2", "", ""));
}

TEST_CASE("Extracting sequence 2 extending past sequence 1",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "TTTTTT", "!!!!!!");
  const fastq expected_2 = fastq("read2", "AAGGGGGG", "12345678");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN().offset(-2), read1, read2);
  REQUIRE(read1 == fastq("read1", "", ""));
  REQUIRE(read2 == fastq("read2", "AA", "12"));
}

TEST_CASE("Extracting both sequences extending past each other",
          "[alignment::extract_adapter_sequences]")
{
  const fastq expected_1 = fastq("read1", "TTTTTTCCC", "ABCDEFGHI");
  const fastq expected_2 = fastq("read2", "AAGGGGGG", "12345678");
  fastq read1 = expected_1;
  fastq read2 = expected_2;
  extract_adapter_sequences(ALN().offset(-2), read1, read2);
  REQUIRE(read1 == fastq("read1", "CCC", "GHI"));
  REQUIRE(read2 == fastq("read2", "AA", "12"));
}

///////////////////////////////////////////////////////////////////////////////
// Brute-force checking of alignment calculations

struct MMNs
{
  size_t mismatches = 0;
  size_t ambiguous = 0;

  bool operator==(const MMNs& other) const
  {
    return mismatches == other.mismatches && ambiguous == other.ambiguous;
  }

  bool operator!=(const MMNs& other) const { return !(*this == other); }
};

char
rotate_nucleotide(char c)
{
  switch (c) {
    case 'A':
      return 'C';
    case 'C':
      return 'G';
    case 'G':
      return 'T';
    case 'T':
      return 'A';

    default:
      AR_FAIL("invalid nucleotide");
  }
}

void
compare(simd::compare_subsequences_func func,
        const std::string& seq_1,
        const std::string& seq_2,
        const size_t length,
        const MMNs& expected)
{
  AR_REQUIRE(seq_1.size() == seq_2.size());

  MMNs alignment{ 0, 0 };
  const auto result = func(alignment.mismatches,
                           alignment.ambiguous,
                           seq_1.c_str(),
                           seq_2.c_str(),
                           length,
                           length * 2);

  AR_REQUIRE(result);

  if (alignment != expected) {
    REQUIRE(alignment == expected);
  }
}

TEST_CASE("Brute-force validation", "[alignment::compare_subsequences]")
{
  // Parameterize tests over supported SIMD instruction sets
  const auto is = PARAMETERIZE_IS;
  const auto func = simd::get_compare_subsequences_func(is);
  const auto padding = simd::padding(is);
  // Randomly generated sequence
  const std::string reference =
    "CTGGTTAAAGATCAGAATCCTTTTATTTGCGGAAATTCGAATTATATCCTGATCAGTCGGTGGCGCTAGTGTCC"
    "AGGGGATCTTGGAATTGGATCCAAAAAGTGCTGGGGAATGCGGATTCCATTATGAGACCTGT";

  // The SECTION identifies the instruction set in failure messages
  SECTION(simd::name(is))
  {
    compare(func, "", "", 0, MMNs{ 0, 0 });

    for (size_t length = 1; length < reference.size(); ++length) {
      SECTION(std::to_string(length))
      {
        std::string seq_1 = reference.substr(0, length);
        seq_1.resize(length + padding, 'N');
        std::string seq_2 = seq_1;

        compare(func, seq_1, seq_2, length, MMNs{ 0, 0 });

        for (size_t i = 0; i < length; ++i) {
          seq_2.at(i) = rotate_nucleotide(seq_2.at(i));
          compare(func, seq_1, seq_2, length, MMNs{ 1, 0 });

          if (i > 0) {
            seq_1.at(i - 1) = rotate_nucleotide(seq_1.at(i - 1));
            compare(func, seq_1, seq_2, length, MMNs{ 2, 0 });

            seq_1.at(i - 1) = 'N';
            compare(func, seq_1, seq_2, length, MMNs{ 1, 1 });

            seq_1.at(i - 1) = reference.at(i - 1);
          }

          seq_2.at(i) = 'N';
          compare(func, seq_1, seq_2, length, MMNs{ 0, 1 });

          if (i > 0) {
            seq_1.at(i - 1) = rotate_nucleotide(seq_1.at(i - 1));
            compare(func, seq_1, seq_2, length, MMNs{ 1, 1 });

            seq_1.at(i - 1) = 'N';
            compare(func, seq_1, seq_2, length, MMNs{ 0, 2 });

            seq_1.at(i - 1) = reference.at(i - 1);
          }

          seq_2.at(i) = reference.at(i);
        }
      }
    }
  }
}

} // namespace adapterremoval

namespace Catch {

using namespace adapterremoval;

template<>
std::string
StringMaker<alignment_info, void>::convert(alignment_info const& value)
{
  std::vector<std::string> labels = { "score",       "offset",
                                      "length",      "n_mismatches",
                                      "n_ambiguous", "adapter_id" };

  std::vector<int64_t> values = {
    static_cast<int64_t>(value.score()),
    static_cast<int64_t>(value.offset),
    static_cast<int64_t>(value.length),
    static_cast<int64_t>(value.n_mismatches),
    static_cast<int64_t>(value.n_ambiguous),
    static_cast<int64_t>(value.adapter_id),
  };

  std::ostringstream stream;
  stream << "alignment_info(";

  bool any_streamed = false;
  for (size_t i = 0; i < labels.size(); ++i) {
    if (values.at(i)) {
      if (any_streamed) {
        stream << ", ";
      }

      stream << labels.at(i) << " = " << values.at(i);
      any_streamed = true;
    }
  }

  stream << ")";

  return stream.str();
}

template<>
std::string
StringMaker<ALN, void>::convert(ALN const& value)
{
  return StringMaker<alignment_info>::convert(value.info);
}

template<>
std::string
StringMaker<adapterremoval::MMNs, void>::convert(
  adapterremoval::MMNs const& value)
{
  std::ostringstream stream;
  stream << "{MMs = " << value.mismatches << ", Ns = " << value.ambiguous
         << "}";
  return stream.str();
}

} // namespace Catch
