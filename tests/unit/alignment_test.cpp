// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "alignment.hpp" // for alignment_info, sequence_merger, extract_...
#include "catch.hpp"
#include "commontypes.hpp"   // for merge_strategy, merge_strategy::determini...
#include "debug.hpp"         // for AR_FAIL
#include "errors.hpp"        // for assert_failed
#include "fastq.hpp"         // for fastq
#include "fastq_enc.hpp"     // for FASTQ_ENCODING_SAM
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for adapter_set
#include "simd.hpp"          // for size_t, instruction_set, supported, get_c...
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <sstream>           // for ostringstream
#include <string>            // for string, basic_string, operator<<

// Ignore nucleotide and quality strings
// spell-checker:ignoreRegExp /"[!-~]+"/g
// Ignore nucleotide comments
// spell-checker:ignoreRegExp /\W[acgtnACGTN]+\W/g

namespace adapterremoval {

struct ALN;

#define TEST_ALIGNMENT_SETTER(TYPE, NAME)                                      \
  ALN& NAME(TYPE value)                                                        \
  {                                                                            \
    info.m_##NAME = value;                                                     \
    return *this;                                                              \
  }

// Parameterize tests over supported SIMD instruction sets
#define PARAMETERIZE_IS GENERATE(from_range(simd::supported()))

// Dummy mismatch threshold
const double DEFAULT_MISMATCH_THRESHOLD = 0.0;

struct ALN
{
  ALN()
    : info()
  {
    // Default to representing an actual alignment
    info.m_adapter_id = 0;
  }

  TEST_ALIGNMENT_SETTER(int, offset);
  TEST_ALIGNMENT_SETTER(size_t, length);
  TEST_ALIGNMENT_SETTER(size_t, n_mismatches);
  TEST_ALIGNMENT_SETTER(size_t, n_ambiguous);
  TEST_ALIGNMENT_SETTER(int, adapter_id);

  ALN& is_good()
  {
    info.m_type = alignment_type::good;
    return *this;
  }

  ALN& can_merge()
  {
    info.m_type = alignment_type::mergeable;
    return *this;
  }

  // NOLINTNEXTLINE(hicpp-explicit-conversions)
  operator const alignment_info&() const { return info; }

  alignment_info info;
};

alignment_info
align_single_ended_sequence(const fastq& read,
                            const adapter_set& adapters,
                            int max_shift)
{
  return sequence_aligner(adapters, PARAMETERIZE_IS, DEFAULT_MISMATCH_THRESHOLD)
    .align_single_end(read, max_shift);
}

alignment_info
align_paired_ended_sequences(const fastq& read1,
                             const fastq& read2,
                             const adapter_set& adapters,
                             int max_shift)
{
  return sequence_aligner(adapters, PARAMETERIZE_IS, DEFAULT_MISMATCH_THRESHOLD)
    .align_paired_end(read1, read2, max_shift);
}

void
REQUIRE_TRUNCATED_PE_IS_UNCHANGED(const alignment_info& alignment,
                                  const fastq& record1,
                                  const fastq& record2)
{
  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(alignment.truncate_paired_end(tmp_record1, tmp_record2) == 0);
  REQUIRE(tmp_record1 == record1);
  REQUIRE(tmp_record2 == record2);
}

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
// Single end alignment using `align_single_ended_sequence`

///////////////////////////////////////////////////////////////////////////////
// Case 1: No overlap between sequences:
//         AAAAAAAAAAA
//                       BBBBBBBBBBBBBBB
// Expected result = suboptimal alignment or no overlap

TEST_CASE("SE: Unalignable sequence yields default alignment",
          "[alignment::single_end]")
{
  const fastq record("Rec", "AAAA", "!!!!");
  const adapter_set adapters = { { "TTTT", "" } };

  REQUIRE(align_single_ended_sequence(record, adapters, 0) == alignment_info());
}

TEST_CASE("SE: Random sequences yields suboptimal alignment",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
  const adapter_set adapters = { { "TGAGACGGT", "" } };

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
  const adapter_set adapters = { { "AGTAAGGT", "" } };
  const alignment_info expected = ALN().offset(4).length(5).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACGT", "1234"));
}

TEST_CASE("SE: Partial alignment with mismatches between ends",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTAA", "123457890");
  const adapter_set adapters = { { "AGGAAGGT", "" } };
  const alignment_info expected = ALN().offset(4).length(5).n_mismatches(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACGT", "1234"));
}

TEST_CASE("SE: Partial alignment with ambiguous between ends",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTAA", "123457890");
  const adapter_set adapters = { { "AGNAAGGT", "" } };
  const alignment_info expected =
    ALN().offset(4).length(5).n_ambiguous(1).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
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
  const adapter_set adapters = { { { record.sequence() }, "" } };
  const alignment_info expected = ALN().length(8).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("SE: Completely overlapping sequences with mismatches",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
  const adapter_set adapters = { { "GCGTAGTA", "" } };
  const alignment_info expected = ALN().length(8).n_mismatches(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("SE: Completely overlapping sequences with mismatches and ambiguous",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "!!!!!!!!");
  const adapter_set adapters = { { "GCGTAGTN", "" } };
  const alignment_info expected =
    ALN().length(8).n_mismatches(1).n_ambiguous(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
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
  const adapter_set adapters = { { "TAGTA", "" } };
  const alignment_info expected = ALN().offset(3).length(5).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACG", "ABC"));
}

TEST_CASE("Complete adapter inside sequence with mismatch",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
  const adapter_set adapters = { { "TATTA", "" } };
  const alignment_info expected = ALN().offset(3).length(5).n_mismatches(1);
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACG", "ABC"));
}

TEST_CASE("Complete adapter inside sequence with ambiguous",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTA", "ABCDEFGH");
  const adapter_set adapters = { { "TAGNA", "" } };
  const alignment_info expected =
    ALN().offset(3).length(5).n_ambiguous(1).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACG", "ABC"));
}

TEST_CASE("Complete sequence inside adapter", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGT", "!!!!");
  const adapter_set adapters = { { "ACGTAGTA", "" } };
  const alignment_info expected = ALN().length(4).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("Complete sequence inside adapter with mismatches",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGT", "!!!!");
  const adapter_set adapters = { { "GCGTAGTA", "" } };
  const alignment_info expected = ALN().length(4).n_mismatches(1);
  const auto result = align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

TEST_CASE("Complete sequence inside adapter with ambiguous",
          "[alignment::single_end]")
{
  const fastq record("Rec", "ACGT", "!!!!");
  const adapter_set adapters = { { "ACGNAGTA", "" } };
  const alignment_info expected = ALN().length(4).n_ambiguous(1).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

///////////////////////////////////////////////////////////////////////////////
// Cases 6 and 7: Sequence A extends past sequence B (and vice versa)
//         AAAAAAAAAAAAAA    AAAA
//               BBBBB     BBBBBBBBBBBBB

TEST_CASE("Sequence extends past adapter", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTATA", "0123456789");
  const adapter_set adapters = { { "AGTA", "" } };
  const alignment_info expected = ALN().offset(4).length(4).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "ACGT", "0123"));
}

TEST_CASE("Sequence extends past adapter, no shift", "[alignment::single_end]")
{
  const fastq record("Rec", "CGTA", "#!%%");
  const adapter_set adapters = { { "ACGTAGTATA", "" } };
  const alignment_info expected = ALN().offset(3).length(1).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "CGT", "#!%"));
}

TEST_CASE("Sequence extends past adapter, shift 1", "[alignment::single_end]")
{
  const fastq record("Rec", "CGTA", "#!%%");
  const adapter_set adapters = { { "ACGTAGTATA", "" } };
  const alignment_info expected = ALN().offset(-1).length(4).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 1);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
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
  const adapter_set adapters = { { "CCGAACGTAGTATA", "" } };
  const alignment_info expected = ALN().offset(-4).length(10).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 4);
  REQUIRE(result == expected);

  fastq tmp_record = record;
  result.truncate_single_end(tmp_record);
  REQUIRE(tmp_record == fastq("Rec", "", ""));
}

///////////////////////////////////////////////////////////////////////////////
// Empty sequence or adapter

TEST_CASE("Empty sequence alignment", "[alignment::single_end]")
{
  const fastq record;
  const adapter_set adapters = { { "CCGAACGTAGTATA", "" } };
  const alignment_info expected;
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Empty adapter alignment", "[alignment::single_end]")
{
  const fastq record("Rec", "ACGTAGTATATAGT", "!!!!!!!!!!!!!!");
  const adapter_set adapters = { { "", "" } };
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
  const adapter_set adapters = { { "AAAAAAAAA", "" } };

  sequence_aligner aligner(adapters,
                           PARAMETERIZE_IS,
                           DEFAULT_MISMATCH_THRESHOLD);
  const auto result = aligner.align_single_end(record, 0);
  const alignment_info expected = ALN().length(9).n_mismatches(1);

  REQUIRE(result == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Misc

TEST_CASE("Lower than possible shift is allowed", "[alignment::single_end]")
{
  const fastq record("Rec", "AAAA", "!!!!");
  const adapter_set adapters = { { "TTTT", "" } };
  const alignment_info expected;
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 10);
  REQUIRE(result == expected);
}

TEST_CASE("Only adapter 1 is used", "[alignment::single_end]")
{
  const adapter_set adapters = { { "AAA", "TTTAAA" } };
  const fastq record("Rec", "CCCCTTTAAA", "0987654321");
  const alignment_info expected = ALN().offset(7).length(3).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Best matching adapter is returned: First", "[alignment::single_end]")
{
  const adapter_set adapters = {
    { "TGCTGC", "" },
    { "TGCTGA", "" },
  };
  const fastq record("Read", "TAGTCGCTATGCTGC", "!!!!!!!!!103459");
  const alignment_info expected = ALN().offset(9).length(6).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Best matching adapter returned: Second", "[alignment::single_end]")
{
  const adapter_set adapters = {
    { "TGCTGA", "" },
    { "TGCTGC", "" },
  };

  const fastq record("Read", "TAGTCGCTATGCTGC", "!!!!!!!!!103459");
  const alignment_info expected =
    ALN().offset(9).length(6).adapter_id(1).is_good();
  const alignment_info result =
    align_single_ended_sequence(record, adapters, 0);
  REQUIRE(result == expected);
}

TEST_CASE("Best matching adapter returned: Neither", "[alignment::single_end]")
{
  const adapter_set barcodes = {
    { "AAAAAA", "" },
    { "CCCCCC", "" },
  };

  const fastq record = fastq("Read", "AACTGTACGTAGTT", "!!!!!!10345923");
  const alignment_info result =
    align_single_ended_sequence(record, barcodes, 0);
  REQUIRE(result == alignment_info());
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
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
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
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
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

TEST_CASE("Partial overlap in sequence pair (1 bp)", "[alignment::paired_end]")
{
  const fastq record1("Rec", "AAAAAAAAA", "!!!!!!!!!");
  const fastq record2("Rec", "ATTTTTTT", "!!!!!!!!");
  const adapter_set adapters = { { "TTTTTT", "TTTTT" } };
  const alignment_info expected = ALN().offset(8).length(1).is_good();
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

TEST_CASE("Partial overlap in sequence pair (2 bp)", "[alignment::paired_end]")
{
  const fastq record1("Rec", "AAAAAAAAA", "!!!!!!!!!");
  const fastq record2("Rec", "AATTTTTT", "!!!!!!!!");
  const adapter_set adapters = { { "TTTTTT", "TTTTT" } };
  const alignment_info expected = ALN().offset(7).length(2).is_good();
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

TEST_CASE("Partial overlap in sequence pair (5 bp)", "[alignment::paired_end]")
{
  const fastq record1("Rec", "ACGTAGTAA", "!!!!!!!!!");
  const fastq record2("Rec", "AGTAAGGT", "!!!!!!!!");
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
  const alignment_info expected = ALN().offset(4).length(5).is_good();
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
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
  const alignment_info expected = ALN().length(8).is_good();
  const alignment_info result =
    align_paired_ended_sequences(record1, record1, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record1);
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
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
  const alignment_info expected = ALN().offset(3).length(5).is_good();
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);
  REQUIRE_TRUNCATED_PE_IS_UNCHANGED(result, record1, record2);
}

TEST_CASE("Sequence B contains sequence A", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "ACGT", "!!!!");
  const fastq record2("Rec2", "ACGTAGTA", "!!!!!!!!");
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
  const alignment_info expected = ALN().length(4).is_good();
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
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
  const alignment_info expected = ALN().offset(4).length(6).is_good();
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(result.truncate_paired_end(tmp_record1, tmp_record2) == 1);
  REQUIRE(tmp_record1 == fastq("Rec1", "ACGTAGTA", "!!!!!!!!"));
  REQUIRE(tmp_record2 == fastq("Rec2", "AGTA", "!!!!"));
}

TEST_CASE("Sequence B extends past sequence A", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "CGTA", "!!!!");
  const fastq record2("Rec2", "ACCGTAGTAT", "!!!!!!!!!!");
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
  const alignment_info expected = ALN().offset(-2).length(6).is_good();
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(result.truncate_paired_end(tmp_record1, tmp_record2) == 1);
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
  const adapter_set adapters = { { "CGCTGA", "GTACA" } };
  const alignment_info expected = ALN().offset(-4).length(18).is_good();
  const alignment_info result =
    align_paired_ended_sequences(record1, record2, adapters, 0);
  REQUIRE(result == expected);

  fastq tmp_record1 = record1;
  fastq tmp_record2 = record2;
  REQUIRE(result.truncate_paired_end(tmp_record1, tmp_record2) == 2);
  REQUIRE(tmp_record1 == fastq("Rec1", "ACGTAGTATA", "!!!!!!!!!!"));
  REQUIRE(tmp_record2 == fastq("Rec2", "ACGTAGTATA", "!!!!!!!!!!"));
}

TEST_CASE("Adapter only sequences", "[alignment::paired_end]")
{
  const fastq record1("Rec1", "CCCGAC", "!!!!!!");
  const fastq record2("Rec2", "ATGCCTT", "!!!!!!!");
  const adapter_set adapters = { { "CCCGACCCGT", "AAGGCATCTT" } };
  const alignment_info expected = ALN().offset(-7).length(13).is_good();
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
  const adapter_set adapters = { { "CCCGACCCGT", "AAGGCATCTT" } };

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
// Alignment flags

TEST_CASE("Positive score is required for good alignment", "[alignment:flags]")
{
  const adapter_set adapters{ { "AT", "GC" } };
  sequence_aligner aligner{ adapters, PARAMETERIZE_IS, 1.0 };

  SECTION("single end (score = 0)")
  {
    const auto result = aligner.align_single_end({ "Read", "AC" }, 0);
    REQUIRE(result == ALN().length(2).n_mismatches(1));
  }

  SECTION("single end (score = 1)")
  {
    const auto result = aligner.align_single_end({ "Read", "AN" }, 0);
    REQUIRE(result == ALN().length(2).n_ambiguous(1).is_good());
  }

  SECTION("paired end (score = 0)")
  {
    const auto result =
      aligner.align_paired_end({ "Read", "CCT" }, { "Read", "CA" }, 0);
    REQUIRE(result == ALN().offset(-1).length(4).n_mismatches(2));
  }

  SECTION("paired end (score = 1)")
  {
    const auto result =
      aligner.align_paired_end({ "Read", "CN" }, { "Read", "CA" }, 0);
    REQUIRE(result == ALN().length(2).n_ambiguous(1).is_good());
  }
}

TEST_CASE("Single end with minimum overlap", "[alignment:flags]")
{
  const adapter_set adapters{ { "ATGC", "CCCC" } };
  sequence_aligner aligner{ adapters, PARAMETERIZE_IS, 1.0 };
  aligner.set_min_se_overlap(2);

  SECTION("does not affect PE")
  {
    const auto result =
      aligner.align_paired_end({ "Read", "AT" }, { "Read", "TA" }, 0);
    REQUIRE(result == ALN().offset(1).length(1).is_good());
  }

  SECTION("SE length 1 too short")
  {
    const auto result = aligner.align_single_end({ "Read", "CGTA" }, 0);
    REQUIRE(result == ALN().offset(3).length(1));
  }

  SECTION("SE length 2 long enough")
  {
    const auto result = aligner.align_single_end({ "Read", "GCAT" }, 0);
    REQUIRE(result == ALN().offset(2).length(2).is_good());
  }

  SECTION("SE length 3 long enough")
  {
    const auto result = aligner.align_single_end({ "Read", "CATG" }, 0);
    REQUIRE(result == ALN().offset(1).length(3).is_good());
  }
}

TEST_CASE("Very short alignments must be perfect", "[alignment:flags]")
{
  const std::string adapter{ "TATGGTTGAT" };
  const adapter_set adapters{ { adapter, "" } };
  sequence_aligner aligner{ adapters, PARAMETERIZE_IS, 1.0 };

  SECTION("less than 6 disallows mismatches (good)")
  {
    std::string seq = adapter.substr(0, GENERATE(2, 3, 4, 5));
    const fastq read{ "Rec", seq };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(seq.length()).is_good());
  }

  SECTION("less than 6 disallows mismatches (not good)")
  {
    std::string seq = adapter.substr(0, GENERATE(2, 3, 4, 5));
    seq.at(0) = 'C'; // single mismatch
    const fastq read{ "Rec", seq };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(seq.length()).n_mismatches(1));
  }

  SECTION("less than 10 allows 1 mismatch (good, mm = 0)")
  {
    std::string seq = adapter.substr(0, GENERATE(6, 7, 8, 9));
    const fastq read{ "Rec", seq };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(seq.length()).is_good());
  }

  SECTION("less than 10 allows 1 mismatch (good, mm = 1)")
  {
    std::string seq = adapter.substr(0, GENERATE(6, 7, 8, 9));
    seq.at(0) = 'C'; // single mismatch
    const fastq read{ "Rec", seq };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(seq.length()).n_mismatches(1).is_good());
  }

  SECTION("less than 10 allows 1 mismatch (bad, mm = 2)")
  {
    std::string seq = adapter.substr(0, GENERATE(6, 7, 8, 9));
    seq.at(0) = seq.at(3) = 'C'; // two mismatch
    const fastq read{ "Rec", seq };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(seq.length()).n_mismatches(2));
  }

  SECTION("10 allows more mismatches")
  {
    std::string seq = adapter;
    seq.at(0) = seq.at(3) = 'C'; // two mismatch
    const fastq read{ "Rec", seq };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(10).n_mismatches(2).is_good());
  }
}

TEST_CASE("is_good depends on allowed mismatch rate", "[alignment:flags]")
{

  const adapter_set adapters{ { "TATGGTTGAT", "" } };
  const fastq read{ "Rec", "TTTGCTCGTT" };

  SECTION("too high error rate")
  {
    sequence_aligner aligner{ adapters,
                              PARAMETERIZE_IS,
                              GENERATE(0.0, 0.1, 0.2, 0.3) };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(10).n_mismatches(4));
  }

  SECTION("acceptable error rate")
  {
    sequence_aligner aligner{ adapters,
                              PARAMETERIZE_IS,
                              GENERATE(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) };
    const auto result = aligner.align_single_end(read, 0);
    REQUIRE(result == ALN().length(10).n_mismatches(4).is_good());
  }
}

TEST_CASE("can_merge requires good alignment", "[alignment:flags]")
{
  const fastq read1("Rec1", "ACCTACTG");
  const fastq read2("Rec2", "CTACTGTT");
  sequence_aligner aligner{ { {} }, PARAMETERIZE_IS, 0.0 };

  SECTION("cannot merge by default")
  {
    const auto result = aligner.align_paired_end(read1, read2, 0);
    REQUIRE(result == ALN().offset(2).length(6).is_good());
  }

  SECTION("can merge when threshold less-than-or-equal")
  {
    // A threshold of 0 is functionally identical to 1
    aligner.set_merge_threshold(GENERATE(0, 1, 2, 3, 4, 5, 6));
    const auto result = aligner.align_paired_end(read1, read2, 0);
    REQUIRE(result == ALN().offset(2).length(6).can_merge());
  }

  SECTION("cannot merge when threshold greater than")
  {
    // A threshold of 0 is functionally identical to 1
    aligner.set_merge_threshold(GENERATE(7, 8, 9, 10));
    const auto result = aligner.align_paired_end(read1, read2, 0);
    REQUIRE(result == ALN().offset(2).length(6).is_good());
  }
}

TEST_CASE("can_merge threshold ignores Ns", "[alignment:flags]")
{
  const fastq read1("Rec1", "ACCTACTG");
  const fastq read2("Rec2", "CTNCTGTT");
  sequence_aligner aligner{ { {} }, PARAMETERIZE_IS, 0.0 };

  SECTION("can merge when threshold less-than-or-equal")
  {
    // A threshold of 0 is functionally identical to 1
    aligner.set_merge_threshold(GENERATE(4, 5));
    const auto result = aligner.align_paired_end(read1, read2, 0);
    REQUIRE(result == ALN().offset(2).length(6).n_ambiguous(1).can_merge());
  }

  SECTION("cannot merge when threshold greater than")
  {
    // A threshold of 0 is functionally identical to 1
    aligner.set_merge_threshold(GENERATE(6, 7));
    const auto result = aligner.align_paired_end(read1, read2, 0);
    REQUIRE(result == ALN().offset(2).length(6).n_ambiguous(1).is_good());
  }
}

///////////////////////////////////////////////////////////////////////////////
// Misc

// Test for bug where SE alignments would occasionally test sequences with fewer
// bp than the current score, that therefore could never be picked.
TEST_CASE("Pointless SE alignments #1", "[alignment::single_end]")
{
  // Subsequent candidates can never be better than the best previous alignment
  const adapter_set adapters{ { "TT", "" }, { "T", "" }, {} };
  sequence_aligner aligner{ adapters,
                            simd::instruction_set::none,
                            DEFAULT_MISMATCH_THRESHOLD };
  const alignment_info expected = ALN().length(2).adapter_id(0).is_good();
  const auto result = aligner.align_single_end({ "Rec", "TTT" }, 10);
  REQUIRE(result == expected);
}

TEST_CASE("Pointless SE alignments #2", "[alignment::single_end]")
{
  // Subsequent candidates can never be better than the best previous alignment
  const adapter_set adapters{ { "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "" },
                              { "TGGAATTCTCGGGTGCCAAGG", "" } };
  const fastq read{ "Rec", "ATCGGAAGAGCACACGTCTGAACTCCAGTCACAAGGGCATCT" };

  sequence_aligner aligner{ adapters, simd::instruction_set::none, 1.0 / 6.0 };
  const alignment_info expected =
    ALN().offset(-2).length(31).adapter_id(0).is_good();
  const auto result = aligner.align_single_end(read, 2);
  REQUIRE(result == expected);
}

TEST_CASE("Adapter-skipping for SE adapters", "[alignment::single_end]")
{
  const adapter_set adapters = {
    { "AAAAA", "" },
    { "AAAA", "" },
    { "AAAAAA", "" },
  };

  auto align = [&adapters](std::string seq) {
    return align_single_ended_sequence({ "R1", seq }, adapters, 0);
  };

  CHECK(align("AAAAA") == ALN().length(5).adapter_id(0).is_good());
  CHECK(align("AAAAAA") == ALN().length(6).adapter_id(2).is_good());
}

TEST_CASE("Adapter-skipping for PE adapters", "[alignment::single_end]")
{
  const adapter_set adapters = {
    { "AAAAA", "AAAAA" },
    { "AAAA", "AAAA" },
    { "AAAAAA", "AAAAAA" },
  };

  auto align = [&adapters](std::string seq1, std::string seq2) {
    return align_paired_ended_sequences({ "R1", seq1 },
                                        { "R2", seq2 },
                                        adapters,
                                        0);
  };

  CHECK(align("AAAAA", "TTTTT") ==
        ALN().length(10).offset(-5).adapter_id(0).is_good());
  CHECK(align("AAAAAA", "TTTTTT") ==
        ALN().length(12).offset(-6).adapter_id(2).is_good());
}

TEST_CASE("Adapter-sorting for SE adapters", "[alignment::single_end]")
{
  const adapter_set adapters = {
    { "AAAAA", "" },
    { "TTTTT", "" },
  };

  sequence_aligner aligner{ adapters,
                            PARAMETERIZE_IS,
                            DEFAULT_MISMATCH_THRESHOLD };
  auto align = [&aligner](std::string seq) {
    return aligner.align_single_end({ "R1", seq }, 0);
  };

  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 0
  REQUIRE(align("TTTTT").adapter_id() == 1); // pos 1; swap
  REQUIRE(align("TTTTT").adapter_id() == 1); // pos 0
  REQUIRE(align("TTTTT").adapter_id() == 1); // pos 0
  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 1
  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 1; swap
  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 0
}

TEST_CASE("Adapter-sorting for PE adapters", "[alignment::single_end]")
{
  const adapter_set adapters = {
    { "AAAAA", "TTTTT" },
    { "TTTTT", "AAAAA" },
  };

  sequence_aligner aligner{ adapters,
                            PARAMETERIZE_IS,
                            DEFAULT_MISMATCH_THRESHOLD };
  auto align = [&aligner](std::string seq) {
    return aligner.align_paired_end({ "R1", seq }, { "R2", seq }, 0);
  };

  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 0
  REQUIRE(align("TTTTT").adapter_id() == 1); // pos 1; swap
  REQUIRE(align("TTTTT").adapter_id() == 1); // pos 0
  REQUIRE(align("TTTTT").adapter_id() == 1); // pos 0
  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 1
  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 1; swap
  REQUIRE(align("AAAAA").adapter_id() == 0); // pos 0
}

///////////////////////////////////////////////////////////////////////////////
//

TEST_CASE("Invalid alignment", "[alignment::paired_end]")
{
  fastq record1("Rec", "", "");
  fastq record2("Rec", "", "");
  const alignment_info alignment = ALN().offset(1);

  REQUIRE_THROWS_AS(alignment.truncate_paired_end(record1, record2),
                    assert_failed);
}

///////////////////////////////////////////////////////////////////////////////
// Merging of reads using the conservative method implemented in AdapterRemoval

TEST_CASE("Merge partial overlap [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "ATATTATA", "01234567");
  fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "ATATTATAA", "JJJJJJJJJ");
  fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
  const alignment_info alignment = ALN().offset(-1);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 2);
  const fastq expected =
    fastq("Rec1", "ATATTATA", "ssssssss", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap for mate 1 [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "ATATTATAG", "JJJJJJJJJ");
  fastq record2("Rec2", "ATATTATA", "JJJJJJJJ");
  const alignment_info alignment;
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 1);
  const fastq expected =
    fastq("Rec1", "ATATTATA", "ssssssss", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap for mate 2 [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "ATATTATA", "JJJJJJJJ");
  fastq record2("Rec2", "AATATTATA", "JJJJJJJJJ");
  const alignment_info alignment = ALN().offset(-1);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 1);
  const fastq expected =
    fastq("Rec1", "ATATTATA", "ssssssss", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 1 shorter [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "ATA", "012");
  fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(3);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATANNNNACGT", "012ABCDEFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 1 shorter, mate 2 extends past "
          "[additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "ATA", "012");
  fastq record2("Rec2", "AANNNNACGT", "90ABCDEFGH");
  const alignment_info alignment = ALN().offset(-2);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 1);
  const fastq expected = fastq("Rec1", "ATANACGT", "012DEFGH");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Unequal sequence lengths, mate 2 shorter [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "ATATTATA", "01234567");
  fastq record2("Rec2", "ACG", "EFG");
  const alignment_info alignment = ALN().offset(8);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ATATTATAACG", "01234567EFG");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Ambiguous sites are filled from mate [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "NNNNNNTATA", "0123456789");
  fastq record2("Rec2", "ACGTNNNNNN", "ABCDEFGHIJ");
  const alignment_info alignment;
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected = fastq("Rec1", "ACGTNNTATA", "ABCD!!6789");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Identical nucleotides gets higher qualities [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "GCATGATATA", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected =
    fastq("Rec1", "GCATGATATATACAAC", "012345(B?_EFGHIJ", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Identical nucleotides gets higher qualities, no more than 41 "
          "[additive]")
{
  sequence_merger merger;
  merger.set_max_recalculated_score(41);
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "GCATGATATA", "0123456789");
  fastq record2("Rec2", "TATATACAAC", "ABCDEFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected =
    fastq("Rec1", "GCATGATATATACAAC", "012345JJJJEFGHIJ", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Higher quality nucleotide is selected [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "GCATGAGCAT", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected =
    fastq("Rec1", "GCATGATAATTACAAC", "012345($5#EFGHIJ", FASTQ_ENCODING_SAM);
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Set conflicts to N/! [additive]")
{
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment;
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);

  const fastq expected = fastq("Rec1", "N", "!");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Offsets past the end throws [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Rec1", "G", "1");
  const fastq record2("Rec2", "T", "1");
  const alignment_info alignment = ALN().offset(2);
  REQUIRE_THROWS_AS(merger.merge(alignment, record1, record2), assert_failed);
}

TEST_CASE("Mate numbering is removed [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Read/1", "ATATTATA", "01234567");
  const fastq record2("Read/2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate numbering removed, meta from 1 kept [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Read/1 Meta1", "ATATTATA", "01234567");
  const fastq record2("Read/2 Meta2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read Meta1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate numbering removed, custom separator [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  merger.set_mate_separator(':');

  fastq record1("Read:1", "ATATTATA", "01234567");
  const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

TEST_CASE("Mate number removed, custom separator, not set [additive]")
{
  sequence_merger merger;
  merger.set_merge_strategy(merge_strategy::additive);
  fastq record1("Read:1", "ATATTATA", "01234567");
  const fastq record2("Read:2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  const fastq expected = fastq("Read:1", "ATATTATAACGT", "01234567EFGH");
  merger.merge(alignment, record1, record2);

  REQUIRE(record1 == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Merging of reads using the additive merging method

TEST_CASE("Merge partial overlap")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATATTATA", "01234567");
  fastq record2("Rec2", "NNNNACGT", "ABCDEFGH");
  const alignment_info alignment = ALN().offset(4);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
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
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 2);
  const fastq expected = fastq("Rec1", "ATATTATA", "JJJJJJJJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Merge complete overlap for mate 1")
{
  sequence_merger merger;
  fastq record1("Rec1", "ATATTATAG", "JJJJJJJJJ");
  fastq record2("Rec2", "ATATTATA", "JJJJJJJJ");
  const alignment_info alignment;
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 1);
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
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 1);
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
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
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
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 1);
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
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
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
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
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
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
  const fastq expected = fastq("Rec1", "GCATGATATATACAAC", "012345(3:AEFGHIJ");
  merger.merge(alignment, record1, record2);
  REQUIRE(record1 == expected);
}

TEST_CASE("Higher quality nucleotide is selected")
{
  sequence_merger merger;
  fastq record1("Rec1", "GCATGAGCAT", "012345!0:A");
  fastq record2("Rec2", "TATATACAAC", "(3&?EFGHIJ");
  const alignment_info alignment = ALN().offset(6);
  REQUIRE(alignment.truncate_paired_end(record1, record2) == 0);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment = ALN().offset(2);
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment = ALN().offset(2);
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment;
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment = ALN().offset(-2);
  alignment.extract_adapter_sequences(read1, read2);
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
  alignment_info alignment = ALN().offset(-2);
  alignment.extract_adapter_sequences(read1, read2);
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

std::ostream&
operator<<(std::ostream& os, const MMNs& value)
{
  os << "MMNs{mismatches=" << value.mismatches
     << ", ambiguous=" << value.ambiguous << "}";

  return os;
}

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

TEST_CASE("rotate nucleotide")
{
  REQUIRE(rotate_nucleotide('A') == 'C');
  REQUIRE(rotate_nucleotide('C') == 'G');
  REQUIRE(rotate_nucleotide('G') == 'T');
  REQUIRE(rotate_nucleotide('T') == 'A');
  REQUIRE(rotate_nucleotide('T') == 'A');
  REQUIRE_THROWS_AS(rotate_nucleotide('X'), assert_failed);
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

    const size_t length = reference.length();
    auto seq_1 = std::string(length, 'A') + std::string(padding, 'N');
    auto seq_2 = std::string(length, 'N') + std::string(padding, 'N');
    compare(func, seq_1, seq_2, length, MMNs{ 0, length });

    seq_2 = std::string(length, 'T') + std::string(padding, 'N');
    compare(func, seq_1, seq_2, length, MMNs{ length, 0 });
  }
}

TEST_CASE("stringmaker for empty alignment_info")
{
  std::ostringstream os;
  os << alignment_info{};

  REQUIRE(os.str() == "alignment_info{score=0, adapter_id=-1, offset=0, "
                      "length=0, n_mismatches=0, n_ambiguous=0, m_type=bad}");
}

TEST_CASE("stringmaker for alignment_info")
{
  alignment_info info = ALN()
                          .offset(2)
                          .length(3)
                          .n_mismatches(4)
                          .n_ambiguous(5)
                          .adapter_id(1)
                          .is_good();

  std::ostringstream os;
  os << info;

  REQUIRE(os.str() == "alignment_info{score=-10, adapter_id=1, offset=2, "
                      "length=3, n_mismatches=4, n_ambiguous=5, m_type=good}");
}

TEST_CASE("stringmaker for ALN")
{
  std::ostringstream os;
  os << alignment_info{};

  REQUIRE(os.str() == "alignment_info{score=0, adapter_id=-1, offset=0, "
                      "length=0, n_mismatches=0, n_ambiguous=0, m_type=bad}");
}

TEST_CASE("stringmaker for MMNs")
{
  std::ostringstream os;
  os << adapterremoval::MMNs{ 1, 2 };

  REQUIRE(os.str() == "MMNs{mismatches=1, ambiguous=2}");
}

} // namespace adapterremoval
