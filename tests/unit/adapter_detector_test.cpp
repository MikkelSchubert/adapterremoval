// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_detector.hpp" // declarations
#include "fastq.hpp"            // for fastq
#include "logging.hpp"          // for log_capture
#include "sequence.hpp"         // for dna_sequence
#include "sequence_sets.hpp"    // for adapter_set
#include "simd.hpp"             // for instruction_set
#include "testing.hpp"          // for TEST_CASE, REQUIRE, ...
#include <initializer_list>     // for initializer_list

namespace adapterremoval {

namespace {
// Parameterize tests over supported SIMD instruction sets
#define PARAMETERIZE_IS GENERATE(from_range(simd::supported()))

#define REQUIRE_CONTAINS(a, b) REQUIRE_THAT((a), Catch::Matchers::Contains(b))

const double DEFAULT_MISMATCH_THRESHOLD = 1.0 / 6.0;

inline adapter_detector
simple_detector(std::initializer_list<string_view_pair> args)
{
  return adapter_detector{ args,
                           PARAMETERIZE_IS,
                           DEFAULT_MISMATCH_THRESHOLD,
                           adapter_detector::include_known::no };
}

using hits_vec = std::vector<adapter_detection_stats::hits>;

} // namespace

////////////////////////////////////////////////////////////////////////////////
// adapter_detection_stats::hits

TEST_CASE("adapter_detection_stats::hits equality operator")
{
  using hits = adapter_detection_stats::hits;

  CHECK(hits{} == hits{});
  CHECK(hits{ 1, 2, 3 } == hits{ 1, 2, 3 });

  CHECK_FALSE(hits{ 1, 2, 3 } == hits{});
  CHECK_FALSE(hits{ 1, 2, 3 } == hits{ 2, 2, 3 });
  CHECK_FALSE(hits{ 1, 2, 3 } == hits{ 1, 3, 3 });
  CHECK_FALSE(hits{ 1, 2, 3 } == hits{ 1, 2, 4 });
}

TEST_CASE("adapter_detection_stats::hits stringify")
{
  using hits = adapter_detection_stats::hits;

  CHECK(Catch::fallbackStringifier(hits{}) ==
        "adapter_detection_stats::hits{hits=0, aligned=0, mismatches=0}");
  CHECK(Catch::fallbackStringifier(hits{ 1, 2, 3 }) ==
        "adapter_detection_stats::hits{hits=1, aligned=2, mismatches=3}");
}

TEST_CASE("adapter_detection_stats::hits merge")
{
  using hits = adapter_detection_stats::hits;

  const hits src{ 1, 2, 3 };
  hits dst{ 4, 5, 6 };
  merge(dst, src);

  REQUIRE(dst == hits{ 5, 7, 9 });
}

////////////////////////////////////////////////////////////////////////////////
// adapter_detection_stats

TEST_CASE("merging adapter detection stats")
{
  auto ad = simple_detector({
    { "GTTATTTA", "ACGTGTTA" },
    { "ACGGACGT", "GGCAGTTA" },
  });

  REQUIRE(ad.sequences() == sequence_vec{ "ACGGACGT"_dna,
                                          "ACGTGTTA"_dna,
                                          "GGCAGTTA"_dna,
                                          "GTTATTTA"_dna });

  adapter_detection_stats stats_1;
  adapter_detection_stats stats_2;

  ad.detect_pe(stats_1,
               fastq{ "read1", "AGTGTTATTTAA" },
               fastq{ "read2", "ACTACGTGTAA" });
  ad.detect_pe(stats_1,
               fastq{ "read1", "ACGGACGTTTA" },
               fastq{ "read2", "ACTACGTGTTT" });
  ad.detect_pe(stats_2,
               fastq{ "read1", "ACGGACGTTTA" },
               fastq{ "read2", "ACTACGTGTTT" });

  REQUIRE(stats_1.mate_1() == hits_vec{ { 1, 8 }, {}, {}, { 1, 8 } });
  REQUIRE(stats_1.mate_2() == hits_vec{ {}, { 2, 16, 2 }, {}, {} });

  REQUIRE(stats_2.mate_1() == hits_vec{ { 1, 8 }, {}, {}, {} });
  REQUIRE(stats_2.mate_2() == hits_vec{ {}, { 1, 8, 1 }, {}, {} });

  stats_1.merge(stats_2);
  REQUIRE(stats_1.mate_1() == hits_vec{ { 2, 16 }, {}, {}, { 1, 8 } });
  REQUIRE(stats_1.mate_2() == hits_vec{ {}, { 3, 24, 3 }, {}, {} });
}

////////////////////////////////////////////////////////////////////////////////
// adapter_detector (common)

TEST_CASE("adapter sets are flattened and sorted")
{
  const sequence_vec expected{ "ACGGACGT"_dna,
                               "ACGTGTTA"_dna,
                               "GGCAGTTA"_dna,
                               "GTTATTTA"_dna };

  SECTION("PE adapters")
  {
    auto ad = simple_detector({
      { "GTTATTTA", "ACGTGTTA" },
      { "ACGGACGT", "GGCAGTTA" },
    });

    REQUIRE(ad.sequences() == expected);
  }

  SECTION("SE adapters")
  {
    auto ad = simple_detector({
      { "ACGTGTTA", {} },
      { "GTTATTTA", {} },
      { "GGCAGTTA", {} },
      { "ACGGACGT", {} },
    });

    REQUIRE(ad.sequences() == expected);
  }
}

TEST_CASE("adapter_detector including known adapters")
{
  std::initializer_list<string_view_pair> seqs{
    { "ACGTGTTA", "GTTATTTA" },
    { "ACGGACGT", "GGCAGTTA" },
  };

  // known adapters only
  adapter_detector ad_1{ {}, PARAMETERIZE_IS, DEFAULT_MISMATCH_THRESHOLD };
  REQUIRE(ad_1.size() > 10);

  // user specified adapters only
  auto ad_2 = simple_detector(seqs);
  REQUIRE(ad_2.size() == 4);

  // known and user specified adapters
  adapter_detector ad_3{ seqs, PARAMETERIZE_IS, DEFAULT_MISMATCH_THRESHOLD };
  REQUIRE(ad_3.size() == ad_1.size() + ad_2.size());
}

TEST_CASE("adapter_detector returns empty sequences if no stats")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "ACGTGTTA", "GTTATTTA" },
    { "ACGGACGT", "GGCAGTTA" },
  });

  REQUIRE(ad.select_best(stats) == sequence_pair{});
}

TEST_CASE("only unique adapters are collected")
{
  SECTION("PE adapters")
  {
    auto ad = simple_detector({
      { "ACGTGTTA", "ACGGACGT" },
      { "ACGGACGT", "ACGTGTTA" },
    });

    REQUIRE(ad.sequences() == sequence_vec{ "ACGGACGT"_dna, "ACGTGTTA"_dna });
  }

  SECTION("SE adapters")
  {
    auto ad = simple_detector({
      { "ACGTGTTA", {} },
      { "GTTATTTA", {} },
      { "ACGGACGT", {} },
      { "GTTATTTA", {} },
    });

    REQUIRE(ad.sequences() ==
            sequence_vec{ "ACGGACGT"_dna, "ACGTGTTA"_dna, "GTTATTTA"_dna });
  }
}

TEST_CASE("very short adapters are skipped")
{
  log::log_capture cap;

  auto ad = simple_detector({
    { "ACGTGTT", "ACGGACGT" },
    { "ACGTGTTA", "ACGGAC" },
  });

  REQUIRE(ad.sequences() == sequence_vec{ "ACGGACGT"_dna, "ACGTGTTA"_dna });
  REQUIRE_CONTAINS(cap.str(), "Adapter sequence 'ACGTGTT' is too short");
  REQUIRE_CONTAINS(cap.str(), "Adapter sequence 'ACGGAC' is too short");
}

TEST_CASE("adapters not detected")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "ACGTGTTA", "GTTATTTA" },
    { "ACGGACGT", "GGCAGTTA" },
  });

  SECTION("single end")
  {
    ad.detect_se(stats, fastq{ "read", "TTTTTTTTTT" });

    CHECK(stats.reads_1() == 1);
    CHECK(stats.reads_2() == 0);
    CHECK(stats.mate_1() == adapter_detection_stats::values{ 4 });
  }

  SECTION("paired end")
  {
    ad.detect_pe(stats,
                 fastq{ "read", "TTTTTTTTTT" },
                 fastq{ "read", "TTTTTTTTTT" });

    CHECK(stats.reads_1() == 1);
    CHECK(stats.reads_2() == 1);
    CHECK(stats.mate_1() == adapter_detection_stats::values{ 4 });
    CHECK(stats.mate_2() == adapter_detection_stats::values{ 4 });
  }
}

////////////////////////////////////////////////////////////////////////////////
// adapter_detector for SE adapters

TEST_CASE("reads with unique adapter matches")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "AGATCGGAAGAGCACACGTCT", {} },
    { "AGATCGGAAGAGCACACGTCTGAAC", {} },
    { "AGATCGGAAGAGCGTCGTGTAGGGA", {} },
  });

  SECTION("best match is partial sequence")
  {
    ad.detect_se(stats, fastq{ "read", "TTAAGATCGGAAGAGCGTCGTGTAG" });
    REQUIRE(stats.mate_1() == hits_vec{ {}, {}, { 1, 22 } });
  }

  SECTION("best match for complete sequence")
  {
    ad.detect_se(stats, fastq{ "read", "TTAAGATCGGAAGAGCACACGTCTGAAC" });
    REQUIRE(stats.mate_1() == hits_vec{ {}, { 1, 25 }, {} });
  }

  SECTION("best match for complete, shorter sequence")
  {
    ad.detect_se(stats, fastq{ "read", "AGATAGATCGGAAGAGCACACGTCTTTTTT" });
    REQUIRE(stats.mate_1() == hits_vec{ { 1, 21 }, {}, {} });
  }

  REQUIRE(stats.mate_2() == hits_vec{});
  REQUIRE(stats.reads_1() == 1);
  REQUIRE(stats.reads_2() == 0);
}

TEST_CASE("reads with shared prefixes")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "AGATCGGAAGAGCACACGTCT", {} },
    { "AGATCGGAAGAGCACACGTCTGAAC", {} },
    { "AGATCGGAAGAGCGTCGTGTAGGGA", {} },
  });

  SECTION("too short sequence")
  {
    ad.detect_se(stats, fastq{ "read", "AGTTTGAAGATCGG" });
    REQUIRE(stats.mate_1() == hits_vec{ {}, {}, {} });
  }

  SECTION("overlapping one")
  {
    ad.detect_se(stats, fastq{ "read", "AGATCGGAAGAGCGTC" });
    REQUIRE(stats.mate_1() == hits_vec{ {}, {}, { 1, 16 } });
  }

  SECTION("overlapping two")
  {
    ad.detect_se(stats, fastq{ "read", "AGTAGATCGGAAGAGCAC" });
    REQUIRE(stats.mate_1() == hits_vec{ { 1, 15 }, { 1, 15 }, {} });
  }

  SECTION("overlapping three")
  {
    ad.detect_se(stats, fastq{ "read", "AGTTTGAAGATCGGA" });
    REQUIRE(stats.mate_1() == hits_vec{ { 1, 8 }, { 1, 8 }, { 1, 8 } });
  }

  REQUIRE(stats.mate_2() == hits_vec{});
  REQUIRE(stats.reads_1() == 1);
  REQUIRE(stats.reads_2() == 0);
}

TEST_CASE("reads with shared mismatches")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({ { "AGATCGGAAGAGCACACGTCT", {} } });

  // 0 mismatches
  ad.detect_se(stats, fastq{ "read", "AGATCGGAAGAGCACACGTCT" });
  REQUIRE(stats.mate_1() == hits_vec{ { 1, 21, 0 } });

  // 1 mismatch
  ad.detect_se(stats, fastq{ "read", "AGATTGGAAGAGCACACGTCT" });
  REQUIRE(stats.mate_1() == hits_vec{ { 2, 42, 1 } });

  // 2 mismatches
  ad.detect_se(stats, fastq{ "read", "AGATTGGAAGAGCTCACGTCT" });
  REQUIRE(stats.mate_1() == hits_vec{ { 3, 63, 3 } });

  // 3 mismatches
  ad.detect_se(stats, fastq{ "read", "TGATTGGAAGAGCTCACGTCT" });
  REQUIRE(stats.mate_1() == hits_vec{ { 4, 84, 6 } });

  // 4 mismatches > 21/6
  ad.detect_se(stats, fastq{ "read", "TGTTTGGAAGAGCTCACGTCT" });
  REQUIRE(stats.mate_1() == hits_vec{ { 4, 84, 6 } });
}

TEST_CASE("reads with ambiguous bases")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({ { "AGATCGGAAGAGCACACGTCT", {} } });

  // 8 bases are required
  ad.detect_se(stats, fastq{ "read", "ANANCNGNAGNGN" });
  REQUIRE(stats.mate_1() == hits_vec{ {} });
  ad.detect_se(stats, fastq{ "read", "AGANCNGNAGNGN" });
  REQUIRE(stats.mate_1() == hits_vec{ { 1, 8 } });
}

////////////////////////////////////////////////////////////////////////////////
// adapter_detector for PE adapters

TEST_CASE("adapter detection for PE adapters")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "AGATCGGAAGAGCACACGTCT", {} },
    { "AGATCGGAAGAGCACACGTCTGAAC", {} },
    { "AGATCGGAAGAGCGTCGTGTAGGGA", {} },
  });

  SECTION("matching different sequences")
  {
    ad.detect_pe(stats,
                 fastq{ "read1", "AGATCGGAAGAGCGTC" },
                 fastq{ "read2", "AGATCGGAAGAGCGTC" });
    REQUIRE(stats.mate_1() == hits_vec{ {}, {}, { 1, 16 } });
    REQUIRE(stats.mate_2() == stats.mate_1());
  }

  SECTION("matching different sequences")
  {
    ad.detect_pe(stats,
                 fastq{ "read1", "AGATCGGAAGAGCGTC" },
                 fastq{ "read2", "AGATCGGAAGAGCACAC" });
    REQUIRE(stats.mate_1() == hits_vec{ {}, {}, { 1, 16 } });
    REQUIRE(stats.mate_2() == hits_vec{ { 1, 17 }, { 1, 17 }, {} });
  }

  REQUIRE(stats.reads_1() == 1);
  REQUIRE(stats.reads_2() == 1);
};

TEST_CASE("mate 1 and mate 2 read counts are independent")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "AGATCGGAAGAGCACACGTCT", {} },
    { "AGATCGGAAGAGCACACGTCTGAAC", {} },
    { "AGATCGGAAGAGCGTCGTGTAGGGA", {} },
  });

  ad.detect_pe(stats,
               fastq{ "read1", "AGATCGGAAGAGCGTC" },
               fastq{ "read2", "AGATCGGAAGAGCACAC" });
  ad.detect_se(stats, fastq{ "read", "AGATCGGAAGAGCACACGTCTGAAC" });

  REQUIRE(stats.mate_1() == hits_vec{ {}, { 1, 25 }, { 1, 16 } });
  REQUIRE(stats.mate_2() == hits_vec{ { 1, 17 }, { 1, 17 }, {} });
  REQUIRE(stats.reads_1() == 2);
  REQUIRE(stats.reads_2() == 1);
};

////////////////////////////////////////////////////////////////////////////////
// adapter_detector selection

TEST_CASE("selection requires 10 hits")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "AGATCGGAAGAGCACACGTCT", {} },
  });

  log::log_capture _;
  REQUIRE(ad.select_best({ 9, { { 9 } } }) == sequence_pair{});
  REQUIRE(ad.select_best({ 9, {}, { { 9 } } }) == sequence_pair{});

  REQUIRE(ad.select_best({ 10, { { 10 } } }) ==
          sequence_pair{ "AGATCGGAAGAGCACACGTCT"_dna, {} });
  REQUIRE(ad.select_best({ 10, {}, { { 10 } } }) ==
          sequence_pair{ {}, "AGATCGGAAGAGCACACGTCT"_dna });

  REQUIRE(
    ad.select_best({ 10, { { 10 } }, { { 10 } } }) ==
    sequence_pair{ "AGATCGGAAGAGCACACGTCT"_dna, "AGATCGGAAGAGCACACGTCT"_dna });
}

TEST_CASE("selection requires 1/0.1 percent sequences")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "AGATCGGAAGAGCACACGTCT", {} },
  });

  log::log_capture _;
  // Common sequences require >= 1% hits
  REQUIRE(ad.select_best({ 10000, { { 99, 9900, 990 } } }) == sequence_pair{});
  REQUIRE(ad.select_best({ 10000, { { 100, 10000, 1000 } } }) ==
          sequence_pair{ "AGATCGGAAGAGCACACGTCT"_dna, {} });

  // Rare sequences rare >= 0.1% sequences and <= 5% errors
  REQUIRE(ad.select_best({ 100000, { { 99, 9900, 990 } } }) == sequence_pair{});
  REQUIRE(ad.select_best({ 100000, { { 99, 9900, 500 } } }) == sequence_pair{});
  REQUIRE(ad.select_best({ 100000, { { 100, 10000, 1000 } } }) ==
          sequence_pair{});
  REQUIRE(ad.select_best({ 100000, { { 100, 10000, 500 } } }) ==
          sequence_pair{ "AGATCGGAAGAGCACACGTCT"_dna, {} });
}

TEST_CASE("selection prefers the most aligned bases for same hits")
{
  adapter_detection_stats stats;
  auto ad = simple_detector({
    { "AGATCGGAAGAGCACACGTCT", {} },
    { "AGATCGGAAGAGCACACGTCTGAAC", {} },
    { "AGATCGGAAGAGCGTCGTGTAGGGA", {} },
  });

  log::log_capture _;
  // Prefer sequence with the most hits
  REQUIRE(ad.select_best({ 100, { { 30, 400 }, { 10, 50 }, { 40, 300 } } }) ==
          sequence_pair{ "AGATCGGAAGAGCGTCGTGTAGGGA"_dna, {} });
  // Prefer sequence with the most bases, for same number of hits
  REQUIRE(ad.select_best({ 100, { { 40, 400 }, { 10, 50 }, { 40, 300 } } }) ==
          sequence_pair{ "AGATCGGAAGAGCACACGTCT"_dna, {} });
}

} // namespace adapterremoval
