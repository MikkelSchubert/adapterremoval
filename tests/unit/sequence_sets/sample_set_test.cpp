// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "commontypes.hpp"   // for barcode_orientation
#include "errors.hpp"        // for parsing_error
#include "linereader.hpp"    // for vec_reader
#include "read_group.hpp"    // for read_group
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for read_group
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <initializer_list>  // for initializer_list
#include <sstream>           // for ostringstream
#include <string>            // for string, operator==
#include <string_view>       // for string_view
#include <vector>            // for vector, operator==

using Contains = Catch::Matchers::StdString::ContainsMatcher;

namespace adapterremoval {

namespace {

constexpr auto CONFIG_SE = barcode_config{}.paired_end_mode(false);
constexpr auto CONFIG_SE_MULTIPLE =
  barcode_config{ CONFIG_SE }.allow_multiple_barcodes();
constexpr auto CONFIG_PE = barcode_config{}.paired_end_mode(true);
constexpr auto CONFIG_PE_MULTIPLE =
  barcode_config{ CONFIG_PE }.allow_multiple_barcodes();
constexpr auto CONFIG_PE_EXPLICIT =
  barcode_config{ CONFIG_PE }.orientation(barcode_table_orientation::explicit_);

sample_set
sample_set_pe(std::initializer_list<std::string_view> lines)
{
  barcode_config config;
  config.paired_end_mode(true);

  return { lines, config };
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- constructor

TEST_CASE("constructor", "[sample_set]")
{
  sample_set ss;

  CHECK(ss.size() == 1);
  CHECK(std::vector<sample>{ ss.begin(), ss.end() } ==
        std::vector<sample>{ sample{} });
  CHECK(ss.adapters() == adapter_set{});
  CHECK(ss.readgroup() == read_group{});

  {
    sample unidentified;
    unidentified.set_read_group(read_group{ "DS:unidentified" });
    CHECK(ss.unidentified() == unidentified);
  }

  CHECK(ss.samples() == std::vector{ sample{} });
}

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- load

TEST_CASE("Loading single barcodes", "[sample_set]")
{
  sample sample_1{ "sample_1",
                   "ACGTA"_dna,
                   ""_dna,
                   barcode_orientation::unspecified };
  sample_1.set_read_group(read_group{ "ID:sample_1\tSM:sample_1\tBC:ACGTA" });
  sample sample_2{ "sample_2",
                   "TGCAT"_dna,
                   ""_dna,
                   barcode_orientation::unspecified };
  sample_2.set_read_group(read_group{ "ID:sample_2\tSM:sample_1\tBC:TGCAT" });

  const sample_set ss{
    "sample_1 ACGTA",
    "sample_2 TGCAT",
  };

  CHECK(ss.size() == 2);
  CHECK(ss.samples() == std::vector{ sample_1, sample_2 });
  CHECK(ss.at(0) == sample_1);
  CHECK(ss.at(1) == sample_2);
  CHECK(ss.at(0).at(0).read_group_ ==
        read_group{ "ID:sample_1\tSM:sample_1\tBC:ACGTA" });
  CHECK(ss.at(1).at(0).read_group_ ==
        read_group{ "ID:sample_2\tSM:sample_2\tBC:TGCAT" });
}

TEST_CASE("Loading two barcodes per sample", "[sample_set]")
{
  sample sample_1{ "sample_1",
                   "ACGTA"_dna,
                   "TTGTC"_dna,
                   barcode_orientation::unspecified };
  sample_1.set_read_group(read_group{ "ID:sample_1\tSM:sample_1" });
  sample sample_2{ "sample_2",
                   "TGCAT"_dna,
                   "CCGAT"_dna,
                   barcode_orientation::unspecified };
  sample_2.set_read_group(read_group{ "ID:sample_2\tSM:sample_1" });

  const sample_set ss{
    "sample_1 ACGTA TTGTC",
    "sample_2 TGCAT CCGAT",
  };

  CHECK(ss.size() == 2);
  CHECK(ss.samples() == std::vector{ sample_1, sample_2 });
  CHECK(ss.at(0) == sample_1);
  CHECK(ss.at(1) == sample_2);
  CHECK(ss.at(0).at(0).read_group_ ==
        read_group{ "ID:sample_1\tSM:sample_1\tBC:ACGTA-TTGTC" });
  CHECK(ss.at(1).at(0).read_group_ ==
        read_group{ "ID:sample_2\tSM:sample_2\tBC:TGCAT-CCGAT" });
}

TEST_CASE("Loading multiple barcodes per sample", "[sample_set]")
{
  sample sample_1{ "sample_1",
                   "ACGTA"_dna,
                   "TTGTC"_dna,
                   barcode_orientation::unspecified };
  sample_1.add_barcodes(dna_sequence{ "GGTTG" },
                        dna_sequence{ "ACGTA" },
                        barcode_orientation::unspecified);
  sample_1.set_read_group(read_group{ "ID:sample_1\tSM:sample_1" });
  sample sample_2{ "sample_2",
                   "TGCAT"_dna,
                   "CCGAT"_dna,
                   barcode_orientation::unspecified };
  sample_2.set_read_group(read_group{ "ID:sample_2\tSM:sample_2" });

  REQUIRE(sample_1.at(0).read_group_ ==
          read_group{ "ID:sample_1.1\tSM:sample_1\tBC:ACGTA-TTGTC" });
  REQUIRE(sample_1.at(1).read_group_ ==
          read_group{ "ID:sample_1.2\tSM:sample_1\tBC:GGTTG-ACGTA" });
  REQUIRE(sample_2.at(0).read_group_ ==
          read_group{ "ID:sample_2\tSM:sample_2\tBC:TGCAT-CCGAT" });

  const sample_set ss{ {
                         "sample_1 ACGTA TTGTC",
                         "sample_2 TGCAT CCGAT",
                         "sample_1 GGTTG ACGTA",
                       },
                       GENERATE(CONFIG_SE_MULTIPLE, CONFIG_PE_MULTIPLE) };

  CHECK(ss.size() == 2);
  CHECK(ss.samples() == std::vector{ sample_1, sample_2 });
  CHECK(ss.at(0) == sample_1);
  CHECK(ss.at(1) == sample_2);
}

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- load -- invalid barcode tables

TEST_CASE("Empty barcode table", "[sample_set]")
{
  CHECK_THROWS_MESSAGE(sample_set({}),
                       parsing_error,
                       "No samples/barcodes found in table");
}

TEST_CASE("Invalid sample name", "[sample_set]")
{
  CHECK_THROWS_MESSAGE(sample_set{ "sample-1 ACGT" },
                       parsing_error,
                       "The sample name 'sample-1' is not a valid sample name; "
                       "only letters ('a' to 'z' and 'A' to 'Z'), numbers (0 "
                       "to 9) and underscores (_) are allowed.");
}

TEST_CASE("Reserved sample name", "[sample_set]")
{
  CHECK_THROWS_MESSAGE(sample_set{ "unidentified ACGT" },
                       parsing_error,
                       "The sample name 'unidentified' is a reserved name, and "
                       "cannot be used!");
}

TEST_CASE("Degenerate bases in barcodes", "[sample_set]")
{
  CHECK_THROWS_MESSAGE(sample_set{ "sample_1 ACNT TGGA" },
                       parsing_error,
                       "Unsupported character found in mate 1 barcode sequence "
                       "'ACNT'. Only bases A, C, G, and T, are supported; "
                       "please fix before continuing");
  CHECK_THROWS_MESSAGE(sample_set{ "sample_1 YCAT TGGA" },
                       parsing_error,
                       "Unsupported character found in mate 1 barcode sequence "
                       "'YCAT'. Only bases A, C, G, and T, are supported; "
                       "please fix before continuing");

  CHECK_THROWS_MESSAGE(sample_set{ "sample_1 ACCT TNGA" },
                       parsing_error,
                       "Unsupported character found in mate 2 barcode sequence "
                       "'TNGA'. Only bases A, C, G, and T, are supported; "
                       "please fix before continuing");
  CHECK_THROWS_MESSAGE(sample_set{ "sample_1 ACCT TAGX" },
                       parsing_error,
                       "Unsupported character found in mate 2 barcode sequence "
                       "'TAGX'. Only bases A, C, G, and T, are supported; "
                       "please fix before continuing");
}

TEST_CASE("Overlapping SE barcodes fail #1", "[sample_set]")
{
  const std::string message = "Sample \'sample1\' (ACGT) and sample "
                              "\'sample2\' (ACGT) have overlapping barcodes. "
                              "Please remove any duplicate entries from the "
                              "barcode table before continuing";

  const std::initializer_list<std::string_view> barcodes{
    "sample1 ACGT",
    "sample2 ACGT",
  };

  CHECK_THROWS_MESSAGE(sample_set(barcodes, CONFIG_SE), parsing_error, message);
  CHECK_THROWS_MESSAGE(sample_set(barcodes, CONFIG_PE), parsing_error, message);
}

TEST_CASE("Overlapping SE barcodes fail #2", "[sample_set]")
{
  const std::initializer_list<std::string_view> barcodes{
    "sample1 ACGT CGTG",
    "sample2 ACGT TATA",
  };

  CHECK_THROWS_MESSAGE(sample_set(barcodes),
                       parsing_error,
                       "Sample \'sample1\' (ACGT-CGTG) and sample \'sample2\' "
                       "(ACGT-TATA) have overlapping barcodes. Note that "
                       "AdapterRemoval cannot distinguish these barcodes in "
                       "single-end mode, even though the second barcodes "
                       "differ. Please remove any duplicate entries from the "
                       "barcode table before continuing");
  CHECK_NOTHROW(sample_set_pe(barcodes));
}

TEST_CASE("Overlapping PE barcodes fail", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(sample_set_pe({
                           "sample1 ACGT CGTG",
                           "sample2 ACGT CGTG",
                         }),
                         parsing_error,
                         "Sample \'sample1\' (ACGT-CGTG) and sample "
                         "\'sample2\' (ACGT-CGTG) have overlapping barcodes. "
                         "Please remove any duplicate entries from the barcode "
                         "table before continuing");
}

TEST_CASE("Variable length SE barcodes fail #1", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(sample_set_pe({
                           "sample1 ACGTA",
                           "sample2 TGCT",
                         }),
                         parsing_error,
                         "Inconsistent mate 1 barcode lengths found: Last "
                         "barcode was 5 base-pairs long, but barcode \'TGCT\' "
                         "is 4 base-pairs long. Variable length barcodes are "
                         "not supported");
}

TEST_CASE("Variable length SE barcodes fail #2", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(sample_set_pe({
                           "sample1 ACGT",
                           "sample2 TGCTA",
                         }),
                         parsing_error,
                         "Inconsistent mate 1 barcode lengths found: Last "
                         "barcode was 4 base-pairs long, but barcode \'TGCTA\' "
                         "is 5 base-pairs long. Variable length barcodes are "
                         "not supported");
}

TEST_CASE("Variable length PE barcodes fail #1", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(sample_set_pe({
                           "sample1 ACGTT CGTG",
                           "sample2 CGTG ACGT",
                         }),
                         parsing_error,
                         "Inconsistent mate 1 barcode lengths found: Last "
                         "barcode was 5 base-pairs long, but barcode \'CGTG\' "
                         "is 4 base-pairs long. Variable length barcodes are "
                         "not supported");
}

TEST_CASE("Variable length PE barcodes fail #2", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(sample_set_pe({
                           "sample1 ACGT CGTGT",
                           "sample2 CGTG ACGT",
                         }),
                         parsing_error,
                         "Inconsistent mate 2 barcode lengths found: Last "
                         "barcode was 5 base-pairs long, but barcode \'ACGT\' "
                         "is 4 base-pairs long. Variable length barcodes are "
                         "not supported");
}

TEST_CASE("Variable length PE barcodes fail #3", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(sample_set_pe({
                           "sample1 ACGT CGTG",
                           "sample2 CGTGA ACGT",
                         }),
                         parsing_error,
                         "Inconsistent mate 1 barcode lengths found: Last "
                         "barcode was 4 base-pairs long, but barcode \'CGTGA\' "
                         "is 5 base-pairs long. Variable length barcodes are "
                         "not supported");
}

TEST_CASE("Variable length PE barcodes fail #4", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(sample_set_pe({
                           "sample1 ACGT CGTG",
                           "sample2 CGTGA ACGTA",
                         }),
                         parsing_error,
                         "Inconsistent mate 1 barcode lengths found: Last "
                         "barcode was 4 base-pairs long, but barcode \'CGTGA\' "
                         "is 5 base-pairs long. Variable length barcodes are "
                         "not supported");
}

TEST_CASE("Duplicate sample names", "[sample_set]")
{
  const std::initializer_list<std::string_view> lines{
    "sample1 ACGT CGTG",
    "sample1 CGTG ACGT",
  };

  SECTION("invalid without multiple")
  {
    CHECK_THROWS_MESSAGE(sample_set(lines, GENERATE(CONFIG_SE, CONFIG_PE)),
                         parsing_error,
                         "Duplicate sample name 'sample1'; multiple barcodes "
                         "per samples is not enabled. Either ensure that all "
                         "sample names are unique or use --multiple-barcodes "
                         "to map multiple barcodes to a single sample");
  }

  // Verify that test data is valid; actual output is tested below
  SECTION("valid with multiple")
  {
    const auto config = GENERATE(CONFIG_SE_MULTIPLE, CONFIG_PE_MULTIPLE);
    CHECK_NOTHROW(sample_set(lines, config));
  }
}

TEST_CASE("Duplicate sample names with different explicit adapters",
          "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 ACCC CTAA forward",
    GENERATE(
      // duplicate barcodes, right orientation
      "sample_1 ACCC CTAA reverse",
      // right barcodes, wrong orientation
      "sample_1 CTAA ACCC forward",
      // wrong barcodes, right orientation
      "sample_1 GCAG TGAC reverse"),
  };

  CHECK_THROWS_MESSAGE(sample_set(lines, CONFIG_PE_EXPLICIT),
                       parsing_error,
                       "Duplicate sample name 'sample_1'; multiple barcodes "
                       "per samples is not enabled. Either ensure that all "
                       "sample names are unique or use --multiple-barcodes to "
                       "map multiple barcodes to a single sample");
}

TEST_CASE("Names for barcodes must be distinct", "[sample_set]")
{
  const auto config =
    GENERATE(CONFIG_SE, CONFIG_PE, CONFIG_SE_MULTIPLE, CONFIG_PE_MULTIPLE);
  const std::initializer_list<std::string_view> lines{
    "sample_1 ACGTA TTGTC",
    "sample_2 TGCAT CCGAT",
    "SAMPLE_1 GGTTG ACGTA",
  };

  CHECK_THROWS_MESSAGE(sample_set(lines, config),
                       parsing_error,
                       "Samples with names 'SAMPLE_1' and 'sample_1' differ "
                       "only by case. Either use the exact same name for both, "
                       "if they the same sample, or give them distinct names");
}

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- load with mixed orientation barcodes

TEST_CASE("Barcodes in forward orientation", "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 CTTGCCCT ACGTTATT",
    "sample_2 CGCCGATG TGCACGGG",
  };

  sample sample_1{ "sample_1",
                   "CTTGCCCT"_dna,
                   "ACGTTATT"_dna,
                   barcode_orientation::forward };
  sample_1.add_barcodes(dna_sequence{ "ACGTTATT" },
                        dna_sequence{ "CTTGCCCT" },
                        barcode_orientation::reverse);
  sample_1.set_read_group(read_group{});

  sample sample_2{ "sample_2",
                   "CGCCGATG"_dna,
                   "TGCACGGG"_dna,
                   barcode_orientation::forward };
  sample_2.add_barcodes(dna_sequence{ "TGCACGGG" },
                        dna_sequence{ "CGCCGATG" },
                        barcode_orientation::reverse);
  sample_2.set_read_group(read_group{});

  std::vector<sample> samples = { sample_1, sample_2 };

  barcode_config config;
  config.paired_end_mode(GENERATE(true, false))
    .allow_multiple_barcodes(GENERATE(true, false))
    .orientation(barcode_table_orientation::forward);

  const sample_set ss{ lines, config };
  REQUIRE(ss.samples() == samples);
}

TEST_CASE("Barcodes in reverse orientation", "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 CTTGCCCT ACGTTATT",
    "sample_2 CGCCGATG TGCACGGG",
  };

  sample sample_1{ "sample_1",
                   "CTTGCCCT"_dna,
                   "ACGTTATT"_dna,
                   barcode_orientation::reverse };
  sample_1.add_barcodes(dna_sequence{ "ACGTTATT" },
                        dna_sequence{ "CTTGCCCT" },
                        barcode_orientation::forward);
  sample_1.set_read_group(read_group{});

  sample sample_2{ "sample_2",
                   "CGCCGATG"_dna,
                   "TGCACGGG"_dna,
                   barcode_orientation::reverse };
  sample_2.add_barcodes(dna_sequence{ "TGCACGGG" },
                        dna_sequence{ "CGCCGATG" },
                        barcode_orientation::forward);
  sample_2.set_read_group(read_group{});

  std::vector<sample> samples = { sample_1, sample_2 };

  barcode_config config;
  config.paired_end_mode(GENERATE(true, false))
    .allow_multiple_barcodes(GENERATE(true, false))
    .orientation(barcode_table_orientation::reverse);

  const sample_set ss{ lines, config };
  REQUIRE(ss.samples() == samples);
}

TEST_CASE("Barcodes in explicit orientation", "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 ACCC CTAA forwArd", "sample_2 GACA ACAG fWd",
    "sample_3 TTCA TCCA +",       "sample_4 CGAG CCAT reversE",
    "sample_5 TAAG GCTT rev",     "sample_6 ATAT CTCT -",
  };

  const auto fwd = barcode_orientation::forward;
  const auto rev = barcode_orientation::reverse;

  std::vector<sample> samples = {
    sample{ "sample_1", "ACCC"_dna, "CTAA"_dna, fwd },
    sample{ "sample_2", "GACA"_dna, "ACAG"_dna, fwd },
    sample{ "sample_3", "TTCA"_dna, "TCCA"_dna, fwd },
    sample{ "sample_4", "CGAG"_dna, "CCAT"_dna, rev },
    sample{ "sample_5", "TAAG"_dna, "GCTT"_dna, rev },
    sample{ "sample_6", "ATAT"_dna, "CTCT"_dna, rev },
  };

  for (auto& sample : samples) {
    sample.set_read_group(read_group{});
  }

  const sample_set ss{ lines, CONFIG_PE_EXPLICIT };
  REQUIRE(ss.samples() == samples);
}

TEST_CASE("Barcodes in explicit orientation for same sample", "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 ACCC CTAA forwArd",
    "sample_1 CTAA ACCC reverse",
  };

  const auto fwd = barcode_orientation::forward;
  const auto rev = barcode_orientation::reverse;
  sample sample_1{ "sample_1",
                   dna_sequence{ "ACCC" },
                   dna_sequence{ "CTAA" },
                   fwd };
  sample_1.add_barcodes(dna_sequence{ "CTAA" }, dna_sequence{ "ACCC" }, rev);
  sample_1.set_read_group(read_group{});

  const sample_set ss{ lines, CONFIG_PE_EXPLICIT };
  REQUIRE(ss.samples() == std::vector{ sample_1 });
}

TEST_CASE("Different names for barcodes in opposite orientation",
          "[sample_set]")
{
  const auto fwd = barcode_orientation::forward;
  const auto rev = barcode_orientation::reverse;
  std::vector<sample> samples = {
    sample{ "sample_1", "ACCC"_dna, "CTAA"_dna, fwd },
    sample{ "sample_2", "CTAA"_dna, "ACCC"_dna, rev },
  };

  for (auto& sample : samples) {
    sample.set_read_group(read_group{});
  }

  const sample_set ss{ {
                         "sample_1 ACCC CTAA forward",
                         "sample_2 CTAA ACCC reverse",
                       },
                       CONFIG_PE_EXPLICIT };
  REQUIRE(ss.samples() == samples);
}

TEST_CASE("Different names for barcodes in same orientation", "[sample_set]")
{
  const auto fwd = barcode_orientation::forward;
  std::vector<sample> samples = {
    sample{ "sample_1", "ACCC"_dna, "CTAA"_dna, fwd },
    sample{ "sample_2", "CTAA"_dna, "ACCC"_dna, fwd },
  };

  for (auto& sample : samples) {
    sample.set_read_group(read_group{});
  }

  const sample_set ss{ {
                         "sample_1 ACCC CTAA forward",
                         "sample_2 CTAA ACCC forward",
                       },
                       CONFIG_PE_EXPLICIT };
  REQUIRE(ss.samples() == samples);
}

TEST_CASE("Multiple barcodes in forward orientation", "[sample_set]")
{
  using seq = dna_sequence;
  const auto fwd = barcode_orientation::forward;
  const auto rev = barcode_orientation::reverse;

  std::initializer_list<std::string_view> lines{
    "sample_1 CTTGCCCT ACGTTATT",
    "sample_1 CGCCGATG TGCACGGG",
  };

  sample sample_1{ "sample_1", seq{ "CTTGCCCT" }, seq{ "ACGTTATT" }, fwd };
  sample_1.add_barcodes(seq{ "ACGTTATT" }, seq{ "CTTGCCCT" }, rev);
  sample_1.add_barcodes(seq{ "CGCCGATG" }, seq{ "TGCACGGG" }, fwd);
  sample_1.add_barcodes(seq{ "TGCACGGG" }, seq{ "CGCCGATG" }, rev);
  sample_1.set_read_group(read_group{});

  std::vector<sample> samples = { sample_1 };

  barcode_config config;
  config.allow_multiple_barcodes(true).orientation(
    barcode_table_orientation::forward);

  const sample_set ss{ lines, config };
  REQUIRE(ss.samples() == samples);
}

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- load with mixed orientation barcodes -- parsing errors

TEST_CASE("Double barcodes are required for barcode orientation",
          "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 CTTGCCCT",
    "sample_2 CGCCGATG",
  };

  barcode_config config;
  config.paired_end_mode(GENERATE(true, false))
    .allow_multiple_barcodes(GENERATE(true, false))
    .orientation(GENERATE(barcode_table_orientation::forward,
                          barcode_table_orientation::reverse));

  CHECK_THROWS_MESSAGE(sample_set(lines, config),
                       parsing_error,
                       "Error at line 1: Expected at least 3 columns, but "
                       "found 2 column(s)");
}

TEST_CASE("Overlap between specified and derived barcodes fails",
          "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 CTTGCCCT ACGTTATT",
    "sample_2 CGCCGATG TGCACGGG",
    "sample_3 ACGTTATT CTTGCCCT",
  };

  barcode_config config;
  config.paired_end_mode(GENERATE(true, false))
    .allow_multiple_barcodes(GENERATE(true, false));

  SECTION("forward")
  {
    config.orientation(barcode_table_orientation::forward);
    CHECK_THROWS_MESSAGE(sample_set(lines, config),
                         parsing_error,
                         "Sample 'sample_1' (ACGTTATT-CTTGCCCT; reverse) and "
                         "sample 'sample_3' (ACGTTATT-CTTGCCCT; forward) have "
                         "overlapping barcodes. Please remove any duplicate "
                         "entries from the barcode table before continuing");
  }

  SECTION("reverse")
  {
    config.orientation(barcode_table_orientation::reverse);
    CHECK_THROWS_MESSAGE(sample_set(lines, config),
                         parsing_error,
                         "Sample 'sample_1' (ACGTTATT-CTTGCCCT; forward) and "
                         "sample 'sample_3' (ACGTTATT-CTTGCCCT; reverse) have "
                         "overlapping barcodes. Please remove any duplicate "
                         "entries from the barcode table before continuing");
  }
}

TEST_CASE("Overlap between explicit barcodes fails", "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 CTTGCCCT ACGTTATT forward",
    GENERATE("sample_1 CTTGCCCT ACGTTATT forward",
             "sample_2 CTTGCCCT ACGTTATT forward",
             "sample_2 CTTGCCCT ACGTTATT reverse"),
  };

  barcode_config config;
  config.allow_multiple_barcodes(true).orientation(
    barcode_table_orientation::explicit_);

  CHECK_THROWS_WITH(sample_set(lines, config),
                    Catch::Contains("have overlapping barcodes"));
}

TEST_CASE("Orientation must be present for explicit barcodes", "[sample_set]")
{
  std::initializer_list<std::string_view> lines{
    "sample_1 CTTG CCCT forward",
    "sample_2 CGCC GATG ",
  };

  CHECK_THROWS_MESSAGE(sample_set(lines, CONFIG_PE_EXPLICIT),
                       parsing_error,
                       "Error at line 2: Expected at least 4 columns, but "
                       "found 3 column(s)");
}

TEST_CASE("Orientation must be a valid value", "[sample_set]")
{
  CHECK_THROWS_MESSAGE(sample_set({ "x CTTG CCCT ." }, CONFIG_PE_EXPLICIT),
                       parsing_error,
                       "Invalid barcode orientation for sample 'x': '.'");
  CHECK_THROWS_MESSAGE(sample_set({ "x CTTG CCCT forwar" }, CONFIG_PE_EXPLICIT),
                       parsing_error,
                       "Invalid barcode orientation for sample 'x': 'forwar'");
}

////////////////////////////////////////////////////////////////////////////////
// Sample set -- set read group

TEST_CASE("set read group for non-demultiplexing", "[sample_set]")
{
  sample_set ss;
  SECTION("without SM")
  {
    ss.set_read_group("DS:demux");
    CHECK(ss.readgroup() == read_group{ "DS:demux" });
    CHECK(ss.at(0).at(0).read_group_ == read_group{ "ID:1\tDS:demux" });
    CHECK(ss.unidentified().at(0).read_group_ ==
          read_group{ "ID:1\tDS:unidentified" });
  }

  SECTION("with SM")
  {
    ss.set_read_group("DS:demux\tSM:my_sample");
    CHECK(ss.readgroup() == read_group{ "DS:demux\tSM:my_sample" });
    CHECK(ss.at(0).at(0).read_group_ ==
          read_group{ "ID:1\tDS:demux\tSM:my_sample" });
    CHECK(ss.unidentified().at(0).read_group_ ==
          read_group{ "ID:1\tDS:unidentified" });
  }

  SECTION("with LB")
  {
    ss.set_read_group("LB:lib1");
    CHECK(ss.readgroup() == read_group{ "LB:lib1" });
    CHECK(ss.at(0).at(0).read_group_ == read_group{ "ID:1\tLB:lib1" });
    CHECK(ss.unidentified().at(0).read_group_ ==
          read_group{ "ID:1\tLB:lib1\tDS:unidentified" });
  }
}

TEST_CASE("set read group for demultiplexing", "[sample_set]")
{
  auto ss = sample_set_pe({
    "sample_1 ACGTA TTGTC",
    "sample_2 TGCAT CCGAT",
  });

  ss.set_read_group(
    GENERATE("DS:demux", "DS:demux\tSM:unused", "ID:unused\tDS:demux"));

  CHECK(ss.at(0).at(0).read_group_ ==
        read_group{ "ID:sample_1\tDS:demux\tSM:sample_1\tBC:ACGTA-TTGTC" });
  CHECK(ss.at(1).at(0).read_group_ ==
        read_group{ "ID:sample_2\tDS:demux\tSM:sample_2\tBC:TGCAT-CCGAT" });
}

TEST_CASE("set read group for demultiplexing before loading", "[sample_set]")
{
  sample_set ss;
  ss.set_read_group(
    GENERATE("DS:demux", "DS:demux\tSM:unused", "ID:unused\tDS:demux"));

  vec_reader reader{ "sample_1 ACGTA TTGTC" };
  ss.load(reader, barcode_config{});

  CHECK(ss.at(0).at(0).read_group_ ==
        read_group{ "ID:sample_1\tDS:demux\tSM:sample_1\tBC:ACGTA-TTGTC" });
}

////////////////////////////////////////////////////////////////////////////////
// Sample set -- set adapters

TEST_CASE("set adapters for non-demultiplexing", "[sample_set]")
{
  sample_set ss;
  const adapter_set as{ { "ACGTA", "TGGAT" } };
  ss.set_adapters(as);
  CHECK(ss.adapters() == as);

  CHECK(ss.at(0).at(0).adapters() == as);
  CHECK(ss.unidentified().at(0).adapters() == as);
}

TEST_CASE("set adapters for demultiplexing", "[sample_set]")
{
  auto ss = sample_set_pe({
    "sample_1 ACGTA TTGTC",
    "sample_2 TGCAT CCGAT",
  });

  adapter_set as{ { "ACGTA", "TGGAT" } };
  ss.set_adapters(as);

  CHECK(ss.adapters() == as);
  CHECK(ss.at(0).at(0).adapters() ==
        adapter_set{ { "GACAAACGTA", "TACGTTGGAT" } });
  CHECK(ss.at(1).at(0).adapters() ==
        adapter_set{ { "ATCGGACGTA", "ATGCATGGAT" } });
}

TEST_CASE("set adapters for demultiplexing before loading", "[sample_set]")
{
  sample_set ss;
  ss.set_adapters({ { "ACGTA", "TGGAT" } });

  vec_reader reader{ { "sample_1 ACGTA TTGTC" } };
  ss.load(reader, barcode_config{});

  CHECK(ss.at(0).at(0).adapters() ==
        adapter_set{ { "GACAAACGTA", "TACGTTGGAT" } });
}

////////////////////////////////////////////////////////////////////////////////
// Sample set -- get sequences

TEST_CASE("get samples", "[sample_set]")
{
  auto ss = sample_set_pe({
    "sample_1 ACGTA TTGTC",
    "sample_2 TGCAT CCGAT",
    "sample_3 GGTTG ACGTA",
  });

  CHECK(ss.samples() == std::vector<sample>{ ss.begin(), ss.end() });
}

////////////////////////////////////////////////////////////////////////////////
// Sample set -- uninitialized adapters

TEST_CASE("set and clear uninitialized adapters in sample set", "[sample_set]")
{
  auto ss = sample_set_pe({
    "sample_1 TT GG",
    "sample_2 AA CC",
  });

  adapter_set as_1{ { "TGGAT", "ACGTA" } };
  ss.set_adapters(as_1);

  REQUIRE(ss.adapters() == as_1);
  REQUIRE(ss.at(0).at(0).adapters() == adapter_set{ { "CCTGGAT", "AAACGTA" } });
  REQUIRE(ss.at(1).at(0).adapters() == adapter_set{ { "GGTGGAT", "TTACGTA" } });
  REQUIRE_THROWS_AS(ss.uninitialized_adapters(), assert_failed);

  ss.flag_uninitialized_adapters();
  REQUIRE_THROWS_AS(ss.adapters(), assert_failed);
  REQUIRE_THROWS_AS(ss.at(0).at(0).adapters(), assert_failed);
  REQUIRE_THROWS_AS(ss.at(1).at(0).adapters(), assert_failed);
  REQUIRE(ss.uninitialized_adapters() == as_1);

  adapter_set as_2{ { "ACGTA", "TGGAT" } };
  ss.set_adapters(as_2);
  REQUIRE(ss.adapters() == as_2);
  REQUIRE(ss.at(0).at(0).adapters() == adapter_set{ { "CCACGTA", "AATGGAT" } });
  REQUIRE(ss.at(1).at(0).adapters() == adapter_set{ { "GGACGTA", "TTTGGAT" } });
  REQUIRE_THROWS_AS(ss.uninitialized_adapters(), assert_failed);
}

TEST_CASE("loading sample set with uninitialized adapters", "[sample_set]")
{
  sample_set ss;

  adapter_set as_1{ { "TGGAT", "ACGTA" } };
  ss.set_adapters(as_1);
  ss.flag_uninitialized_adapters();

  vec_reader reader{ { "sample_1 TT GG" } };
  ss.load(reader, barcode_config{});

  REQUIRE_THROWS_AS(ss.at(0).at(0).adapters(), assert_failed);
  REQUIRE(ss.uninitialized_adapters() == as_1);

  adapter_set as_2{ { "ACGTA", "TGGAT" } };
  ss.set_adapters(as_2);
  REQUIRE(ss.adapters() == as_2);
  REQUIRE(ss.at(0).at(0).adapters() == adapter_set{ { "CCACGTA", "AATGGAT" } });
  REQUIRE_THROWS_AS(ss.uninitialized_adapters(), assert_failed);
}

////////////////////////////////////////////////////////////////////////////////
// Sample set -- debug string

TEST_CASE("sample set to string", "[sample_set]")
{
  std::ostringstream os;
  os << sample_set{};

  CHECK(
    os.str() ==
    "sample_set{samples=[sample{name='', "
    "barcodes=[sample_sequences{has_read_group=false, "
    "read_group=read_group{id='1', header='@RG\\tID:1'}, "
    "barcode_1=dna_sequence{''}, barcode_2=dna_sequence{''}, "
    "orientation=barcode_orientation::unspecified, adapters=adapter_set{[]}, "
    "uninitialized_adapters=false}]}], unidentified=sample{name='', "
    "barcodes=[sample_sequences{has_read_group=true, "
    "read_group=read_group{id='1', header='@RG\\tID:1\\tDS:unidentified'}, "
    "barcode_1=dna_sequence{''}, barcode_2=dna_sequence{''}, "
    "orientation=barcode_orientation::unspecified, adapters=adapter_set{[]}, "
    "uninitialized_adapters=false}]}, read_group=read_group{id='1', "
    "header='@RG\\tID:1'}, adapters=adapter_set{[]}}");
}

} // namespace adapterremoval
