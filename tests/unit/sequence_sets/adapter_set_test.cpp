// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "commontypes.hpp"
#include "errors.hpp"        // for parsing_error
#include "linereader.hpp"    // for vec_reader
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for read_group
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <string>            // for string==, string
#include <vector>            // for vector

using Contains = Catch::Matchers::StdString::ContainsMatcher;

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Adapter set

TEST_CASE("default constructor", "[adapter_set]")
{
  adapter_set set;

  REQUIRE(set.size() == 0);
  REQUIRE(set.empty());
  REQUIRE(sequence_pair_vec{ set.begin(), set.end() } == sequence_pair_vec{});
}

TEST_CASE("initializer list", "[adapter_set]")
{
  // read orientation
  adapter_set set{ { "ACGTA", "TGGAT" }, { "CCGAT", "AGAGT" } };

  // alignment orientation
  const sequence_pair pair_1{ "ACGTA"_dna, "ATCCA"_dna };
  const sequence_pair pair_2{ "CCGAT"_dna, "ACTCT"_dna };

  REQUIRE(set.size() == 2);
  REQUIRE(!set.empty());
  REQUIRE(set.at(0) == pair_1);
  REQUIRE(set.at(1) == pair_2);
  REQUIRE(sequence_pair_vec{ set.begin(), set.end() } ==
          sequence_pair_vec{ pair_1, pair_2 });
}

///////////////////////////////////////////////////////////////////////////////
// Adapter set -- load

namespace {

adapter_set
load_adapters(std::vector<std::string> lines, bool paired_end)
{
  adapter_set adapters;
  vec_reader reader{ lines };
  adapters.load(reader, paired_end);
  return adapters;
}

} // namespace

TEST_CASE("load table", "[adapter_set]")
{
  SECTION("one adapter")
  {
    auto adapters = load_adapters({ "ACACA", "GTTAGA" }, false);
    CHECK(adapters == adapter_set{ { "ACACA", {} }, { "GTTAGA", {} } });
  }

  SECTION("two adapters")
  {
    const auto paired_end = GENERATE(true, false);
    auto adapters =
      load_adapters({ "ACACA GTCGT", "GTTAGA CCTGAG" }, paired_end);
    CHECK(adapters ==
          adapter_set{ { "ACACA", "GTCGT" }, { "GTTAGA", "CCTGAG" } });
  }

  SECTION("whitespace and comments")
  {
    const auto paired_end = GENERATE(true, false);
    auto adapters = load_adapters({ "",
                                    "# comment",
                                    " ",
                                    "  ACACA  GTCGT  ",
                                    "GTTAGA CCTGAG",
                                    "",
                                    " # second comment" },
                                  paired_end);
    CHECK(adapters ==
          adapter_set{ { "ACACA", "GTCGT" }, { "GTTAGA", "CCTGAG" } });
  }
}

TEST_CASE("load empty table", "[adapter_set]")
{
  const auto paired_end = GENERATE(true, false);

  CHECK_THROWS_MESSAGE(load_adapters({}, paired_end),
                       parsing_error,
                       "No adapter sequences in table");
}

TEST_CASE("load table with wrong number of columns", "[adapter_set]")
{
  // it only makes sense to test PE, since SE with too few values is an empty
  CHECK_THROWS_MESSAGE(load_adapters({ "ACGTTA" }, true),
                       parsing_error,
                       "Error at line 1: Expected at least 2 columns, but "
                       "found 1 column(s)");

  CHECK_THROWS_MESSAGE(load_adapters({ "GTTAGA CCTGAG", "ACGTTA" }, false),
                       parsing_error,
                       "Error at line 2: Inconsistent number of columns; "
                       "expected 2 column(s) but found 1");
}

TEST_CASE("load table with too many columns", "[adapter_set]")
{
  const auto paired_end = GENERATE(true, false);

  CHECK_THROWS_MESSAGE(load_adapters({ "ACGTTA ACCGTA bad" }, paired_end),
                       parsing_error,
                       "Error at line 1: Expected at most 2 columns, but found "
                       "3 column(s)");
}

///////////////////////////////////////////////////////////////////////////////
// Adapter set -- other functions

TEST_CASE("add squences", "[adapter_set]")
{
  const dna_sequence set_1a{ "ACT" };
  const dna_sequence set_1b{ "CTA" };
  const dna_sequence set_2a{ "ATT" };
  const dna_sequence set_2b{ "CTT" };

  adapter_set set;
  CHECK(std::vector(set.begin(), set.end()) == sequence_pair_vec{});

  set.add("ACT"_dna, "CTA"_dna);
  CHECK(std::vector(set.begin(), set.end()) ==
        sequence_pair_vec{ { "ACT"_dna, "TAG"_dna } });

  set.add("ATT"_dna, "CTT"_dna);
  CHECK(
    std::vector(set.begin(), set.end()) ==
    sequence_pair_vec{ { "ACT"_dna, "TAG"_dna }, { "ATT"_dna, "AAG"_dna } });
}

TEST_CASE("add barcodes to adapters", "[adapter_set]")
{
  // adapter_set add_barcodes(const dna_sequence& barcode1, const dna_sequence&
  // barcode2) const;
  const string_pair set_1{ "ACT", "CTA" };
  const string_pair set_2{ "ATT", "CTT" };
  const dna_sequence barcode_1{ "AAAA" };
  const dna_sequence barcode_2{ "GGGG" };

  SECTION("empty set")
  {
    adapter_set set;
    CHECK(set.add_barcodes(barcode_1, barcode_2) == adapter_set{});
    CHECK(set == adapter_set{});
  }

  SECTION("empty adapters")
  {
    adapter_set set{ { "", "" } };
    CHECK(set.add_barcodes(barcode_1, barcode_2) ==
          adapter_set{ { "CCCC", "TTTT" } });
    CHECK(set == adapter_set{ { "", "" } });
  }

  SECTION("single adapter pair")
  {
    adapter_set set{ set_1 };
    CHECK(set.add_barcodes(barcode_1, barcode_2) ==
          adapter_set{ { "CCCCACT", "TTTTCTA" } });
    CHECK(set == adapter_set{ set_1 });
  }

  SECTION("multiple adapter pairs")
  {
    adapter_set set{ set_1, set_2 };
    CHECK(set.add_barcodes(barcode_1, barcode_2) ==
          adapter_set{ { "CCCCACT", "TTTTCTA" }, { "CCCCATT", "TTTTCTT" } });
    CHECK(set == adapter_set{ set_1, set_2 });
  }
}

TEST_CASE("to read orientation", "[adapter_set]")
{
  const auto set_1 = string_pair{ "ACT", "CTA" };
  const auto set_2 = string_pair{ "ATT", "CTT" };

  CHECK(adapter_set{}.to_read_orientation() == sequence_pair_vec{});
  CHECK(adapter_set{ set_1 }.to_read_orientation() ==
        sequence_pair_vec{ sequence_pair{ "ACT", "CTA" } });
  CHECK(adapter_set{ set_1, set_2 }.to_read_orientation() ==
        sequence_pair_vec{ sequence_pair{ "ACT", "CTA" },
                           sequence_pair{ "ATT", "CTT" } });
}

TEST_CASE("equality operator", "[adapter_set]")
{
  auto empty = adapter_set{};
  auto set_1a = adapter_set{ { "ACT", "CTA" } };
  auto set_1b = set_1a;
  auto set_2 = adapter_set{ { "ATT", "CTT" } };

  CHECK(empty == empty);
  CHECK(empty == adapter_set{});
  CHECK(set_1a == set_1a);
  CHECK(set_1a == set_1b);
  CHECK_FALSE(set_1a == empty);
  CHECK_FALSE(set_1a == set_2);
  CHECK_FALSE(set_2 == empty);
}

///////////////////////////////////////////////////////////////////////////////
// Adapter set -- to debug string

TEST_CASE("adapter set to debug string", "[adapter_set]")
{
  using ::Catch::fallbackStringifier;

  CHECK(fallbackStringifier(adapter_set{}) == "adapter_set{[]}");

  CHECK(fallbackStringifier(adapter_set{ { "ACGTA", "GTTCA" } }) ==
        "adapter_set{[pair{first=dna_sequence{'ACGTA'}, "
        "second=dna_sequence{'TGAAC'}}]}");

  CHECK(fallbackStringifier(
          adapter_set{ { "ACGTA", "GTTCA" }, { "TTTGG", "CCGAT" } }) ==
        "adapter_set{[pair{first=dna_sequence{'ACGTA'}, "
        "second=dna_sequence{'TGAAC'}}, pair{first=dna_sequence{'TTTGG'}, "
        "second=dna_sequence{'ATCGG'}}]}");
}

} // namespace adapterremoval
