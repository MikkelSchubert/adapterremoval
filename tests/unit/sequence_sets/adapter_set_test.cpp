// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"        // for parsing_error
#include "linereader.hpp"    // for vec_reader
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for read_group
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <string>            // for string==, string
#include <string_view>       // for basic_string_view
#include <utility>           // for pair, operator==
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
  adapter_set set{ { dna_sequence{ "ACGTA" }, dna_sequence{ "TGGAT" } },
                   { dna_sequence{ "CCGAT" }, dna_sequence{ "AGAGT" } } };

  // alignment orientation
  const auto pair_1 =
    sequence_pair{ dna_sequence{ "ACGTA" }, dna_sequence{ "ATCCA" } };
  const auto pair_2 =
    sequence_pair{ dna_sequence{ "CCGAT" }, dna_sequence{ "ACTCT" } };

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
  const auto set_1a = dna_sequence{ "ACT" };
  const auto set_1b = dna_sequence{ "CTA" };
  const auto set_2a = dna_sequence{ "ATT" };
  const auto set_2b = dna_sequence{ "CTT" };

  adapter_set set;
  CHECK(std::vector(set.begin(), set.end()) == sequence_pair_vec{});

  set.add(dna_sequence{ "ACT" }, dna_sequence{ "CTA" });
  CHECK(std::vector(set.begin(), set.end()) ==
        sequence_pair_vec{ { dna_sequence{ "ACT" }, dna_sequence{ "TAG" } } });

  set.add(dna_sequence{ "ATT" }, dna_sequence{ "CTT" });
  CHECK(std::vector(set.begin(), set.end()) ==
        sequence_pair_vec{ { dna_sequence{ "ACT" }, dna_sequence{ "TAG" } },
                           { dna_sequence{ "ATT" }, dna_sequence{ "AAG" } } });
}

TEST_CASE("add barcodes to adapters", "[adapter_set]")
{
  // adapter_set add_barcodes(const dna_sequence& barcode1, const dna_sequence&
  // barcode2) const;
  const sequence_pair set_1{ dna_sequence{ "ACT" }, dna_sequence{ "CTA" } };
  const sequence_pair set_2{ dna_sequence{ "ATT" }, dna_sequence{ "CTT" } };
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
  const auto set_1 =
    sequence_pair{ dna_sequence{ "ACT" }, dna_sequence{ "CTA" } };
  const auto set_2 =
    sequence_pair{ dna_sequence{ "ATT" }, dna_sequence{ "CTT" } };

  CHECK(adapter_set{}.to_read_orientation() == sequence_pair_vec{});
  CHECK(adapter_set{ set_1 }.to_read_orientation() ==
        sequence_pair_vec{ set_1 });
  CHECK(adapter_set{ set_1, set_2 }.to_read_orientation() ==
        sequence_pair_vec{ set_1, set_2 });
}

TEST_CASE("equality operator", "[adapter_set]")
{
  auto empty = adapter_set{};
  auto set_1a = adapter_set{ { dna_sequence{ "ACT" }, dna_sequence{ "CTA" } } };
  auto set_1b = set_1a;
  auto set_2 = adapter_set{ { dna_sequence{ "ATT" }, dna_sequence{ "CTT" } } };

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
