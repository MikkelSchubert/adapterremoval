// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"        // for parsing_error
#include "sequence_sets.hpp" // for read_group
#include "strutils.hpp"      // for log_escape
#include "testing.hpp"       // for catch.hpp, StringMaker
#include <stdexcept>         // for invalid_argument

using Contains = Catch::Matchers::StdString::ContainsMatcher;

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Read groups

TEST_CASE("default read group")
{
  SECTION("implicit")
  {
    read_group rg;
    REQUIRE(rg.id() == "1");
    REQUIRE(rg.header() == "@RG\tID:1");
  }

  SECTION("explicit")
  {
    read_group rg{ "" };
    REQUIRE(rg.id() == "1");
    REQUIRE(rg.header() == "@RG\tID:1");
  }
}

TEST_CASE("minimal read group with PG")
{
  std::string_view header = GENERATE("PG:foo", "@RG\tPG:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "1");
  REQUIRE(rg.header() == "@RG\tID:1\tPG:foo");
}

TEST_CASE("minimal read group with ID")
{
  std::string_view header = GENERATE("ID:foo", "@RG\tID:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "foo");
  REQUIRE(rg.header() == "@RG\tID:foo");

  std::string id = GENERATE("2", "longer");

  rg.set_id(id);
  REQUIRE(rg.id() == id);
  REQUIRE(rg.header() == std::string("@RG\tID:") + id + "");
}

TEST_CASE("minimal read group with SM")
{
  std::string_view header = GENERATE("SM:foo", "@RG\tSM:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "1");
  REQUIRE(rg.header() == "@RG\tID:1\tSM:foo");

  std::string name = GENERATE("2", "longer");

  rg.set_sample(name);
  REQUIRE(rg.header() == std::string("@RG\tID:1\tSM:") + name + "");
}

TEST_CASE("unsetting read group fields")
{
  read_group rg;
  REQUIRE(rg.header() == "@RG\tID:1");
  rg.set_sample("foo");
  REQUIRE(rg.header() == "@RG\tID:1\tSM:foo");
  rg.set_comment("comment");
  REQUIRE(rg.header() == "@RG\tID:1\tSM:foo\tCO:comment");

  SECTION("SM first")
  {
    rg.set_sample("");
    REQUIRE(rg.header() == "@RG\tID:1\tCO:comment");
    rg.set_comment("");
  }

  SECTION("CO first")
  {
    rg.set_comment("");
    REQUIRE(rg.header() == "@RG\tID:1\tSM:foo");
    rg.set_sample("");
  }

  REQUIRE(rg.header() == "@RG\tID:1");
}

TEST_CASE("invalid read groups")
{
  REQUIRE_THROWS_MESSAGE(read_group{ "ID:" },
                         std::invalid_argument,
                         "tags must be at least 4 characters long");

  REQUIRE_THROWS_MESSAGE(read_group{ "1D:1" },
                         std::invalid_argument,
                         "first character in tag name must be a letter");
  REQUIRE_THROWS_MESSAGE(read_group{ "!D:1" },
                         std::invalid_argument,
                         "first character in tag name must be a letter");

  REQUIRE_THROWS_MESSAGE(read_group{ "D!:1" },
                         std::invalid_argument,
                         "second character in tag name must be a letter or "
                         "number");

  REQUIRE_THROWS_MESSAGE(read_group{ "ID:1\tID:2" },
                         std::invalid_argument,
                         "multiple ID tags found");

  REQUIRE_THROWS_MESSAGE(read_group{ "ID?1" },
                         std::invalid_argument,
                         "third character in tag must be a colon");
  REQUIRE_THROWS_MESSAGE(read_group{ "ID:1\nSM:foo" },
                         std::invalid_argument,
                         "only characters in the range ' ' to '~' are allowed");
}

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- initializer list

TEST_CASE("Overlapping SE barcodes fail", "[sample_set::initializer_list]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{} };
  sample sample2{ "sample2", dna_sequence{ "ACGT" }, dna_sequence{} };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Duplicate barcode pairs found in \'initializer_list\' with barcodes "
    "\'ACGT\' and \'\'. please verify correctness of the barcode table and "
    "remove any duplicate entries!");
}

TEST_CASE("Overlapping PE barcodes fail", "[sample_set::initializer_list]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Duplicate barcode pairs found in \'initializer_list\' with barcodes "
    "\'ACGT\' and \'CGTG\'. please verify correctness of the barcode table and "
    "remove any duplicate entries!");
}

TEST_CASE("Variable length SE barcodes fail #1",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGTA" }, dna_sequence{} };
  sample sample2{ "sample2", dna_sequence{ "TGCT" }, dna_sequence{} };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 5 base-pairs "
    "long, but barcode \'TGCT\' is 4 base-pairs long. Variable length barcodes "
    "are not supported");
}

TEST_CASE("Variable length SE barcodes fail #2",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{} };
  sample sample2{ "sample2", dna_sequence{ "TGCTA" }, dna_sequence{} };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 4 base-pairs "
    "long, but barcode \'TGCTA\' is 5 base-pairs long. Variable length "
    "barcodes are not supported");
}

TEST_CASE("Variable length PE barcodes fail #1",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGTT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 5 base-pairs "
    "long, but barcode \'CGTG\' is 4 base-pairs long. Variable length barcodes "
    "are not supported");
}

TEST_CASE("Variable length PE barcodes fail #2",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTGT" } };
  sample sample2{ "sample2", dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Inconsistent mate 2 barcode lengths found: Last barcode was 5 base-pairs "
    "long, but barcode \'ACGT\' is 4 base-pairs long. Variable length barcodes "
    "are not supported");
}

TEST_CASE("Variable length PE barcodes fail #3",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "CGTGA" }, dna_sequence{ "ACGT" } };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 4 base-pairs "
    "long, but barcode \'CGTGA\' is 5 base-pairs long. Variable length "
    "barcodes are not supported");
}

TEST_CASE("Variable length PE barcodes fail #4",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "CGTGA" }, dna_sequence{ "ACGTA" } };
  REQUIRE_THROWS_MESSAGE(
    sample_set({ sample1, sample2 }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 4 base-pairs "
    "long, but barcode \'CGTGA\' is 5 base-pairs long. Variable length "
    "barcodes are not supported");
}

} // namespace adapterremoval

namespace Catch {

using namespace adapterremoval;

template<>
std::string
StringMaker<parsing_error, void>::convert(parsing_error const& value)
{
  return log_escape(value.what());
}

} // namespace Catch