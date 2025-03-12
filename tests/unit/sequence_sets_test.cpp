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
#include "errors.hpp"        // for parsing_error
#include "sequence_sets.hpp" // for read_group
#include "testing.hpp"       // for catch.hpp, StringMaker

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
  std::string_view header = GENERATE("PG:foo", "@RG\tPG:foo", "@RG\\tPG:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "1");
  REQUIRE(rg.header() == "@RG\tID:1\tPG:foo");
}

TEST_CASE("black-slash in read-group")
{
  read_group rg{ "ID:foo\\\\bar" };
  REQUIRE(rg.id() == "foo\\bar");
  REQUIRE(rg.header() == "@RG\tID:foo\\bar");
}

TEST_CASE("minimal read group with ID")
{
  std::string_view header = GENERATE("ID:foo", "@RG\tID:foo", "@RG\\tID:foo");

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
  std::string_view header = GENERATE("SM:foo", "@RG\tSM:foo", "@RG\\tSM:foo");

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
  REQUIRE_THROWS_WITH(read_group{ "ID:" },
                      "tags must be at least 4 characters long");

  REQUIRE_THROWS_WITH(read_group{ "1D:1" },
                      "first character in tag name must be a letter");
  REQUIRE_THROWS_WITH(read_group{ "!D:1" },
                      "first character in tag name must be a letter");

  REQUIRE_THROWS_WITH(
    read_group{ "D!:1" },
    "second character in tag name must be a letter or number");

  REQUIRE_THROWS_WITH(read_group{ "ID:1\tID:2" }, "multiple ID tags found");

  REQUIRE_THROWS_WITH(read_group{ "ID?1" },
                      "third character in tag must be a colon");
  REQUIRE_THROWS_WITH(read_group{ "ID:1\nSM:foo" },
                      "only characters in the range ' ' to '~' are allowed");

  REQUIRE_THROWS_WITH(read_group{ "SM:foo\\x" },
                      "invalid escape sequence '\\x'");
  REQUIRE_THROWS_WITH(read_group{ "SM:foo\\" },
                      "incomplete escape sequence at end of string");
}

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- initializer list

TEST_CASE("Overlapping SE barcodes fail", "[sample_set::initializer_list]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{} };
  sample sample2{ "sample2", dna_sequence{ "ACGT" }, dna_sequence{} };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

TEST_CASE("Overlapping PE barcodes fail", "[sample_set::initializer_list]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

TEST_CASE("Variable length SE barcodes fail #1",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGTA" }, dna_sequence{} };
  sample sample2{ "sample2", dna_sequence{ "TGCT" }, dna_sequence{} };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

TEST_CASE("Variable length SE barcodes fail #2",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{} };
  sample sample2{ "sample2", dna_sequence{ "TGCTA" }, dna_sequence{} };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #1",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGTT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #2",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTGT" } };
  sample sample2{ "sample2", dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #3",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "CGTGA" }, dna_sequence{ "ACGT" } };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #4",
          "[sample_set::initializer_liset]")
{
  sample sample1{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } };
  sample sample2{ "sample2", dna_sequence{ "CGTGA" }, dna_sequence{ "ACGTA" } };
  REQUIRE_THROWS_AS(sample_set({ sample1, sample2 }), parsing_error);
}

} // namespace adapterremoval
