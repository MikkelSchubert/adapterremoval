// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "commontypes.hpp"
#include "read_group.hpp" // for read_group
#include "testing.hpp"    // for TEST_CASE, REQUIRE, ...
#include <stdexcept>      // for invalid_argument

using Contains = Catch::Matchers::StdString::ContainsMatcher;

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Read groups

TEST_CASE("default read group", "[read_group]")
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

TEST_CASE("minimal read group with PG", "[read_group]")
{
  std::string_view header = GENERATE("PG:foo", "@RG\tPG:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "1");
  REQUIRE(rg.header() == "@RG\tID:1\tPG:foo");
}

TEST_CASE("minimal read group with ID", "[read_group]")
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

TEST_CASE("minimal read group with SM", "[read_group]")
{
  std::string_view header = GENERATE("SM:foo", "@RG\tSM:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "1");
  REQUIRE(rg.header() == "@RG\tID:1\tSM:foo");

  std::string name = GENERATE("2", "longer");

  rg.set_sample(name);
  REQUIRE(rg.header() == std::string("@RG\tID:1\tSM:") + name + "");
}

TEST_CASE("minimal read group with barcodes", "[read_group]")
{
  {
    read_group rg;
    rg.set_barcodes("ACGT");
    CHECK(rg.header() == "@RG\tID:1\tBC:ACGT");
  }

  {
    read_group rg;
    rg.set_barcodes("ACGT-TTGA");
    CHECK(rg.header() == "@RG\tID:1\tBC:ACGT-TTGA");
  }
}

TEST_CASE("unsetting read group fields", "[read_group]")
{
  read_group rg;
  REQUIRE(rg.header() == "@RG\tID:1");
  rg.set_sample("foo");
  REQUIRE(rg.header() == "@RG\tID:1\tSM:foo");
  rg.set_description("comment");
  REQUIRE(rg.header() == "@RG\tID:1\tSM:foo\tDS:comment");

  SECTION("SM first")
  {
    rg.set_sample("");
    REQUIRE(rg.header() == "@RG\tID:1\tDS:comment");
    rg.set_description("");
  }

  SECTION("DS first")
  {
    rg.set_description("");
    REQUIRE(rg.header() == "@RG\tID:1\tSM:foo");
    rg.set_sample("");
  }

  REQUIRE(rg.header() == "@RG\tID:1");
}

TEST_CASE("read group equality operator", "[read_group]")
{
  CHECK(read_group{} == read_group{});
  CHECK(read_group{} == read_group{ "@RG\tID:1" });
  CHECK(read_group{} == read_group{ "ID:1" });
  CHECK_FALSE(read_group{} == read_group{ "@RG\tID:sample" });
  CHECK_FALSE(read_group{} == read_group{ "ID:sample" });
  CHECK_FALSE(read_group{} == read_group{ "@RG\tID:1\tSM:sample" });
  CHECK_FALSE(read_group{} == read_group{ "ID:1\tSM:sample" });
}

TEST_CASE("invalid read groups", "[read_group]")
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

} // namespace adapterremoval
