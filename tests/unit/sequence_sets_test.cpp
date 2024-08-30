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
    REQUIRE(rg.header() == "@RG\tID:1\tPG:adapterremoval");
  }

  SECTION("explicit")
  {
    read_group rg{ "" };
    REQUIRE(rg.id() == "1");
    REQUIRE(rg.header() == "@RG\tID:1\tPG:adapterremoval");
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
  REQUIRE(rg.header() == "@RG\tID:foo\\bar\tPG:adapterremoval");
}

TEST_CASE("minimal read group with ID")
{
  std::string_view header = GENERATE("ID:foo", "@RG\tID:foo", "@RG\\tID:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "foo");
  REQUIRE(rg.header() == "@RG\tID:foo\tPG:adapterremoval");

  std::string id = GENERATE("2", "longer");

  rg.set_id(id);
  REQUIRE(rg.id() == id);
  REQUIRE(rg.header() == std::string("@RG\tID:") + id + "\tPG:adapterremoval");
}

TEST_CASE("minimal read group with SM")
{
  std::string_view header = GENERATE("SM:foo", "@RG\tSM:foo", "@RG\\tSM:foo");

  read_group rg{ header };
  REQUIRE(rg.id() == "1");
  REQUIRE(rg.header() == "@RG\tID:1\tSM:foo\tPG:adapterremoval");

  std::string name = GENERATE("2", "longer");

  rg.set_sample(name);
  REQUIRE(rg.header() ==
          std::string("@RG\tID:1\tSM:") + name + "\tPG:adapterremoval");
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

} // namespace adapterremoval
