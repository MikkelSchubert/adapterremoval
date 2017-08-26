/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <limits>
#include <stdexcept>

#include "testing.hpp"
#include "strutils.hpp"

namespace ar
{

///////////////////////////////////////////////////////////////////////////////
// Tests for 'toupper'

TEST_CASE("ASCII letters are uppercased", "[strutils::toupper]")
{
    REQUIRE(toupper("") == "");
    REQUIRE(toupper("a1{2BZ`zAdeK") == "A1{2BZ`ZADEK");
};


///////////////////////////////////////////////////////////////////////////////
// Tests for 'indent_lines'

TEST_CASE("Empty lines are not indented", "[strutils::indent_lines]")
{
    REQUIRE(indent_lines("") == "");
    REQUIRE(indent_lines("\n") == "\n");
    REQUIRE(indent_lines("\n\n") == "\n\n");
    REQUIRE(indent_lines("\n\n\n") == "\n\n\n");
}


TEST_CASE("Singular lines are indented", "[strutils::indent_lines]")
{
    REQUIRE(indent_lines("this is a test") == "    this is a test");
    REQUIRE(indent_lines("this is a test\n") == "    this is a test\n");
}


TEST_CASE("Multiple lines are indented", "[strutils::indent_lines]")
{
    REQUIRE(indent_lines("this is a test\nline #2")
            == "    this is a test\n    line #2");
    REQUIRE(indent_lines("this is a test\nline #2\n")
            == "    this is a test\n    line #2\n");
}


TEST_CASE("Mixed lines are partially indented", "[strutils::indent_lines]")
{
    REQUIRE(indent_lines("this is a test\nline #2\n\n")
            == "    this is a test\n    line #2\n\n");
    REQUIRE(indent_lines("this is a test\n\nline #2\n")
            == "    this is a test\n\n    line #2\n");
    REQUIRE(indent_lines("this is a test\nline #2\n\n\n")
            == "    this is a test\n    line #2\n\n\n");
}


///////////////////////////////////////////////////////////////////////////////
// Tests for 'columnize_text'

TEST_CASE("Whitespace only is stripped", "[strutils::columnize_text]")
{
    REQUIRE(columnize_text("") == "");
    REQUIRE(columnize_text(" ") == "");
    REQUIRE(columnize_text("\n") == "");
    REQUIRE(columnize_text("\n\n") == "");
    REQUIRE(columnize_text("\n \n") == "");
    REQUIRE(columnize_text("\n ") == "");
    REQUIRE(columnize_text(" \n") == "");
    REQUIRE(columnize_text(" \n ") == "");
}


TEST_CASE("Whitespace between words is stripped", "[strutils::columnize_text]")
{
    REQUIRE(columnize_text("foo bar") == "foo bar");
    REQUIRE(columnize_text("foo  bar") == "foo bar");
    REQUIRE(columnize_text("foo\nbar") == "foo bar");
    REQUIRE(columnize_text("foo\nbar\n") == "foo bar");
    REQUIRE(columnize_text("\nfoo\nbar") == "foo bar");
    REQUIRE(columnize_text("foo\nbar\n\n") == "foo bar");
    REQUIRE(columnize_text("\nfoo\nbar\n") == "foo bar");
    REQUIRE(columnize_text("\n \n foo \n bar \n") == "foo bar");
}


TEST_CASE("Linebreaks are added by max width", "[strutils::columnize_text]")
{
    REQUIRE(columnize_text("foo bar\nzood", 12) == "foo bar zood");
    REQUIRE(columnize_text("foo bar\nzood", 11) == "foo bar\nzood");
    REQUIRE(columnize_text("foo bar\nzood", 7) == "foo bar\nzood");
    REQUIRE(columnize_text("foo bar\nzood", 6) == "foo\nbar\nzood");
    REQUIRE(columnize_text("foo bar\nzood", 3) == "foo\nbar\nzood");
    REQUIRE(columnize_text("foo bar\nzood", 1) == "foo\nbar\nzood");
    REQUIRE(columnize_text("foo bar\nzood", 0) == "foo\nbar\nzood");
}


TEST_CASE("Subsequent lines are indented", "[strutils::columnize_text]")
{
    REQUIRE(columnize_text("foo bar\nzood", 12, 2) == "foo bar zood");
    REQUIRE(columnize_text("foo bar\nzood", 11, 0) == "foo bar\nzood");
    REQUIRE(columnize_text("foo bar\nzood", 11, 1) == "foo bar\n zood");
    REQUIRE(columnize_text("foo bar\nzood", 11, 2) == "foo bar\n  zood");
    REQUIRE(columnize_text("foo bar\nzood", 7, 2) == "foo bar\n  zood");
    REQUIRE(columnize_text("foo bar\nzood", 6, 2) == "foo\n  bar\n  zood");
    REQUIRE(columnize_text("foo bar\nzood", 3, 2) == "foo\n  bar\n  zood");
    REQUIRE(columnize_text("foo bar\nzood", 1, 2) == "foo\n  bar\n  zood");
    REQUIRE(columnize_text("foo bar\nzood", 0, 2) == "foo\n  bar\n  zood");
}


TEST_CASE("Proper numbers are converted", "[strutils::str_to_unsigned]")
{
    REQUIRE(str_to_unsigned("0") == 0);
    REQUIRE(str_to_unsigned("1") == 1);
    REQUIRE(str_to_unsigned("2") == 2);
    REQUIRE(str_to_unsigned("2147483647") == 2147483647);
    REQUIRE(str_to_unsigned("2147483648") == 2147483648);
    REQUIRE(str_to_unsigned("2147483649") == 2147483649);
    REQUIRE(str_to_unsigned("4294967293") == 4294967293);
    REQUIRE(str_to_unsigned("4294967294") == 4294967294);
    REQUIRE(str_to_unsigned("4294967295") == 4294967295);
}


TEST_CASE("Whitespace is ignored", "[strutils::str_to_unsigned]")
{
    REQUIRE(str_to_unsigned(" 123") == 123);
    REQUIRE(str_to_unsigned("123 ") == 123);
    REQUIRE(str_to_unsigned(" 123 ") == 123);
    REQUIRE(str_to_unsigned("\n123\n") == 123);
}


TEST_CASE("Numbers outside range are rejected", "[strutils::str_to_unsigned]")
{
    REQUIRE_THROWS_AS(str_to_unsigned("-2"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("-1"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("4294967296"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("4294967297"), std::invalid_argument);
}


TEST_CASE("Non numeric strings are rejected", "[strutils::str_to_unsigned]")
{
    REQUIRE_THROWS_AS(str_to_unsigned(" "), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("x"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("a b c"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("\n"), std::invalid_argument);
}


TEST_CASE("Mixed strings are rejected", "[strutils::str_to_unsigned]")
{
    REQUIRE_THROWS_AS(str_to_unsigned("1a"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("123x"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("123 x"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("x123"), std::invalid_argument);
    REQUIRE_THROWS_AS(str_to_unsigned("x 123"), std::invalid_argument);
}


TEST_CASE("Base 10 is assumed", "[strutils::str_to_unsigned]")
{
    REQUIRE(str_to_unsigned("010") == 10);
}


} // namespace ar
