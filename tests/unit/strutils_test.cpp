/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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

#include "debug.hpp"
#include "strutils.hpp"
#include "testing.hpp"

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Tests for 'levenshtein'

TEST_CASE("Levenshtein distance", "[strutils::levenshtein]")
{
  REQUIRE(levenshtein("", "") == 0);
  REQUIRE(levenshtein("a", "") == 1);
  REQUIRE(levenshtein("", "a") == 1);
  REQUIRE(levenshtein("a", "a") == 0);
  REQUIRE(levenshtein("ab", "") == 2);
  REQUIRE(levenshtein("", "ab") == 2);
  REQUIRE(levenshtein("ab", "c") == 2);
  REQUIRE(levenshtein("c", "ab") == 2);
  REQUIRE(levenshtein("b", "abc") == 2);
  REQUIRE(levenshtein("abc", "b") == 2);
  REQUIRE(levenshtein("abc", "") == 3);
  REQUIRE(levenshtein("", "abc") == 3);
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'tolower'

TEST_CASE("ASCII letters are lowercased", "[strutils::tolower]")
{
  REQUIRE(tolower("") == "");
  REQUIRE(tolower("a1{2BZ`zAdeK") == "a1{2bz`zadek");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'toupper'

TEST_CASE("ASCII letters are uppercased", "[strutils::toupper]")
{
  REQUIRE(toupper("") == "");
  REQUIRE(toupper("a1{2BZ`zAdeK") == "A1{2BZ`ZADEK");
}

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
  REQUIRE(indent_lines("this is a test\nline #2") ==
          "    this is a test\n    line #2");
  REQUIRE(indent_lines("this is a test\nline #2\n") ==
          "    this is a test\n    line #2\n");
}

TEST_CASE("Mixed lines are partially indented", "[strutils::indent_lines]")
{
  REQUIRE(indent_lines("this is a test\nline #2\n\n") ==
          "    this is a test\n    line #2\n\n");
  REQUIRE(indent_lines("this is a test\n\nline #2\n") ==
          "    this is a test\n\n    line #2\n");
  REQUIRE(indent_lines("this is a test\nline #2\n\n\n") ==
          "    this is a test\n    line #2\n\n\n");
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

///////////////////////////////////////////////////////////////////////////////
// Tests for 'template_replace'

TEST_CASE("Empty target", "[strutils::template_replace]")
{
  REQUIRE(template_replace("", "tmpl", "abc") == "");
  REQUIRE(template_replace("", "tmpl", "") == "");
}

TEST_CASE("Simple replace", "[strutils::template_replace]")
{
  REQUIRE(template_replace("{tmpl}", "tmpl", "abc") == "abc");
  REQUIRE(template_replace("{.tmpl}", "tmpl", "abc") == ".abc");
  REQUIRE(template_replace("{tmpl.}", "tmpl", "abc") == "abc.");
  REQUIRE(template_replace("{.tmpl.}", "tmpl", "abc") == ".abc.");
  REQUIRE(template_replace("{tmpl}", "tmpl", "") == "");
}

TEST_CASE("Replace once in string", "[strutils::template_replace]")
{
  REQUIRE(template_replace("foo.{tmpl}.gz", "tmpl", "abc") == "foo.abc.gz");
  REQUIRE(template_replace("foo{.tmpl}.gz", "tmpl", "abc") == "foo.abc.gz");
  REQUIRE(template_replace("foo.{tmpl.}gz", "tmpl", "abc") == "foo.abc.gz");
  REQUIRE(template_replace("foo{.tmpl.}gz", "tmpl", "abc") == "foo.abc.gz");
  REQUIRE(template_replace("foo{.tmpl}.gz", "tmpl", "") == "foo.gz");
}

TEST_CASE("Replace twice in string", "[strutils::template_replace]")
{
  REQUIRE(template_replace("{tmpl}/foo.{tmpl}.gz", "tmpl", "abc") ==
          "abc/foo.abc.gz");
  REQUIRE(template_replace("{tmpl}/foo{.tmpl}.gz", "tmpl", "abc") ==
          "abc/foo.abc.gz");
  REQUIRE(template_replace("{tmpl}/foo.{tmpl.}gz", "tmpl", "abc") ==
          "abc/foo.abc.gz");
  REQUIRE(template_replace("{tmpl}/foo{.tmpl.}gz", "tmpl", "abc") ==
          "abc/foo.abc.gz");
  REQUIRE(template_replace("{tmpl}/foo{.tmpl}.gz", "tmpl", "") == "/foo.gz");
}

TEST_CASE("Partials ignored", "[strutils::template_replace]")
{
  REQUIRE(template_replace("tmpl}", "tmpl", "abc") == "tmpl}");
  REQUIRE(template_replace(".tmpl}", "tmpl", "abc") == ".tmpl}");
  REQUIRE(template_replace("{tmpl", "tmpl", "abc") == "{tmpl");
  REQUIRE(template_replace("{tmpl.", "tmpl", "abc") == "{tmpl.");
  REQUIRE(template_replace("tmpl", "tmpl", "") == "tmpl");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'shell_escape'

TEST_CASE("shell_escape preserves empty values")
{
  REQUIRE(shell_escape("") == "''");
}

TEST_CASE("shell_escape passes simple values unchanged")
{
  REQUIRE(shell_escape("my_file.txt") == "my_file.txt");
}

TEST_CASE("shell_escape escapes whitespace")
{
  REQUIRE(shell_escape("my file.txt") == "'my file.txt'");
}

TEST_CASE("shell_escape escapes single quotes")
{
  REQUIRE(shell_escape("bob's_file.txt") == "'bob\\'s_file.txt'");
}

TEST_CASE("shell_escape escapes multiple values")
{
  REQUIRE(shell_escape("bob's file.txt") == "'bob\\'s file.txt'");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'shell_escape_command'

TEST_CASE("shell_escape_command with empty list")
{
  REQUIRE(shell_escape_command({}) == "");
}

TEST_CASE("shell_escape_command with mixed command")
{
  REQUIRE(shell_escape_command({ "echo", "a", "m$x", "of", "$(things)" }) ==
          "echo a 'm$x' of '$(things)'");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'format_thousand_sep'

TEST_CASE("format_thousand_sep")
{
  REQUIRE(format_thousand_sep(0) == "0");
  REQUIRE(format_thousand_sep(999) == "999");
  REQUIRE(format_thousand_sep(1000) == "1,000");
  REQUIRE(format_thousand_sep(999999) == "999,999");
  REQUIRE(format_thousand_sep(1000000) == "1,000,000");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'format_rough_number'

TEST_CASE("format_rough_number")
{
  REQUIRE(format_rough_number(0) == "0");
  REQUIRE(format_rough_number(999) == "999");
  REQUIRE(format_rough_number(1230) == "1.23 K");
  REQUIRE(format_rough_number(12300) == "12.3 K");
  REQUIRE(format_rough_number(999100) == "999 K");
  REQUIRE(format_rough_number(999999) == "1.00 M");
  REQUIRE(format_rough_number(1234000) == "1.23 M");
  REQUIRE(format_rough_number(12340000) == "12.3 M");
  REQUIRE(format_rough_number(999100000) == "999 M");
  REQUIRE(format_rough_number(999999999) == "1.00 G");
  REQUIRE(format_rough_number(1234000000) == "1.23 G");
  REQUIRE(format_rough_number(12340000000) == "12.3 G");
  REQUIRE(format_rough_number(999100000000) == "999 G");
  REQUIRE(format_rough_number(999999999999) == "1.00 T");
  REQUIRE(format_rough_number(1234000000000) == "1.23 T");
  REQUIRE(format_rough_number(12340000000000) == "12.3 T");
  REQUIRE(format_rough_number(999100000000000) == "999 T");
  REQUIRE(format_rough_number(999999999999999) == "1.00 P");
  REQUIRE(format_rough_number(1234000000000000) == "1.23 P");
  REQUIRE(format_rough_number(12340000000000000) == "12.3 P");
  REQUIRE(format_rough_number(1234000000000000000) == "1230 P");
}

TEST_CASE("format_rough_number with precision")
{
  REQUIRE(format_rough_number(0, 1) == "0");
  REQUIRE(format_rough_number(0, 2) == "0");

  REQUIRE(format_rough_number(129, 1) == "100");
  REQUIRE(format_rough_number(129, 2) == "130");
  REQUIRE(format_rough_number(129, 3) == "129");
  REQUIRE(format_rough_number(129, 4) == "129");

  REQUIRE(format_rough_number(1239, 1) == "1 K");
  REQUIRE(format_rough_number(1239, 2) == "1.2 K");
  REQUIRE(format_rough_number(1239, 3) == "1.24 K");
  REQUIRE(format_rough_number(1239, 4) == "1239");

  REQUIRE(format_rough_number(12349, 1) == "10 K");
  REQUIRE(format_rough_number(12349, 2) == "12 K");
  REQUIRE(format_rough_number(12349, 3) == "12.3 K");
  REQUIRE(format_rough_number(12349, 4) == "12.35 K");
  REQUIRE(format_rough_number(12349, 5) == "12349");
}

TEST_CASE("format_rough_number with no precision")
{
  REQUIRE_THROWS_AS(format_rough_number(12349, 0), assert_failed);
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'format_fraction'

TEST_CASE("format_fraction")
{
  REQUIRE(format_fraction(0, 0) == "NA");
  REQUIRE(format_fraction(1, 0) == "NA");
  REQUIRE(format_fraction(55, 300) == "0.18");
  REQUIRE(format_fraction(55, 300, 0) == "0");
  REQUIRE(format_fraction(55, 300, 1) == "0.2");
  REQUIRE(format_fraction(55, 300, 2) == "0.18");
  REQUIRE(format_fraction(55, 300, 3) == "0.183");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'format_fraction'

TEST_CASE("format_percentage")
{
  REQUIRE(format_percentage(0, 0) == "NA");
  REQUIRE(format_percentage(1, 0) == "NA");
  REQUIRE(format_percentage(55, 300) == "18.3");
  REQUIRE(format_percentage(55, 300, 0) == "18");
  REQUIRE(format_percentage(55, 300, 1) == "18.3");
  REQUIRE(format_percentage(55, 300, 2) == "18.33");
}

} // namespace adapterremoval
