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
#include "errors.hpp"   // for assert_failed
#include "strutils.hpp" // for format_rough_number, wrap_text, str_to_unsigned
#include "testing.hpp"  // for catch.hpp, StringMaker
#include <cstdint>      // for INT64_MAX, INTPTR_MAX
#include <stdexcept>    // for invalid_argument
#include <string>       // for basic_string, operator==, string
#include <vector>       // for operator==, vector

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
// Tests for 'to_lower'

TEST_CASE("ASCII letters are lowercased", "[strutils::to_lower]")
{
  REQUIRE(to_lower("") == "");
  REQUIRE(to_lower("a1{2BZ`zAdeK") == "a1{2bz`zadek");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'to_upper'

TEST_CASE("ASCII letters are uppercased", "[strutils::to_upper]")
{
  REQUIRE(to_upper("") == "");
  REQUIRE(to_upper("a1{2BZ`zAdeK") == "A1{2BZ`ZADEK");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'ends_with'

TEST_CASE("ends_with")
{
  REQUIRE(ends_with("", ""));
  REQUIRE(ends_with("morb", ""));
  REQUIRE(ends_with("morb", "b"));
  REQUIRE(ends_with("morb", "rb"));
  REQUIRE(ends_with("morb", "orb"));
  REQUIRE(ends_with("morb", "morb"));

  REQUIRE_FALSE(ends_with("", "x"));
  REQUIRE_FALSE(ends_with("morb", "x"));
  REQUIRE_FALSE(ends_with("morb", "xb"));
  REQUIRE_FALSE(ends_with("morb", "xrb"));
  REQUIRE_FALSE(ends_with("morb", "xorb"));
  REQUIRE_FALSE(ends_with("morb", "xmorb"));
  REQUIRE_FALSE(ends_with("morb", "MORB"));
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
// Tests for 'wrap_text'

using vec = std::vector<std::string>;

TEST_CASE("Whitespace only is stripped", "[strutils::wrap_text]")
{
  REQUIRE(wrap_text("") == vec{});
  REQUIRE(wrap_text(" ") == vec{});
  REQUIRE(wrap_text("\n") == vec{});
  REQUIRE(wrap_text("\n\n") == vec{});
  REQUIRE(wrap_text("\n \n") == vec{});
  REQUIRE(wrap_text("\n ") == vec{});
  REQUIRE(wrap_text(" \n") == vec{});
  REQUIRE(wrap_text(" \n ") == vec{});
}

TEST_CASE("Whitespace between words is stripped", "[strutils::wrap_text]")
{
  REQUIRE(wrap_text("foo bar") == vec{ "foo bar" });
  REQUIRE(wrap_text("foo  bar") == vec{ "foo bar" });
  REQUIRE(wrap_text("foo\nbar") == vec{ "foo bar" });
  REQUIRE(wrap_text("foo\nbar\n") == vec{ "foo bar" });
  REQUIRE(wrap_text("\nfoo\nbar") == vec{ "foo bar" });
  REQUIRE(wrap_text("foo\nbar\n\n") == vec{ "foo bar" });
  REQUIRE(wrap_text("\nfoo\nbar\n") == vec{ "foo bar" });
  REQUIRE(wrap_text("\n \n foo \n bar \n") == vec{ "foo bar" });
}

TEST_CASE("Line-breaks are added by max width", "[strutils::wrap_text]")
{
  REQUIRE(wrap_text("foo bar\nzood", 12) == vec{ "foo bar zood" });
  REQUIRE(wrap_text("foo bar\nzood", 11) == vec{ "foo bar", "zood" });
  REQUIRE(wrap_text("foo bar\nzood", 7) == vec{ "foo bar", "zood" });
  REQUIRE(wrap_text("foo bar\nzood", 6) == vec{ "foo", "bar", "zood" });
  REQUIRE(wrap_text("foo bar\nzood", 3) == vec{ "foo", "bar", "zood" });
  REQUIRE(wrap_text("foo bar\nzood", 1) == vec{ "foo", "bar", "zood" });
  REQUIRE(wrap_text("foo bar\nzood", 0) == vec{ "foo", "bar", "zood" });
}

TEST_CASE("Subsequent lines are indented", "[strutils::wrap_text]")
{
  REQUIRE(wrap_text("foo bar\nzood", 12, 2) == vec{ "foo bar zood" });
  REQUIRE(wrap_text("foo bar\nzood", 11, 0) == vec{ "foo bar", "zood" });
  REQUIRE(wrap_text("foo bar\nzood", 11, 1) == vec{ "foo bar", " zood" });
  REQUIRE(wrap_text("foo bar\nzood", 11, 2) == vec{ "foo bar", "  zood" });
  REQUIRE(wrap_text("foo bar\nzood", 7, 2) == vec{ "foo bar", "  zood" });
  REQUIRE(wrap_text("foo bar\nzood", 6, 2) == vec{ "foo", "  bar", "  zood" });
  REQUIRE(wrap_text("foo bar\nzood", 3, 2) == vec{ "foo", "  bar", "  zood" });
  REQUIRE(wrap_text("foo bar\nzood", 1, 2) == vec{ "foo", "  bar", "  zood" });
  REQUIRE(wrap_text("foo bar\nzood", 0, 2) == vec{ "foo", "  bar", "  zood" });
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
// Tests for 'shell_escape' / 'log_escape'

TEST_CASE("shell_escape tests")
{
  // empty values must be quoted
  REQUIRE(shell_escape("") == "''");
  // passes simple values unchanged
  REQUIRE(shell_escape("my_file.txt") == "my_file.txt");
  // whitespace in a value must be quoted
  REQUIRE(shell_escape("my file.txt") == "'my file.txt'");
  // various characters must be escaped (using standard escapes)
  REQUIRE(shell_escape("bob's_file.txt") == "'bob\\'s_file.txt'");
  REQUIRE(shell_escape("my\\file.txt") == "'my\\\\file.txt'");
  REQUIRE(shell_escape("my\nfile.txt") == "'my\\nfile.txt'");
  REQUIRE(shell_escape("my\rfile.txt") == "'my\\rfile.txt'");
  REQUIRE(shell_escape("my\tfile.txt") == "'my\\tfile.txt'");
  REQUIRE(shell_escape("my\bfile.txt") == "'my\\bfile.txt'");
  REQUIRE(shell_escape("my\ffile.txt") == "'my\\ffile.txt'");
  // mixed strings must be quoted and escaped
  REQUIRE(shell_escape("bob's file.txt") == "'bob\\'s file.txt'");
  // unprintable characters must be hex encoded
  REQUIRE(shell_escape("\1") == "'\\x1'");
  REQUIRE(shell_escape("\x7f") == "'\\x7f'");
}

TEST_CASE("log_escape tests")
{
  // empty values must be quoted
  REQUIRE(log_escape("") == "''");
  // passes simple values unchanged
  REQUIRE(log_escape("my_file.txt") == "'my_file.txt'");
  // whitespace in a value must be quoted
  REQUIRE(log_escape("my file.txt") == "'my file.txt'");
  // various characters must be escaped (using standard escapes)
  REQUIRE(log_escape("bob's_file.txt") == "'bob\\'s_file.txt'");
  REQUIRE(log_escape("my\\file.txt") == "'my\\\\file.txt'");
  REQUIRE(log_escape("my\bfile.txt") == "'my\\bfile.txt'");
  REQUIRE(log_escape("my\ffile.txt") == "'my\\ffile.txt'");
  REQUIRE(log_escape("my\nfile.txt") == "'my\\nfile.txt'");
  REQUIRE(log_escape("my\rfile.txt") == "'my\\rfile.txt'");
  REQUIRE(log_escape("my\tfile.txt") == "'my\\tfile.txt'");
  // mixed strings must be quoted and escaped
  REQUIRE(log_escape("bob's file.txt") == "'bob\\'s file.txt'");
  // unprintable characters must be hex encoded
  REQUIRE(log_escape("\1") == "'\\x1'");
  REQUIRE(log_escape("\x7f") == "'\\x7f'");
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

  // 64 bit platforms only
#if INTPTR_MAX >= INT64_MAX
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
#endif
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
