// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"   // for assert_failed
#include "strutils.hpp" // for format_rough_number, wrap_text, str_to_u32
#include "testing.hpp"  // for TEST_CASE, REQUIRE, ...
#include <algorithm>    // for find
#include <cstdint>      // for INT64_MAX, INTPTR_MAX
#include <limits>       // for numeric_limits
#include <stdexcept>    // for invalid_argument
#include <string>       // for basic_string, operator==, string
#include <string_view>  // for string_view
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
  REQUIRE(to_lower("").empty());
  REQUIRE(to_lower("a1{2BZ`zAdeK") == "a1{2bz`zadek");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'to_upper'

TEST_CASE("ASCII letters are uppercased", "[strutils::to_upper]")
{
  REQUIRE(to_upper("").empty());
  REQUIRE(to_upper("a1{2BZ`zAdeK") == "A1{2BZ`ZADEK");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'starts_with'

TEST_CASE("starts_with")
{
  REQUIRE(starts_with("", ""));
  REQUIRE(starts_with("morb", ""));
  REQUIRE(starts_with("morb", "m"));
  REQUIRE(starts_with("morb", "mo"));
  REQUIRE(starts_with("morb", "mor")); // codespell:ignore mor
  REQUIRE(starts_with("morb", "morb"));

  REQUIRE_FALSE(starts_with("", "x"));
  REQUIRE_FALSE(starts_with("morb", "x"));
  REQUIRE_FALSE(starts_with("morb", "mx"));
  REQUIRE_FALSE(starts_with("morb", "mox"));
  REQUIRE_FALSE(starts_with("morb", "morx"));
  REQUIRE_FALSE(starts_with("morb", "morbx"));
  REQUIRE_FALSE(starts_with("morb", "MORB"));
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
// Tests for 'split_text' / 'split_lines'

TEST_CASE("split_text")
{
  REQUIRE(split_text("", ':') == string_vec{ "" });
  REQUIRE(split_text("x", ':') == string_vec{ "x" });
  REQUIRE(split_text(":", ':') == string_vec{ "", "" });
  REQUIRE(split_text("x:", ':') == string_vec{ "x", "" });
  REQUIRE(split_text(":x", ':') == string_vec{ "", "x" });
  REQUIRE(split_text("::", ':') == string_vec{ "", "", "" });
  REQUIRE(split_text("a:bc:d", ':') == string_vec{ "a", "bc", "d" });
}

TEST_CASE("split_lines")
{
  REQUIRE(split_lines("") == string_vec{ "" });
  REQUIRE(split_lines("x") == string_vec{ "x" });
  REQUIRE(split_lines("\n") == string_vec{ "", "" });
  REQUIRE(split_lines("x\n") == string_vec{ "x", "" });
  REQUIRE(split_lines("\nx") == string_vec{ "", "x" });
  REQUIRE(split_lines("\n\n") == string_vec{ "", "", "" });
  REQUIRE(split_lines("a\nbc\nd") == string_vec{ "a", "bc", "d" });
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'indent_lines'

TEST_CASE("Empty lines are not indented", "[strutils::indent_lines]")
{
  REQUIRE(indent_lines("").empty());
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
  REQUIRE(wrap_text("").empty());
  REQUIRE(wrap_text(" ").empty());
  REQUIRE(wrap_text("\n").empty());
  REQUIRE(wrap_text("\n\n").empty());
  REQUIRE(wrap_text("\n \n").empty());
  REQUIRE(wrap_text("\n ").empty());
  REQUIRE(wrap_text(" \n").empty());
  REQUIRE(wrap_text(" \n ").empty());
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

///////////////////////////////////////////////////////////////////////////////
// Tests for 'is_ascii_*'

TEST_CASE("is_ascii_letter is true for lower/uppercase letters")
{
  const unsigned min = std::numeric_limits<unsigned char>::min();
  const unsigned max = std::numeric_limits<unsigned char>::max();
  const std::string_view letters =
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

  for (auto i = min; i <= max; ++i) {
    const auto c = static_cast<char>(i);
    const bool is_letter = letters.find(c) != std::string_view::npos;

    REQUIRE(is_ascii_letter(c) == is_letter);
  }
}

TEST_CASE("is_ascii_letter_or_digit is true for letters and digits")
{
  const unsigned min = std::numeric_limits<unsigned char>::min();
  const unsigned max = std::numeric_limits<unsigned char>::max();
  const std::string_view letters =
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

  for (auto i = min; i <= max; ++i) {
    const auto c = static_cast<char>(i);
    const bool is_letter_or_digit = letters.find(c) != std::string_view::npos;

    REQUIRE(is_ascii_letter_or_digit(c) == is_letter_or_digit);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'str_to_u32'

TEST_CASE("Proper numbers are converted", "[strutils::str_to_u32]")
{
  REQUIRE(str_to_u32("0") == 0);
  REQUIRE(str_to_u32("1") == 1);
  REQUIRE(str_to_u32("2") == 2);
  REQUIRE(str_to_u32("2147483647") == 2147483647);
  REQUIRE(str_to_u32("2147483648") == 2147483648);
  REQUIRE(str_to_u32("2147483649") == 2147483649);
  REQUIRE(str_to_u32("4294967293") == 4294967293);
  REQUIRE(str_to_u32("4294967294") == 4294967294);
  REQUIRE(str_to_u32("4294967295") == 4294967295);
}

TEST_CASE("Whitespace is ignored", "[strutils::str_to_u32]")
{
  REQUIRE(str_to_u32(" 123") == 123);
  REQUIRE(str_to_u32("123 ") == 123);
  REQUIRE(str_to_u32(" 123 ") == 123);
  REQUIRE(str_to_u32("\n123\n") == 123);
}

TEST_CASE("Numbers outside range are rejected", "[strutils::str_to_u32]")
{
  REQUIRE_THROWS_AS(str_to_u32("-2"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("-1"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("4294967296"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("4294967297"), std::invalid_argument);
}

TEST_CASE("Non numeric strings are rejected", "[strutils::str_to_u32]")
{
  REQUIRE_THROWS_AS(str_to_u32(" "), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("x"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("a b c"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("\n"), std::invalid_argument);
}

TEST_CASE("Mixed strings are rejected", "[strutils::str_to_u32]")
{
  REQUIRE_THROWS_AS(str_to_u32("1a"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("123x"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("123 x"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("x123"), std::invalid_argument);
  REQUIRE_THROWS_AS(str_to_u32("x 123"), std::invalid_argument);
}

TEST_CASE("Base 10 is assumed", "[strutils::str_to_u32]")
{
  REQUIRE(str_to_u32("010") == 10);
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
  REQUIRE(shell_escape("\x01") == "'\\x01'");
  REQUIRE(shell_escape("\x7f") == "'\\x7f'");
  REQUIRE(shell_escape("\xff") == "'\\xff'");
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
  REQUIRE(log_escape("\x01") == "'\\x01'");
  REQUIRE(log_escape("\x7f") == "'\\x7f'");
  REQUIRE(log_escape("\xff") == "'\\xff'");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'shell_escape_command'

TEST_CASE("shell_escape_command with empty list")
{
  REQUIRE(shell_escape_command({}).empty());
}

TEST_CASE("shell_escape_command with mixed command")
{
  REQUIRE(shell_escape_command({ "echo", "a", "m$x", "of", "$(things)" }) ==
          "echo a 'm$x' of '$(things)'");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'html_escape'

TEST_CASE("html_escape")
{
  REQUIRE(html_escape("") == "");
  REQUIRE(html_escape(" abc1-9 ") == " abc1-9 ");
  REQUIRE(html_escape(" '&<>' ") == " &#39;&amp;&lt;&gt;&#39; ");
  REQUIRE(html_escape("\"quote\"") == "&quot;quote&quot;");
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
  REQUIRE(format_percentage(55, 300) == "18.3 %");
  REQUIRE(format_percentage(55, 300, 0) == "18 %");
  REQUIRE(format_percentage(55, 300, 1) == "18.3 %");
  REQUIRE(format_percentage(55, 300, 2) == "18.33 %");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'join_text'

TEST_CASE("join_text on list of ints")
{

  using int_vec = std::vector<int>;

  REQUIRE(join_text(int_vec{}, ", ").empty());
  REQUIRE(join_text(int_vec{}, ", ", ", ").empty());

  REQUIRE(join_text(int_vec{ 1 }, ", ") == "1");

  REQUIRE(join_text(int_vec{ 1, 2 }, ", ") == "1, 2");
  REQUIRE(join_text(int_vec{ 1, 2, 3 }, ", ") == "1, 2, 3");

  REQUIRE(join_text(int_vec{ 1 }, ", ", ", or ") == "1");
  REQUIRE(join_text(int_vec{ 1, 2 }, ", ", ", or ") == "1, or 2");
  REQUIRE(join_text(int_vec{ 1, 2, 3 }, ", ", ", or ") == "1, 2, or 3");

  REQUIRE(join_text(int_vec{ 1, 2, 3 }, ". ", ", or ") == "1. 2, or 3");
}

TEST_CASE("join_text on list of strings")
{
  using str_vec = std::vector<std::string>;
  REQUIRE(join_text(str_vec{ "foo" }, ", ") == "foo");
  REQUIRE(join_text(str_vec{ "foo", "bar" }, ", ") == "foo, bar");
  REQUIRE(join_text(str_vec{ "foo", "bar", "zod" }, ", ") == "foo, bar, zod");
  REQUIRE(join_text(str_vec{ "foo", "bar", "zod" }, ", ", ", and ") ==
          "foo, bar, and zod");
}

TEST_CASE("join_text on list of strings views")
{
  using str_vec = std::vector<std::string_view>;
  REQUIRE(join_text(str_vec{ "foo" }, ", ") == "foo");
  REQUIRE(join_text(str_vec{ "foo", "bar" }, ", ") == "foo, bar");
  REQUIRE(join_text(str_vec{ "foo", "bar", "zod" }, ", ") == "foo, bar, zod");
  REQUIRE(join_text(str_vec{ "foo", "bar", "zod" }, ", ", ", and ") ==
          "foo, bar, and zod");
}

///////////////////////////////////////////////////////////////////////////////
// Tests for 'trim_ascii_whitespace'

TEST_CASE("trim_ascii_whitespace")
{
  REQUIRE(trim_ascii_whitespace("") == "");
  REQUIRE(trim_ascii_whitespace("\f\n\r\t\v") == "");
  REQUIRE(trim_ascii_whitespace("\f\n\r\t\vabc") == "abc");
  REQUIRE(trim_ascii_whitespace("abc\f\n\r\t\v") == "abc");
  REQUIRE(trim_ascii_whitespace("\f\n\r\t\vabc\f\n\r\t\v") == "abc");
  REQUIRE(trim_ascii_whitespace("\f\n\r\t\va b c\f\n\r\t\v") == "a b c");
  REQUIRE(trim_ascii_whitespace("abc") == "abc");
  REQUIRE(trim_ascii_whitespace("a b c") == "a b c");
}

} // namespace adapterremoval
