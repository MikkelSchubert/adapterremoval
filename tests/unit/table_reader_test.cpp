// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"       // for parsing_error, io_error, assert_failed
#include "linereader.hpp"   // for vec_reader, line_reader_base
#include "table_reader.hpp" // declarations
#include "testing.hpp"      // for TEST_CASE, REQUIRE, ...
#include <sstream>          // for ostringstream
#include <stdexcept>        // for std::out_of_range
#include <string>           // for string
#include <string_view>      // for string_view
#include <vector>           // for vector

namespace adapterremoval {

namespace {

class throwing_reader : public line_reader_base
{
public:
  bool getline(std::string& /* dst */) override
  {
    throw io_error("simulated io failure");
  }
};

} // namespace

////////////////////////////////////////////////////////////////////////////////
// table_row

TEST_CASE("table_row constructing empty row")
{
  SECTION("default constructor")
  {
    const table_row row;
    REQUIRE(row.line_num() == 0);
    REQUIRE(row.size() == 0);
    REQUIRE_THROWS_AS(row.at(0), std::out_of_range);
  }

  SECTION("explicitly empty")
  {
    const table_row row{ 0, {} };
    REQUIRE(row.line_num() == 0);
    REQUIRE(row.size() == 0);
    REQUIRE_THROWS_AS(row.at(0), std::out_of_range);
  }
}

TEST_CASE("table_row constructing row with values")
{
  const table_row row{ 1234, { "alpha", "beta", "gamma" } };
  REQUIRE(row.line_num() == 1234);
  REQUIRE(row.size() == 3);
  REQUIRE(row.at(0) == "alpha");
  REQUIRE(row.at(1) == "beta");
  REQUIRE(row.at(2) == "gamma");
  REQUIRE_THROWS_AS(row.at(4), std::out_of_range);
}

TEST_CASE("table_row equality operator")
{
  CHECK(table_row{} == table_row{});
  CHECK_FALSE(table_row{} != table_row{});

  CHECK(table_row{ 1, { "x", "y" } } == table_row{ 1, { "x", "y" } });
  CHECK_FALSE(table_row{ 1, { "x", "y" } } != table_row{ 1, { "x", "y" } });

  CHECK_FALSE(table_row{ 2, { "x", "y" } } == table_row{ 1, { "x", "y" } });
  CHECK(table_row{ 2, { "x", "y" } } != table_row{ 1, { "x", "y" } });

  CHECK_FALSE(table_row{ 1, { "x", "y" } } == table_row{ 1, { "y", "x" } });
  CHECK(table_row{ 1, { "x", "y" } } != table_row{ 1, { "y", "x" } });
}

TEST_CASE("table_row to string")
{
  std::ostringstream os;

  SECTION("empty row")
  {
    os << table_row{};
    CHECK(os.str() == "table_row{line_num=NA, values=[]}");
  }

  SECTION("non-empty row")
  {
    os << table_row{ 0, { "foo", "bar" } };
    CHECK(os.str() == "table_row{line_num=NA, values=['foo', 'bar']}");
  }

  SECTION("non-empty row with line number")
  {
    os << table_row{ 3, { "foo", "bar" } };
    CHECK(os.str() == "table_row{line_num=3, values=['foo', 'bar']}");
  }
}

////////////////////////////////////////////////////////////////////////////////
// table_reader builder methods

TEST_CASE("table_reader returns self for chaining")
{
  table_reader reader;
  REQUIRE(&reader.with_name("my table") == &reader);
  REQUIRE(&reader.with_comment_char('#') == &reader);
  REQUIRE(&reader.with_min_columns(1) == &reader);
  REQUIRE(&reader.with_max_columns(5) == &reader);
}

TEST_CASE("table_reader columns must be non-zero")
{
  REQUIRE_THROWS_AS(table_reader{}.with_min_columns(0), assert_failed);
  REQUIRE_THROWS_AS(table_reader{}.with_max_columns(0), assert_failed);
}

TEST_CASE("table_reader min columns must be no more than max columns")
{
  table_reader reader;
  reader.with_max_columns(2);
  REQUIRE_THROWS_AS(reader.with_min_columns(3), assert_failed);
  REQUIRE_NOTHROW(reader.with_min_columns(2));
  REQUIRE_NOTHROW(reader.with_min_columns(1));
}

TEST_CASE("table_reader max columns must be at least min columns")
{
  table_reader reader;
  reader.with_min_columns(3);
  REQUIRE_THROWS_AS(reader.with_max_columns(2), assert_failed);
  REQUIRE_NOTHROW(reader.with_max_columns(3));
  REQUIRE_NOTHROW(reader.with_max_columns(4));
}

////////////////////////////////////////////////////////////////////////////////
// table_reader - parsing

TEST_CASE("table_reader parse on empty input returns empty table #1")
{
  vec_reader lr({});
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.empty());
}

TEST_CASE("table_reader parse on empty input returns empty table #2")
{
  vec_reader lr({ "" });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.empty());
}

TEST_CASE("table_reader parse single row single column")
{
  vec_reader lr({ "hello" });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 1, { "hello" } });
}

TEST_CASE("table_reader parse single row multiple columns")
{
  vec_reader lr({ "alpha beta gamma" });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 1, { "alpha", "beta", "gamma" } });
}

TEST_CASE("table_reader parse multiple rows")
{
  vec_reader lr({ "a b", "c d", "e f" });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.size() == 3);
  REQUIRE(rows.at(0) == table_row{ 1, { "a", "b" } });
  REQUIRE(rows.at(1) == table_row{ 2, { "c", "d" } });
  REQUIRE(rows.at(2) == table_row{ 3, { "e", "f" } });
}

TEST_CASE("table_reader empty lines are ignored, but reflect in line numbers")
{
  vec_reader lr({ "", "a b", " ", "c d", "\t" });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.size() == 2);
  REQUIRE(rows.at(0) == table_row{ 2, { "a", "b" } });
  REQUIRE(rows.at(1) == table_row{ 4, { "c", "d" } });
}

TEST_CASE("table_reader all whitespace is ignored")
{
  vec_reader lr({ "  a \t b   c \r\n" });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 1, { "a", "b", "c" } });
}

TEST_CASE("table_reader parse inconsistent rows throws")
{
  vec_reader lr({ "a b", "c d e" });
  REQUIRE_THROWS_MESSAGE(table_reader{}.parse(lr),
                         parsing_error,
                         "Error at line 2: Inconsistent number of columns; "
                         "expected 2 column(s), but found 3 column(s)");
}

////////////////////////////////////////////////////////////////////////////////
// table_reader - comments

TEST_CASE("table_reader does not handle comments by default (#)")
{
  vec_reader lr({ "a b # c" });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 1, { "a", "b", "#", "c" } });
}

TEST_CASE("table_reader does not handle comments by default (\0)")
{
  using namespace std::string_literals;

  // '\0' is the sentinel value and should not be treated as a comment char
  vec_reader lr({ "a b \0 c"s });
  const auto rows = table_reader{}.parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 1, { "a", "b", "\0"s, "c" } });
}

TEST_CASE("table_reader parse strips inline comment")
{
  vec_reader lr({ "a b # this is ignored", "c d" });
  const auto rows = table_reader{}.with_comment_char('#').parse(lr);
  REQUIRE(rows.size() == 2);
  REQUIRE(rows.at(0) == table_row{ 1, { "a", "b" } });
  REQUIRE(rows.at(1) == table_row{ 2, { "c", "d" } });
}

TEST_CASE("table_reader parse skips full-line comment")
{
  vec_reader lr({ "#", "# this is a comment", "a b" });
  const auto rows = table_reader{}.with_comment_char('#').parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 3, { "a", "b" } });
}

TEST_CASE("table_reader parse comment char terminates line")
{
  vec_reader lr({ "a#b" });
  const auto rows = table_reader{}.with_comment_char('#').parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 1, { "a" } });
}

TEST_CASE("table_reader parse all comment lines returns empty table")
{
  vec_reader lr({ "# first", "# second", "# third" });
  const auto rows = table_reader{}.with_comment_char('#').parse(lr);
  REQUIRE(rows.empty());
}

TEST_CASE("table_reader parse hash comment with whitespace before it")
{
  vec_reader lr({ "a  # comment" });
  const auto rows = table_reader{}.with_comment_char('#').parse(lr);
  REQUIRE(rows.size() == 1);
  REQUIRE(rows.at(0) == table_row{ 1, { "a" } });
}

////////////////////////////////////////////////////////////////////////////////
// table_reader - column count enforcement

TEST_CASE("table_reader parse accepts exact min_columns")
{
  vec_reader lr({ "a b" });
  REQUIRE_NOTHROW(table_reader{}.with_min_columns(2).parse(lr));
}

TEST_CASE("table_reader parse accepts exact max_columns")
{
  vec_reader lr({ "a b c" });
  REQUIRE_NOTHROW(table_reader{}.with_max_columns(3).parse(lr));
}

TEST_CASE("table_reader parse accepts column count within min-max range")
{
  vec_reader lr({ "a b c" });
  REQUIRE_NOTHROW(
    table_reader{}.with_min_columns(2).with_max_columns(4).parse(lr));
}

TEST_CASE("table_reader parse min columns is enforced")
{
  vec_reader lr({ "a b", "c" });
  REQUIRE_THROWS_MESSAGE(
    table_reader{}.with_min_columns(2).parse(lr),
    parsing_error,
    "Error at line 2: Expected at least 2 column(s), but found 1 column(s)");
}

TEST_CASE("table_reader parse max columns is enforced")
{
  vec_reader lr({ "a b", "c d e" });
  REQUIRE_THROWS_MESSAGE(
    table_reader{}.with_max_columns(2).parse(lr),
    parsing_error,
    "Error at line 2: Expected at most 2 column(s), but found 3 column(s)");
}

TEST_CASE(
  "table_reader parse rows within min-max but inconsistent still throws")
{
  vec_reader lr({ "a b", "c d e" });
  REQUIRE_THROWS_MESSAGE(
    table_reader{}.with_min_columns(2).with_max_columns(4).parse(lr),
    parsing_error,
    "Error at line 2: Inconsistent number of columns; expected 2 column(s), "
    "but found 3 column(s)");
}

////////////////////////////////////////////////////////////////////////////////
// table_reader - with named table

TEST_CASE("table_reader name appears in min_columns error")
{
  vec_reader lr({ "a" });
  REQUIRE_THROWS_WITH(
    table_reader{}.with_name("myfile.txt").with_min_columns(2).parse(lr),
    Catch::Contains("in table myfile.txt"));
}

TEST_CASE("table_reader name appears in max_columns error")
{
  vec_reader lr({ "a b c" });
  REQUIRE_THROWS_WITH(
    table_reader{}.with_name("data.tsv").with_max_columns(2).parse(lr),
    Catch::Contains("in table data.tsv"));
}

TEST_CASE("table_reader name appears in inconsistent columns error")
{
  vec_reader lr({ "a b", "c d e" });
  REQUIRE_THROWS_WITH(table_reader{}.with_name("test.tsv").parse(lr),
                      Catch::Contains("in table test.tsv"));
}

TEST_CASE("table_reader parse without name omits 'in table' from error")
{
  vec_reader lr({ "a" });
  REQUIRE_THROWS_WITH(table_reader{}.with_min_columns(2).parse(lr),
                      !Catch::Contains("in table"));
}

TEST_CASE("table_reader parse error includes both line number and table name")
{
  vec_reader lr({ "a" });
  REQUIRE_THROWS_WITH(
    table_reader{}.with_name("src.txt").with_min_columns(2).parse(lr),
    Catch::Contains("at line 1"));
}

////////////////////////////////////////////////////////////////////////////////
// table_reader::parse - io_error handling

TEST_CASE("table_reader parse rethrows io_error as io_error")
{
  throwing_reader lr;
  REQUIRE_THROWS_AS(table_reader{}.parse(lr), io_error);
}

TEST_CASE("table_reader parse io_error message is preserved")
{
  throwing_reader lr;
  REQUIRE_THROWS_WITH(table_reader{}.parse(lr),
                      Catch::Contains("simulated io failure"));
}

TEST_CASE("table_reader parse io_error with name includes table name")
{
  throwing_reader lr;
  REQUIRE_THROWS_WITH(table_reader{}.with_name("mydata").parse(lr),
                      Catch::Contains("in table mydata"));
}

TEST_CASE("table_reader parse io_error without name omits 'in table'")
{
  throwing_reader lr;
  REQUIRE_THROWS_WITH(table_reader{}.parse(lr), !Catch::Contains("in table"));
}

} // namespace adapterremoval
