// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "argparse.hpp" // for argument, parser, str_sink, parse_result
#include "errors.hpp"   // for assert_failed
#include "logging.hpp"  // for log_capture
#include "strutils.hpp" // for string_vec
#include "testing.hpp"  // for TEST_CASE, REQUIRE, ...
#include <cstddef>      // for size_t
#include <cstdint>      // for uint32_t
#include <limits>       // for numeric_limits
#include <stdexcept>    // for invalid_argument
#include <string>       // for basic_string, operator==, string, allocator
#include <string_view>  // for string_view
#include <vector>       // for vector, operator==

namespace adapterremoval {

using argparse::argument_ptr;

const size_t invalid_choice = static_cast<size_t>(-2);

using Catch::Matchers::Contains;

#define REQUIRE_POSTFIX(a, b) REQUIRE_THAT((a), Catch::Matchers::EndsWith(b))
#define REQUIRE_CONTAINS(a, b) REQUIRE_THAT((a), Catch::Matchers::Contains(b))

///////////////////////////////////////////////////////////////////////////////
// boolean sink

TEST_CASE("bool sink is initialized", "[argparse::bool_sink]")
{
  bool value = true;
  argparse::bool_sink sink(&value);
  REQUIRE_FALSE(value);
}

TEST_CASE("bool sink value", "[argparse::bool_sink]")
{
  bool value = false;
  argparse::bool_sink sink(&value);
  REQUIRE(sink.value() == "off");
  value = true;
  REQUIRE(sink.value() == "on");
}

TEST_CASE("bool sink has_default", "[argparse::bool_sink]")
{
  bool value = false;
  argparse::bool_sink sink(&value);
  REQUIRE_FALSE(sink.has_default());
}

TEST_CASE("bool sink does not require values", "[argparse::bool_sink]")
{
  bool value = false;
  argparse::bool_sink sink(&value);

  string_vec values;
  REQUIRE(sink.consume(values.begin(), values.end()) == 0);
  REQUIRE(sink.value() == "on");
  REQUIRE(value);
}

TEST_CASE("bool sink rejects values", "[argparse::bool_sink]")
{
  string_vec values{ "0" };
  bool value = false;
  argparse::bool_sink sink(&value);

  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("bool does not support default values", "[argparse::bool_sink]")
{
  bool value = false;
  argparse::bool_sink sink(&value);

  REQUIRE_FALSE(sink.has_default());
  REQUIRE_THROWS_AS(sink.default_value(), assert_failed);
}

///////////////////////////////////////////////////////////////////////////////
// uint sink

TEST_CASE("uint sink is required", "[argparse::u32_sink]")
{
  REQUIRE_THROWS_AS(argparse::u32_sink(nullptr), assert_failed);
}

TEST_CASE("uint sink is initialized", "[argparse::u32_sink]")
{
  uint32_t value = 12345;
  argparse::u32_sink sink(&value);
  REQUIRE(value == 0);
}

TEST_CASE("uint sink value", "[argparse::u32_sink]")
{
  uint32_t value = 1234567;
  argparse::u32_sink sink(&value);
  REQUIRE(sink.value() == "0");
  value = 1234567;
  REQUIRE(sink.value() == "1234567");
}

TEST_CASE("uint sink with_default", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  REQUIRE_FALSE(sink.has_default());
  REQUIRE_THROWS_AS(sink.default_value(), assert_failed);
  sink.with_default(1234567);

  REQUIRE(sink.has_default());
  REQUIRE(sink.default_value() == "1234567");
  REQUIRE(value == 1234567);
}

TEST_CASE("uint sink requires value", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  string_vec values;
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("uint sink consumes single value", "[argparse::u32_sink]")
{
  string_vec values{ "1234567" };

  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == 1234567);
}

TEST_CASE("uint sink consumes single value only", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  const string_vec values{ "1234567", "8901234" };
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("uint requires valid integer value", "[argparse::u32_sink]")
{
  string_vec values{ "abc" };

  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()),
                    std::invalid_argument);
  REQUIRE(value == 0);
}

TEST_CASE("uint requires positive integer value", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  string_vec values{ "-123" };
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()),
                    std::invalid_argument);
  REQUIRE(value == 0);
}

TEST_CASE("uint accepts uint32_t lower bound", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);
  sink.with_default(123456);

  string_vec values{ "0" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == 0);
}

TEST_CASE("uint accepts uint32_t upper bound", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  string_vec values{ "4294967295" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == 4294967295);
}

TEST_CASE("uint rejects past uint32_t upper bound", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  string_vec values{ "4294967296" };
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()),
                    std::invalid_argument);
  REQUIRE(value == 0);
}

TEST_CASE("uint disallows trailing garbage", "[argparse::u32_sink]")
{
  string_vec values{ "123abc" };

  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()),
                    std::invalid_argument);
  REQUIRE(value == 0);
}

TEST_CASE("uint with_minimum requires valid default")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);

  SECTION("fail")
  {
    sink.with_minimum(0);
    REQUIRE_THROWS_AS(sink.with_minimum(1), assert_failed);
  }

  SECTION("pass")
  {
    sink.with_default(101).with_minimum(100);
    sink.with_default(100).with_minimum(100);
  }
}

TEST_CASE("with uint minimum", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);
  sink.with_default(100);
  sink.with_minimum(100);

  SECTION("pass")
  {
    const auto expected = GENERATE(100LU, 101LU, 4294967295LU);
    string_vec values{ std::to_string(expected) };
    REQUIRE(sink.consume(values.begin(), values.end()) == 1);
    REQUIRE(value == expected);
  }

  SECTION("fail")
  {
    string_vec values{ "99" };
    REQUIRE_THROWS_MESSAGE(sink.consume(values.begin(), values.end()),
                           std::invalid_argument,
                           "value must be at least 100");
    REQUIRE(value != 0);
  }
}

TEST_CASE("uint with_maximum requires valid default")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);
  sink.with_default(100);

  SECTION("fail")
  {
    sink.with_maximum(100);
    REQUIRE_THROWS_AS(sink.with_maximum(99), assert_failed);
  }

  SECTION("pass")
  {
    sink.with_maximum(100);
    sink.with_maximum(101);
    sink.with_maximum(std::numeric_limits<uint32_t>::max());
  }
}

TEST_CASE("with uint maximum", "[argparse::u32_sink]")
{
  uint32_t value = 0;
  argparse::u32_sink sink(&value);
  sink.with_maximum(100);

  SECTION("pass")
  {
    const auto expected = GENERATE(0LU, 99LU, 100LU);
    string_vec values{ std::to_string(expected) };
    REQUIRE(sink.consume(values.begin(), values.end()) == 1);
    REQUIRE(value == expected);
  }

  SECTION("fail")
  {
    string_vec values{ "101" };
    REQUIRE_THROWS_MESSAGE(sink.consume(values.begin(), values.end()),
                           std::invalid_argument,
                           "value must be at most 100");
    REQUIRE(value != 101);
  }
}

///////////////////////////////////////////////////////////////////////////////
// double sink

TEST_CASE("double sink is required", "[argparse::double_sink]")
{
  REQUIRE_THROWS_AS(argparse::double_sink(nullptr), assert_failed);
}

TEST_CASE("double sink is initialized", "[argparse::double_sink]")
{
  double value = 12345;
  argparse::double_sink sink(&value);
  REQUIRE(value == 0);
}

TEST_CASE("double sink value", "[argparse::double_sink]")
{
  double value = 12345.67;
  argparse::double_sink sink(&value);
  REQUIRE(sink.value() == "0");
  value = 12345.67;
  REQUIRE(sink.value() == "12345.67");
  value = 1234567;
  REQUIRE(sink.value() == "1234567");
}

TEST_CASE("double sink with_default", "[argparse::double_sink]")
{
  double value = 0;
  argparse::double_sink sink(&value);
  sink.with_default(12345.67);

  REQUIRE(value == 12345.67);
}

TEST_CASE("double sink has_default", "[argparse::double_sink]")
{
  double value = 0;
  argparse::double_sink sink(&value);
  REQUIRE_FALSE(sink.has_default());
  REQUIRE_THROWS_AS(sink.default_value(), assert_failed);
  sink.with_default(12345.67);
  REQUIRE(sink.has_default());
  REQUIRE(sink.default_value() == "12345.67");
}

TEST_CASE("double sink trims trailing zeroes in default value",
          "[argparse::double_sink]")
{
  double value = 0;
  argparse::double_sink sink(&value);

  REQUIRE_FALSE(sink.has_default());
  sink.with_default(1.0);
  REQUIRE(sink.default_value() == "1");
}

TEST_CASE("double sink require value", "[argparse::double_sink]")
{
  double value = 0;
  argparse::double_sink sink(&value);

  string_vec values;
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("double sink consumes single value", "[argparse::double_sink]")
{
  string_vec values{ "12345.67" };

  double value = 0;
  argparse::double_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == 12345.67);
}

TEST_CASE("double sink consumes single value only", "[argparse::double_sink]")
{

  double value = 0;
  argparse::double_sink sink(&value);

  const string_vec values{ "12345.67", "8901234" };
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("double requires valid double #1", "[argparse::double_sink]")
{
  string_vec values{ "abc" };

  double value = 0;
  argparse::double_sink sink(&value);

  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()),
                    std::invalid_argument);
  REQUIRE(value == 0);
}

TEST_CASE("double requires valid double #2", "[argparse::double_sink]")
{
  string_vec values{ "123.0abc" };

  double value = 0;
  argparse::double_sink sink(&value);

  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()),
                    std::invalid_argument);
  REQUIRE(value == 0);
}

TEST_CASE("double with min/max value", "[argparse::double_sink]")
{
  double value = 0;
  argparse::double_sink sink(&value);

  string_vec values{ "0" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == Approx(0));

  SECTION("minimum")
  {
    sink.with_minimum(-5.0);

    values = string_vec{ "-5.0" };
    REQUIRE(sink.consume(values.begin(), values.end()) == 1);
    REQUIRE(value == Approx(-5.0));

    values = string_vec{ "-5.1" };
    REQUIRE_THROWS_MESSAGE(sink.consume(values.begin(), values.end()),
                           std::invalid_argument,
                           "value must be at least -5");
  }

  SECTION("maximum")
  {
    sink.with_maximum(5.0);

    values = string_vec{ "5.0" };
    REQUIRE(sink.consume(values.begin(), values.end()) == 1);
    REQUIRE(value == Approx(5.0));

    values = string_vec{ "5.1" };
    REQUIRE_THROWS_MESSAGE(sink.consume(values.begin(), values.end()),
                           std::invalid_argument,
                           "value must be at most 5");
  }
}

///////////////////////////////////////////////////////////////////////////////
// str sink

TEST_CASE("str sink is initialized", "[argparse::str_sink]")
{
  std::string value = "foo";
  argparse::str_sink sink(&value);
  REQUIRE(value.empty());
}

TEST_CASE("str sink value", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  value = "foobar";
  REQUIRE(sink.value() == "foobar");
}

TEST_CASE("str sink value escapes", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  value = "foo bar";
  REQUIRE(sink.value() == "foo bar");
}

TEST_CASE("str sink with_default (char*)", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  REQUIRE(value == "");
  REQUIRE_THROWS_AS(sink.default_value(), assert_failed);
  sink.with_default("foobar");
  REQUIRE(value == "foobar");
  REQUIRE(sink.default_value() == "foobar");
}

TEST_CASE("str sink has_default (char*)", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  REQUIRE_FALSE(sink.has_default());
  sink.with_default("foobar");
  REQUIRE(sink.has_default());
}

TEST_CASE("str sink with_default (std::string)", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_default(std::string("foobar"));
  REQUIRE(value == "foobar");
}

TEST_CASE("str sink has_default", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  REQUIRE_FALSE(sink.has_default());
  sink.with_default(std::string("foobar"));
  REQUIRE(sink.has_default());
}

TEST_CASE("str sink require value", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);

  const string_vec values;
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("str sink consumes single value", "[argparse::str_sink]")
{
  string_vec values{ "foobar" };

  std::string value;
  argparse::str_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == "foobar");
}

TEST_CASE("str sink consumes single value only", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);

  const string_vec values{ "foo", "bar" };
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("str sink consumes empty string", "[argparse::str_sink]")
{
  string_vec values{ "" };

  std::string value;
  argparse::str_sink sink(&value);
  sink.with_default("foo");

  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value.empty());
}

TEST_CASE("str sink accepts any value if no choices", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_choices({});

  string_vec values{ "ghi" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == "ghi");
}

TEST_CASE("str sink accepts value in choices", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_choices({ "abc", "def", "ghi" });

  string_vec values{ GENERATE("ghi", "GHI") };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == "ghi");
}

TEST_CASE("str sink rejects values not in choices", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_choices({ "abc", "def", "ghi" });

  string_vec values{ "foo" };
  REQUIRE(sink.consume(values.begin(), values.end()) == invalid_choice);
  REQUIRE(value.empty());
}

TEST_CASE("str sink case-insensitive, returns choice", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_choices({ "Abc", "dEf", "ghI" });

  string_vec values{ "Ghi" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == "ghI");
}

TEST_CASE("str sink choices are set", "[argparse::str_sink]")
{
  const string_vec choices{ "Abc", "dEf", "ghI" };

  std::string value;
  argparse::str_sink sink(&value);
  sink.with_choices(choices);

  REQUIRE(sink.choices() == choices);
}

TEST_CASE("str sink with implicit argument", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_implicit_argument("foo");
  REQUIRE(value == "");

  string_vec values{ "bar" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == "bar");

  REQUIRE(sink.consume(values.end(), values.end()) == 0);
  REQUIRE(value == "foo");
}

TEST_CASE("str with pre-processor returns modified value", "[argparse::parser]")
{
  std::string sink;
  argparse::parser p;
  auto& arg = p.add("--foo").bind_str(&sink).with_preprocessor(
    [](std::string& value) { value = to_lower(value); });

  REQUIRE(p.parse_args({ "exe", "--foo", "ABC" }) ==
          argparse::parse_result::ok);
  REQUIRE(arg.value() == "abc");
}

///////////////////////////////////////////////////////////////////////////////
// vec sink

TEST_CASE("vec sink is required", "[argparse::vec_sink]")
{
  REQUIRE_THROWS_AS(argparse::vec_sink(nullptr), assert_failed);
}

TEST_CASE("vec sink is initialized", "[argparse::vec_sink]")
{
  string_vec value{ "foo" };
  argparse::vec_sink sink(&value);
  REQUIRE(value.empty());
}

TEST_CASE("vec sink to_vec #1", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);
  value.push_back("foobar");
  REQUIRE(sink.value() == "foobar");
}

TEST_CASE("vec sink to_vec #2", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);
  value.push_back("foo");
  value.push_back("bar");
  REQUIRE(sink.value() == "foo;bar");
}

TEST_CASE("vec sink to_vec escapes", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);
  value.push_back("foo");
  value.push_back("1 2");
  value.push_back("bar");
  REQUIRE(sink.value() == "foo;'1 2';bar");
}

TEST_CASE("vec sink require value", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);

  const string_vec values;
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("vec sink consumes single value", "[argparse::vec_sink]")
{
  string_vec values{ "foobar" };

  string_vec value;
  argparse::vec_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value.size() == 1);
  REQUIRE(value.front() == "foobar");
}

TEST_CASE("vec sink consumes multiple values", "[argparse::vec_sink]")
{
  string_vec values{ "foo", "bar" };

  string_vec value;
  argparse::vec_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == 2);
  REQUIRE(value.size() == 2);
  REQUIRE(value.front() == "foo");
  REQUIRE(value.back() == "bar");
}

TEST_CASE("vec sink consumes empty string", "[argparse::vec_sink]")
{
  string_vec values{ "" };

  string_vec value;
  argparse::vec_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value.size() == 1);
  REQUIRE(value.front().empty());
}

TEST_CASE("vec sink with minimum n values", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);
  sink.with_min_values(2);
  REQUIRE(sink.min_values() == 2);
}

TEST_CASE("vec sink with maximum n values", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);

  sink.with_max_values(2);
  REQUIRE(sink.max_values() == 2);
}

TEST_CASE("vec does not support default values", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);

  REQUIRE_FALSE(sink.has_default());
  REQUIRE_THROWS_AS(sink.default_value(), assert_failed);
}

///////////////////////////////////////////////////////////////////////////////
// argument

TEST_CASE("argument properties", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");

  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.key() == "--12345");
  REQUIRE(arg.metavar() == "67890");
}

TEST_CASE("help without default", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");
  REQUIRE(arg.help().empty());

  arg.help("pointless gibberish");
  REQUIRE(arg.help() == "pointless gibberish");
}

TEST_CASE("help with default", "[argparse::argument]")
{
  uint32_t value = 0;
  argparse::argument arg("--12345", "67890");
  arg.help("pointless gibberish").bind_u32(&value).with_default(7913);

  REQUIRE(arg.help() == "pointless gibberish [default: 7913]");
}

TEST_CASE("help with custom default", "[argparse::argument]")
{
  uint32_t value = 0;
  argparse::argument arg("--12345", "67890");
  arg.help("pointless gibberish [default: foo]")
    .bind_u32(&value)
    .with_default(7913);

  REQUIRE(arg.help() == "pointless gibberish [default: foo]");
}

TEST_CASE("help with choices", "[argparse::argument]")
{
  std::string value;
  argparse::argument arg("--12345", "67890");
  auto& sink = arg.help("gibberish").bind_str(&value);

  SECTION("single choice")
  {
    sink.with_choices({ "foo" });
    REQUIRE(arg.help() == "gibberish. Possible values are foo");
  }

  SECTION("two choices")
  {
    sink.with_choices({ "foo", "bar" });
    REQUIRE(arg.help() == "gibberish. Possible values are foo, and bar");
  }

  SECTION("three choices")
  {
    sink.with_choices({ "foo", "bar", "x" });
    REQUIRE(arg.help() == "gibberish. Possible values are foo, bar, and x");
  }
}

TEST_CASE("deprecated argument", "[argparse::argument]")
{
  argparse::argument arg("--12345");

  REQUIRE_FALSE(arg.is_deprecated());
  REQUIRE_FALSE(arg.is_hidden());
  arg.deprecated();
  REQUIRE(arg.is_deprecated());
  REQUIRE(arg.is_hidden());

  log::log_capture ss;

  string_vec values{ "--12345" };
  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE_POSTFIX(ss.str(),
                  "[WARNING] Option --12345 is deprecated and will be "
                  "removed in the future.\n");
}

TEST_CASE("removed argument without help", "[argparse::argument]")
{
  argparse::argument arg("--12345");

  REQUIRE_FALSE(arg.is_removed());
  REQUIRE_FALSE(arg.is_hidden());
  arg.removed();
  REQUIRE(arg.is_removed());
  REQUIRE(arg.is_hidden());

  log::log_capture ss;

  string_vec values{ "--12345" };
  REQUIRE(arg.parse(values.begin(), values.end()) ==
          argparse::argument::parsing_failed);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Option --12345 has been removed in AdapterRemoval v3. For more "
    "information, see "
    "https://adapterremoval.readthedocs.io/en/stable/migrating.html\n");
}

TEST_CASE("removed argument with help", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  arg.help("additional help text").removed();

  log::log_capture ss;

  string_vec values{ "--12345" };
  REQUIRE(arg.parse(values.begin(), values.end()) ==
          argparse::argument::parsing_failed);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Option --12345 has been removed in AdapterRemoval v3: additional "
    "help text. For more information, see "
    "https://adapterremoval.readthedocs.io/en/stable/migrating.html\n");
}

TEST_CASE("hidden argument", "[argparse::argument]")
{
  argparse::argument arg("--12345");

  REQUIRE_FALSE(arg.is_deprecated());
  REQUIRE_FALSE(arg.is_hidden());
  arg.hidden();
  REQUIRE_FALSE(arg.is_deprecated());
  REQUIRE(arg.is_hidden());
}

TEST_CASE("default argument sink", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");

  REQUIRE(arg.value() == "off");
}

TEST_CASE("canonical key", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");

  REQUIRE(arg.keys() == string_vec{ "--12345" });
  REQUIRE_FALSE(arg.is_deprecated_alias("--12345"));
}

TEST_CASE("short argument alias", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");
  arg.abbreviation('a');

  REQUIRE(arg.key() == "--12345");
  REQUIRE(arg.short_key() == "-a");
  REQUIRE(arg.keys() == string_vec{ "--12345", "-a" });
  REQUIRE_FALSE(arg.is_deprecated_alias("--12345"));
  REQUIRE_FALSE(arg.is_deprecated_alias("-1"));
}

TEST_CASE("deprecated argument alias", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");
  arg.deprecated_alias("--foo");

  REQUIRE(arg.key() == "--12345");
  REQUIRE(arg.short_key().empty());
  REQUIRE(arg.keys() == string_vec{ "--foo", "--12345" });
  REQUIRE_FALSE(arg.is_deprecated_alias("--12345"));
  REQUIRE(arg.is_deprecated_alias("--foo"));
}

TEST_CASE("argument requires", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");
  arg.depends_on("--bar");
  arg.depends_on("--zod");

  REQUIRE(arg.depends_on() == string_vec{ "--bar", "--zod" });
}

TEST_CASE("required option does not exist", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").depends_on("--bar");

  REQUIRE_THROWS_AS(p.parse_args({ "exe", "--foo" }), assert_failed);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Option --foo requires unknown command-line option --bar\n");
}

TEST_CASE("argument conflicts with", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");
  arg.conflicts_with("--bar");
  arg.conflicts_with("--zod");

  REQUIRE(arg.conflicts_with() == string_vec{ "--bar", "--zod" });
}

TEST_CASE("default bind", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345" };

  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE(arg.is_set());
}

TEST_CASE("bind bool", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345" };

  bool sink = false;
  arg.bind_bool(&sink);

  REQUIRE_FALSE(sink);
  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE(arg.is_set());
  REQUIRE(sink);
}

TEST_CASE("bind uint", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345", "7913" };

  uint32_t sink = 0;
  arg.bind_u32(&sink);

  REQUIRE(sink == 0);
  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.parse(values.begin(), values.end()) == 2);
  REQUIRE(arg.is_set());
  REQUIRE(sink == 7913);
}

TEST_CASE("bind double", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345", "7.913" };

  double sink = 0;
  arg.bind_double(&sink);

  REQUIRE(sink == 0);
  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.parse(values.begin(), values.end()) == 2);
  REQUIRE(arg.is_set());
  REQUIRE_THAT(sink, Catch::WithinAbs(7.913, 1e-6));
}

TEST_CASE("bound double takes negative values", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345", "-123" };

  double sink = 0;
  arg.bind_double(&sink);

  REQUIRE(sink == 0);
  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.parse(values.begin(), values.end()) == 2);
  REQUIRE(arg.is_set());
  REQUIRE(sink == -123);
}

TEST_CASE("bind str", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345", "abcdef" };

  std::string sink;
  arg.bind_str(&sink);

  REQUIRE(sink.empty());
  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.parse(values.begin(), values.end()) == 2);
  REQUIRE(arg.is_set());
  REQUIRE(sink == "abcdef");
}

TEST_CASE("bind vec", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345", "abcdef", "7913" };

  string_vec sink;
  arg.bind_vec(&sink);

  REQUIRE(sink.empty());
  REQUIRE_FALSE(arg.is_set());
  REQUIRE(arg.parse(values.begin(), values.end()) == 3);
  REQUIRE(arg.is_set());
  REQUIRE(sink == string_vec{ "abcdef", "7913" });
}

TEST_CASE("bind vec clears previous values", "[argparse::argument]")
{
  string_vec sink;
  argparse::argument arg("--12345");
  arg.bind_vec(&sink);

  string_vec values_1{ "--12345", "abcdef", "7913" };
  REQUIRE(arg.parse(values_1.begin(), values_1.end()) == 3);
  string_vec values_2{ "--12345", "foo", "bar" };
  REQUIRE(arg.parse(values_2.begin(), values_2.end()) == 3);
  REQUIRE(sink == string_vec{ "foo", "bar" });
}

TEST_CASE("parse wrong argument asserts", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--54321" };

  REQUIRE_THROWS_AS(arg.parse(values.begin(), values.end()), assert_failed);
}

TEST_CASE("warning on first duplicate argument", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  string_vec values{ "--12345" };

  log::log_capture ss;

  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE(ss.str().empty());
  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE_POSTFIX(ss.str(),
                  "[WARNING] Command-line option --12345 has been specified "
                  "more than once.\n");
  // Only display the warning once per argument
  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE_POSTFIX(ss.str(),
                  "[WARNING] Command-line option --12345 has been specified "
                  "more than once.\n");
}

TEST_CASE("no warning on main alias", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  arg.deprecated_alias("--foo");
  string_vec values{ "--12345" };

  log::log_capture ss;

  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE(ss.str().empty());
}

TEST_CASE("warning on deprecated alias", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  arg.deprecated_alias("--foo");
  string_vec values{ "--foo" };

  log::log_capture ss;

  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE_POSTFIX(ss.str(),
                  "[WARNING] Option --foo has been renamed to --12345. "
                  "Support for the old name will be removed in the future.\n");
}

///////////////////////////////////////////////////////////////////////////////
// parser

namespace {

const std::string_view HELP_HEADER =
  "My App v1234\n\n"
  "basic help\n\n"
  "OPTIONS:\n"
  "   -h, --help      Display this message\n"
  "   -v, --version   Print the version string\n"
  "   --licenses      Print licenses for this software\n\n";

} // namespace

TEST_CASE("--version", "[argparse::parser]")
{
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");

  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", GENERATE("-v", "--version") }) ==
          argparse::parse_result::exit);
  REQUIRE(ss.str() == "My App v1234\n");
}

TEST_CASE("--licenses", "[argparse::parser]")
{
  argparse::parser p;
  p.set_name("My App");
  p.set_licenses("line 1\nline 2\nline 3");

  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "--licenses" }) ==
          argparse::parse_result::exit);
  REQUIRE(ss.str() == "line 1\nline 2\nline 3\n");
}

TEST_CASE("--help", "[argparse::parser]")
{
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.set_preamble("basic help");

  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", GENERATE("-h", "--help") }) ==
          argparse::parse_result::exit);
  REQUIRE(ss.str() == HELP_HEADER);
}

TEST_CASE("--help with deprecated argument", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--foo").deprecated();
  const string_vec args = { "exe", "--help" };

  log::log_capture ss;

  REQUIRE(p.parse_args(args) == argparse::parse_result::exit);
  REQUIRE_THAT(ss.str(), !Contains("--foo"));
}

TEST_CASE("--help with hidden argument", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--foo").hidden();

  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "--help" }) == argparse::parse_result::exit);
  REQUIRE_THAT(ss.str(), !Contains("--foo"));
}

TEST_CASE("--help must not end in dot", "[argparse::parser]")
{
  argparse::parser p;
  auto& arg = p.add("--foo");

  REQUIRE_THROWS_AS(arg.help("invalid help text."), assert_failed);
}

TEST_CASE("unexpected positional argument", "[argparse::parser]")
{
  argparse::parser p;
  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "foo" }) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(), "[ERROR] Unexpected positional argument 'foo'\n");
}

TEST_CASE("unexpected argument", "[argparse::parser]")
{
  argparse::parser p;
  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "--foo" }) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(), "[ERROR] Unknown argument '--foo'\n");
}

TEST_CASE("typo in argument", "[argparse::parser]")
{
  argparse::parser p;
  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "--halp" }) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--halp'. "
                      "Did you mean --help?\n");
}

TEST_CASE("partial argument", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--partofalongargument");

  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "--part" }) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--part'. "
                      "Did you mean --partofalongargument?\n");
}

TEST_CASE("two possible arguments", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--arg1");
  p.add("--arg2");

  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "--arg" }) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--arg'. "
                      "Did you mean --arg1, or --arg2?\n");
}

TEST_CASE("three possible arguments", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--arg1");
  p.add("--arg2");
  p.add("--arg3");

  log::log_capture ss;

  REQUIRE(p.parse_args({ "exe", "--arg" }) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--arg'. "
                      "Did you mean --arg1, --arg2, or --arg3?\n");
}

TEST_CASE("parse multiple arguments", "[argparse::parser]")
{
  uint32_t sink = 0;
  argparse::parser p;
  p.add("--arg1").bind_u32(&sink);
  p.add("--arg2");
  p.add("--arg3");

  REQUIRE(p.parse_args({ "exe", "--arg1", "1234", "--arg3" }) ==
          argparse::parse_result::ok);
  REQUIRE(p.is_set("--arg1"));
  REQUIRE_FALSE(p.is_set("--arg2"));
  REQUIRE(p.is_set("--arg3"));
  REQUIRE(sink == 1234);
}

TEST_CASE("value returns value as str", "[argparse::parser]")
{
  uint32_t usink = 0;
  std::string ssink;
  argparse::parser p;
  p.add("--arg1").bind_u32(&usink);
  p.add("--arg2").bind_str(&ssink);

  REQUIRE(p.parse_args({ "exe", "--arg1", "1234", "--arg2", "foo bar*" }) ==
          argparse::parse_result::ok);
  REQUIRE(p.value("--arg1") == "1234");
  REQUIRE(p.value("--arg2") == "foo bar*");
}

TEST_CASE("user supplied argument", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.set_preamble("basic help");
  p.add("--test");

  p.print_help();

  REQUIRE(ss.str() == std::string(HELP_HEADER) + "   --test\n");
}

TEST_CASE("user supplied argument with meta-var", "[argparse::parser]")
{
  uint32_t sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.set_preamble("basic help");
  p.add("--test", "META").bind_u32(&sink);

  p.print_help();

  REQUIRE(ss.str() == "My App v1234\n\n"
                      "basic help\n\n"
                      "OPTIONS:\n"
                      "   -h, --help      Display this message\n"
                      "   -v, --version   Print the version string\n"
                      "   --licenses      Print licenses for this software\n\n"
                      "   --test <META>\n");
}

TEST_CASE("user supplied argument with meta-var and help", "[argparse::parser]")
{
  uint32_t sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.set_preamble("basic help");
  p.add("--test", "META")
    .help("A long help message that exceeds the limit of 60 characters by some "
          "amount in order to test the line break functionality")
    .bind_u32(&sink);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() ==
          "My App v1234\n\n"
          "basic help\n\n"
          "OPTIONS:\n"
          "   -h, --help      Display this message\n"
          "   -v, --version   Print the version string\n"
          "   --licenses      Print licenses for this software\n\n"
          "   --test <META>   A long help message that exceeds the\n"
          "                   limit of 60 characters by some amount in\n"
          "                   order to test the line break\n"
          "                   functionality\n");
}

TEST_CASE("help with default value", "[argparse::parser]")
{
  uint32_t sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.set_preamble("basic help");
  p.add("--test", "META")
    .help("A long help message that exceeds the limit of 60 characters by some "
          "amount in order to test the line break functionality")
    .bind_u32(&sink)
    .with_default(1234);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() ==
          "My App v1234\n\n"
          "basic help\n\n"
          "OPTIONS:\n"
          "   -h, --help      Display this message\n"
          "   -v, --version   Print the version string\n"
          "   --licenses      Print licenses for this software\n\n"
          "   --test <META>   A long help message that exceeds the\n"
          "                   limit of 60 characters by some amount in\n"
          "                   order to test the line break\n"
          "                   functionality [default: 1234]\n");
}

TEST_CASE("required option missing", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").depends_on("--bar");
  p.add("--bar");

  REQUIRE(p.parse_args({ "exe", "--foo" }) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(),
                  "[ERROR] Option --bar is required when using option --foo\n");
}

TEST_CASE("required option supplied", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--foo").depends_on("--bar");
  p.add("--bar");

  REQUIRE(p.parse_args({ "exe", "--foo", "--bar" }) ==
          argparse::parse_result::ok);
}

TEST_CASE("conflicting option missing", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--foo").conflicts_with("--bar");
  p.add("--bar");

  REQUIRE(p.parse_args({ "exe", "--foo" }) == argparse::parse_result::ok);
}

TEST_CASE("conflicting option supplied", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").conflicts_with("--bar");
  p.add("--bar");

  REQUIRE(p.parse_args({ "exe", "--foo", "--bar" }) ==
          argparse::parse_result::error);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Option --bar cannot be used together with option --foo\n");
}

TEST_CASE("conflicting option does not exist", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").conflicts_with("--bar");

  REQUIRE_THROWS_AS(p.parse_args({ "exe", "--foo" }), assert_failed);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Option --foo conflicts with unknown command-line option --bar\n");
}

TEST_CASE("missing value", "[argparse::parser]")
{
  uint32_t sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").bind_u32(&sink);

  REQUIRE(p.parse_args({ "exe", "--foo" }) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(),
                  "[ERROR] Command-line argument --foo takes 1 value, but 0 "
                  "values were provided!\n");
}

TEST_CASE("missing values, singular", "[argparse::parser]")
{
  string_vec sink;
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").bind_vec(&sink).with_min_values(1).with_max_values(2);

  REQUIRE(p.parse_args({ "exe", "--foo" }) == argparse::parse_result::error);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Command-line argument --foo takes at least 1 value, but 0 "
    "values were provided!\n");
}

TEST_CASE("missing values, plural", "[argparse::parser]")
{
  string_vec sink;
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").bind_vec(&sink).with_min_values(2).with_max_values(3);

  REQUIRE(p.parse_args({ "exe", "--foo", "bar" }) ==
          argparse::parse_result::error);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Command-line argument --foo takes at least 2 values, but 1 "
    "value was provided!\n");
}

TEST_CASE("excessive values, takes zero", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo");

  SECTION("singular")
  {
    REQUIRE(p.parse_args({ "exe", "--foo", "foo" }) ==
            argparse::parse_result::error);
    REQUIRE_POSTFIX(ss.str(),
                    "[ERROR] Command-line argument --foo takes 0 values, but 1 "
                    "value was provided!\n");
  }

  SECTION("plural")
  {
    REQUIRE(p.parse_args({ "exe", "--foo", "foo", "bar" }) ==
            argparse::parse_result::error);
    REQUIRE_POSTFIX(ss.str(),
                    "[ERROR] Command-line argument --foo takes 0 values, but 2 "
                    "values were provided!\n");
  }
}

TEST_CASE("excessive values, takes non-zero", "[argparse::parser]")
{
  string_vec vec;
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").bind_vec(&vec).with_max_values(2);

  REQUIRE(p.parse_args({ "exe", "--foo", "foo", "bar", "zod" }) ==
          argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(),
                  "[ERROR] Command-line argument --foo takes at most 2 values, "
                  "but 3 values were provided!\n");
}

TEST_CASE("invalid value", "[argparse::parser]")
{
  uint32_t sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").bind_u32(&sink);

  REQUIRE(p.parse_args({ "exe", "--foo", "one" }) ==
          argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(),
                  "[ERROR] Invalid command-line argument --foo one: value is "
                  "not a valid number\n");
}

TEST_CASE("invalid choice", "[argparse::parser]")
{
  std::string sink;
  log::log_capture cap;
  argparse::parser p;
  p.add("--foo").bind_str(&sink).with_choices({ "abc", "def", "ghi" });

  REQUIRE(p.parse_args({ "exe", "--foo", "one" }) ==
          argparse::parse_result::error);
  REQUIRE_POSTFIX(cap.str(),
                  "[ERROR] Invalid command-line argument --foo one: Valid "
                  "values for --foo are abc, def, and ghi\n");
}

TEST_CASE("help with finite max number of values", "[argparse::parser]")
{
  string_vec sink;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.add("--test", "X")
    .help("Help text")
    .bind_vec(&sink)
    .with_min_values(0)
    .with_max_values(2);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() == "My App v1234\n\n"
                      "OPTIONS:\n"
                      "   -h, --help       Display this message\n"
                      "   -v, --version    Print the version string\n"
                      "   --licenses       Print licenses for this software\n\n"
                      "   --test [X] [X]   Help text\n");
}

TEST_CASE("help with infinite number of values", "[argparse::parser]")
{
  string_vec sink;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.add("--test", "X").help("Help text").bind_vec(&sink).with_min_values(0);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() ==
          "My App v1234\n\n"
          "OPTIONS:\n"
          "   -h, --help        Display this message\n"
          "   -v, --version     Print the version string\n"
          "   --licenses        Print licenses for this software\n\n"
          "   --test [X, ...]   Help text\n");
}

TEST_CASE("help with lower bound of values", "[argparse::parser]")
{
  string_vec sink;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.add("--test", "X").help("Help text").bind_vec(&sink).with_min_values(1);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() ==
          "My App v1234\n\n"
          "OPTIONS:\n"
          "   -h, --help            Display this message\n"
          "   -v, --version         Print the version string\n"
          "   --licenses            Print licenses for this software\n\n"
          "   --test <X> [X, ...]   Help text\n");
}

TEST_CASE("help with lower and upper bound of values", "[argparse::parser]")
{
  string_vec sink;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("1234");
  p.add("--test", "X")
    .help("Help text")
    .bind_vec(&sink)
    .with_min_values(2)
    .with_max_values(3);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() ==
          "My App v1234\n\n"
          "OPTIONS:\n"
          "   -h, --help           Display this message\n"
          "   -v, --version        Print the version string\n"
          "   --licenses           Print licenses for this software\n\n"
          "   --test <X> <X> [X]   Help text\n");
}

TEST_CASE("parse_result to string", "[argparse]")
{
  using Catch::fallbackStringifier;

  CHECK(fallbackStringifier(argparse::parse_result::exit) ==
        "parse_result::exit");
  CHECK(fallbackStringifier(argparse::parse_result::error) ==
        "parse_result::error");
  CHECK(fallbackStringifier(argparse::parse_result::ok) == "parse_result::ok");

  CHECK_THROWS_AS(fallbackStringifier(static_cast<argparse::parse_result>(-1)),
                  assert_failed);
}

} // namespace adapterremoval
