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
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "argparse.hpp"
#include "debug.hpp"
#include "logging.hpp" // for log
#include "testing.hpp"

namespace adapterremoval {

using argparse::argument;
using argparse::argument_ptr;

const size_t parsing_failed = static_cast<size_t>(-1);

using Catch::Matchers::Contains;

#define REQUIRE_POSTFIX(a, b) REQUIRE_THAT((a), Catch::Matchers::EndsWith(b))
#define REQUIRE_CONTAINS(a, b) REQUIRE_THAT((a), Catch::Matchers::Contains(b))

size_t
consume(argparse::sink& dst, const string_vec& values)
{
  return dst.consume(values.begin(), values.end());
}

///////////////////////////////////////////////////////////////////////////////
// boolean sink

TEST_CASE("bool sink is required", "[argparse::bool_sink]")
{
  REQUIRE_THROWS_AS(argparse::bool_sink(nullptr), assert_failed);
}

TEST_CASE("bool sink is initialized", "[argparse::bool_sink]")
{
  bool value = true;
  argparse::bool_sink sink(&value);
  REQUIRE_FALSE(value);
}

TEST_CASE("bool sink to_str", "[argparse::bool_sink]")
{
  bool value = false;
  argparse::bool_sink sink(&value);
  REQUIRE(sink.to_str() == "off");
  value = true;
  REQUIRE(sink.to_str() == "on");
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
  REQUIRE(sink.to_str() == "on");
  REQUIRE(value);
}

TEST_CASE("bool sink rejects values", "[argparse::bool_sink]")
{
  string_vec values{ "0" };
  bool value = false;
  argparse::bool_sink sink(&value);

  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

///////////////////////////////////////////////////////////////////////////////
// uint sink

TEST_CASE("uint sink is required", "[argparse::uint_sink]")
{
  REQUIRE_THROWS_AS(argparse::uint_sink(nullptr), assert_failed);
}

TEST_CASE("uint sink is initialized", "[argparse::uint_sink]")
{
  unsigned value = 12345;
  argparse::uint_sink sink(&value);
  REQUIRE(value == 0);
}

TEST_CASE("uint sink to_str", "[argparse::uint_sink]")
{
  unsigned value = 1234567;
  argparse::uint_sink sink(&value);
  REQUIRE(sink.to_str() == "0");
  value = 1234567;
  REQUIRE(sink.to_str() == "1234567");
}

TEST_CASE("uint sink with_default", "[argparse::uint_sink]")
{
  unsigned value = 0;
  argparse::uint_sink sink(&value);

  REQUIRE_FALSE(sink.has_default());
  sink.with_default(1234567);

  REQUIRE(sink.has_default());
  REQUIRE(value == 1234567);
}

TEST_CASE("uint sink requires value", "[argparse::uint_sink]")
{
  unsigned value = 0;
  argparse::uint_sink sink(&value);

  string_vec values;
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("uint sink consumes single value", "[argparse::uint_sink]")
{
  string_vec values{ "1234567" };

  unsigned value = 0;
  argparse::uint_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == 1234567);
}

TEST_CASE("uint sink consumes single value only", "[argparse::uint_sink]")
{
  unsigned value = 0;
  argparse::uint_sink sink(&value);

  const string_vec values{ "1234567", "8901234" };
  REQUIRE_THROWS_AS(sink.consume(values.begin(), values.end()), assert_failed);
}

TEST_CASE("uint requires valid integer value", "[argparse::uint_sink]")
{
  string_vec values{ "abc" };

  unsigned value = 0;
  argparse::uint_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == parsing_failed);
  REQUIRE(value == 0);
}

TEST_CASE("uint requires positive integer value", "[argparse::uint_sink]")
{
  unsigned value = 0;
  argparse::uint_sink sink(&value);

  string_vec values{ "-123" };
  REQUIRE(sink.consume(values.begin(), values.end()) == parsing_failed);
  REQUIRE(value == 0);
}

TEST_CASE("uint accepts unsigned lower bound", "[argparse::uint_sink]")
{
  unsigned value = 0;
  argparse::uint_sink sink(&value);
  sink.with_default(123456);

  string_vec values{ "0" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == 0);
}

TEST_CASE("uint accepts unsigned upper bound", "[argparse::uint_sink]")
{
  unsigned value = 0;
  argparse::uint_sink sink(&value);

  string_vec values{ "4294967295" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == 4294967295);
}

TEST_CASE("uint rejects past unsigned upper bound", "[argparse::uint_sink]")
{
  unsigned value = 0;
  argparse::uint_sink sink(&value);

  string_vec values{ "4294967296" };
  REQUIRE(sink.consume(values.begin(), values.end()) == parsing_failed);
  REQUIRE(value == 0);
}

TEST_CASE("uint disallows trailing garbage", "[argparse::uint_sink]")
{
  string_vec values{ "123abc" };

  unsigned value = 0;
  argparse::uint_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == parsing_failed);
  REQUIRE(value == 0);
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

TEST_CASE("double sink to_str", "[argparse::double_sink]")
{
  double value = 12345.67;
  argparse::double_sink sink(&value);
  REQUIRE(sink.to_str() == "0");
  value = 12345.67;
  REQUIRE(sink.to_str() == "12345.67");
  value = 1234567;
  REQUIRE(sink.to_str() == "1234567");
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
  sink.with_default(12345.67);
  REQUIRE(sink.has_default());
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

  REQUIRE(sink.consume(values.begin(), values.end()) == parsing_failed);
  REQUIRE(value == 0);
}

TEST_CASE("double requires valid double #2", "[argparse::double_sink]")
{
  string_vec values{ "123.0abc" };

  double value = 0;
  argparse::double_sink sink(&value);

  REQUIRE(sink.consume(values.begin(), values.end()) == parsing_failed);
  REQUIRE(value == 0);
}

///////////////////////////////////////////////////////////////////////////////
// str sink

TEST_CASE("str sink is required", "[argparse::str_sink]")
{
  REQUIRE_THROWS_AS(argparse::str_sink(nullptr), assert_failed);
}

TEST_CASE("str sink is initialized", "[argparse::str_sink]")
{
  std::string value = "foo";
  argparse::str_sink sink(&value);
  REQUIRE(value.empty());
}

TEST_CASE("str sink to_str", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  value = "foobar";
  REQUIRE(sink.to_str() == "foobar");
}

TEST_CASE("str sink to_str escapes", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  value = "foo bar";
  REQUIRE(sink.to_str() == "foo bar");
}

TEST_CASE("str sink with_default (char*)", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_default("foobar");
  REQUIRE(value == "foobar");
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
  REQUIRE(value == "");
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

  string_vec values{ "ghi" };
  REQUIRE(sink.consume(values.begin(), values.end()) == 1);
  REQUIRE(value == "ghi");
}

TEST_CASE("str sink rejectes values not in choices", "[argparse::str_sink]")
{
  std::string value;
  argparse::str_sink sink(&value);
  sink.with_choices({ "abc", "def", "ghi" });

  string_vec values{ "foo" };
  REQUIRE(sink.consume(values.begin(), values.end()) == parsing_failed);
  REQUIRE(value == "");
}

TEST_CASE("str sink case-insenstive, returns choice", "[argparse::str_sink]")
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
  REQUIRE(sink.to_str() == "foobar");
}

TEST_CASE("vec sink to_vec #2", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);
  value.push_back("foo");
  value.push_back("bar");
  REQUIRE(sink.to_str() == "foo;bar");
}

TEST_CASE("vec sink to_vec escapes", "[argparse::vec_sink]")
{
  string_vec value;
  argparse::vec_sink sink(&value);
  value.push_back("foo");
  value.push_back("1 2");
  value.push_back("bar");
  REQUIRE(sink.to_str() == "foo;'1 2';bar");
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
  REQUIRE(value.front() == "");
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
  REQUIRE(arg.help() == "");

  arg.help("pointless gibberish");
  REQUIRE(arg.help() == "pointless gibberish");
}

TEST_CASE("help with default", "[argparse::argument]")
{
  unsigned value = 0;
  argparse::argument arg("--12345", "67890");
  arg.help("pointless gibberish").bind_uint(&value).with_default(7913);

  REQUIRE(arg.help() == "pointless gibberish [default: 7913]");
}

TEST_CASE("help with custom default", "[argparse::argument]")
{
  unsigned value = 0;
  argparse::argument arg("--12345", "67890");
  arg.help("pointless gibberish [default: foo]")
    .bind_uint(&value)
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
  arg.deprecated();
  REQUIRE(arg.is_deprecated());

  log::log_capture ss;

  string_vec values{ "--12345" };
  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE_POSTFIX(ss.str(),
                  "[WARNING] Option --12345 is deprecated and will be "
                  "removed in the future.\n");
}

TEST_CASE("default argument sink", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");

  REQUIRE(arg.to_str() == "off");
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
  arg.abbreviation('1');

  REQUIRE(arg.key() == "--12345");
  REQUIRE(arg.short_key() == "-1");
  REQUIRE(arg.keys() == string_vec{ "--12345", "-1" });
  REQUIRE_FALSE(arg.is_deprecated_alias("--12345"));
  REQUIRE_FALSE(arg.is_deprecated_alias("-1"));
}

TEST_CASE("deprecated argument alias", "[argparse::argument]")
{
  argparse::argument arg("--12345", "67890");
  arg.deprecated_alias("--foo");

  REQUIRE(arg.key() == "--12345");
  REQUIRE(arg.short_key() == "");
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

  unsigned sink = 0;
  arg.bind_uint(&sink);

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
  REQUIRE(ss.str() == "");
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
  REQUIRE(ss.str() == "");
}

TEST_CASE("warning on deprecated alias", "[argparse::argument]")
{
  argparse::argument arg("--12345");
  arg.deprecated_alias("--foo");
  string_vec values{ "--foo" };

  log::log_capture ss;

  REQUIRE(arg.parse(values.begin(), values.end()) == 1);
  REQUIRE_POSTFIX(ss.str(),
                  "[WARNING] Option --foo is deprecated and will be removed "
                  "in the future. Please use --12345 instead.\n");
}

///////////////////////////////////////////////////////////////////////////////
// parser

const std::string HELP_HEADER =
  "My App v1234\n\n"
  "basic help\n"
  "OPTIONS:\n"
  "   -h, --help      Display this message.\n"
  "   -v, --version   Print the version string.\n"
  "   --licenses      Print licences for this software.\n\n";

TEST_CASE("--version", "[argparse::parser]")
{
  argparse::parser p;
  p.set_name("My App");
  p.set_version("v1234");

  const auto arg = GENERATE("-v", "--version");
  const char* args[] = { "exe", arg };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::exit);
  REQUIRE(ss.str() == "My App v1234\n");
}

TEST_CASE("--help", "[argparse::parser]")
{
  argparse::parser p;
  p.set_name("My App");
  p.set_version("v1234");
  p.set_preamble("basic help");
  const auto arg = GENERATE("-h", "--help");
  const char* args[] = { "exe", arg };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::exit);
  REQUIRE(ss.str() == HELP_HEADER);
}

TEST_CASE("--help with deprecated argument", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--foo").deprecated();
  const char* args[] = { "exe", "--help" };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::exit);
  REQUIRE_THAT(ss.str(), !Contains("--foo"));
}

TEST_CASE("unexpected positional argument", "[argparse::parser]")
{
  argparse::parser p;
  const char* args[] = { "exe", "foo" };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(), "[ERROR] Unexpected positional argument 'foo'\n");
}

TEST_CASE("unexpected argument", "[argparse::parser]")
{
  argparse::parser p;
  const char* args[] = { "exe", "--foo" };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(), "[ERROR] Unknown argument '--foo'\n");
}

TEST_CASE("typo in argument", "[argparse::parser]")
{
  argparse::parser p;
  const char* args[] = { "exe", "--halp" };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--halp'. "
                      "Did you mean --help?\n");
}

TEST_CASE("partial argument", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--partofalongargument");
  const char* args[] = { "exe", "--part" };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--part'. "
                      "Did you mean --partofalongargument?\n");
}

TEST_CASE("two possible arguments", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--arg1");
  p.add("--arg2");
  const char* args[] = { "exe", "--arg" };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--arg'. "
                      "Did you mean --arg1 or --arg2?\n");
}

TEST_CASE("three possible arguments", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--arg1");
  p.add("--arg2");
  p.add("--arg3");
  const char* args[] = { "exe", "--arg" };

  log::log_capture ss;

  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE(ss.str() == "[ERROR] Unknown argument '--arg'. "
                      "Did you mean --arg1, --arg2 or --arg3?\n");
}

TEST_CASE("parse multiple arguments", "[argparse::parser]")
{
  unsigned sink = 0;
  argparse::parser p;
  p.add("--arg1").bind_uint(&sink);
  p.add("--arg2");
  p.add("--arg3");
  const char* args[] = { "exe", "--arg1", "1234", "--arg3" };

  REQUIRE(p.parse_args(4, args) == argparse::parse_result::ok);
  REQUIRE(p.is_set("--arg1"));
  REQUIRE_FALSE(p.is_set("--arg2"));
  REQUIRE(p.is_set("--arg3"));
  REQUIRE(sink == 1234);
}

TEST_CASE("to_str returns value as str", "[argparse::parser]")
{
  unsigned usink = 0;
  std::string ssink;
  argparse::parser p;
  p.add("--arg1").bind_uint(&usink);
  p.add("--arg2").bind_str(&ssink);
  const char* args[] = { "exe", "--arg1", "1234", "--arg2", "foo bar*" };

  REQUIRE(p.parse_args(5, args) == argparse::parse_result::ok);
  REQUIRE(p.to_str("--arg1") == "1234");
  REQUIRE(p.to_str("--arg2") == "foo bar*");
}

TEST_CASE("user supplied argument", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("v1234");
  p.set_preamble("basic help");
  p.add("--test");

  p.print_help();

  REQUIRE(ss.str() == HELP_HEADER + "   --test\n");
}

TEST_CASE("user supplied argument with meta-var", "[argparse::parser]")
{
  unsigned sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("v1234");
  p.set_preamble("basic help");
  p.add("--test", "META").bind_uint(&sink);

  p.print_help();

  REQUIRE(ss.str() == "My App v1234\n\n"
                      "basic help\n"
                      "OPTIONS:\n"
                      "   -h, --help      Display this message.\n"
                      "   -v, --version   Print the version string.\n"
                      "   --licenses      Print licences for this software.\n\n"
                      "   --test <META>\n");
}

TEST_CASE("user supplied argument with meta-var and help", "[argparse::parser]")
{
  unsigned sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("v1234");
  p.set_preamble("basic help");
  p.add("--test", "META")
    .help("A long help message that exceeds the limit of 60 characters by some "
          "amount in order to test the line break functionality")
    .bind_uint(&sink);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() ==
          "My App v1234\n\n"
          "basic help\n"
          "OPTIONS:\n"
          "   -h, --help      Display this message.\n"
          "   -v, --version   Print the version string.\n"
          "   --licenses      Print licences for this software.\n\n"
          "   --test <META>   A long help message that exceeds the\n"
          "                   limit of 60 characters by some amount in\n"
          "                   order to test the line break\n"
          "                   functionality\n");
}

TEST_CASE("help with default value", "[argparse::parser]")
{
  unsigned sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.set_name("My App");
  p.set_version("v1234");
  p.set_preamble("basic help");
  p.add("--test", "META")
    .help("A long help message that exceeds the limit of 60 characters by some "
          "amount in order to test the line break functionality")
    .bind_uint(&sink)
    .with_default(1234);

  p.set_terminal_width(60);
  p.print_help();

  REQUIRE(ss.str() ==
          "My App v1234\n\n"
          "basic help\n"
          "OPTIONS:\n"
          "   -h, --help      Display this message.\n"
          "   -v, --version   Print the version string.\n"
          "   --licenses      Print licences for this software.\n\n"
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

  const char* args[] = { "exe", "--foo" };
  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(),
                  "[ERROR] Option --bar is required when using option --foo\n");
}

TEST_CASE("required option supplied", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--foo").depends_on("--bar");
  p.add("--bar");

  const char* args[] = { "exe", "--foo", "--bar" };
  REQUIRE(p.parse_args(3, args) == argparse::parse_result::ok);
}

TEST_CASE("conflicting option missing", "[argparse::parser]")
{
  argparse::parser p;
  p.add("--foo").conflicts_with("--bar");
  p.add("--bar");

  const char* args[] = { "exe", "--foo" };
  REQUIRE(p.parse_args(2, args) == argparse::parse_result::ok);
}

TEST_CASE("conflicting option supplied", "[argparse::parser]")
{
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").conflicts_with("--bar");
  p.add("--bar");

  const char* args[] = { "exe", "--foo", "--bar" };
  REQUIRE(p.parse_args(3, args) == argparse::parse_result::error);
  REQUIRE_POSTFIX(
    ss.str(),
    "[ERROR] Option --bar cannot be used together with option --foo\n");
}

TEST_CASE("missing value", "[argparse::parser]")
{
  unsigned sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").bind_uint(&sink);

  const char* args[] = { "exe", "--foo" };
  REQUIRE(p.parse_args(2, args) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(),
                  "[ERROR] Command-line argument --foo takes 1 value, but 0 "
                  "values were provided!\n");
}

TEST_CASE("invalid value", "[argparse::parser]")
{
  unsigned sink = 0;
  log::log_capture ss;
  argparse::parser p;
  p.add("--foo").bind_uint(&sink);

  const char* args[] = { "exe", "--foo", "one" };
  REQUIRE(p.parse_args(3, args) == argparse::parse_result::error);
  REQUIRE_POSTFIX(ss.str(), "[ERROR] Invalid value for --foo: one\n");
}

} // namespace adapterremoval
