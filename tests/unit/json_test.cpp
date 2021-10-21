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
#include <sstream>

#include "counts.hpp"
#include "debug.hpp"
#include "json.hpp"
#include "testing.hpp"

TEST_CASE("empty_json", "[json::json]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
  }

  REQUIRE(ss.str() == "{}");
}

TEST_CASE("single_value", "[json::json]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_int("test", 17);
  }

  REQUIRE(ss.str() == "{\n  \"test\": 17\n}");
}

TEST_CASE("write", "[json::write]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write("value", "a\rb\"c");
  }

  REQUIRE(ss.str() == "{\n  \"value\": \"a\\rb\\\"c\"\n}");
}

TEST_CASE("write_vector([])", "[json::write_vector]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    std::vector<std::string> vector;
    json.write("value", vector);
  }

  REQUIRE(ss.str() == "{\n  \"value\": []\n}");
}

TEST_CASE("write_vector([1])", "[json::write_vector]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    std::vector<std::string> vector;
    vector.push_back("foo");
    json.write("value", vector);
  }

  REQUIRE(ss.str() == "{\n  \"value\": [\"foo\"]\n}");
}

TEST_CASE("write_vector([2])", "[json::write_vector]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    std::vector<std::string> vector;
    vector.push_back("foo");
    vector.push_back("bar");
    json.write("value", vector);
  }

  REQUIRE(ss.str() == "{\n  \"value\": [\"foo\", \"bar\"]\n}");
}

TEST_CASE("write_vector([3])", "[json::write_vector]")
{
  rates r;
  r.inc(0, 0);
  r.inc(1, std::nan(""));
  r.inc(2, 0.5);

  std::stringstream ss;
  {
    json_writer json(ss);
    json.write("value", r);
  }

  REQUIRE(ss.str() == "{\n  \"value\": [0.000, null, 0.500]\n}");
}

TEST_CASE("write_int", "[json::write_int]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_int("value", -12345);
  }

  REQUIRE(ss.str() == "{\n  \"value\": -12345\n}");
}

TEST_CASE("write_int_vector", "[json::write_int_vector]")
{
  std::stringstream ss;
  {
    json_writer json(ss);

    counts values;
    values.inc(0, 7);
    values.inc(1, 9);
    values.inc(2, 13);

    json.write("value", values);
  }

  REQUIRE(ss.str() == "{\n  \"value\": [7, 9, 13]\n}");
}

TEST_CASE("write_float", "[json::write_float]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_float("value", 1.3);
  }

  REQUIRE(ss.str() == "{\n  \"value\": 1.300\n}");
}

TEST_CASE("write_nan", "[json::write_nan]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_float("value", std::nan(""));
  }

  REQUIRE(ss.str() == "{\n  \"value\": null\n}");
}

TEST_CASE("write_bool(true)", "[json::write_bool]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_bool("value", true);
  }

  REQUIRE(ss.str() == "{\n  \"value\": true\n}");
}

TEST_CASE("write_bool(false)", "[json::write_bool]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_bool("value", false);
  }

  REQUIRE(ss.str() == "{\n  \"value\": false\n}");
}

TEST_CASE("write_null", "[json::write_null]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_null("value");
  }

  REQUIRE(ss.str() == "{\n  \"value\": null\n}");
}

TEST_CASE("multiple_values", "[json::json]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write_bool("value1", true);
    json.write_float("value2", 13.7);
  }

  REQUIRE(ss.str() == "{\n  \"value1\": true,\n  \"value2\": 13.700\n}");
}

TEST_CASE("empty_nested_dict", "[json::json]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.start("foo");
  }

  REQUIRE(ss.str() == "{\n  \"foo\": {}\n}");
}

TEST_CASE("nested_dict", "[json::json]")
{
  std::stringstream ss;
  {
    json_writer json(ss);
    json.write("value1", "foo");

    if (auto _ = json.start("value2")) {
      json.write_int("bar", 1234);
    }
  }

  REQUIRE(
    ss.str() ==
    "{\n  \"value1\": \"foo\",\n  \"value2\": {\n    \"bar\": 1234\n  }\n}");
}
