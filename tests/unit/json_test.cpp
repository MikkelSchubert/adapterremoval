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

template<typename T>
std::string
_write_json(const T& json)
{
  std::stringstream ss;
  json.write(ss);

  return ss.str();
}

////////////////////////////////////////////////////////////////////////////////
// json_token::constructor

TEST_CASE("empty json string to string")
{
  const auto s = "";
  const json_token t(s);

  REQUIRE(t.to_string() == s);
}

TEST_CASE("simple json string")
{
  const auto s = "simple test 123";
  const json_token t(s);

  REQUIRE(t.to_string() == s);
}

TEST_CASE("special characters are not encoded in json token")
{
  const auto s = "\\simple\ttest\n123";
  const json_token t(s);

  REQUIRE(t.to_string() == s);
}

TEST_CASE("json string to_string and write produces identical results")
{
  const auto s = "\\simple\ttest\n123";
  const json_token t(s);

  REQUIRE(t.to_string() == _write_json(t));
}

////////////////////////////////////////////////////////////////////////////////
// json_token::from_str

TEST_CASE("empty json_token::from_str")
{
  auto s = json_token::from_str("");

  REQUIRE(s->to_string() == "\"\"");
}

TEST_CASE("simple json_token::from_str")
{
  auto s = json_token::from_str("simple test 123");

  REQUIRE(s->to_string() == "\"simple test 123\"");
}

TEST_CASE("special characters are encoded by json_token::from_str")
{
  auto s = json_token::from_str("\\simple\ttest\n123");

  REQUIRE(s->to_string() == "\"\\\\simple\\ttest\\n123\"");
}

TEST_CASE("identical to_string and write for json_token::from_str")
{
  auto s = json_token::from_str("\\simple\ttest\n123");
  REQUIRE(s->to_string() == _write_json(*s));
}

////////////////////////////////////////////////////////////////////////////////
// json_token::from_i64

TEST_CASE("json_token::from_i64 for 0")
{
  REQUIRE(json_token::from_i64(0)->to_string() == "0");
}

TEST_CASE("json_token::from_i64 for postive")
{
  REQUIRE(json_token::from_i64(1324)->to_string() == "1324");
}

TEST_CASE("json_token::from_i64 for negative")
{
  REQUIRE(json_token::from_i64(-8671)->to_string() == "-8671");
}

TEST_CASE("identical to_string and write for json_token::from_i64")
{
  auto s = json_token::from_i64(-83754);
  REQUIRE(s->to_string() == _write_json(*s));
}

////////////////////////////////////////////////////////////////////////////////
// json_token::from_double

TEST_CASE("json_token::from_f64 increases digits to 3")
{
  auto s = json_token::from_f64(6);

  REQUIRE(s->to_string() == "6.000");
}

TEST_CASE("json_token::from_f64 rounds to 3 digits")
{
  auto s = json_token::from_f64(6.4568);

  REQUIRE(s->to_string() == "6.457");
}

TEST_CASE("json_token::from_f64 writes NaN as null")
{
  auto s = json_token::from_f64(std::nan(""));

  REQUIRE(s->to_string() == "null");
}

TEST_CASE("identical to_string and write for json_token::from_double")
{
  auto s = json_token::from_f64(-6.4568);
  REQUIRE(s->to_string() == _write_json(*s));
}

////////////////////////////////////////////////////////////////////////////////
// json_token::from_boolean

TEST_CASE("json_token::from_boolean is true")
{
  auto s = json_token::from_boolean(true);

  REQUIRE(s->to_string() == "true");
}

TEST_CASE("json_token::from_boolean is false")
{
  auto s = json_token::from_boolean(false);

  REQUIRE(s->to_string() == "false");
}

TEST_CASE("identical to_string and write for json_token::from_boolean")
{
  auto s = json_token::from_boolean(false);
  REQUIRE(s->to_string() == _write_json(*s));
}

////////////////////////////////////////////////////////////////////////////////
// json_token::from_null

TEST_CASE("json_token::from_null is null")
{
  auto s = json_token::from_null();

  REQUIRE(s->to_string() == "null");
}

TEST_CASE("identical to_string and write for json_token::from_null")
{
  auto s = json_token::from_null();
  REQUIRE(s->to_string() == _write_json(*s));
}

////////////////////////////////////////////////////////////////////////////////
// json_dict::from_i64_vec

TEST_CASE("json_token::from_i64_vec is null")
{
  counts c(3);
  c.inc(0, 3);
  c.inc(1, 2);
  c.inc(2, 1);
  auto s = json_token::from_i64_vec(c);

  REQUIRE(s->to_string() == "[3, 2, 1]");
}

TEST_CASE("identical to_string and write for json_token::from_i64_vec")
{
  counts c(3);
  c.inc(0, 3);
  c.inc(1, 2);
  c.inc(2, 1);

  auto s = json_token::from_i64_vec(c);
  REQUIRE(s->to_string() == _write_json(*s));
}

////////////////////////////////////////////////////////////////////////////////
// json_dict::from_f64_vec

TEST_CASE("json_token::from_f64_vec")
{
  rates c(3);
  c.inc(0, 1);
  c.inc(1, 2.71);
  c.inc(2, 3.1415);
  auto s = json_token::from_f64_vec(c);

  REQUIRE(s->to_string() == "[1.000, 2.710, 3.142]");
}

TEST_CASE("identical to_string and write for json_token::from_f64_vec")
{
  rates c(3);
  c.inc(0, 1);
  c.inc(1, 2.71);
  c.inc(2, 3.1415);

  auto s = json_token::from_f64_vec(c);
  REQUIRE(s->to_string() == _write_json(*s));
}

TEST_CASE("json_token::from_f64_vec writes NaN as null")
{
  rates c(3);
  c.inc(0, 1);
  c.inc(1, std::nan(""));
  c.inc(2, 3.1415);
  auto s = json_token::from_f64_vec(c);

  REQUIRE(s->to_string() == "[1.000, null, 3.142]");
}

////////////////////////////////////////////////////////////////////////////////
// json_dict::from_str_vec

TEST_CASE("json_token::from_str_vec")
{
  string_vec v{ "foo", "b\tr" };
  auto s = json_token::from_str_vec(v);

  REQUIRE(s->to_string() == R"(["foo", "b\tr"])");
}

TEST_CASE("identical to_string and write for json_token::from_str_vec")
{
  string_vec v{ "foo", "b\tr" };

  auto s = json_token::from_str_vec(v);
  REQUIRE(s->to_string() == _write_json(*s));
}

////////////////////////////////////////////////////////////////////////////////
// json_dict

TEST_CASE("empty json_dict is single line")
{
  json_dict s;

  REQUIRE(s.to_string() == "{}");
}

TEST_CASE("json_dict with int")
{
  json_dict s;
  s.i64("foo", -1235);

  REQUIRE(s.to_string() == "{\n  \"foo\": -1235\n}");
}

TEST_CASE("json_dict with float")
{
  json_dict s;
  s.f64("foo", 1234.4567);

  REQUIRE(s.to_string() == R"({
  "foo": 1234.457
})");
}

TEST_CASE("json_dict with string")
{
  json_dict s;
  s.str("foo", "\"bar\"");

  REQUIRE(s.to_string() == R"({
  "foo": "\"bar\""
})");
}

TEST_CASE("json_dict with bool")
{
  json_dict s;
  s.boolean("foo", false);

  REQUIRE(s.to_string() == R"({
  "foo": false
})");
}

TEST_CASE("json_dict with null")
{
  json_dict s;
  s.null("foo");

  REQUIRE(s.to_string() == R"({
  "foo": null
})");
}

TEST_CASE("json_dict with counts")
{
  counts c(3);
  c.inc(0, 3);
  c.inc(1, 2);
  c.inc(2, 1);

  json_dict s;
  s.i64_vec("foo", c);

  REQUIRE(s.to_string() == R"({
  "foo": [3, 2, 1]
})");
}

TEST_CASE("json_dict with rates")
{
  rates c(3);
  c.inc(0, 1);
  c.inc(1, 2.71);
  c.inc(2, 3.1415);

  json_dict s;
  s.f64_vec("foo", c);

  REQUIRE(s.to_string() == R"({
  "foo": [1.000, 2.710, 3.142]
})");
}

TEST_CASE("json_dict with vec")
{
  json_dict s;
  s.str_vec("foo", string_vec{ "foo", "b\tr" });

  REQUIRE(s.to_string() == R"({
  "foo": ["foo", "b\tr"]
})");
}

TEST_CASE("json_dict multiple values")
{
  json_dict s;
  s.boolean("foo", false);
  s.i64("some int", 17);

  REQUIRE(s.to_string() == R"({
  "foo": false,
  "some int": 17
})");
}

TEST_CASE("json_dict with empty list")
{
  json_dict s;
  auto l = s.list("foo");

  REQUIRE(s.to_string() == R"({
  "foo": []
})");
}

TEST_CASE("json_dict with list")
{
  json_dict s;
  auto l = s.list("foo");
  auto d = l->dict();
  d->f64("pi", 3);

  REQUIRE(s.to_string() == R"({
  "foo": [
    {
      "pi": 3.000
    }
  ]
})");
}

TEST_CASE("json_dict with empty sub-dict")
{
  json_dict s;
  s.dict("foo");

  REQUIRE(s.to_string() == R"({
  "foo": {}
})");
}

TEST_CASE("complex json_dict")
{
  json_dict s;
  s.i64("foo", 1234);
  auto ss = s.dict("bar");
  ss->str("name", "test");
  ss->null("value");
  s.i64("count", 1);

  REQUIRE(s.to_string() == R"({
  "foo": 1234,
  "bar": {
    "name": "test",
    "value": null
  },
  "count": 1
})");
}

////////////////////////////////////////////////////////////////////////////////
// json_list

TEST_CASE("empty json_list")
{
  json_list l;

  REQUIRE(l.to_string() == "[]");
}

TEST_CASE("json_list with empty dict")
{
  json_list l;
  l.dict();

  REQUIRE(l.to_string() == R"([
  {}
])");
}

TEST_CASE("complex json_list ")
{
  json_list l;
  l.dict();
  auto d = l.dict();
  d->str("foo", "bar");

  REQUIRE(l.to_string() == R"([
  {},
  {
    "foo": "bar"
  }
])");
}

////////////////////////////////////////////////////////////////////////////////
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
  auto r = rates(3);
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

    auto values = counts(3);
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
