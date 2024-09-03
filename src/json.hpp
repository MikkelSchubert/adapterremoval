/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#pragma once

#include "commontypes.hpp" // for string_vec
#include "counts.hpp"      // for counts, rates
#include <cstddef>         // for size_t
#include <cstdint>         // for int64_t
#include <functional>      // for less
#include <map>             // for map
#include <memory>          // for shared_ptr
#include <ostream>         // for ostream
#include <string>          // for string, basic_string
#include <vector>          // for vector

namespace adapterremoval {

class json_list;
class json_dict;
class json_value;

using json_dict_ptr = std::shared_ptr<json_dict>;
using json_list_ptr = std::shared_ptr<json_list>;
using json_ptr = std::shared_ptr<json_value>;

/** Base class for json values */
class json_value
{
public:
  json_value() = default;
  virtual ~json_value() = default;

  /** Returns the value as a valid JSON string */
  virtual std::string to_string() const;

  /**
   * Writes the object to a stream. The first line is not indented, but
   * subsequent lines are indented by `indent` spaces.
   */
  virtual void write(std::ostream& out, size_t indent = 0) const = 0;

  json_value(const json_value&) = delete;
  json_value(json_value&&) = delete;
  json_value& operator=(const json_value&) = delete;
  json_value& operator=(json_value&&) = delete;
};

class json_token : public json_value
{
public:
  /** Constructor takes a single assumed-to-be-valid JSON token/string */
  explicit json_token(std::string value);

  /** Create a JSON token from a string; special characters are escaped */
  static json_ptr from_str(std::string_view value);
  /** Create a single-line JSON token representing a list strings */
  static json_ptr from_str_vec(const string_vec& values);

  /** Create a JSON token from an integer value. */
  static json_ptr from_i64(int64_t value);
  /** Create a single-line JSON token representing a list of count values */
  static json_ptr from_i64_vec(const counts& values);

  /** Create a JSON token from an integer value. */
  static json_ptr from_u64(uint64_t value);

  /** Create a JSON token from a floating point value; formatted as '%.3f' */
  static json_ptr from_f64(double value);
  /** Create a single-line JSON token representing a list of '%.3f' values */
  static json_ptr from_f64_vec(const rates& values);

  /** Create a JSON token from a boolean value */
  static json_ptr from_boolean(bool value);
  /** Create a JSON token representing null */
  static json_ptr from_null();

  /** See json_value::to_string */
  std::string to_string() const override { return m_value; }

  /** See json_value::write */
  void write(std::ostream& out, size_t indent = 0) const override;

private:
  /** Create a single-line JSON list from a set of JSON encoded values */
  static json_ptr from_raw_vec(const string_vec& values);

  //! A valid JSON token/string representing or more objects
  const std::string m_value;
};

/** Represents a list of JSON objects written across multiple lines */
class json_list : public json_value
{
public:
  /** Creates empty list */
  json_list() = default;

  /** See json_value::write */
  void write(std::ostream& out, size_t indent = 0) const override;

  /** Appends an empty JSON dictionary to the list and returns it */
  json_dict_ptr dict();

  /** Appends an empty inline JSON dictionary to the list and returns it */
  json_dict_ptr inline_dict();

private:
  //! Items in the list
  std::vector<json_ptr> m_values{};
};

/**
 * Simple write-only JSON dictionary serializer.
 *
 * Keys are written in the order they are (first) assigned.
 */
class json_dict : public json_value
{
public:
  /** Creates empty dictionary */
  json_dict() = default;

  /** See json_value::write */
  void write(std::ostream& out, size_t indent = 0) const override;

  /** Add and return a sub-dict with the specified key */
  json_dict_ptr dict(std::string_view key);
  /** Add and return an inline sub-dict with the specified key */
  json_dict_ptr inline_dict(std::string_view key);
  /** Add and return a sub-list with the specified key */
  json_list_ptr list(std::string_view key);

  /** Set key with string value. */
  void str(std::string_view key, std::string_view value);
  /** Assign string value array. */
  void str_vec(std::string_view key, const string_vec& value);

  /** Assign integer value. */
  void i64(std::string_view key, int64_t value);
  /** Assign integer value array. */
  void i64_vec(std::string_view key, const counts& value);

  /** Assign integer value. */
  void u64(std::string_view key, uint64_t value);

  /** Assign floating point value. */
  void f64(std::string_view key, double value);
  /** Assign float value array. */
  void f64_vec(std::string_view key, const rates& value);

  /** Assign boolean value (true/false). */
  void boolean(std::string_view key, bool value);
  /** Assign null value. */
  void null(std::string_view key);

private:
  friend class json_list;

  /** Assigns a JSON object to a given key and updates the list of keys */
  void _set(std::string_view key, const json_ptr& ptr);

  //! List of keys in insertion order
  std::vector<std::pair<std::string, json_ptr>> m_items{};
  //! Multi-line or single (inline) dictionary
  bool m_multi_line = true;
};

} // namespace adapterremoval
