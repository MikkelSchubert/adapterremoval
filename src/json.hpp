/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#pragma once

#include <map>      // for map
#include <memory>   // for shared_ptr
#include <ostream>  // for ostream
#include <stddef.h> // for size_t
#include <stdint.h> // for int64_t
#include <string>   // for string, char_traits

#include "commontypes.hpp" // for string_vec
#include "counts.hpp"      // for counts, rates

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
  virtual ~json_value();

  /** Returns the value as a valid JSON string */
  virtual std::string to_string() const;

  /**
   * Writes the object to a stream. The first line is not indented, but
   * subsequent lines are indented by `indent` spaces.
   */
  virtual void write(std::ostream& out, size_t indent = 0) const = 0;
};

class json_token : public json_value
{
public:
  /** Constructor takes a single assumed-to-be-valid JSON token/string */
  json_token(const std::string& value);

  /** Create a JSON token from a string; special characters are escaped */
  static json_ptr from_str(const std::string& value);
  /** Create a single-line JSON token representing a list strings */
  static json_ptr from_str_vec(const string_vec& values);

  /** Create a JSON token from an integer value. */
  static json_ptr from_i64(const int64_t value);
  /** Create a single-line JSON token representing a list of count values */
  static json_ptr from_i64_vec(const counts& values);

  /** Create a JSON token from a floating point value; formatted as '%.3f' */
  static json_ptr from_f64(const double value);
  /** Create a single-line JSON token representing a list of '%.3f' values */
  static json_ptr from_f64_vec(const rates& values);

  /** Create a JSON token from a boolean value */
  static json_ptr from_boolean(const bool value);
  /** Create a JSON token representing null */
  static json_ptr from_null();

  /** See json_value::to_string */
  virtual std::string to_string() const override;
  /** See json_value::write */
  virtual void write(std::ostream& out, size_t indent = 0) const override;

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
  json_list();

  /** See json_value::write */
  virtual void write(std::ostream& out, size_t indent = 0) const override;

  /** Appends an empty JSON dictionary to the list and returns it */
  json_dict_ptr dict();

private:
  //! Items in the list
  std::vector<json_ptr> m_values;
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
  json_dict();

  /** See json_value::write */
  virtual void write(std::ostream& out, size_t indent = 0) const override;

  /** Add and return a sub-dict with the specified key */
  json_dict_ptr dict(const std::string& key);
  /** Add and return a sub-list with the specified key */
  json_list_ptr list(const std::string& key);

  /** Set key with string value. */
  void str(const std::string& key, const std::string& value);
  /** Assign string value array. */
  void str_vec(const std::string& key, const string_vec& value);

  /** Assign integer value. */
  void i64(const std::string& key, const int64_t value);
  /** Assign integer value array. */
  void i64_vec(const std::string& key, const counts& value);

  /** Assign floating point value. */
  void f64(const std::string& key, const double value);
  /** Assign float value array. */
  void f64_vec(const std::string& key, const rates& value);

  /** Assign boolean value (true/false). */
  void boolean(const std::string& key, const bool value);
  /** Assign null value. */
  void null(const std::string& key);

private:
  /** Assigns a JSON object to a given key and updates the list of keys */
  void _set(const std::string& key, const json_ptr& ptr);

  //! List of keys in insertion order
  std::vector<std::string> m_keys;
  //! Map of unencoded keys to JSON objects
  std::map<std::string, json_ptr> m_values;
};
