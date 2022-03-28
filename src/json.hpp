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

#include <cmath>    // for isnan
#include <sstream>  // for operator<<, stringstream, basic_ostream, fixed
#include <stddef.h> // for size_t
#include <stdint.h> // for int64_t
#include <string>   // for string, char_traits
#include <vector>   // for vector

class json_writer;
template<typename T>
class counts_tmpl;

/** Helper class for creating sub-sections in JSON output. */
class json_section
{
public:
  /** Move constructor. **/
  json_section(json_section&& other) noexcept;

  /** Dictionary as value; for use with json_writer::start_list. */
  explicit json_section(json_writer* parent);
  /** Dictionary as key/value pair; . */
  explicit json_section(json_writer* parent, const std::string& key);

  /** Ends section. */
  ~json_section();

  inline operator bool() const { return true; }

  //! Copy construction not supported
  json_section(const json_section&) = delete;
  //! Assignment not supported
  json_section& operator=(const json_section&) = delete;

private:
  json_writer* m_parent;
};

/**
 * Writes pretty-printed JSON to a stream.
 *
 * Values are written as they are pushed and no checks for duplicate values are
 * performed. This class is meant as a simple alternative to writing simple JSON
 * output by hand, not as a full-fledged JSON serializer.
 */
class json_writer
{
public:
  json_writer(std::ostream& out);
  ~json_writer();

  /** Start new sub-section. */
  json_section start();
  json_section start(const std::string& key);

  /** Start new list. */
  void start_list(const std::string& key);
  /** Terminate a list. */
  void end_list();

  /** Write key with string value. */
  void write(const std::string& key, const std::string& value);
  /** Write key with integer value array. */
  template<typename T>
  void write(const std::string& key, const counts_tmpl<T>& value);

  /** Write key with string value array. */
  void write(const std::string& key, const std::vector<std::string>& value);

  /** Write key with integer value. */
  void write_int(const std::string& key, const int64_t value);
  /** Write key with floating point value. */
  void write_float(const std::string& key, const double value);
  /** Write key with boolean value (true/false). */
  void write_bool(const std::string& key, const bool value);
  /** Write key with string value. */
  void write_null(const std::string& key);

  //! Copy construction not supported
  json_writer(const json_writer&) = delete;
  //! Assignment not supported
  json_writer& operator=(const json_writer&) = delete;

private:
  /** End current sub-section. */
  void end(char c);
  /** Writes a (assumed to be) valid JSON value. */
  void _write(const std::string& value);
  /** Writes a key and value. The value is assumed to be a valid JSON value. */
  void _write(const std::string& key, const std::string& value);

  std::ostream& m_stream;
  size_t m_indent;
  bool m_values;

  friend class json_section;
};

template<typename T>
void
json_writer::write(const std::string& key, const counts_tmpl<T>& value)
{
  std::stringstream ss;
  ss.precision(3);
  ss << std::fixed << "[";
  for (size_t i = 0; i < value.size(); ++i) {
    if (i) {
      ss << ", ";
    }

    if (std::isnan(value.get(i))) {
      ss << "null";
    } else {
      ss << value.get(i);
    }
  }

  ss << "]";

  _write(key, ss.str());
}
