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
#ifndef JSON_H
#define JSON_H

#include <ostream>
#include <string>
#include <vector>

#include "debug.hpp"

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
  void start(const std::string& key);
  /** End current sub-section. */
  void end();

  /** Write key with string value. */
  void write(const std::string& key, const std::string& value);
  /** Write key with string value array. */
  void write_vector(const std::string& key,
                    const std::vector<std::string>& value);

  /** Write key with integer value. */
  void write_int(const std::string& key, const int64_t value);
  /** Write key with integer value array. */
  void write_int_vector(const std::string& key,
                        const std::vector<int64_t>& value);

  /** Write key with floating point value. */
  void write_float(const std::string& key, const double value);
  /** Write key with boolean value (true/false). */
  void write_bool(const std::string& key, const bool value);
  /** Write key with string value. */
  void write_null(const std::string& key);

private:
  void _write(const std::string& key, const std::string& value);

  std::ostream& m_stream;
  size_t m_indent;
  bool m_values;
};

#endif
