/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "table_reader.hpp" // declarations
#include "debug.hpp"
#include "errors.hpp"     // for parsing_error
#include "linereader.hpp" // for line_reader
#include "strutils.hpp"   // for trim_ascii_whitespace
#include <sstream>
#include <string_view> // for string_view

namespace adapterremoval {

namespace {

template<typename T>
void
throw_table_error(std::string_view name, size_t linenum, std::string_view error)
{
  std::ostringstream ss;
  if (!name.empty() || linenum != static_cast<size_t>(-1)) {
    ss << "Error ";

    if (linenum != static_cast<size_t>(-1)) {
      ss << " at line " << linenum;
    }

    if (!name.empty()) {
      ss << " in table " << name;
    }

    ss << ": ";
  }

  ss << error;

  throw T(ss.str());
}

} // namespace

table_reader&
table_reader::with_comment_char(char value)
{
  m_comment_char = value;
  return *this;
}

table_reader&
table_reader::with_min_columns(size_t value)
{
  AR_REQUIRE(value > 0 && value <= m_max_columns);
  m_min_columns = value;
  return *this;
}

table_reader&
table_reader::with_max_columns(size_t value)
{
  AR_REQUIRE(value > 0 && value >= m_min_columns);
  m_max_columns = value;
  return *this;
}

table_reader&
table_reader::with_name(std::string_view value)
{
  m_name = value;
  return *this;
}

table_reader::value_type
table_reader::parse(line_reader_base& lr) const
{
  std::string field;
  std::string buffer;
  value_type table;

  try {
    for (size_t linenum = 1; lr.getline(buffer); ++linenum) {
      const size_t index = buffer.find('#');
      if (index != std::string::npos) {
        buffer.resize(index);
      }

      std::string_view line = trim_ascii_whitespace(buffer);
      if (line.empty()) {
        continue;
      }

      table_row row{ linenum, {} };
      std::istringstream instream{ buffer };
      for (; instream >> field; field.clear()) {
        row.values.push_back(field);
      }

      if (row.size() < m_min_columns) {
        std::ostringstream message;
        message << "Expected at least " << m_min_columns
                << " columns, but found " << row.size() << " column(s)";

        throw_table_error<parsing_error>(m_name, linenum, message.str());
      } else if (row.size() > m_max_columns) {
        std::ostringstream message;
        message << "Expected at most " << m_max_columns
                << " columns, but found " << row.size() << " column(s)";

        throw_table_error<parsing_error>(m_name, linenum, message.str());
      } else if (!table.empty() && table.back().size() != row.size()) {
        std::ostringstream message;
        message << "Inconsistent number of columns; expected "
                << table.back().size() << " column(s) but found " << row.size();

        throw_table_error<parsing_error>(m_name, linenum, message.str());
      }

      table.emplace_back(std::move(row));
    }
  } catch (const io_error& error) {
    throw_table_error<io_error>(m_name, static_cast<size_t>(-1), error.what());
  }

  return table;
}

} // namespace adapterremoval
