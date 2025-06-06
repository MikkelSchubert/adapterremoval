// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "table_reader.hpp" // declarations
#include "debug.hpp"        // for AR_REQUIRE
#include "errors.hpp"       // for io_error, parsing_error
#include "linereader.hpp"   // for line_reader_base
#include "strutils.hpp"     // for trim_ascii_whitespace
#include <sstream>          // for ostringstream
#include <string_view>      // for operator<<, string_view
#include <utility>          // for move

namespace adapterremoval {

namespace {

template<typename T>
void
throw_table_error(std::string_view name, size_t linenum, std::string_view error)
{
  std::ostringstream ss;
  if (!name.empty() || linenum != static_cast<size_t>(-1)) {
    ss << "Error";

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
