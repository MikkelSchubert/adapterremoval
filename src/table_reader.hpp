// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef>     // for size_t
#include <limits>      // for numeric_limits
#include <ostream>     // for ostream
#include <string>      // for string
#include <string_view> // for string_view
#include <vector>      // for vector

namespace adapterremoval {

class line_reader_base;

/** Represents a row in a fixed-size table */
class table_row
{
public:
  table_row() = default;
  table_row(size_t line_num, std::vector<std::string> values);

  /** Returns the 1-based line number of the row in the input file */
  [[nodiscard]] size_t line_num() const { return m_line_num; }

  /** Returns the number of columns in the row */
  [[nodiscard]] size_t size() const { return m_values.size(); }

  /** Returns the cell of the 0-based column */
  [[nodiscard]] const std::string& at(size_t idx) const
  {
    return m_values.at(idx);
  }

  /** Returns true if both line number and values match */
  [[nodiscard]] bool operator==(const table_row& other) const noexcept;

  /** Returns true if line number or values do not match */
  [[nodiscard]] bool operator!=(const table_row& other) const noexcept;

  friend std::ostream& operator<<(std::ostream& os, const table_row& row);

private:
  //! 1-based line number; 0 for no line number
  size_t m_line_num = 0;
  //! One or more white-space separated values
  std::vector<std::string> m_values{};
};

/** Simple parser for whitespace separated tables. Empty rows are ignored */
class table_reader
{
public:
  using value_type = std::vector<table_row>;

  table_reader() = default;

  /** Ignore all text on a line past this character. Set to `\0` to disable */
  table_reader& with_comment_char(char value);
  /** Table must contain at least this many columns; must be a non-zero value */
  table_reader& with_min_columns(size_t value);
  /** Table must contain at most this many columns; must be a non-zero value */
  table_reader& with_max_columns(size_t value);
  /** Set filename/name of table for use in error messages */
  table_reader& with_name(std::string_view value);

  /** Reads table from line reader, throwing either parsing_error or io_error */
  value_type parse(line_reader_base& lr) const;

private:
  //! Optional comment char; NUL signifies no comment char
  char m_comment_char = '\0';
  //! Minimum number of columns
  size_t m_min_columns = 1;
  //! Maximum number of columns
  size_t m_max_columns = std::numeric_limits<size_t>::max();
  //! Optional (file)name used in error messages
  std::string m_name{};
};

} // namespace adapterremoval
