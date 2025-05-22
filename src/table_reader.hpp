// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef>     // for size_t
#include <string>      // for string
#include <string_view> // for string_view
#include <vector>      // for vector

namespace adapterremoval {

class line_reader_base;

struct table_row
{
  //! 1-based line number
  size_t line_num = 0;
  //! One or more white-space separated values
  std::vector<std::string> values;

  [[nodiscard]] size_t size() const { return this->values.size(); }

  [[nodiscard]] const std::string& at(size_t idx) const
  {
    return this->values.at(idx);
  }
};

/** Simple parser for whitespace separated tables. Empty rows are ignored */
class table_reader
{
public:
  using value_type = std::vector<table_row>;

  table_reader() = default;

  /** Skip lines starting with this character; leading whitespace allowed */
  table_reader& with_comment_char(char value);
  /** Table must contain at least this many columns; most be a non-zero value */
  table_reader& with_min_columns(size_t value);
  /** Table must contain at most this many columns; most be a non-zero value */
  table_reader& with_max_columns(size_t value);
  /** Set filename/name of table for use in error messages */
  table_reader& with_name(std::string_view value);

  /** Reads table from line reader, throwing either parsing_error or io_error */
  value_type parse(line_reader_base& lr) const;

private:
  //! Optional comment char; NUL signifies no comment char
  char m_comment_char = '\0';
  //! Minimum number of columns; must be at least one
  size_t m_min_columns = 1;
  //! Maximum number of columns; must be greater than m_min_columns
  size_t m_max_columns = static_cast<size_t>(-1);
  //! Optional (file)name used in error messages
  std::string m_name{};
};

} // namespace adapterremoval
