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

#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t
#include <ostream>     // for ostream
#include <sstream>     // for ostringstream
#include <string>      // for string
#include <string_view> // for string_view
#include <vector>      // for vector

namespace adapterremoval {

using string_vec = std::vector<std::string>;

const size_t DEFAULT_MAX_COLUMNS = 78;
const size_t DEFAULT_INDENTATION = 4;

/** Computes the Levenshtein distance between two strings. */
size_t
levenshtein(std::string_view s, std::string_view t);

/**
 * Returns a timestamp in the specified format, see
 * https://en.cppreference.com/w/cpp/io/manip/put_time
 */
std::string
timestamp(const char* format, bool milliseconds = false);

/**
 * Convert a string to a uint32_t.
 *
 * Throws std::invalid_argument if the string does not contain a proper number,
 * if it contains more than just a number, or if the number overflows a
 * 32 bit unsigned integer. Whitespace is ignored.
 */
uint32_t
str_to_u32(std::string_view s);

/** Convert a string to a double using the same rules as str_to_u32  */
double
str_to_double(const std::string& s);

/** Lowercases letters in the range a-z */
constexpr char
to_lower(char c)
{
  return (c >= 'A' && c <= 'Z') ? static_cast<char>(c + 32) : c;
}

/** Lowercases letters in the range a-z */
std::string
to_lower(std::string str);

/** Uppercase letters in the range a-z */
constexpr char
to_upper(char c)
{
  return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 32) : c;
}

/** Uppercase letters in the range a-z */
std::string
to_upper(std::string str);

/** Returns true if str1 ends with str2 (case sensitive) */
bool
starts_with(std::string_view str1, std::string_view str2);

/** Returns true if str1 ends with str2 (case sensitive) */
bool
ends_with(std::string_view str1, std::string_view str2);

string_vec
split_text(std::string_view text, char separator);

/** Split text by newlines */
inline string_vec
split_lines(std::string_view text)
{
  return split_text(text, '\n');
}

/** Split text by newlines and add fixed indentation following newlines. */
std::string
indent_lines(std::string_view lines, size_t indentation = DEFAULT_INDENTATION);

/**
 * Formats text into fixed-width columns.
 *
 * @param value Text representing a single paragraph to be formatted.
 * @param max_width Maximum width of output lines in characters.
 * @param ljust Indent lines after the first line by this amount of characters.
 *
 * Note that all whitespace in the input string is consumed; output words are
 * separated by a single space, and the terminal line does not end with a
 * newline.
 */
string_vec
wrap_text(const std::string& value,
          size_t max_width = DEFAULT_MAX_COLUMNS,
          size_t ljust = 0);

template<typename T>
std::string
join_text(const std::vector<T>& values,
          std::string_view sep,
          std::string_view final_sep)
{
  std::ostringstream stream;
  for (auto it = values.begin(); it != values.end(); ++it) {
    if (it != values.begin()) {
      if (it + 1 == values.end()) {
        stream << final_sep;
      } else {
        stream << sep;
      }
    }

    stream << *it;
  }

  return stream.str();
}

template<typename T>
std::string
join_text(const std::vector<T>& values, std::string_view sep)
{
  return join_text(values, sep, sep);
}

/** Trims whitespace from both ends of a string view */
std::string_view
trim_ascii_whitespace(std::string_view s);

/**
 * Wrapper around 'indent_lines' and 'wrap_text'.
 *
 * Defaults to 78 character columns, with 4 character indentation, and no ljust,
 * corresponding to function defaults. Unlike simply combining 'indent_lines'
 * and 'wrap_text', the 'format' function preserves newlines, allowing
 * paragraphs to be printed easily.
 */
class cli_formatter
{
public:
  /** Creates formatter using default parameters (see above). */
  cli_formatter();

  /** Include (or not) indentation on the first line of output. */
  cli_formatter& set_indent_first_line(bool value);
  /** Maximum number of columns in each line in output. */
  cli_formatter& set_column_width(size_t value);
  /** Number of spaces to indent line 2+ in each paragraph. */
  cli_formatter& set_ljust(size_t value);
  /** Fixed indentation for each line (see also set_indent_first_line). */
  cli_formatter& set_indent(size_t value);

  /** Formats string using the current settings. */
  std::string format(std::string_view value) const;

private:
  //! Specifies whether or not to indent the first line of output.
  bool m_indent_first;
  //! Number of spaces to indent the 2+ lines in each paragraph.
  size_t m_ljust;
  //! The maximum number of columns (in characters).
  size_t m_columns;
  //! The number of spaces to indent each line (see also m_indent_first_line)
  size_t m_indentation;
};

/** Escapes a string to make it safe to use in copy/pasted commands */
std::string
shell_escape(std::string_view s);

/**
 * Escapes a string for use in logging; this is equivalent to shell_escape
 * except that strings are always quoted for readability.
 */
std::string
log_escape(std::string_view s);

/** Escapes a full shell command using shell_escape */
std::string
shell_escape_command(const string_vec& v);

/** Performs basic HTML escape of [<>'"&] */
std::string
html_escape(std::string_view s);

/** Adds thousand separators to a number */
std::string
format_thousand_sep(size_t count);

/** Rounds a number using K, M, etc. units. */
std::string
format_rough_number(size_t value, size_t out_digits = 3);

/** Formats a fraction, returning "NA" if denominator is 0 */
std::string
format_fraction(uint64_t num, uint64_t denom, size_t precision = 2);

/** Formats percentage and adds trailing " %". Returns "NA" if denom is 0 */
std::string
format_percentage(uint64_t num, uint64_t denom, size_t precision = 1);

} // namespace adapterremoval
