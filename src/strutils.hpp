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

#include <stddef.h> // for size_t
#include <string>   // for string
#include <vector>   // for vector

namespace adapterremoval {

using string_vec = std::vector<std::string>;

const size_t DEFAULT_MAX_COLUMNS = 78;
const size_t DEFAULT_INDENTATION = 4;

/** Computes the Levenshtein distance between two strings. */
size_t
levenshtein(const std::string& s, const std::string& t);

/**
 * Convert a string to an unsigned integer.
 *
 * Throws std::invalid_argument if the string does not contain a proper number,
 * if it contains more than just a number, or if the number overflows a
 * unsigned integer.
 */
unsigned
str_to_unsigned(const std::string& s);

/** Lowercases letters in the range a-z */
char
to_lower(char c);

/** Lowercases letters in the range a-z */
std::string
to_lower(const std::string& str);

/** Lowercases letters in the range a-z */
char
to_upper(char c);

/** Uppercases letters in the range a-z */
std::string
to_upper(const std::string& str);

/** Returns true if str1 ends with str2 (case sensitive) */
bool
ends_with(const std::string& str1, const std::string& str2);

/** Split text by newlines */
string_vec
split_lines(const std::string& text);

/** Split text by newlines and add fixed indentation following newlines. */
std::string
indent_lines(const std::string& lines, size_t identation = DEFAULT_INDENTATION);

/**
 * Formats text into fixed-width columns.
 *
 * @param value Text representing a single paragraph to be formatted.
 * @param max_width Maximum width of output lines in characters.
 * @param ljust Indent lines after the first line by this amount of characters.
 *
 * Note that all whitespace in the input string is consumed; output words are
 * seperated by a single space, and the terminal line does not end with a
 * newline.
 */
string_vec
wrap_text(const std::string& value,
          size_t max_width = DEFAULT_MAX_COLUMNS,
          size_t ljust = 0);

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
  std::string format(const std::string& value) const;

private:
  //! Specifices whether or not to indent the first line of output.
  bool m_indent_first;
  //! Number of spaces to indent the 2+ lines in each paragraph.
  size_t m_ljust;
  //! The maximum number of columns (in characters).
  size_t m_columns;
  //! The number of spaces to indent each line (see also m_indent_first_line)
  size_t m_indentation;
};

/** Escapes a string to make safe to use in copy/pasted commands */
std::string
shell_escape(const std::string& s);

/** Escapes a full shell command using shell_escape */
std::string
shell_escape_command(const string_vec& v);

/** Adds thousand separators to a number */
std::string
format_thousand_sep(size_t count);

/** Rounds a number using K, M, etc. units. */
std::string
format_rough_number(size_t count, size_t digits = 3);

/** Formats a fraction, returning "NA" if denominator is 0 */
std::string
format_fraction(uint64_t num, uint64_t denom, size_t precision = 2);

/** Formats a percentage, returning "NA" if denominator is 0 */
std::string
format_percentage(uint64_t num, uint64_t denom, size_t precision = 1);

} // namespace adapterremoval
