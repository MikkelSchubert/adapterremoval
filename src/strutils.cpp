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
#include "strutils.hpp" // declarations
#include "debug.hpp"    // for AR_REQUIRE
#include <algorithm>    // for min, reverse, max
#include <cctype>       // for isprint, isalnum, tolower, toupper
#include <chrono>       // for system_clock
#include <cmath>        // for log10, pow, round
#include <cstdint>      // for uint64_t, int64_t
#include <iomanip>      // for operator<<, setprecision
#include <limits>       // for numeric_limits
#include <sstream>      // for ostringstream, operator<<, basic_ostream, bas...
#include <stdexcept>    // for invalid_argument
#include <unistd.h>     // for STDOUT_FILENO
#include <vector>       // for vector, swap

namespace adapterremoval {

namespace {

string_vec
indent(string_vec lines, size_t n_indent, bool indent_first)
{
  const std::string indentation(n_indent, ' ');
  for (auto it = lines.begin(); it != lines.end(); ++it) {
    if (!it->empty() && (indent_first || it != lines.begin())) {
      it->insert(0, indentation);
    }
  }

  return lines;
}

} // namespace

size_t
levenshtein(const std::string& s, const std::string& t)
{
  std::vector<size_t> v0(t.size() + 1, 0);
  std::vector<size_t> v1(t.size() + 1, 0);

  for (size_t i = 0; i < v0.size(); ++i) {
    v0.at(i) = i;
  }

  for (size_t i = 0; i < s.size(); ++i) {
    v1.at(0) = i + 1;

    for (size_t j = 0; j < t.size(); ++j) {
      const auto del = v0.at(j + 1) + 1;
      const auto ins = v1.at(j) + 1;
      const auto sub = s.at(i) == t.at(j) ? v0.at(j) : v0.at(j) + 1;

      v1.at(j + 1) = std::min(del, std::min(ins, sub));
    }

    std::swap(v0, v1);
  }

  return v0.back();
}

std::string
timestamp(const char* format, const bool milliseconds)
{
  AR_REQUIRE(format);
  using namespace std::chrono;

  const auto now = system_clock::now();
  const auto in_time_t = system_clock::to_time_t(now);

  tm in_localtime{};
  std::ostringstream ss;
  ss << std::put_time(localtime_r(&in_time_t, &in_localtime), format);

  if (milliseconds) {
    const auto ms =
      duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
    ss << '.' << std::setfill('0') << std::setw(3) << (ms.count() % 1000);
  }

  return ss.str();
}

unsigned
str_to_unsigned(const std::string& s)
{
  std::istringstream stream(s);
  int64_t temp = 0;

  if (!(stream >> temp)) {
    throw std::invalid_argument("value is not a valid number");
  }

  // Failing on trailing, non-numerical values
  std::string trailing;
  if (stream >> trailing) {
    throw std::invalid_argument("value contains trailing text");
  }

  if (temp < 0 || temp > std::numeric_limits<unsigned>::max()) {
    throw std::invalid_argument("numerical value overflows");
  }

  return static_cast<unsigned>(temp);
}

std::string
to_lower(const std::string& str)
{
  std::string lowercase = str;
  for (auto& current : lowercase) {
    current = to_lower(current);
  }

  return lowercase;
}

std::string
to_upper(const std::string& str)
{
  std::string uppercase = str;
  for (auto& current : uppercase) {
    current = to_upper(current);
  }

  return uppercase;
}

/** Returns true if str1 ends with str2 (case sensitive) */
bool
ends_with(const std::string& str1, const std::string& str2)
{
  if (str1.size() < str2.size()) {
    return false;
  }

  auto it1 = str1.rbegin();
  auto it2 = str2.rbegin();
  for (; it2 != str2.rend(); ++it1, ++it2) {
    if (*it1 != *it2) {
      return false;
    }
  }

  return true;
}

string_vec
split_lines(const std::string& text)
{
  string_vec lines;

  size_t start = 0;
  size_t end = std::string::npos;
  do {
    end = text.find('\n', start);

    lines.push_back(text.substr(start, end - start));

    start = end + 1;
  } while (end != std::string::npos);

  return lines;
}

std::string
indent_lines(const std::string& lines, size_t indentation)
{
  return join_text(indent(split_lines(lines), indentation, true), "\n");
}

string_vec
wrap_text(const std::string& value, size_t max_width, size_t ljust)
{
  size_t current_width = 0;
  size_t current_ljust = 0;
  std::istringstream lines_in(value);
  string_vec lines_out;

  std::string substr;
  while (lines_in >> substr) {
    if (current_width) {
      if (current_ljust + current_width + 1 + substr.length() > max_width) {
        current_ljust = ljust;
        lines_out.emplace_back(current_ljust, ' ');
        lines_out.back().append(substr);
        current_width = substr.length();
      } else {
        lines_out.back().push_back(' ');
        lines_out.back().append(substr);
        current_width += substr.length() + 1;
      }
    } else {
      lines_out.push_back(substr);
      current_width += substr.length();
    }
  }

  return lines_out;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'cli_formatter'

cli_formatter::cli_formatter()
  : m_indent_first(true)
  , m_ljust(0)
  , m_columns(DEFAULT_MAX_COLUMNS)
  , m_indentation(4)
{
}

cli_formatter&
cli_formatter::set_column_width(size_t value)
{
  m_columns = value;

  return *this;
}

cli_formatter&
cli_formatter::set_ljust(size_t value)
{
  m_ljust = value;

  return *this;
}

cli_formatter&
cli_formatter::set_indent(size_t value)
{
  m_indentation = value;

  return *this;
}

cli_formatter&
cli_formatter::set_indent_first_line(bool value)
{
  m_indent_first = value;

  return *this;
}

std::string
cli_formatter::format(const std::string& value) const
{
  std::ostringstream lines_out;

  for (const auto& line : split_lines(value)) {
    const auto block = wrap_text(line, m_columns, m_ljust);

    lines_out << join_text(indent(block, m_indentation, m_indent_first), "\n");
  }

  return lines_out.str();
}

std::string
shell_escape(const std::string& s)
{
  if (s.empty()) {
    return "''";
  }

  for (const auto c : s) {
    // Conservative list of safe values; better safe than sorry
    if (!isalnum(c) && c != '_' && c != '.' && c != '/' && c != '-') {
      return log_escape(s);
    }
  }

  return s;
}

std::string
log_escape(const std::string& s)
{
  std::string out;
  out.push_back('\'');

  for (const auto c : s) {
    switch (c) {
      case '\'':
        out.append("\\'");
        break;
      case '\\':
        out.append("\\\\");
        break;
      case '\b':
        out.append("\\b");
        break;
      case '\f':
        out.append("\\f");
        break;
      case '\n':
        out.append("\\n");
        break;
      case '\r':
        out.append("\\r");
        break;
      case '\t':
        out.append("\\t");
        break;
      default:
        if (!std::isprint(c)) {
          std::ostringstream ss;
          ss << "\\x" << std::hex << static_cast<int>(c);

          out.append(ss.str());

        } else {
          out.push_back(c);
        }
    }
  }

  out.push_back('\'');

  return out;
}

std::string
shell_escape_command(const string_vec& values)
{
  std::ostringstream ss;
  for (size_t i = 0; i < values.size(); ++i) {
    if (i) {
      ss << ' ';
    }

    ss << shell_escape(values.at(i));
  }

  return ss.str();
}

std::string
format_thousand_sep(size_t count)
{
  if (!count) {
    return "0";
  }

  std::string ss;
  for (size_t i = 0; count; ++i) {
    ss.push_back('0' + (count % 10));
    count /= 10;
    if (count && i && i % 3 == 2) {
      ss.push_back(',');
    }
  }

  std::reverse(ss.begin(), ss.end());

  return ss;
}

std::string
format_rough_number(size_t value, size_t out_digits)
{
  AR_REQUIRE(out_digits > 0);
  if (value == 0) {
    return "0";
  }

  auto rounded = static_cast<double>(value);
  auto in_digits = static_cast<size_t>(std::log10(rounded));
  if (out_digits > in_digits) {
    return std::to_string(value);
  }

  // Round to desired number of significant digits
  const auto tmp = std::pow(10, in_digits - out_digits + 1);
  rounded = std::round(rounded / tmp) * tmp;

  // Rounding up may result in the number of digits increasing
  in_digits = static_cast<size_t>(std::log10(rounded));

  const std::string units = "KMGTP";
  const size_t unit = std::min<size_t>(units.size(), in_digits / 3);
  const double scaled = rounded / std::pow(10.0, unit * 3);
  const size_t precision =
    out_digits - std::min<size_t>(out_digits, in_digits - unit * 3 + 1);

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(precision) << scaled;

  if (unit) {
    ss << " " << units.at(unit - 1);
  }

  return ss.str();
}

std::string
format_fraction(uint64_t num, uint64_t denom, size_t precision)
{
  if (denom) {
    const double fraction = static_cast<double>(num) / denom;

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << fraction;

    return ss.str();
  } else {
    return "NA";
  }
}

std::string
format_percentage(uint64_t num, uint64_t denom, size_t precision)
{
  return format_fraction(num * 100, denom, precision);
}

} // namespace adapterremoval
