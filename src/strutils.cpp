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
#include <algorithm>   // for reverse
#include <cmath>       // for floor
#include <iomanip>     // for operator<<, setw
#include <limits>      // for numeric_limits
#include <sstream>     // for istringstream, ostringstream
#include <stdexcept>   // for invalid_argument
#include <stdint.h>    // for int64_t
#include <sys/ioctl.h> // for ioctl, winsize, TIOCGWINSZ
#include <unistd.h>    // for STDOUT_FILENO
#include <vector>      // for vector

#include "debug.hpp" // for AR_REQUIRE
#include "strutils.hpp"

namespace adapterremoval {

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
tolower(const std::string& str)
{
  std::string uppercased = str;
  for (auto& current : uppercased) {
    if (current >= 'A' && current <= 'Z') {
      current += 32;
    }
  }

  return uppercased;
}

std::string
toupper(const std::string& str)
{
  std::string uppercased = str;
  for (auto& current : uppercased) {
    if (current >= 'a' && current <= 'z') {
      current -= 32;
    }
  }

  return uppercased;
}

std::string
indent_lines(const std::string& lines, size_t n_indent)
{
  std::string line;
  std::ostringstream lines_out;
  const std::string indentation = std::string(n_indent, ' ');

  size_t last_pos = 0;
  while (true) {
    const size_t next_pos = lines.find('\n', last_pos);
    if (next_pos == std::string::npos) {
      if (lines.size() - last_pos) {
        lines_out << indentation;
      }

      lines_out << lines.substr(last_pos);
      break;
    } else {
      if (next_pos - last_pos) {
        lines_out << indentation;
      }

      lines_out << lines.substr(last_pos, next_pos - last_pos + 1);
      last_pos = next_pos + 1;
    }
  }

  return lines_out.str();
}

std::string
template_replace(const std::string& haystack,
                 const std::string& needle,
                 const std::string& value)
{
  AR_REQUIRE(needle.size());

  std::string result;
  std::size_t last = 0;
  std::size_t pos = haystack.find(needle, 1);
  while (pos != std::string::npos) {
    bool prefix = false;
    if (pos >= 2 && haystack.at(pos - 2) == '{') {
      prefix = true;
    } else if (haystack.at(pos - 1) != '{') {
      pos = haystack.find(needle, pos + 1);
      continue;
    }

    bool postfix = false;
    const size_t end = pos + needle.size();
    if (end + 1 < haystack.size() && haystack.at(end + 1) == '}') {
      postfix = true;
    } else if (end >= haystack.size() || haystack.at(end) != '}') {
      pos = haystack.find(needle, pos + 1);
      continue;
    }

    result.append(haystack.substr(last, pos - last - 1 - prefix));

    if (value.size()) {
      if (prefix) {
        result.push_back(haystack.at(pos - 1));
      }

      result.append(value);

      if (postfix) {
        result.push_back(haystack.at(pos + needle.length()));
      }
    }

    last = pos + needle.size() + 1 + postfix;
    pos = haystack.find(needle, last);
  }

  result.append(haystack.substr(last));

  return result;
}

std::string
columnize_text(const std::string& value, size_t max_width, size_t ljust)
{
  size_t current_width = 0;
  size_t current_ljust = 0;
  std::ostringstream lines_out;
  std::istringstream lines_in(value);
  std::string substr;

  while (lines_in >> substr) {
    if (current_width) {
      if (current_ljust + current_width + 1 + substr.length() > max_width) {
        current_ljust = ljust;
        lines_out << "\n"
                  << std::left << std::setw(current_ljust) << "" << substr;
        current_width = substr.length();
      } else {
        lines_out << " " << substr;
        current_width += substr.length() + 1;
      }
    } else {
      lines_out << substr;
      current_width += substr.length();
    }
  }

  return lines_out.str();
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'cli_formatter'

cli_formatter::cli_formatter()
  : m_indent_first(true)
  , m_ljust(0)
  , m_columns(DEFAULT_MAX_COLUMNS)
  , m_indentation(4)
{
  struct winsize size;
  if (!ioctl(STDOUT_FILENO, TIOCGWINSZ, &size)) {
    m_columns = size.ws_col;
  }
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
cli_formatter::format(const std::string& lines) const
{
  std::string line;
  std::ostringstream lines_out;

  size_t last_pos = 0;
  while (true) {
    const size_t next_pos = lines.find('\n', last_pos);

    // Current line, excluding terminal newline
    const size_t end_pos =
      (next_pos == std::string::npos) ? next_pos : (next_pos - last_pos);
    line = lines.substr(last_pos, end_pos);

    // Format into fixed width columns, indenting by ljust after first line
    line = columnize_text(line, m_columns, m_ljust);

    // Add fixed width indentation
    if (m_indentation) {
      line = indent_lines(line, m_indentation);
    }

    lines_out << line;

    if (next_pos == std::string::npos) {
      break;
    } else if (lines.at(next_pos) == '\n') {
      lines_out << '\n';
      last_pos = next_pos + 1;
    }
  }

  line = lines_out.str();
  if (!m_indent_first) {
    line = line.substr(m_indentation);
  }

  return line;
}

std::string
cli_formatter::fmt(const std::string& value)
{
  cli_formatter formatter;

  return formatter.format(value);
}

std::string
cli_formatter::fmt(const std::string& prefix, const std::string& value)
{
  cli_formatter formatter;
  formatter.set_indent(prefix.size());
  formatter.set_indent_first_line(false);

  return prefix + value;
}

std::string
shell_escape(const std::string& s)
{
  if (s.empty()) {
    return "''";
  }

  bool must_escape = false;
  for (const auto c : s) {
    // Conservative list of safe values; better safe than sorry
    if (!isalnum(c) && c != '_' && c != '.' && c != '/' && c != '-') {
      must_escape = true;
      break;
    }
  }

  if (!must_escape) {
    return s;
  }

  std::string out;
  out.push_back('\'');

  for (const auto c : s) {
    if (c == '\'') {
      out.push_back('\\');
    }

    out.push_back(c);
  }

  out.push_back('\'');

  return out;
}

std::string
shell_escape_command(const std::vector<std::string>& values)
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

static const std::vector<std::pair<double, char>> rough_units = {
  { 1e3, 'K' }, { 1e6, 'M' }, { 1e9, 'G' }, { 1e12, 'T' }, { 1e15, 'P' },
};

std::string
format_rough_number(size_t count, size_t precision)
{
  if (count < 1000) {
    return std::to_string(count);
  }

  size_t scale_i = 0;
  for (size_t i = 0; i < rough_units.size(); ++i) {
    if (rough_units.at(i).first >= count) {
      break;
    }

    scale_i = i;
  }

  const double shift = std::pow(10, precision);
  const double scale = rough_units.at(scale_i).first;
  double tmp = std::round((count / scale) * shift);

  if (tmp >= shift * 1e3 && scale_i + 1 < rough_units.size()) {
    // count was rounded up so that it falls into the next unit
    scale_i++;
    tmp /= 1e3;
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(precision) << (tmp / shift) << ' '
     << rough_units.at(scale_i).second;

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
