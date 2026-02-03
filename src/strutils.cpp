// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "strutils.hpp" // declarations
#include "debug.hpp"    // for AR_REQUIRE
#include <algorithm>    // for min, reverse, max
#include <cctype>       // for isprint, isalnum, tolower, toupper
#include <charconv>     // for from_chars
#include <cmath>        // for log10, pow, round
#include <cstdint>      // for uint64_t, int64_t
#include <iomanip>      // for operator<<, setprecision
#include <sstream>      // for ostringstream, operator<<, basic_ostream, bas...
#include <stdexcept>    // for invalid_argument
#include <string_view>  // for string_view
#include <system_error> // for errc
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
levenshtein(const std::string_view s, const std::string_view t)
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

namespace {

template<typename T>
inline std::string
stringify_int(T value)
{
  // WORKAROUND: to_chars(int) requires GCC >= 8 and LLVM >= 14?
  std::ostringstream os;
  os << value;

  return os.str();
}

} // namespace

std::string
stringify(int value)
{
  return stringify_int(static_cast<long long>(value));
}

std::string
stringify(long value)
{
  return stringify_int(static_cast<long long>(value));
}

std::string
stringify(long long value)
{
  return stringify_int(value);
}

std::string
stringify(unsigned value)
{
  return stringify_int(static_cast<unsigned long long>(value));
}

std::string
stringify(unsigned long value)
{
  return stringify_int(static_cast<unsigned long long>(value));
}

std::string
stringify(unsigned long long value)
{
  return stringify_int(value);
}

std::string
stringify(double value)
{
  // WORKAROUND: to_chars(double) requires GCC >= 11 and LLVM >= 14?
  std::ostringstream os;
  os << std::fixed << std::setprecision(6) << value;

  return os.str();
}

namespace {

template<typename T>
inline T
str_to(std::string_view s)
{
  T value{};

  s = trim_ascii_whitespace(s);
  const auto* begin = s.data();
  const auto* end = begin + s.size();
  const auto result = std::from_chars(begin, end, value);

  if (result.ec != std::errc{}) {
    throw std::invalid_argument("value is not a valid number");
  } else if (result.ptr != end) {
    throw std::invalid_argument("number contains trailing text");
  }

  return value;
}

} // namespace

uint32_t
str_to_u32(std::string_view s)
{
  return str_to<uint32_t>(s);
}

double
str_to_double(const std::string& s)
{
  // WORKAROUND: Should use `str_to`, but that is not supported by older Clang
  double value = 0;
  std::istringstream stream(s);
  if (!(stream >> value)) {
    throw std::invalid_argument("value is not a valid number");
  }

  char trailing = 0;
  if (stream >> trailing) {
    throw std::invalid_argument("number contains trailing text");
  }

  return value;
}

std::string
to_lower(std::string str)
{
  for (auto& current : str) {
    current = to_lower(current);
  }

  return str;
}

std::string
to_upper(std::string str)
{
  for (auto& current : str) {
    current = to_upper(current);
  }

  return str;
}

bool
starts_with(const std::string_view str1, const std::string_view str2)
{
  if (str1.size() < str2.size()) {
    return false;
  }

  return str1.substr(0, str2.size()) == str2;
}

bool
ends_with(const std::string_view str1, const std::string_view str2)
{
  if (str1.size() < str2.size()) {
    return false;
  }

  return str1.substr(str1.size() - str2.size()) == str2;
}

string_vec
split_text(const std::string_view text, const char separator)
{
  string_vec lines;

  size_t start = 0;
  size_t end = std::string::npos;
  do {
    end = text.find(separator, start);

    lines.push_back(std::string{ text.substr(start, end - start) });

    start = end + 1;
  } while (end != std::string::npos);

  return lines;
}

std::string
indent_lines(const std::string_view lines, size_t indentation)
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

std::string_view
trim_ascii_whitespace(std::string_view s)
{
  constexpr std::string_view whitespace = " \f\n\r\t\v";
  auto pos = s.find_first_not_of(whitespace);
  if (pos == std::string_view::npos) {
    return std::string_view{};
  }

  s.remove_prefix(pos);
  pos = s.find_last_not_of(whitespace);
  // Assumes `pos != npos` since `s` is guaranteed to contain non-whitespace
  s.remove_suffix(s.size() - pos - 1);

  return s;
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
cli_formatter::format(const std::string_view value) const
{
  std::ostringstream lines_out;

  for (const auto& line : split_lines(value)) {
    const auto block = wrap_text(line, m_columns, m_ljust);

    lines_out << join_text(indent(block, m_indentation, m_indent_first), "\n");
  }

  return lines_out.str();
}

std::string
shell_escape(const std::string_view s)
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

  return std::string{ s };
}

std::string
log_escape(const std::string_view s)
{
  std::string out;
  out.push_back('\'');

  constexpr std::string_view HEX{ "0123456789abcdef" };

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

          out.append("\\x");
          out.push_back(HEX.at((static_cast<uint8_t>(c) & 0xf0) >> 4));
          out.push_back(HEX.at((static_cast<uint8_t>(c) & 0xf)));

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
html_escape(std::string_view s)
{
  std::string out;

  for (const auto c : s) {
    switch (c) {
      case '\'':
        out.append("&#39;");
        break;
      case '"':
        out.append("&quot;");
        break;
      case '&':
        out.append("&amp;");
        break;
      case '<':
        out.append("&lt;");
        break;
      case '>':
        out.append("&gt;");
        break;
      default:
        out.push_back(c);
    }
  }

  return out;
}

std::string
format_thousand_sep(size_t count)
{

  std::string ss;
  if (count) {
    for (size_t i = 0; count; ++i) {
      ss.push_back('0' + (count % 10));
      count /= 10;
      if (count && i && i % 3 == 2) {
        ss.push_back(',');
      }
    }

    std::reverse(ss.begin(), ss.end());
  } else {
    ss.push_back('0');
  }

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
    return stringify(value);
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
  if (denom) {
    return format_fraction(num * 100, denom, precision) + " %";
  } else {
    return "NA";
  }
}

} // namespace adapterremoval
