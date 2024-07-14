/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "logging.hpp"
#include "debug.hpp"    // for AR_REQUIRE, AR_FAIL
#include "main.hpp"     // for NAME, VERSION
#include "strutils.hpp" // for split_lines, cli_formatter, string_vec
#include <algorithm>    // for max, min
#include <iomanip>      // for operator<<, put_time
#include <iostream>     // for cerr
#include <limits>       // for numeric_limits
#include <mutex>        // for mutex, unique_lock
#include <sys/ioctl.h>  // for ioctl, winsize, TIOCGWINSZ
#include <unistd.h>     // for size_t, STDERR_FILENO
#include <vector>       // for vector

namespace adapterremoval {

namespace log {

namespace {

//! Shared mutex for STDOUT / STDERR
std::recursive_mutex g_log_mutex;
//! Pointer to logging stream; normally this will be std::cerr;
std::ostream* g_log_out = &std::cerr;

//! Minimum log-level to print
level g_log_level = level::debug;
//! Include timestamps when writing log lines
bool g_log_timestamps = true;
//! Use ANSI colors when writing log lines
bool g_log_colors = false;
//! Indicates if the last message was transient
bool g_log_transient = false;
//! Indicates that the preamble has been logged
bool g_preamble_logged = false;

//! ANSI color codes: https://en.wikipedia.org/wiki/ANSI_escape_code#Colors
enum class color
{
  black = 30,
  red,
  green,
  yellow,
  blue,
  magenta,
  cyan,
  white,
  reset
};

/** Returns the color associated with a given log level */
color
level_color(level l)
{
  switch (l) {
    case level::debug:
      return color::cyan;
    case level::info:
      return color::green;
    case level::warning:
      return color::yellow;
    case level::error:
      return color::red;
    default:
      AR_FAIL("invalid logging level");
  }
}

/** Writes the name of a given log level */
std::ostream&
operator<<(std::ostream& out, level l)
{
  switch (l) {
    case level::debug:
      return out << "DEBUG";
    case level::info:
      return out << "INFO";
    case level::warning:
      return out << "WARNING";
    case level::error:
      return out << "ERROR";
    default:
      AR_FAIL("invalid logging level");
  }
}

/** Writes a color control code */
std::ostream&
operator<<(std::ostream& out, color c)
{
  if (c == color::reset) {
    out << "\033[0m";
  } else {
    out << "\033[0;" << static_cast<int>(c) << "m";
  }

  return out;
}

/** Create timestamp+level prefix for logs */
std::string
log_header(level l, bool colors = false)
{
  std::ostringstream header;
  if (g_log_timestamps) {
    header << timestamp("%Y-%m-%d %X", true) << " ";
  }

  if (colors) {
    header << "[" << level_color(l) << l << color::reset << "] ";
  } else {
    header << "[" << l << "] ";
  }

  return header.str();
}

size_t
log_linewidth(const std::ostream& out)
{
  // Piped logs are not pretty-printed, to make analyses easier
  if (&out == &std::cerr) {
    return get_terminal_width();
  }

  return std::numeric_limits<size_t>::max();
}

void
log_print_lines(std::ostream& out,
                const std::string& head,
                const std::string& msg)
{
  // Account for unprinted color codes:
  // Size of color ("\033[0;XXm" = 7) + reset ("\033[0m" = 4)
  const int color_width = g_log_colors ? 11 : 0;
  const int head_width = head.size() - color_width;
  const int line_width = log_linewidth(out);

  for (const auto& line : split_lines(msg)) {
    const auto indent = line.find_first_not_of(' ');
    if (indent == std::string::npos) {
      out << head << "\n";
      continue;
    }

    cli_formatter fmt;
    fmt.set_ljust(2);
    fmt.set_indent(indent);
    fmt.set_column_width(line_width - indent - head_width);

    bool first_line = true;
    for (const auto& fragment : split_lines(fmt.format(line))) {
      if (first_line) {
        out << head + fragment << "\n";
        first_line = false;
      } else {
        out << std::string(head_width, ' ') << fragment << "\n";
      }
    }
  }
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

void
log_preamble()
{
  std::unique_lock<std::recursive_mutex> lock(g_log_mutex);
  if (!g_preamble_logged) {
    g_preamble_logged = true;
    log::info() << NAME << " " << VERSION;
  }
}

void
set_level(level l)
{
  std::unique_lock<std::recursive_mutex> lock(g_log_mutex);

  g_log_level = l;
}

void
set_colors(bool enabled)
{
  std::unique_lock<std::recursive_mutex> lock(g_log_mutex);

  g_log_colors = enabled;
}

void
set_timestamps(bool enabled)
{
  std::unique_lock<std::recursive_mutex> lock(g_log_mutex);

  g_log_timestamps = enabled;
}

size_t
get_terminal_width()
{
  struct winsize params = {};
  // Attempt to retrieve the number of columns in the terminal
  if (ioctl(STDERR_FILENO, TIOCGWINSZ, &params) == 0) {
    return std::min<size_t>(120, std::max<size_t>(60, params.ws_col));
  }

  return std::numeric_limits<size_t>::max();
}

////////////////////////////////////////////////////////////////////////////////

log_stream::log_stream(level lvl)
  : m_level(lvl)
  , m_transient()
  , m_stream()
{
  AR_REQUIRE(lvl < level::none);
}

log_stream::~log_stream()
{
  std::unique_lock<std::recursive_mutex> lock(g_log_mutex);
  log_preamble();

  if (g_log_transient) {
    *g_log_out << "\r\033[K";
    g_log_transient = false;
  }

  if (m_level == level::cerr) {
    *g_log_out << m_stream.str();
    g_log_out->flush();
    g_log_transient = m_transient;
  } else if (m_level >= g_log_level) {
    const auto head = log_header(m_level, g_log_colors);

    // Trailing newlines are permitted for ease of use, e.g. when writing
    // complex, multi-line log messages
    auto msg = m_stream.str();
    while (msg.size() && msg.back() == '\n') {
      msg.pop_back();
    }

    log_print_lines(*g_log_out, head, msg);

    g_log_out->flush();
  }
}

log_stream&
log_stream::transient()
{
  AR_REQUIRE(m_level == level::cerr);
  m_transient = true;

  return *this;
}

////////////////////////////////////////////////////////////////////////////////

log_capture::log_capture()
  : m_level()
  , m_colors()
  , m_timestamps()
  , m_stream()
{
  std::unique_lock<std::recursive_mutex> lock(g_log_mutex);

  AR_REQUIRE(g_log_out == &std::cerr);
  g_log_out = &m_stream;

  m_level = g_log_level;
  m_timestamps = g_log_timestamps;
  m_colors = g_log_colors;
}

log_capture::~log_capture()
{
  std::unique_lock<std::recursive_mutex> lock(g_log_mutex);

  AR_REQUIRE(g_log_out == &m_stream);
  g_log_out = &std::cerr;
  g_log_level = m_level;
  g_log_timestamps = m_colors;
  g_log_colors = m_colors;
}

std::string
log_capture::str() const
{
  return m_stream.str();
}

} // namespace log

} // namespace adapterremoval
