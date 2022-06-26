/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert, et al. (2016). AdapterRemoval v2: rapid adapter trimming,   *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <chrono>   // for system_clock
#include <iomanip>  // for put_time
#include <iostream> // for cerr, cout, endl
#include <mutex>    // for mutex, unique_lock

#include "debug.hpp"   // for AR_FAIL
#include "logging.hpp" // declarations

namespace adapterremoval {

namespace log {

namespace {

//! Shared mutex for STDOUT / STDERR
std::mutex g_log_mutex;
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

/** Returns the color assosiated with a given log level */
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
  auto now = std::chrono::system_clock::now();
  auto in_time_t = std::chrono::system_clock::to_time_t(now);

  std::ostringstream prefixss;
  if (g_log_timestamps) {
    prefixss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X") << " ";
  }

  if (colors) {
    prefixss << "[" << level_color(l) << l << color::reset << "] ";
  } else {
    prefixss << "[" << l << "] ";
  }

  return prefixss.str();
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

void
set_level(level l)
{
  std::unique_lock<std::mutex> lock(g_log_mutex);

  g_log_level = l;
}

void
set_colors(bool enabled)
{
  std::unique_lock<std::mutex> lock(g_log_mutex);

  g_log_colors = enabled;
}

void
set_timestamps(bool enabled)
{
  std::unique_lock<std::mutex> lock(g_log_mutex);

  g_log_timestamps = enabled;
}

////////////////////////////////////////////////////////////////////////////////

log_stream::log_stream(level lvl)
  : m_level(lvl)
  , m_transient()
  , m_stream()
{
}

log_stream::~log_stream()
{
  std::unique_lock<std::mutex> lock(g_log_mutex);

  if (g_log_transient) {
    *g_log_out << "\r\033[K";
    g_log_transient = false;
  }

  if (m_level == level::none) {
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

    size_t start = 0;
    size_t end = std::string::npos;
    do {
      end = msg.find('\n', start);

      *g_log_out << head << msg.substr(start, end - start) << "\n";
      start = end + 1;
    } while (end != std::string::npos);

    g_log_out->flush();
  }
}

log_stream&
log_stream::transient()
{
  AR_REQUIRE(m_level == level::none);
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
  std::unique_lock<std::mutex> lock(g_log_mutex);

  AR_REQUIRE(g_log_out == &std::cerr);
  g_log_out = &m_stream;

  m_level = g_log_level;
  m_timestamps = g_log_timestamps;
  m_colors = g_log_colors;
}

log_capture::~log_capture()
{
  std::unique_lock<std::mutex> lock(g_log_mutex);

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
