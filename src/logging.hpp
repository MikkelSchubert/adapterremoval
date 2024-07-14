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
#pragma once

#include <sstream> // for ostringstream
#include <string>  // for string

namespace adapterremoval {

namespace log {

/** Log levels */
enum class level
{
  // Debugging messages; used only during development
  debug,
  // General information messages; feature flags, progress, etc.
  info,
  // Warning messages for mostly harmless things; duplicate arguments, etc.
  warning,
  // Error messages for bugs, bad data, likely mistakes, and the like. Errors
  // should always be followed by the program terminating.
  error,
  // Non-log messages; this is for raw text written to stderr
  cerr,
  // Disables log messages when used with log::set_level; not usable or messages
  none
};

/**
 * Helper class for writing complex log messages.
 *
 * The resulting message is written on destruction. Multi-line messages are
 * split across multiple log-lines and any trailing newlines are stripped.
 */
class log_stream
{
public:
  /* Crates a log stream with the specified log level */
  explicit log_stream(level lvl);
  /** Writes the message if the corresponding log-level is enabled */
  ~log_stream();

  /**
   * Indicates that a message is transient.
   *
   * The message is assumed to contain no new-lines. In addition, this message
   * will be cleared before any further messages are printed. This can only be
   * used with level::none.
   */
  log_stream& transient();

  /** Writes value to a cache */
  template<typename T>
  log_stream& operator<<(const T& value)
  {
    m_stream << value;

    return *this;
  }

  /** Not implemented */
  log_stream(log_stream&&) noexcept = default;
  /** Not implemented */
  log_stream(const log_stream&) = delete;
  /** Not implemented */
  log_stream& operator=(const log_stream&) = delete;

private:
  //! Log level of the current message
  const level m_level;
  //! Indicates if the message is transient
  bool m_transient;
  //! Stream for caching log output prior to writing
  std::ostringstream m_stream;
};

/**
 * Helper class that captures any output that would otherwise have been written
 * to STDERR. Only one instance of `log_capture` is allowed to exist at any one
 * time. For use during testing.
 */
class log_capture
{
public:
  /** Starts capturing log messages; aborts if another instance exists. */
  log_capture();
  /** Restores output to stderr and restores previous settings. */
  ~log_capture();

  /** Returns all output captured log messages / output written to cerr. */
  std::string str() const;

  /** Not implemented */
  log_capture(log_capture&&) = default;
  /** Not implemented */
  log_capture(const log_capture&) = delete;
  /** Not implemented */
  log_capture& operator=(const log_capture&) = delete;

private:
  //! Original log level
  level m_level;
  //! Original colors setting
  bool m_colors;
  //! Original timestamps setting
  bool m_timestamps;
  //! Stream containing text written using `log_stream`
  std::ostringstream m_stream;
};

/**
 * Writes to stderr (or log_capture); unlike other logging streams, `cerr`
 *does not write a timestamp or a log-level, nor does it add a trailing
 *newline.
 **/
inline log_stream
cerr()
{
  return log_stream(level::cerr);
}

/** Logs a debug-level message; only output when verbose logging is enabled */
inline log_stream
debug()
{
  return log_stream(level::debug);
}

/** Logs a info-level message; not output if quiet logging is enabled */
inline log_stream
info()
{
  return log_stream(level::info);
}

/** Logs a warning message; is always output */
inline log_stream
warn()
{
  return log_stream(level::warning);
}

/** Logs a error message; is always output */
inline log_stream
error()
{
  return log_stream(level::error);
}

/** Sets the minimum log level to print */
void
set_level(level l);

/** Enable/disable colors while logging (off by default) */
void
set_colors(bool enabled);

/** Enable/disable timestamps while logging (on by default) */
void
set_timestamps(bool enabled);

/** Returns the terminal column width if available, otherwise (size_t)-1 */
size_t
get_terminal_width();

} // namespace log

} // namespace adapterremoval
