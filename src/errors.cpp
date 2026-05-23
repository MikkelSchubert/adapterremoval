// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"   // declarations
#include "strutils.hpp" // for log_escape, stringify
#include <array>        // for array
#include <cstring>      // for strerror_r
#include <exception>    // for exception
#include <ostream>      // for ostream
#include <sstream>      // for ostringstream
#include <stdexcept>    // for logic_error, runtime_error
#include <string>       // for string
#include <string_view>  // for string_view

namespace adapterremoval {

// The function signature of strerror_r depends on the libc used. Rather than
// attempting to match defines, simply select the appropriate wrapper based on
// the signature of the function
using glibc_strerror_r = char* (*)(int, char*, size_t);
using posix_strerror_r = int (*)(int, char*, size_t);
using posix_strerror_s = int (*)(char*, size_t, int);

//! Recommended buffer size according to `strerror` man page
using strerror_buf = std::array<char, 1024>;

namespace {

[[maybe_unused]] std::string
strerror_r_wrapper(glibc_strerror_r f, int error_number)
{
  strerror_buf buf{};
  return f(error_number, buf.data(), buf.size());
}

[[maybe_unused]] std::string
strerror_r_wrapper(posix_strerror_r f, int error_number)
{
  strerror_buf buf{};
  if (f(error_number, buf.data(), buf.size())) {
    return "unknown error " + stringify(error_number);
  }

  return buf.data();
}

[[maybe_unused]] std::string
strerror_s_wrapper(posix_strerror_s f, int error_number)
{
  strerror_buf buf{};
  if (f(buf.data(), buf.size(), error_number)) {
    return "unknown error " + stringify(error_number);
  }

  return buf.data();
}

std::ostream&
format_exception(std::ostream& os,
                 std::string_view name,
                 const std::exception& value)
{
  return os << name << "{" << log_escape(value.what()) << "}";
}

} // namespace

std::string
format_io_error(std::string_view message, const int error_number)
{
  if (error_number) {
    std::ostringstream stream;
#ifdef _WIN32
    stream << message << ": " << strerror_s_wrapper(::strerror_s, error_number);
#else
    stream << message << ": " << strerror_r_wrapper(::strerror_r, error_number);
#endif

    return stream.str();
  } else {
    return std::string{ message };
  }
}

assert_failed::assert_failed(const std::string& what)
  : std::logic_error(what)
{
}

program_failure::program_failure(const std::string& message)
  : std::runtime_error(message)
{
}

io_error::io_error(const std::string& message)
  : program_failure(message)
{
}

io_error::io_error(const std::string& message, int error_number)
  : program_failure(format_io_error(message, error_number))
{
}

std::ostream&
operator<<(std::ostream& os, const assert_failed& value)
{
  return format_exception(os, "assert_failed", value);
}

std::ostream&
operator<<(std::ostream& os, const program_failure& value)
{
  return format_exception(os, value.kind(), value);
}

} // namespace adapterremoval
