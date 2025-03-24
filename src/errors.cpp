// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"   // declarations
#include "strutils.hpp" // for log_escape
#include <array>        // for array
#include <cerrno>       // for errno
#include <cstdio>       // for sys_errlist, sys_nerr
#include <cstring>      // for strerror_r
#include <exception>    // for exception
#include <ostream>      // for ostream
#include <string>       // for string
#include <string_view>  // for string_view

namespace adapterremoval {

// The function signature of strerror_r depends on the libc used. Rather than
// attempting to match defines, simply select the appropriate wrapper based on
// the signature of the function
using glibc_strerror_r = char* (*)(int, char*, size_t);
using posix_strerror_r = int (*)(int, char*, size_t);
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
    return "unknown error " + std::to_string(error_number);
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
format_io_error(const std::string& message, const int error_number)
{
  if (error_number) {
    std::ostringstream stream;
    stream << message << " ('" << strerror_r_wrapper(::strerror_r, error_number)
           << "')";

    return stream.str();
  } else {
    return message;
  }
}

assert_failed::assert_failed(const std::string& what)
  : std::logic_error(what)
{
}

io_error::io_error(const std::string& message, int error_number)
  : std::ios_base::failure(
      message,
      std::error_code(error_number, std::generic_category()))
{
}

gzip_error::gzip_error(const std::string& message)
  : io_error(message)
{
}

parsing_error::parsing_error(const std::string& message)
  : std::runtime_error(message)
{
}

fastq_error::fastq_error(const std::string& message)
  : parsing_error(message)
{
}

std::ostream&
operator<<(std::ostream& os, const assert_failed& value)
{
  return format_exception(os, "assert_failed", value);
}

std::ostream&
operator<<(std::ostream& os, const io_error& value)
{
  return format_exception(os, "io_error", value);
}

std::ostream&
operator<<(std::ostream& os, const gzip_error& value)
{
  return format_exception(os, "gzip_error", value);
}

std::ostream&
operator<<(std::ostream& os, const parsing_error& value)
{
  return format_exception(os, "parsing_error", value);
}

std::ostream&
operator<<(std::ostream& os, const fastq_error& value)
{
  return format_exception(os, "fastq_error", value);
}

} // namespace adapterremoval
