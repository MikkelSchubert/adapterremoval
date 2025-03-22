// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

#include "errors.hpp" // declarations
#include <array>      // for array
#include <cerrno>     // for errno
#include <cstdio>     // for sys_errlist, sys_nerr
#include <cstring>    // for strerror_r
#include <sstream>    // for operator<<, basic_ostream
#include <string>     // for string

namespace adapterremoval {

// The function signature of strerror_r depends on the libc used. Rather than
// attempting to match defines, simply select the appropriate wrapper based on
// the signature of the function
using glibc_strerror_r = char* (*)(int, char*, size_t);
using posix_strerror_r = int (*)(int, char*, size_t);
//! Recommended buffer size according to `strerror` man page
using strerror_buf = std::array<char, 1024>;

std::string
strerror_r_wrapper(glibc_strerror_r f, int error_number)
{
  strerror_buf buf{};
  return f(error_number, buf.data(), buf.size());
}

std::string
strerror_r_wrapper(posix_strerror_r f, int error_number)
{
  strerror_buf buf{};
  if (f(error_number, buf.data(), buf.size())) {
    return "unknown error " + std::to_string(error_number);
  }

  return buf.data();
}

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

} // namespace adapterremoval
