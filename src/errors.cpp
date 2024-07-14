/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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
