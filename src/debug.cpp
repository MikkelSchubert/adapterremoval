/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
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
#include <cstdlib> // for abort, size_t
#include <sstream> // for ostringstream
#include <string>  // for string

#include "debug.hpp"

namespace adapterremoval {

[[noreturn]] void
terminate(const std::string& message);

assert_failed::assert_failed(const assert_failed& errror)
  : m_what(errror.m_what)
{
}

assert_failed::assert_failed(const std::string& what)
  : m_what(what)
{
}

assert_failed::~assert_failed() {}

const char*
assert_failed::what() const noexcept
{
  return m_what.c_str();
}

void
debug_raise_assert(const char* filename,
                   size_t lineno,
                   const std::string& test,
                   const std::string& msg)
{
  std::ostringstream message;
  message << "Assertion failed at " << filename << ":" << lineno << ": ";
  if (test.empty() || msg.empty()) {
    message << test << msg;
  } else {
    message << msg << ": " << test;
  }

  terminate(message.str());
}

} // namespace adapterremoval
