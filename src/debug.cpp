/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "debug.hpp"
#include <sstream> // for operator<<, basic_ostream, ostringstream
#include <string>  // for char_traits, operator<<, string

namespace adapterremoval {

[[noreturn]] void
terminate(const std::string& message);

assert_failed::assert_failed(const std::string& what)
  : std::logic_error(what)
{
}

void
debug_raise_assert(const char* funcname,
                   const char* filename,
                   size_t lineno,
                   const std::string& test,
                   const std::string& msg)
{
  std::ostringstream message;

  message << "Assertion ";
  if (!test.empty()) {
    message << "'" << test << "' ";
  }

  message << "failed in " << funcname << " at " << filename << ":" << lineno;
  if (!msg.empty()) {
    message << ": " << msg;
  }

  terminate(message.str());
}

} // namespace adapterremoval
