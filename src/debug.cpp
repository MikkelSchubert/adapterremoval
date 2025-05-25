// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "debug.hpp" // declarations
#include <sstream>   // for ostringstream
#include <string>    // for string

namespace adapterremoval {

[[noreturn]] void
terminate(const std::string& message);

void
debug_raise_assert(std::string_view funcname,
                   std::string_view filename,
                   unsigned lineno,
                   std::string_view test,
                   std::string_view msg)
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
