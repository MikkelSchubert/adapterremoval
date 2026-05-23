// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "debug.hpp" // declarations
#include <sstream>   // for ostringstream

namespace adapterremoval {

[[noreturn]] void
terminate_on_assert(std::string_view funcname,
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

  terminate_on_bug(message.str());
}

} // namespace adapterremoval
