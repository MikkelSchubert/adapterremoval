// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2017 Mikkel Schubert <mikkelsch@gmail.com>
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "debug.hpp"  // for assert_failed
#include "errors.hpp" // for assert_failed
#include <string>     // for string

namespace adapterremoval {

[[noreturn]] void
terminate(const std::string& message)
{
  throw assert_failed(message);
}

} // namespace adapterremoval
