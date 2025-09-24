// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2017 Mikkel Schubert <mikkelsch@gmail.com>
#define CATCH_CONFIG_MAIN

#include "errors.hpp"  // for assert_failed
#include "testing.hpp" // for TEST_CASE, REQUIRE, ...
#include <string>      // for string

namespace adapterremoval {

[[noreturn]] void
terminate(const std::string& message)
{
  throw assert_failed(message);
}

TEST_CASE("terminate")
{
  REQUIRE_THROWS_MESSAGE(terminate("foobar"), assert_failed, "foobar");
}

} // namespace adapterremoval
