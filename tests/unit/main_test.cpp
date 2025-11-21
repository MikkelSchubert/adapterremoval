// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2017 Mikkel Schubert <mikkelsch@gmail.com>
#define CATCH_CONFIG_MAIN

#include "errors.hpp"  // for assert_failed
#include "testing.hpp" // for TEST_CASE, REQUIRE, ...
#include <string_view> // for string_view

namespace adapterremoval {

[[noreturn]] void
terminate(std::string_view message) // NOLINT(misc-use-internal-linkage)
{
  throw assert_failed(std::string{ message });
}

TEST_CASE("terminate")
{
  REQUIRE_THROWS_MESSAGE(terminate("foobar"), assert_failed, "foobar");
}

} // namespace adapterremoval
