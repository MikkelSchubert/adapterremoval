// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2017 Mikkel Schubert <mikkelsch@gmail.com>
#define CATCH_CONFIG_MAIN

#include "debug.hpp"   // for terminate_on_assert
#include "errors.hpp"  // for assert_failed
#include "logging.hpp" // for level, set_level
#include "testing.hpp" // for TEST_CASE, REQUIRE, ...
#include <string_view> // for string_view

namespace adapterremoval {

namespace {

/** Resets log settings, setting log level to `none` before each test is run */
struct ConfigureLogListener : Catch::TestEventListenerBase
{
  using TestEventListenerBase::TestEventListenerBase;

  void testCaseStarting(const Catch::TestCaseInfo& /* testInfo */) override
  {
    log::set_level(log::level::none);
    log::set_preamble({});
    log::set_colors(false);
    log::set_timestamps(true);
  }
};

} // namespace

CATCH_REGISTER_LISTENER(ConfigureLogListener)

[[noreturn]] void
terminate_on_bug(std::string_view message)
{
  throw assert_failed(std::string{ message });
}

[[noreturn]] void
terminate_on_failure(std::string_view message)
{
  throw assert_failed(std::string{ message });
}

TEST_CASE("terminate")
{
  REQUIRE_THROWS_MESSAGE(terminate_on_bug("foobar"), assert_failed, "foobar");
  REQUIRE_THROWS_MESSAGE(terminate_on_failure("foobar"),
                         assert_failed,
                         "foobar");
}

} // namespace adapterremoval
