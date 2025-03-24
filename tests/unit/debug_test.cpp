// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "debug.hpp"   // for  AR_REQUIRE, AR_REQUIRE_SINGLE_TH...
#include "errors.hpp"  // for assert_failed
#include "testing.hpp" // for TEST_CASE, REQUIRE, ...
#include <mutex>       // for unique_lock, mutex
#include <sstream>     // for operator<<, basic_ostream, char_traits, ostring...
#include <string>      // for basic_string, operator==, string

namespace adapterremoval {

TEST_CASE("assert on true value")
{
  AR_REQUIRE(1 == 1);
}

TEST_CASE("assert on false value")
{
  std::string what;
  try {
    AR_REQUIRE(1 == 2);
  } catch (const assert_failed& error) {
    what = error.what();
  }

  std::ostringstream message;
  message << "Assertion '1 == 2' failed in "
          << static_cast<const char*>(__PRETTY_FUNCTION__) << " at " << __FILE__
          << ":" << __LINE__ - 8;

  REQUIRE(what == message.str());
}

TEST_CASE("assert on true value with message")
{
  AR_REQUIRE(!false, "this is dummy message");
}

TEST_CASE("assert on false value with message")
{
  std::string what;
  try {
    AR_REQUIRE(!!false, "message goes here");
  } catch (const assert_failed& error) {
    what = error.what();
  }

  std::ostringstream message;
  message << "Assertion '!!false' failed in "
          << static_cast<const char*>(__PRETTY_FUNCTION__) << " at " << __FILE__
          << ":" << __LINE__ - 8 << ": message goes here";

  REQUIRE(what == message.str());
}

TEST_CASE("assert fail")
{
  std::string what;
  try {
    AR_FAIL("big fail");
  } catch (const assert_failed& error) {
    what = error.what();
  }

  std::ostringstream message;
  message << "Assertion failed in " << static_cast<const char*>(__FUNCTION__)
          << " at " << __FILE__ << ":" << __LINE__ - 7 << ": big fail";

  REQUIRE(what == message.str());
}

TEST_CASE("assert single thread")
{
  std::mutex lock;

  {
    AR_REQUIRE_SINGLE_THREAD(lock);
  }

  {
    AR_REQUIRE_SINGLE_THREAD(lock);
  }
}

TEST_CASE("assert single thread failed")
{
  std::mutex lock;

  AR_REQUIRE_SINGLE_THREAD(lock);

  bool caught_exception = false;
  try {
    AR_REQUIRE_SINGLE_THREAD(lock);
  } catch (const assert_failed&) {
    caught_exception = true;
  }

  REQUIRE(caught_exception);
}

} // namespace adapterremoval
