/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "debug.hpp"   // for  AR_REQUIRE, AR_REQUIRE_SINGLE_TH...
#include "errors.hpp"  // for assert_failed
#include "testing.hpp" // for catch.hpp, StringMaker
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
  message << "Assertion '1 == 2' failed in " << __FUNCTION__ << " at "
          << __FILE__ << ":" << __LINE__ - 7;

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
  message << "Assertion '!!false' failed in " << __FUNCTION__ << " at "
          << __FILE__ << ":" << __LINE__ - 7 << ": message goes here";

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
  message << "Assertion failed in " << __FUNCTION__ << " at " << __FILE__ << ":"
          << __LINE__ - 7 << ": big fail";

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
