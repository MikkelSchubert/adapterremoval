/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <mutex>
#include <sstream>

#include "debug.hpp"
#include "testing.hpp"

TEST_CASE("assert on true value")
{
  AR_ASSERT(1 == 1);
}

TEST_CASE("assert on false value")
{
  std::string what;
  try {
    AR_ASSERT(1 == 2);
  } catch (const assert_failed& error) {
    what = error.what();
  }

  std::stringstream message;
  message << "Assertion failed at " << __FILE__ << ":" << __LINE__ - 6
          << ": 1 == 2";

  REQUIRE(what == message.str());
}

TEST_CASE("assert on true value with message")
{
  AR_ASSERT(!false, "this is dummy message");
}

TEST_CASE("assert on false value with message")
{
  std::string what;
  try {
    AR_ASSERT(!!false, "message goes here");
  } catch (const assert_failed& error) {
    what = error.what();
  }

  std::stringstream message;
  message << "Assertion failed at " << __FILE__ << ":" << __LINE__ - 6
          << ": message goes here: !!false";

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

  std::stringstream message;
  message << "Assertion failed at " << __FILE__ << ":" << __LINE__ - 6
          << ": big fail";

  REQUIRE(what == message.str());
}

TEST_CASE("assert single thread")
{
  std::mutex lock;

  {
    AR_ASSERT_SINGLE_THREAD(lock);
  }

  {
    AR_ASSERT_SINGLE_THREAD(lock);
  }
}

TEST_CASE("assert single thread failed")
{
  std::mutex lock;

  AR_ASSERT_SINGLE_THREAD(lock);

  bool caught_exception = false;
  try {
    AR_ASSERT_SINGLE_THREAD(lock);
  } catch (const assert_failed&) {
    caught_exception = true;
  }

#ifdef AR_DEBUG_BUILD
  REQUIRE(caught_exception);
#else
  REQUIRE_FALSE(caught_exception);
#endif
}
