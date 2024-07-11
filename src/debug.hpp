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
#pragma once

#include <cstddef>   // for size_t
#include <stdexcept> // for logic_error
#include <string>    // for string

namespace adapterremoval {

/**
 * Aborts after printing the filename, line-number, and message, plus
 * instructions for how to report the problem.
 */
[[noreturn]] void
debug_raise_assert(const char* funcname,
                   const char* filename,
                   size_t lineno,
                   const std::string& test,
                   const std::string& msg);

/** Custom assert which prints various information on failure; always enabled.
 */
#define AR_REQUIRE_2_(test, msg)                                               \
  do {                                                                         \
    if (!(test)) {                                                             \
      debug_raise_assert(__FUNCTION__, __FILE__, __LINE__, #test, msg);        \
    }                                                                          \
  } while (0)

#define AR_REQUIRE_1_(test) AR_REQUIRE_2_(test, std::string())

#define AR_REQUIRE_GET_(_1, _2, NAME, ...) NAME
#define AR_REQUIRE(...)                                                        \
  AR_REQUIRE_GET_(__VA_ARGS__, AR_REQUIRE_2_, AR_REQUIRE_1_, )(__VA_ARGS__)

/** Raise an assert failure with a user-specified message. */
#define AR_FAIL(msg)                                                           \
  adapterremoval::debug_raise_assert(                                          \
    __FUNCTION__, __FILE__, __LINE__, std::string(), msg)

#define AR_MERGE1_(a, b) a##b
#define AR_MERGE_(a, b) AR_MERGE1_(a, b)

/** Raise a failure if a scope is accessed more than once at the same time. */
#define AR_REQUIRE_SINGLE_THREAD(lock)                                         \
  std::unique_lock<std::mutex> AR_MERGE_(locker, __LINE__)(lock,               \
                                                           std::defer_lock);   \
  do {                                                                         \
    if (!AR_MERGE_(locker, __LINE__).try_lock()) {                             \
      AR_FAIL("race condition detected");                                      \
    }                                                                          \
  } while (0)

} // namespace adapterremoval
