/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#pragma once

#include <stddef.h> // for size_t

#include <string> // for string

/** Exception explaining 'abort' calls when running unit-tests. */
class assert_failed : public std::exception
{
public:
  /** Copy constructor. */
  assert_failed(const assert_failed& errror);
  /** Creates exception with the specified error message. */
  assert_failed(const std::string& what);

  /** Does nothing. */
  virtual ~assert_failed() override;

  /** Returns user supplied error message; owned by object. */
  virtual const char* what() const noexcept override;

private:
  //! User supplied error message
  const std::string m_what;
};

/**
 * Aborts after printing the filename, line-number, and message, plus
 * instructions for how to report the problem.
 */
[[noreturn]] void
debug_raise_assert(const char* filename,
                   size_t lineno,
                   const std::string& test,
                   const std::string& msg);

/** Custom assert which prints various information on failure; always enabled.
 */
#define AR_REQUIRE_2_(test, msg)                                               \
  do {                                                                         \
    if (!(test)) {                                                             \
      debug_raise_assert(__FILE__, __LINE__, #test, msg);                      \
    }                                                                          \
  } while (0)

#define AR_REQUIRE_1_(test) AR_REQUIRE_2_(test, std::string())

#define AR_REQUIRE_GET_(_1, _2, NAME, ...) NAME
#define AR_REQUIRE(...)                                                        \
  AR_REQUIRE_GET_(__VA_ARGS__, AR_REQUIRE_2_, AR_REQUIRE_1_, )(__VA_ARGS__)

/** Raise an assert failure with a user-specified message. */
#define AR_FAIL(msg) debug_raise_assert(__FILE__, __LINE__, std::string(), msg)

#ifdef DEBUG

#define AR_MERGE1_(a, b) a##b
#define AR_MERGE_(a, b) AR_MERGE1_(a, b)

/** Raise a failure if a scope is accessed more than once at the same time. */
#define AR_ASSERT_SINGLE_THREAD(lock)                                          \
  std::unique_lock<std::mutex> AR_MERGE_(locker, __LINE__)(lock,               \
                                                           std::defer_lock);   \
  do {                                                                         \
    if (!AR_MERGE_(locker, __LINE__).try_lock()) {                             \
      AR_FAIL("race condition detected");                                      \
    }                                                                          \
  } while (0)

#define AR_ASSERT(...)                                                         \
  AR_REQUIRE_GET_(__VA_ARGS__, AR_REQUIRE_2_, AR_REQUIRE_1_, )(__VA_ARGS__)

#else

#define AR_ASSERT_SINGLE_THREAD(lock)
#define AR_ASSERT(...)

#endif
