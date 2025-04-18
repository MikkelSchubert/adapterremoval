// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <string_view> // for string_view

namespace adapterremoval {

/**
 * Aborts after printing the filename, line-number, and message, plus
 * instructions for how to report the problem.
 */
[[noreturn]] void
debug_raise_assert(std::string_view funcname,
                   std::string_view filename,
                   unsigned lineno,
                   std::string_view test,
                   std::string_view msg);

#ifdef __has_builtin
#if __has_builtin(__builtin_expect)
#define AR_LIKELY(condition) __builtin_expect(!!(condition), 1)
#define AR_UNLIKELY(condition) __builtin_expect(!!(condition), 0)
#endif
#endif

#ifndef AR_LIKELY
#define AR_LIKELY(condition) (!!(condition))
#define AR_UNLIKELY(condition) (!!(condition))
#endif

/** Custom assert which prints various information on failure; always enabled */
#define AR_REQUIRE_2_(test, msg)                                               \
  do {                                                                         \
    /* NOLINTNEXTLINE(readability-simplify-boolean-expr) */                    \
    if (AR_UNLIKELY(!(test))) {                                                \
      debug_raise_assert(__PRETTY_FUNCTION__, __FILE__, __LINE__, #test, msg); \
    }                                                                          \
  } while (0)

#define AR_REQUIRE_1_(test) AR_REQUIRE_2_(test, {})

#define AR_REQUIRE_GET_(_1, _2, NAME, ...) NAME
#define AR_REQUIRE(...)                                                        \
  AR_REQUIRE_GET_(__VA_ARGS__, AR_REQUIRE_2_, AR_REQUIRE_1_, )(__VA_ARGS__)

/** Raise an assert failure with a user-specified message. */
#define AR_FAIL(msg)                                                           \
  adapterremoval::debug_raise_assert(__FUNCTION__, __FILE__, __LINE__, {}, msg)

#define AR_MERGE1_(a, b) a##b
#define AR_MERGE_(a, b) AR_MERGE1_(a, b)

/** Raise a failure if a scope is accessed more than once at the same time. */
#define AR_REQUIRE_SINGLE_THREAD(lock)                                         \
  std::unique_lock<std::mutex> AR_MERGE_(locker, __LINE__)(lock,               \
                                                           std::defer_lock);   \
  do {                                                                         \
    if (AR_UNLIKELY(!AR_MERGE_(locker, __LINE__).try_lock())) {                \
      AR_FAIL("race condition detected");                                      \
    }                                                                          \
  } while (0)

} // namespace adapterremoval
