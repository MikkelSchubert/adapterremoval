// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <mutex>       // for recursive_mutex
#include <string_view> // for string_view

namespace adapterremoval {

/** Terminate on internal error; prints how to file issue on github */
[[noreturn]] void
terminate_on_bug(std::string_view message);

/** Terminate on external error; prints warning not to use program output */
[[noreturn]] void
terminate_on_failure(std::string_view message);

/** Terminate on unhandled exception */
[[noreturn]] void
terminate_on_exception();

/** Terminate on failed assert; prints location and how to file issue */
[[noreturn]] void
terminate_on_assert(std::string_view funcname,
                    std::string_view filename,
                    unsigned lineno,
                    std::string_view test,
                    std::string_view msg);

#ifdef __has_builtin
#if __has_builtin(__builtin_expect)
#define AR_LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define AR_UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)
#endif
#endif

#ifndef AR_LIKELY
#define AR_LIKELY(condition) (!!(condition))
#define AR_UNLIKELY(condition) (!!(condition))
#endif

/** Custom assert which prints various information on failure; always enabled */
#define AR_REQUIRE_2_(test, msg)                                               \
  do {                                                                         \
    if (!AR_LIKELY(test)) {                                                    \
      ::adapterremoval::terminate_on_assert(__PRETTY_FUNCTION__,               \
                                            __FILE__,                          \
                                            __LINE__,                          \
                                            #test,                             \
                                            msg);                              \
    }                                                                          \
  } while (0)

#define AR_REQUIRE_1_(test) AR_REQUIRE_2_(test, {})

#define AR_REQUIRE_GET_(_1, _2, NAME, ...) NAME
#define AR_REQUIRE(...)                                                        \
  AR_REQUIRE_GET_(__VA_ARGS__, AR_REQUIRE_2_, AR_REQUIRE_1_, )(__VA_ARGS__)

/** Raise an assert failure with a user-specified message. */
#define AR_FAIL(msg)                                                           \
  ::adapterremoval::terminate_on_assert(__PRETTY_FUNCTION__,                   \
                                        __FILE__,                              \
                                        __LINE__,                              \
                                        {},                                    \
                                        msg)

#define AR_MERGE1_(a, b) a##b
#define AR_MERGE_(a, b) AR_MERGE1_(a, b)

/** Raise a failure if scope is accessed by multiple threads simultaneously */
#define AR_REQUIRE_SINGLE_THREAD(lock)                                         \
  std::unique_lock<std::recursive_mutex> AR_MERGE_(locker,                     \
                                                   __LINE__)(lock,             \
                                                             std::defer_lock); \
  do {                                                                         \
    if (!AR_LIKELY(AR_MERGE_(locker, __LINE__).try_lock())) {                  \
      AR_FAIL("race condition detected");                                      \
    }                                                                          \
  } while (0);                                                                 \
  /* Ensure that this macro cannot be placed after if or while without {} */   \
  [[maybe_unused]] int AR_MERGE_(_braces_are_required, __LINE__)

} // namespace adapterremoval
