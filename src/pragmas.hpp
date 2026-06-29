// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#ifdef __has_builtin
#if __has_builtin(__builtin_expect)
//! Indicate that a condition is likely to happen
#define AR_LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
//! Indicate that a condition is unlikely to happen
#define AR_UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)
#endif
#endif

#ifndef AR_LIKELY
#define AR_LIKELY(condition) (!!(condition))
#define AR_UNLIKELY(condition) (!!(condition))
#endif

//! Helper macro to allow application of pragmas
#define _APPLY_PRAGMA(x) _Pragma(#x)

//! Unroll a loop n times
#ifdef __clang__
#define AR_UNROLL(n) _APPLY_PRAGMA(unroll n)
#elif defined(__GNUC__)
#define AR_UNROLL(n) _APPLY_PRAGMA(GCC unroll n)
#else
#define AR_UNROLL(n) // no-op fallback
#endif
