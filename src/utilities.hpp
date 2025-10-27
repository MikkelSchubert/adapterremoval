// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2021 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <array>       // for array
#include <cstddef>     // for size_t
#include <memory>      // for allocator, unique_ptr
#include <type_traits> // for enable_if_t, is_floating_point, is_integral
#include <utility>     // for ignore

namespace adapterremoval {

/** Returns a seed value for a PRNG; not intended to be strongly random. */
unsigned int
prng_seed();

template<typename A>
typename std::enable_if_t<std::is_integral<A>::value ||
                          std::is_floating_point<A>::value>
merge(A& dst, const A& src)
{
  dst += src;
}

template<typename T, std::size_t N>
void
merge(std::array<T, N>& dst, const std::array<T, N>& src)
{
  auto dst_it = std::begin(dst);
  auto src_it = std::begin(src);
  while (src_it != std::end(src)) {
    merge(*dst_it++, *src_it++);
  }
}

template<template<typename, typename> class C,
         typename T,
         typename A = std::allocator<T>>
void
merge(C<T, A>& dst, const C<T, A>& src)
{
  if (dst.size() < src.size()) {
    dst.resize(src.size());
  }

  auto dst_it = std::begin(dst);
  auto src_it = std::begin(src);
  while (src_it != std::end(src)) {
    merge(*dst_it++, *src_it++);
  }
}

/** Returns the underlying value for an enum/enum class member */
template<typename T>
constexpr auto
underlying_value(T value)
{
  return static_cast<std::underlying_type_t<T>>(value);
}

/**
 * Perform dynamic cast on pointer stored in unique_ptr. If the cast succeeds,
 * overship of the pointer transferred is the new unique_ptr, otherwise the
 * source pointer remains unchanged.
 */
template<typename T, typename U>
std::unique_ptr<T, std::default_delete<T>>
dynamic_cast_unique(std::unique_ptr<U, std::default_delete<U>>& src)
{
  if (auto* dst = dynamic_cast<T*>(src.get())) {
    std::ignore = src.release();      // noexcept, silence warning
    return std::unique_ptr<T>{ dst }; // noexcept
  }

  return {};
}

#ifdef __clang__
#define NO_OPTIMIZE_CLANG __attribute__((optnone))
#define NO_OPTIMIZE_GCC
#else
#define NO_OPTIMIZE_CLANG
#define NO_OPTIMIZE_GCC __attribute__((optimize("O0")))
#endif

/** Unoptimized to prevent calculations from being elided by the compiler */
template<typename T>
void NO_OPTIMIZE_GCC
blackbox(T& /* unused */) NO_OPTIMIZE_CLANG
{
}

} // namespace adapterremoval
