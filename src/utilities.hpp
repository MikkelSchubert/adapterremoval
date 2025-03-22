// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2021 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <array>       // for array
#include <cstddef>     // for size_t
#include <cstdint>     // for uint32_t
#include <memory>      // for allocator
#include <type_traits> // for enable_if_t, is_floating_point, is_integral

namespace adapterremoval {

/** Returns a seed value for a PRNG; not intended to be strongly random. */
uint32_t
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

} // namespace adapterremoval
