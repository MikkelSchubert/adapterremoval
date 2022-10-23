/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2021 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include <memory>   // for allocator
#include <stdint.h> // for uint32_t

namespace adapterremoval {

/** Returns a seed value for a PRNG; not intended to be strongly random. */
uint32_t
prng_seed();

template<typename A>
typename std::enable_if<std::is_integral<A>::value ||
                        std::is_floating_point<A>::value>::type
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
