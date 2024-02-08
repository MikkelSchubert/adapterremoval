/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include "debug.hpp"
#include <cstdint> // for int64_t, uint64_t
#include <memory>  // for allocator
#include <numeric> // for accumulate
#include <vector>  // for vector

namespace adapterremoval {

/** Calculates the arithmetic mean */
template<template<typename, typename> class C,
         typename T,
         typename A = std::allocator<T>>
double
arithmetic_mean(const C<T, A>& values)
{
  AR_REQUIRE(values.size());

  // Summation using https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  double sum = 0.0;
  double acc = 0.0;
  for (const double v : values) {
    double y = v - acc;
    double t = sum + y;
    acc = (t - sum) - y;
    sum = t;
  }

  return sum / static_cast<double>(values.size());
}

/** Calculates the sample standard deviation */
template<template<typename, typename> class C,
         typename T,
         typename A = std::allocator<T>>
double
standard_deviation(const C<T, A>& values)
{
  AR_REQUIRE(values.size() > 1);

  double error = 0;
  const auto m = arithmetic_mean(values);
  for (const auto it : values) {
    error += (it - m) * (it - m);
  }

  return std::sqrt(error / (values.size() - 1));
}

} // namespace adapterremoval
