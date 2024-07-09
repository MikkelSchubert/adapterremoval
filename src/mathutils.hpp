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

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t
#include <vector>  // for vector

namespace adapterremoval {

/** Calculates the arithmetic mean */
double
arithmetic_mean(const std::vector<uint64_t>& values);

/** Calculates the sample standard deviation */
double
standard_deviation(const std::vector<uint64_t>& values);

/** Returns critical value for student's t distribution at 99.5% */
double
students_t_critical_value(size_t df);

/**
 * Prune extreme values from a vector using Grubb's test:
 *   https://en.wikipedia.org/wiki/Grubbs%27s_test
 *
 * Returns true if values were pruned.
 */
bool
grubbs_test_prune(std::vector<uint64_t>& values);

} // namespace adapterremoval
