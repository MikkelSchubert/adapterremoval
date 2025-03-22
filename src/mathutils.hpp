// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
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
