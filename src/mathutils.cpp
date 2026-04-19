// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "mathutils.hpp" // declarations
#include "debug.hpp"     // for AR_REQUIRES
#include <algorithm>     // for min
#include <array>         // for array
#include <cmath>         // for sqrt
#include <numeric>       // for accumulate

namespace adapterremoval {

namespace {

//! Critical Values for Grubb's test calculated in R:
//!   > n = 0:100
//!   > t = qt(1 - 0.01/(2 * n), n - 2)
//!   > round(((n - 1) / sqrt(n)) * sqrt((t^2 / (n - 2 + t^2))), 4)
constexpr std::array<double, 101> GRUBBS_TEST_CRITICAL_VALUES = {
  NAN,    NAN,    NAN,    1.1547, 1.4963, 1.7637, 1.9728, 2.1391, 2.2744,
  2.3868, 2.4821, 2.5641, 2.6357, 2.6990, 2.7554, 2.8061, 2.8521, 2.8940,
  2.9325, 2.9680, 3.0008, 3.0314, 3.0599, 3.0866, 3.1117, 3.1353, 3.1577,
  3.1788, 3.1989, 3.2179, 3.2361, 3.2534, 3.2700, 3.2858, 3.3010, 3.3156,
  3.3296, 3.3431, 3.3561, 3.3686, 3.3807, 3.3924, 3.4037, 3.4146, 3.4252,
  3.4354, 3.4454, 3.4551, 3.4645, 3.4736, 3.4825, 3.4911, 3.4995, 3.5077,
  3.5157, 3.5235, 3.5311, 3.5386, 3.5458, 3.5529, 3.5598, 3.5666, 3.5733,
  3.5798, 3.5861, 3.5924, 3.5985, 3.6044, 3.6103, 3.6161, 3.6217, 3.6272,
  3.6327, 3.6380, 3.6433, 3.6484, 3.6535, 3.6585, 3.6633, 3.6682, 3.6729,
  3.6775, 3.6821, 3.6866, 3.6911, 3.6955, 3.6998, 3.7040, 3.7082, 3.7123,
  3.7164, 3.7204, 3.7243, 3.7282, 3.7320, 3.7358, 3.7396, 3.7432, 3.7469,
  3.7505, 3.7540,
};

double
grubbs_test_critical_values(size_t n)
{
  n = std::min(n, GRUBBS_TEST_CRITICAL_VALUES.size() - 1);

  return GRUBBS_TEST_CRITICAL_VALUES.at(n);
}

} // namespace

double
arithmetic_mean(const std::vector<uint64_t>& values)
{
  AR_REQUIRE(!values.empty());

  return std::accumulate(values.begin(), values.end(), uint64_t()) /
         static_cast<double>(values.size());
}

double
standard_deviation(const std::vector<uint64_t>& values)
{
  AR_REQUIRE(values.size() > 1);

  double error = 0;
  const auto m = arithmetic_mean(values);
  for (const auto it : values) {
    error += (it - m) * (it - m);
  }

  return std::sqrt(error / (values.size() - 1));
}

bool
grubbs_test_prune(std::vector<uint64_t>& values)
{
  const size_t original_size = values.size();

  while (values.size() >= 3) {
    const double m = arithmetic_mean(values);
    const double sd = standard_deviation(values);
    const double cutoff = grubbs_test_critical_values(values.size());

    size_t greatest_outlier_i = 0;
    auto greatest_outlier = std::abs(values.at(0) - m) / sd;
    for (size_t i = 1; i < values.size(); ++i) {
      const auto outlier = std::abs(values.at(i) - m) / sd;

      if (outlier > greatest_outlier) {
        greatest_outlier = outlier;
        greatest_outlier_i = i;
      }
    }

    if (greatest_outlier > cutoff) {
      values.at(greatest_outlier_i) = values.back();
      values.pop_back();
    } else {
      break;
    }
  }

  return original_size != values.size();
}

} // namespace adapterremoval
