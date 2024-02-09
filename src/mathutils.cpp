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
#include "mathutils.hpp"
#include "debug.hpp" // for AR_REQUIRES
#include <cmath>     // for sqrt
#include <numeric>   // for accumulate

namespace adapterremoval {

namespace {

//! Critical Values of the Student's t Distribution at 0.995 for 0 to 100 df,
//! calculated via `qt(0.995, 1:100)` in R
const std::vector<double> STUDENTS_T_CRITICAL_VALUES = {
  NAN,    63.657, 9.9248, 5.8409, 4.6041, 4.0321, 3.7074, 3.4995, 3.3554,
  3.2498, 3.1693, 3.1058, 3.0545, 3.0123, 2.9768, 2.9467, 2.9208, 2.8982,
  2.8784, 2.8609, 2.8453, 2.8314, 2.8188, 2.8073, 2.7969, 2.7874, 2.7787,
  2.7707, 2.7633, 2.7564, 2.7500, 2.7440, 2.7385, 2.7333, 2.7284, 2.7238,
  2.7195, 2.7154, 2.7116, 2.7079, 2.7045, 2.7012, 2.6981, 2.6951, 2.6923,
  2.6896, 2.6870, 2.6846, 2.6822, 2.6800, 2.6778, 2.6757, 2.6737, 2.6718,
  2.6700, 2.6682, 2.6665, 2.6649, 2.6633, 2.6618, 2.6603, 2.6589, 2.6575,
  2.6561, 2.6549, 2.6536, 2.6524, 2.6512, 2.6501, 2.6490, 2.6479, 2.6469,
  2.6459, 2.6449, 2.6439, 2.6430, 2.6421, 2.6412, 2.6403, 2.6395, 2.6387,
  2.6379, 2.6371, 2.6364, 2.6356, 2.6349, 2.6342, 2.6335, 2.6329, 2.6322,
  2.6316, 2.6309, 2.6303, 2.6297, 2.6291, 2.6286, 2.6280, 2.6275, 2.6269,
  2.6264, 2.6259
};

} // namespace

double
arithmetic_mean(const std::vector<uint64_t>& values)
{
  AR_REQUIRE(values.size());

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

double
students_t_critical_value(size_t df)
{
  return STUDENTS_T_CRITICAL_VALUES.at(
    std::min<size_t>(df, STUDENTS_T_CRITICAL_VALUES.size() - 1));
}

bool
grubbs_test_prune(std::vector<uint64_t>& values)
{
  const size_t original_size = values.size();
  bool found_any = false;

  do {
    const double N = values.size();
    const double t = students_t_critical_value(N - 2);
    const auto t2 = t * t;

    const double m = arithmetic_mean(values);
    const double sd = standard_deviation(values);

    const double cutoff = ((N - 1) / std::sqrt(N)) * sqrt(t2 / (N - 2 + t2));

    found_any = false;
    for (size_t i = 0; i < values.size(); ++i) {
      const auto G = (values.at(i) - m) / sd;
      if (G < -cutoff || G > cutoff) {
        values.at(i) = values.back();
        values.pop_back();
        found_any = true;
        break;
      }
    }
  } while (found_any);

  return original_size != values.size();
}

} // namespace adapterremoval
