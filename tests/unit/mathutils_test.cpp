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
#include "errors.hpp"    // for assert_failed
#include "mathutils.hpp" // for arithmetic_mean, standard_deviation
#include "testing.hpp"   // for catch.hpp, StringMaker
#include <vector>        // for vector

namespace adapterremoval {

namespace {

using u64vec = std::vector<uint64_t>;

} // namespace

TEST_CASE("arithmetic mean")
{
  REQUIRE_THROWS_AS(arithmetic_mean(u64vec()), assert_failed);
  REQUIRE(arithmetic_mean(u64vec{ 13 }) == Approx(13));
  REQUIRE(arithmetic_mean(u64vec{ 10, 15 }) == Approx(12.5));
  REQUIRE(arithmetic_mean(u64vec{ 119509, 48979, 38789 }) == Approx(69092.33));
}

TEST_CASE("sample standard deviation")
{
  REQUIRE_THROWS_AS(standard_deviation(u64vec{}), assert_failed);
  REQUIRE_THROWS_AS(standard_deviation(u64vec{ 1 }), assert_failed);
  REQUIRE(standard_deviation(u64vec{ 13, 13 }) == Approx(0.0));
  REQUIRE(standard_deviation(u64vec{ 27405, 118434, 94082 }) ==
          Approx(47125.93));
  REQUIRE(standard_deviation(u64vec{ 6541, 57558, 19385 }) == Approx(26535.76));
  REQUIRE(standard_deviation(u64vec{ 119509, 48979, 38789 }) ==
          Approx(43958.38));
}

} // namespace adapterremoval
