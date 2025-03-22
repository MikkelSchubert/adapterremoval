// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
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
