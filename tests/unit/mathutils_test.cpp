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
#include "catch.hpp"     // for SourceLineInfo, operator""_catch_sr, Assert...
#include "mathutils.hpp" // for arithmetic_mean, standard_deviation
#include <vector>        // for vector

namespace adapterremoval {

namespace {

using i32vec = std::vector<int32_t>;
using i64vec = std::vector<int64_t>;
using u32vec = std::vector<uint32_t>;
using u64vec = std::vector<uint64_t>;
using dvec = std::vector<double>;

} // namespace

TEST_CASE("arithmetic mean")
{
  SECTION("i32")
  {
    REQUIRE(arithmetic_mean(i32vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(i32vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(i32vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(i32vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(i32vec{ 2147483647, 124 }) == Approx(1073741885.5));
    REQUIRE(arithmetic_mean(i32vec{ -2147483647, -124 }) ==
            Approx(-1073741885.5));
  }

  SECTION("i64")
  {
    REQUIRE(arithmetic_mean(i64vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(i64vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(i64vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(i64vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(i64vec{ 2147483647, 124 }) == Approx(1073741885.5));
    REQUIRE(arithmetic_mean(i64vec{ -2147483647, -124 }) ==
            Approx(-1073741885.5));
  }

  SECTION("u32")
  {
    REQUIRE_THROWS_AS(arithmetic_mean(u32vec()), assert_failed);
    REQUIRE(arithmetic_mean(u32vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(u32vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(u32vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(u32vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(u32vec{ 2147483647, 124 }) == Approx(1073741885.5));
    REQUIRE(arithmetic_mean(u32vec{ 4294967295, 124 }) == Approx(2147483709.5));
  }

  SECTION("u64")
  {
    REQUIRE_THROWS_AS(arithmetic_mean(u64vec()), assert_failed);
    REQUIRE(arithmetic_mean(u64vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(u64vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(u64vec{ 13 }) == Approx(13));
    REQUIRE(arithmetic_mean(u64vec{ 10, 15 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(u64vec{ 9223372036854775807llu, 124 }) ==
            Approx(4611686018427387965.5));
    REQUIRE(arithmetic_mean(u64vec{ 18446744073709551615llu, 124 }) ==
            Approx(9223372036854775869.5));
  }

  SECTION("d64")
  {
    REQUIRE_THROWS_AS(arithmetic_mean(dvec()), assert_failed);
    REQUIRE(arithmetic_mean(dvec{ 13.0 }) == Approx(13));
    REQUIRE(arithmetic_mean(dvec{ 10.0, 15.0 }) == Approx(12.5));
    REQUIRE(arithmetic_mean(dvec{ 2147483647.0, 124.0 }) ==
            Approx(1073741885.5));
    REQUIRE(arithmetic_mean(dvec{ 4294967295.0, 124.0 }) ==
            Approx(2147483709.5));
    REQUIRE(arithmetic_mean(
              dvec{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }) == 0.1);
  }
}

TEST_CASE("sample standard deviation")
{
  REQUIRE_THROWS_AS(standard_deviation(dvec()), assert_failed);
  REQUIRE_THROWS_AS(standard_deviation(dvec(1.0)), assert_failed);
  REQUIRE(standard_deviation(u32vec{ 13, 13 }) == Approx(0.0));
  REQUIRE(standard_deviation(dvec{ 0.2236275, 0.2638758, -0.8673293 }) ==
          Approx(0.6417985));
  REQUIRE(standard_deviation(dvec{ -1.4489741, -0.1122757, -0.4124156 }) ==
          Approx(0.701344));
  REQUIRE(standard_deviation(dvec{ 0.08877328, -0.59031026, 0.85331033 }) ==
          Approx(0.7222317));
}

} // namespace adapterremoval
