/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <array>
#include <vector>

#include "testing.hpp"
#include "utilities.hpp"

namespace adapterremoval {

using arr = std::array<size_t, 3>;
using vec = std::vector<size_t>;
using arrvec = std::vector<arr>;
using vecvec = std::vector<vec>;

TEST_CASE("merging POD")
{
  SECTION("integers")
  {
    int dst = 17;
    merge(dst, 4);
    REQUIRE(dst == 21);
  }

  SECTION("float")
  {
    double dst = 3.141;
    merge(dst, 2.718);
    REQUIRE_THAT(dst, Catch::WithinAbs(5.859, 1e-6));
  }
}

TEST_CASE("merging arrays")
{
  arr dst = { 1, 2, 3 };
  merge(dst, arr{ 7, 9, 13 });
  REQUIRE(dst == arr{ 8, 11, 16 });
}

TEST_CASE("merging vectors")
{
  SECTION("empty into empty")
  {
    vec dst;
    merge(dst, {});
    REQUIRE(dst == vec{});
  }

  SECTION("values into empty")
  {
    vec dst;
    merge(dst, { 4, 5, 6 });
    REQUIRE(dst == vec{ 4, 5, 6 });
  }

  SECTION("values into values")
  {
    vec dst = { 1, 2, 3 };
    merge(dst, vec{ 7, 9, 13 });
    REQUIRE(dst == vec{ 8, 11, 16 });
  }

  SECTION("values into fewer values")
  {
    vec dst = { 1, 2 };
    merge(dst, vec{ 7, 9, 13 });
    REQUIRE(dst == vec{ 8, 11, 13 });
  }

  SECTION("values into more values")
  {
    vec dst = { 1, 2, 3, 4 };
    merge(dst, vec{ 7, 9, 13 });
    REQUIRE(dst == vec{ 8, 11, 16, 4 });
  }
}

TEST_CASE("merging vector of vectors")
{
  vecvec dst{ { 10, 20 } };
  merge(dst, { { 1, 2 }, { 3 }, { 4, 5, 6 } });
  REQUIRE(dst == vecvec{ { 11, 22 }, { 3 }, { 4, 5, 6 } });
}

TEST_CASE("merging vector of arrays")
{
  arrvec dst{ { 10, 20 } };
  merge(dst, { { 1, 2 }, { 3 }, { 4, 5, 6 } });
  REQUIRE(dst == arrvec{ { 11, 22, 0 }, { 3, 0, 0 }, { 4, 5, 6 } });
}

} // namespace adapterremoval
