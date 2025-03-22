// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "testing.hpp"   // for catch.hpp, StringMaker
#include "utilities.hpp" // for merge
#include <algorithm>     // for max
#include <array>         // for array, operator==
#include <cstddef>       // for size_t
#include <vector>        // for vector, allocator, operator==

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
    REQUIRE(dst.empty());
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
