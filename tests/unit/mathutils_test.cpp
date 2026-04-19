// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"    // for assert_failed
#include "mathutils.hpp" // for arithmetic_mean, standard_deviation
#include "testing.hpp"   // for TEST_CASE, REQUIRE, ...
#include <cstdint>       // for uint64_t
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

TEST_CASE("grubb's test pruning -- too few values")
{
  std::vector<uint64_t> values;

  CHECK_FALSE(grubbs_test_prune(values));
  REQUIRE(values == std::vector<uint64_t>{});

  values.emplace_back(1);
  CHECK_FALSE(grubbs_test_prune(values));
  REQUIRE(values == std::vector<uint64_t>{ 1 });

  values.emplace_back(2);
  CHECK_FALSE(grubbs_test_prune(values));
  REQUIRE(values == std::vector<uint64_t>{ 1, 2 });
}

TEST_CASE("grubb's test pruning -- no outliers")
{
  SECTION("minimal sample")
  {
    std::vector<uint64_t> values = { 1, 2, 3 };
    CHECK_FALSE(grubbs_test_prune(values));
    REQUIRE(values == std::vector<uint64_t>{ 1, 2, 3 });
  }

  SECTION("grubb's test pruning -- many values")
  {
    const std::vector<uint64_t> expected = {
      756, 1608, 281, 239,  649, 968, 180, 1717, 811,  726, 1386, 342,
      287, 2748, 221, 2352, 813, 353, 328, 128,  1501, 537, 306,
    };

    auto observed = expected;
    CHECK_FALSE(grubbs_test_prune(observed));
    REQUIRE(observed == expected);
  }
}

TEST_CASE("grubb's test pruning -- with outliers")
{
  std::vector<uint64_t> values = {
    756, 1608, 281, 239,  649, 968, 180, 1717, 811,  726, 8800, 1386, 342,
    287, 2748, 221, 2352, 813, 353, 328, 128,  1501, 537, 7922, 306,
  };
  const std::vector<uint64_t> expected = {
    756, 1608, 281,  239, 649,  968, 180, 1717, 811, 726,  306, 1386,
    342, 287,  2748, 221, 2352, 813, 353, 328,  128, 1501, 537,
  };

  CHECK(grubbs_test_prune(values));
  REQUIRE(values == expected);
}

} // namespace adapterremoval
