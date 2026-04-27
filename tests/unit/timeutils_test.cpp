// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"    // for assert_failed
#include "testing.hpp"   // for TEST_CASE, REQUIRE, ...
#include "timeutils.hpp" // for seconds_to_duration

namespace adapterremoval {

TEST_CASE("seconds_to_duration", "[timeutils]")
{
  CHECK_THROWS_AS(seconds_to_duration(-1.0), assert_failed);
  CHECK(seconds_to_duration(0) == "0.0s");
  CHECK(seconds_to_duration(0.09) == "0.1s");
  CHECK(seconds_to_duration(0.1) == "0.1s");
  CHECK(seconds_to_duration(12.3) == "12.3s");
  CHECK(seconds_to_duration(60) == "1:00.0s");
  CHECK(seconds_to_duration(83.4) == "1:23.4s");
  CHECK(seconds_to_duration(1234.5) == "20:34.5s");
  CHECK(seconds_to_duration(3600) == "1:00:00.0s");
  CHECK(seconds_to_duration(12345.6) == "3:25:45.6s");
}

} // namespace adapterremoval
