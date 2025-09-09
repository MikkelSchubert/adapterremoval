// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"    // for assertion_failed
#include "scheduler.hpp" // for threadstate
#include "testing.hpp"   // for TEST_CASE, REQUIRE, ...
#include <memory>        // for unique_ptr

namespace Catch {

template<>
std::string
fallbackStringifier(const std::unique_ptr<int>& value)
{
  ReusableStringStream ss;

  if (value) {
    ss << "unique_ptr({" << *value << "})";
  } else {
    ss << "unique_ptr(null)";
  }

  return ss.str();
}

} // namespace Catch

namespace adapterremoval {

TEST_CASE("empty threadstate", "[scheduler::threadstate]")
{
  threadstate<int> state;

  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
  REQUIRE(state.try_acquire() == threadstate<int>::pointer());
}

TEST_CASE("threadstate emplace back", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back(1234);

  SECTION("acquire")
  {
    auto value = state.acquire();
    REQUIRE(value);
    REQUIRE(*value == 1234);
    REQUIRE_THROWS_AS(state.acquire(), assert_failed);
  }

  SECTION("try acquire")
  {
    auto value = state.try_acquire();
    REQUIRE(value);
    REQUIRE(*value == 1234);
    REQUIRE(state.try_acquire() == threadstate<int>::pointer());
  }
}

TEST_CASE("threadstate emplace back n", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back_n(3, 1234);

  for (size_t i = 0; i < 3; ++i) {
    auto value = state.acquire();
    REQUIRE(value);
    REQUIRE(*value == 1234);
  }

  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate release", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back(1234);

  auto value = state.acquire();
  REQUIRE(value);
  REQUIRE(*value == 1234);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);

  state.release(std::move(value));

  value = state.acquire();
  REQUIRE(value);
  REQUIRE(*value == 1234);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate release lvalue", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.release(std::make_unique<int>(1234));

  auto value = state.acquire();
  REQUIRE(value);
  REQUIRE(*value == 1234);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate merg into", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back(1);
  state.emplace_back(3);
  state.emplace_back(5);

  int destination = 10;
  state.merge_into(destination);
  REQUIRE(destination == 19);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

} // namespace adapterremoval
