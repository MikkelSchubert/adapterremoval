// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "scheduler.hpp" // for step_id
#include "testing.hpp"   // for CHECK, CHECK_FALSE, TEST_CASE
#include <sstream>       // for ostringstream

namespace adapterremoval {

step_id
step_id_test_helper(size_t i)
{
  return step_id{ i };
}

namespace {

TEST_CASE("step_id equality", "[step_id]")
{
  CHECK(step_id::UNDEFINED == step_id::UNDEFINED);
  CHECK_FALSE(step_id::UNDEFINED != step_id::UNDEFINED);

  CHECK_FALSE(step_id_test_helper(1) == step_id::UNDEFINED);
  CHECK(step_id_test_helper(1) != step_id::UNDEFINED);

  CHECK_FALSE(step_id::UNDEFINED == step_id_test_helper(1));
  CHECK(step_id::UNDEFINED != step_id_test_helper(1));

  CHECK_FALSE(step_id_test_helper(2) == step_id_test_helper(1));
  CHECK(step_id_test_helper(2) != step_id_test_helper(1));

  const auto s = step_id_test_helper(1);
  CHECK(s == s);
  CHECK_FALSE(s != s);
}

TEST_CASE("step_id get", "[step_id]")
{
  CHECK(step_id_test_helper(0).get() == 0);
  CHECK(step_id_test_helper(12345).get() == 12345);
  CHECK(step_id_test_helper(9223372036854775808LLU).get() ==
        9223372036854775808LLU);
  CHECK(step_id::UNDEFINED.get() == static_cast<size_t>(-1));
}

TEST_CASE("step_id id to string", "[step_id]")
{
  std::ostringstream os;
  os << step_id_test_helper(12345);

  CHECK(os.str() == "step_id{12345}");
}

TEST_CASE("step_id undefined to string", "[step_id]")
{
  std::ostringstream os;
  os << step_id::UNDEFINED;

  CHECK(os.str() == "step_id::UNDEFINED");
}

} // namespace

} // namespace adapterremoval
