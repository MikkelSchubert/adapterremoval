// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"  // for assert_failed
#include "testing.hpp" // for declarations
#include <array>       // for array
#include <ostream>     // for ostream
#include <stdexcept>   // for invalid_argument
#include <utility>     // for pair
#include <vector>      // for vector

////////////////////////////////////////////////////////////////////////////////
// fallbackStringifier

namespace {

class stringifier_test
{
public:
  explicit stringifier_test(int i)
    : m_i(i)
  {
  }

  [[nodiscard]] int get_i() const { return m_i; }

  friend std::ostream& operator<<(std::ostream& os,
                                  const stringifier_test& value)
  {
    return os << "stringifier_test{ " << value.get_i() << " }";
  }

private:
  int m_i{};
};

} // namespace

TEST_CASE("fallbackStringifier", "[testing:fallbackStringifier]")
{
  REQUIRE(Catch::fallbackStringifier(stringifier_test{ 17 }) ==
          "stringifier_test{ 17 }");
}

////////////////////////////////////////////////////////////////////////////////
// custom StringMaker for containers and exceptions

TEST_CASE("StringMaker for array", "[StringMaker]")
{
  REQUIRE(Catch::StringMaker<std::array<int, 3>>::convert({ 1, 2, 3 }) ==
          "array{ 1, 2, 3 }");
}

TEST_CASE("StringMaker for vector", "[StringMaker]")
{
  REQUIRE(Catch::StringMaker<std::vector<int>>::convert({ 1, 2, 3 }) ==
          "vector{ 1, 2, 3 }");
}

TEST_CASE("StringMaker for invalid_argument", "[StringMaker]")
{
  REQUIRE(Catch::StringMaker<std::invalid_argument>::convert(
            std::invalid_argument{ "some argument" }) ==
          "invalid_argument{ 'some argument' }");
}

TEST_CASE("StringMaker for pair", "[StringMaker]")
{
  auto value = std::pair<stringifier_test, int>(stringifier_test{ 17 }, 3);

  REQUIRE(Catch::StringMaker<decltype(value)>::convert(value) ==
          "pair{ first=stringifier_test{ 17 }, second=3 }");
}

////////////////////////////////////////////////////////////////////////////////
// ExceptionContainsMatcher

TEST_CASE("ExceptionContainsMatcher", "[ExceptionContainsMatcher]")
{
  using adapterremoval::assert_failed;
  using Catch::ExceptionContainsMatcher;

  auto throw_error = [](std::string message) { throw assert_failed(message); };

  REQUIRE_THROWS_MATCHES(throw_error("exact message"),
                         assert_failed,
                         ExceptionContainsMatcher("exact message"));
  REQUIRE_THROWS_MATCHES(throw_error("wrong message"),
                         assert_failed,
                         !ExceptionContainsMatcher("exact message"));

  REQUIRE_THROWS_MATCHES(throw_error("this is a long message"),
                         assert_failed,
                         ExceptionContainsMatcher("long message"));
  REQUIRE_THROWS_MATCHES(throw_error("this is a short message"),
                         assert_failed,
                         !ExceptionContainsMatcher("long message"));
}
