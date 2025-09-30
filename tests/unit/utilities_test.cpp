// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "testing.hpp"   // for TEST_CASE, REQUIRE, ...
#include "utilities.hpp" // for merge
#include <array>         // for array, operator==
#include <cstddef>       // for size_t
#include <memory>        // for unique_ptr
#include <string>        // for string
#include <type_traits>   // for is_same_v
#include <vector>        // for vector, allocator, operator==

namespace adapterremoval {

namespace merging {

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

} // namespace merging

////////////////////////////////////////////////////////////////////////////////
// underlying_value

TEST_CASE("underlying_value enum value", "[underlying_value]")
{
  enum my_enum
  {
    a = -1,
    b = 0,
    c = 23,
  };

  CHECK(underlying_value(my_enum::a) == -1);
  CHECK(underlying_value(my_enum::b) == 0);
  CHECK(underlying_value(my_enum::c) == 23);
}

TEST_CASE("underlying_value enum class value", "[underlying_value]")
{
  enum class my_enum : unsigned
  {
    a = 0,
    b = 4,
    c = 8,
  };

  CHECK(underlying_value(my_enum::a) == 0);
  CHECK(underlying_value(my_enum::b) == 4);
  CHECK(underlying_value(my_enum::c) == 8);
}

TEST_CASE("underlying_value returns same type", "[underlying_value]")
{
  // These tests do not do anything at runtime
  {
    enum my_enum : char;
    auto uvalue = underlying_value(my_enum{});
    static_assert(std::is_same_v<decltype(uvalue), char>);
  }
  {
    enum my_enum : unsigned char;
    auto uvalue = underlying_value(my_enum{});
    static_assert(std::is_same_v<decltype(uvalue), unsigned char>);
  }

  {
    enum my_enum : int;
    auto uvalue = underlying_value(my_enum{});
    static_assert(std::is_same_v<decltype(uvalue), int>);
  }
  {
    enum my_enum : unsigned int;
    auto uvalue = underlying_value(my_enum{});
    static_assert(std::is_same_v<decltype(uvalue), unsigned int>);
  }

  {
    enum class my_enum : long;
    auto uvalue = underlying_value(my_enum{});
    static_assert(std::is_same_v<decltype(uvalue), long>);
  }
  {
    enum class my_enum : unsigned long;
    auto uvalue = underlying_value(my_enum{});
    static_assert(std::is_same_v<decltype(uvalue), unsigned long>);
  }
}

////////////////////////////////////////////////////////////////////////////////
// dynamic_cast_unique

namespace {

class test_base_class
{
public:
  test_base_class() = default;
  test_base_class(const test_base_class&) = delete;
  test_base_class(test_base_class&&) = delete;
  virtual ~test_base_class() = default;

  test_base_class& operator=(const test_base_class&) = delete;
  test_base_class& operator=(test_base_class&&) = delete;
};

class test_sub_class_1 : public test_base_class
{
public:
  explicit test_sub_class_1(int value_)
    : value(value_) {};

  test_sub_class_1(const test_sub_class_1&) = delete;
  test_sub_class_1(test_sub_class_1&&) = delete;
  ~test_sub_class_1() override = default;

  test_sub_class_1& operator=(const test_sub_class_1&) = delete;
  test_sub_class_1& operator=(test_sub_class_1&&) = delete;

  int value = 0;
};

class test_sub_class_2 : public test_base_class
{
public:
  explicit test_sub_class_2() = default;
  test_sub_class_2(const test_sub_class_2&) = delete;
  test_sub_class_2(test_sub_class_2&&) = delete;
  ~test_sub_class_2() override = default;

  test_sub_class_2& operator=(const test_sub_class_2&) = delete;
  test_sub_class_2& operator=(test_sub_class_2&&) = delete;
};

} // namespace

TEST_CASE("dynamic_cast_unique on nulllptr returns nullptr")
{
  std::unique_ptr<test_base_class> ptr;
  REQUIRE_FALSE(ptr.get());
  REQUIRE_FALSE(dynamic_cast_unique<test_sub_class_1>(ptr).get());
  REQUIRE_FALSE(ptr.get());
  REQUIRE_FALSE(dynamic_cast_unique<test_sub_class_2>(ptr).get());
  REQUIRE_FALSE(ptr.get());
}

TEST_CASE("dynamic_cast_unique on wrong type returns nullptr")
{
  auto ptr = std::make_unique<test_sub_class_1>(12345);
  REQUIRE_FALSE(dynamic_cast_unique<test_sub_class_2>(ptr).get());
  REQUIRE(ptr.get());
}

TEST_CASE("dynamic_cast_unique on correct type casts")
{
  auto ptr = std::make_unique<test_sub_class_1>(12345);
  auto data = dynamic_cast_unique<test_sub_class_1>(ptr);
  REQUIRE(data.get());
  REQUIRE(data->value == 12345);
  REQUIRE_FALSE(ptr.get());
}

} // namespace adapterremoval
