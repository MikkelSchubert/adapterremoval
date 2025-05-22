// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <array>     // for array
#include <cstddef>   // for size_t
#include <stdexcept> // for invalid_argument
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector

#define CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FALLBACK_STRINGIFIER ::Catch::fallbackStringifier

namespace Catch {

/** Helper macro for typed catching of thrown messages */
#define REQUIRE_THROWS_MESSAGE(expr, exceptionType, message)                   \
  REQUIRE_THROWS_MATCHES(expr, exceptionType, Catch::Message(message))

/** Helper macro for typed catching of thrown messages */
#define CHECK_THROWS_MESSAGE(expr, exceptionType, message)                     \
  CHECK_THROWS_MATCHES(expr, exceptionType, Catch::Message(message))

/** This function needs to be declared before "catch.hpp" is included */
template<typename T>
std::string
fallbackStringifier(const T& value);

} // namespace Catch

#include "catch.hpp" // for Catch, StringMaker, stringify

namespace Catch {

template<typename T>
std::string
fallbackStringifier(const T& value)
{
  // A stringstream is used instead of convertUnstreamable, since the latter
  // will simply omit the contents of unknown object types. This way, the
  // appropriate stream operator must be defined up front
  ReusableStringStream ss;
  ss << value;
  return ss.str();
}

template<typename T, size_t N>
struct StringMaker<std::array<T, N>>
{
  static std::string convert(const std::array<T, N>& range)
  {
    return "array" + rangeToString(range);
  }
};

template<typename T, typename U>
struct StringMaker<std::pair<T, U>>
{
  static std::string convert(const std::pair<T, U>& value)
  {
    ReusableStringStream os;
    os << "pair{ first=" << ::Catch::Detail::stringify(value.first)
       << ", second=" << ::Catch::Detail::stringify(value.second) << " }";
    return os.str();
  }
};

template<typename T>
struct StringMaker<std::vector<T>>
{
  static std::string convert(const std::vector<T>& range)
  {
    return "vector" + rangeToString(range);
  }
};

template<>
struct StringMaker<std::invalid_argument>
{
  static std::string convert(const std::invalid_argument& value);
};

} // namespace Catch
