// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "catch.hpp"
#include <utility>

namespace Catch {

// This template is designed to match all types not handled natively by Catch,
// in order to ensure that the StringMaker is implemented for all tested types
template<typename T>
struct StringMaker<
  T,
  typename std::enable_if<std::is_class<T>::value && !is_range<T>::value>::type>
{
  static std::string convert(T const& value);
};

/** Helper macro for typed catching of thrown messages */
#define REQUIRE_THROWS_MESSAGE(expr, exceptionType, message)                   \
  REQUIRE_THROWS_MATCHES(expr, exceptionType, Catch::Message(message))

/** Helper macro for typed catching of thrown messages */
#define CHECK_THROWS_MESSAGE(expr, exceptionType, message)                     \
  CHECK_THROWS_MATCHES(expr, exceptionType, Catch::Message(message))

} // namespace Catch