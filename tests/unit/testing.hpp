/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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