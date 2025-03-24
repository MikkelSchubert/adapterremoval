// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"  // declarations
#include "testing.hpp" // for TEST_CASE, REQUIRE, ...
#include <sstream>     // for ostringstream

namespace adapterremoval {

TEST_CASE("assert_failed to string", "[errors]")
{
  std::ostringstream ss;
  ss << assert_failed("error message");
  CHECK(ss.str() == "assert_failed{'error message'}");
}

TEST_CASE("parsing_error to string", "[errors]")
{
  std::ostringstream ss;
  ss << parsing_error("error message");
  CHECK(ss.str() == "parsing_error{'error message'}");
}

TEST_CASE("fastq_error to string", "[errors]")
{
  std::ostringstream ss;
  ss << fastq_error("error message");
  CHECK(ss.str() == "fastq_error{'error message'}");
}

} // namespace adapterremoval
