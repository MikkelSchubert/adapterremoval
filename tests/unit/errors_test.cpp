// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp" // declarations
#include "strutils.hpp"
#include "testing.hpp" // for TEST_CASE, REQUIRE, ...
#include <sstream>     // for ostringstream

namespace adapterremoval {

TEST_CASE("assert_failed to string", "[errors]")
{
  std::ostringstream ss;
  ss << assert_failed("66031296-fe26-4378-85e7-c7af0f473f9a");
  CHECK(ss.str() == "assert_failed{'66031296-fe26-4378-85e7-c7af0f473f9a'}");
}

TEST_CASE("io_error to string", "[errors]")
{
  std::ostringstream ss;
  ss << io_error("02b1544a-4ab2-40cc-8ad2-4dc64a918137");
  CHECK(ss.str() == "io_error{'02b1544a-4ab2-40cc-8ad2-4dc64a918137'}");
}

TEST_CASE("io_error to string with errno", "[errors]")
{
  SECTION("errno zero")
  {
    std::ostringstream ss;
    ss << io_error("c4da2513-5d4b-4f58-8cea-2693e08f7ec0", 0);
    CHECK(ss.str() == "io_error{'c4da2513-5d4b-4f58-8cea-2693e08f7ec0'}");
  }

  SECTION("no-zero errno")
  {
    std::ostringstream expected;
    expected << "io_error{"
             << log_escape(
                  format_io_error("a7bb599c-2ef4-47ff-9190-d960f00779af", 1))
             << "}";

    std::ostringstream ss;
    ss << io_error("a7bb599c-2ef4-47ff-9190-d960f00779af", 1);
    CHECK(ss.str() == expected.str());
  }
}

TEST_CASE("gzip_error to string", "[errors]")
{
  SECTION("as gzip_error")
  {
    std::ostringstream ss;
    ss << gzip_error("74e5165e-0417-47a3-a64e-df766f7c57f4");
    CHECK(ss.str() == "gzip_error{'74e5165e-0417-47a3-a64e-df766f7c57f4'}");
  }

  SECTION("as io_error")
  {
    const gzip_error error{ "87a2fa19-8476-435d-8715-794f13004a9d" };
    std::ostringstream ss;
    ss << dynamic_cast<const io_error&>(error);
    CHECK(ss.str() == "gzip_error{'87a2fa19-8476-435d-8715-794f13004a9d'}");
  }
}

TEST_CASE("parsing_error to string", "[errors]")
{
  std::ostringstream ss;
  ss << parsing_error("cea794d8-d3dd-4d88-a113-2f8fe1d95f52");
  CHECK(ss.str() == "parsing_error{'cea794d8-d3dd-4d88-a113-2f8fe1d95f52'}");
}

TEST_CASE("serializing_error to string", "[errors]")
{
  std::ostringstream ss;
  ss << serializing_error("bfeb14b5-817c-47bf-ade1-30cdb0cb3c55");
  CHECK(ss.str() ==
        "serializing_error{'bfeb14b5-817c-47bf-ade1-30cdb0cb3c55'}");
}

TEST_CASE("fastq_error to string", "[errors]")
{
  SECTION("as fastq_error")
  {
    std::ostringstream ss;
    ss << fastq_error("d611b2a5-6501-45b0-94dd-5bec7971f84c");
    CHECK(ss.str() == "fastq_error{'d611b2a5-6501-45b0-94dd-5bec7971f84c'}");
  }

  SECTION("as parsing_error")
  {
    const fastq_error error{ "7f61b662-eea0-4451-85b1-300272d13abc" };
    std::ostringstream ss;
    ss << dynamic_cast<const parsing_error&>(error);
    CHECK(ss.str() == "fastq_error{'7f61b662-eea0-4451-85b1-300272d13abc'}");
  }
}

TEST_CASE("fatal_error to string", "[errors]")
{
  std::ostringstream ss;
  ss << fatal_error("d611b2a5-6501-45b0-94dd-5bec7971f84c");
  CHECK(ss.str() == "fatal_error{'d611b2a5-6501-45b0-94dd-5bec7971f84c'}");
}

} // namespace adapterremoval
