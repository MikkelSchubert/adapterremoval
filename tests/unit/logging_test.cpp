// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"  // for assert_failed
#include "logging.hpp" // for log_capture, info, log_stream, cerr, error, warn
#include "testing.hpp" // for catch.hpp, StringMaker
#include <string>      // for basic_string, operator==, string

namespace adapterremoval {

TEST_CASE("log capture timestamps")
{
  log::set_timestamps(false);
  // Force logging of preamble prior to tests
  log::log_preamble();

  SECTION("timestamps disabled")
  {
    log::log_capture cap;
    log::info() << "test 1";
    REQUIRE(cap.str() == "[INFO] test 1\n");
  }

  SECTION("timestamps enabled")
  {
    log::log_capture cap;
    log::set_timestamps(true);
    log::info() << "test 2";
    REQUIRE_THAT(cap.str(), Catch::StartsWith("20"));
    REQUIRE_THAT(cap.str(), Catch::EndsWith("[INFO] test 2\n"));
  }

  SECTION("timestamps restored")
  {
    log::log_capture cap;
    log::info() << "test 3";
    REQUIRE(cap.str() == "[INFO] test 3\n");
  }
}

TEST_CASE("log capture colors")
{
  log::set_colors(false);
  // Force logging of preamble prior to tests
  log::log_preamble();

  SECTION("colors enabled")
  {
    log::log_capture cap;
    log::set_colors(true);
    log::info() << "test 1";
    REQUIRE(cap.str() == "[\033[0;32mINFO\033[0m] test 1\n");
  }

  SECTION("colors restored")
  {
    log::log_capture cap;
    log::info() << "test 2";
    REQUIRE(cap.str() == "[INFO] test 2\n");
  }
}

TEST_CASE("log capture levels")
{
  log::set_level(log::level::debug);
  // Force logging of preamble prior to tests
  log::log_preamble();

  SECTION("default levels")
  {
    log::log_capture cap;
    log::debug() << "test 1";
    log::info() << "test 2";
    log::warn() << "test 3";
    log::error() << "test 4";
    REQUIRE(cap.str() == "[DEBUG] test 1\n"
                         "[INFO] test 2\n"
                         "[WARNING] test 3\n"
                         "[ERROR] test 4\n");
  }

  SECTION("filter levels")
  {
    log::log_capture cap;
    log::set_level(log::level::warning);
    log::debug() << "test 1";
    log::info() << "test 2";
    log::warn() << "test 3";
    log::error() << "test 4";
    REQUIRE(cap.str() == "[WARNING] test 3\n"
                         "[ERROR] test 4\n");
  }

  SECTION("levels restored")
  {
    log::log_capture cap;
    log::info() << "test 2";
    log::warn() << "test 3";
    log::error() << "test 4";
    REQUIRE(cap.str() == "[INFO] test 2\n"
                         "[WARNING] test 3\n"
                         "[ERROR] test 4\n");
  }
}

TEST_CASE("log::cerr is written as is")
{
  SECTION("no newlines")
  {
    log::log_capture cap;
    log::cerr() << "message without newlines";
    REQUIRE(cap.str() == "message without newlines");
  }

  SECTION("one newline")
  {
    log::log_capture cap;
    log::cerr() << "message with newline\n";
    REQUIRE(cap.str() == "message with newline\n");
  }

  SECTION("multiple newlines")
  {
    log::log_capture cap;
    log::cerr() << "message with multiple newlines\n\n";
    REQUIRE(cap.str() == "message with multiple newlines\n\n");
  }
}

TEST_CASE("transient messages are cleared by log messages")
{
  log::log_capture cap;
  log::cerr().transient() << "Transient message\n";
  log::info() << "Standard message 1";
  log::info() << "Standard message 2";

  REQUIRE(cap.str() == "Transient message\n"
                       "\r\033[K"
                       "[INFO] Standard message 1\n"
                       "[INFO] Standard message 2\n");
}

TEST_CASE("transient messages are cleared by transient messages")
{
  log::log_capture cap;
  log::cerr().transient() << "Transient message 1\n";
  log::cerr().transient() << "Transient message 2\n";
  log::info() << "Standard message 1";

  REQUIRE(cap.str() == "Transient message 1\n"
                       "\r\033[K"
                       "Transient message 2\n"
                       "\r\033[K"
                       "[INFO] Standard message 1\n");
}

TEST_CASE("transient logging not allowed for log messages")
{
  log::log_capture cap;
  REQUIRE_THROWS_AS(log::debug().transient(), assert_failed);
  REQUIRE_THROWS_AS(log::info().transient(), assert_failed);
  REQUIRE_THROWS_AS(log::warn().transient(), assert_failed);
  REQUIRE_THROWS_AS(log::error().transient(), assert_failed);
}

} // namespace adapterremoval
