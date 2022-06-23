/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <limits>
#include <stdexcept>

#include "logging.hpp"
#include "testing.hpp"

namespace adapterremoval {

TEST_CASE("log capture timestamps")
{
  log::set_timestamps(false);

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

} // namespace adapterremoval
