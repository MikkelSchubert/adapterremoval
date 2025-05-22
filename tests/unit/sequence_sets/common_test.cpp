// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "commontypes.hpp"   // for barcode_orientation, ...
#include "sequence_sets.hpp" // for read_group
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <stdexcept>         // for invalid_argument
#include <string>            // for string==

using Contains = Catch::Matchers::StdString::ContainsMatcher;

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// parse_table_orientation

TEST_CASE("parse_table_orientation for valid values", "[table_orientation]")
{
  CHECK(parse_table_orientation("unspecified") ==
        barcode_table_orientation::unspecified);
  CHECK(parse_table_orientation("forward") ==
        barcode_table_orientation::forward);
  CHECK(parse_table_orientation("reverse") ==
        barcode_table_orientation::reverse);
  CHECK(parse_table_orientation("explicit") ==
        barcode_table_orientation::explicit_);
}

TEST_CASE("parse_table_orientation ignores case and whitespace",
          "[table_orientation]")
{
  CHECK(parse_table_orientation(" \tReVeRsE\n") ==
        barcode_table_orientation::reverse);
}

TEST_CASE("parse_table_orientation throws on invalid value",
          "[table_orientation]")
{
  CHECK_THROWS_AS(parse_table_orientation(""), std::invalid_argument);
  CHECK_THROWS_AS(parse_table_orientation("revers"), std::invalid_argument);
  // The trailing underline should not be accepted/used anywhere
  CHECK_THROWS_AS(parse_table_orientation("implicit"), std::invalid_argument);
  // The trailing underline should not be accepted/used anywhere
  CHECK_THROWS_AS(parse_table_orientation("explicit_"), std::invalid_argument);
}

///////////////////////////////////////////////////////////////////////////////
// barcode_orientation* to debug strings

TEST_CASE("barcode_orientation to debug string", "[barcode_orientation]")
{
  using ::Catch::fallbackStringifier;

  CHECK(fallbackStringifier(barcode_orientation::unspecified) ==
        "barcode_orientation::unspecified");
  CHECK(fallbackStringifier(barcode_orientation::forward) ==
        "barcode_orientation::forward");
  CHECK(fallbackStringifier(barcode_orientation::reverse) ==
        "barcode_orientation::reverse");
}

TEST_CASE("barcode_table_orientation to debug string", "[barcode_orientation]")
{
  using ::Catch::fallbackStringifier;

  CHECK(fallbackStringifier(barcode_table_orientation::unspecified) ==
        "barcode_table_orientation::unspecified");
  CHECK(fallbackStringifier(barcode_table_orientation::forward) ==
        "barcode_table_orientation::forward");
  CHECK(fallbackStringifier(barcode_table_orientation::reverse) ==
        "barcode_table_orientation::reverse");
  CHECK(fallbackStringifier(barcode_table_orientation::explicit_) ==
        "barcode_table_orientation::explicit_");
}

} // namespace adapterremoval
