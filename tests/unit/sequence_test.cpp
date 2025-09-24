// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "sequence.hpp" // for dna_sequence
#include "testing.hpp"  // for TEST_CASE, REQUIRE, ...

namespace adapterremoval {

using Catch::fallbackStringifier;

TEST_CASE("empty sequence", "[dna_sequence]")
{
  dna_sequence seq;

  CHECK(seq.empty());
  CHECK(seq.length() == 0);
  CHECK(seq == seq);
  CHECK(seq.reverse_complement() == seq);
  CHECK(seq.as_string() == "");
  CHECK(fallbackStringifier(seq) == "dna_sequence{''}");
}

TEST_CASE("non-empty dna_sequence", "[dna_sequence]")
{
  dna_sequence seq{ "ACCNTAT" };

  CHECK(!seq.empty());
  CHECK(seq.length() == 7);
  CHECK(seq.reverse_complement() == dna_sequence{ "ATANGGT" });
  CHECK(seq.as_string() == "ACCNTAT");
  CHECK(fallbackStringifier(seq) == "dna_sequence{'ACCNTAT'}");
}

TEST_CASE("dna_equence addition", "[dna_sequence]")
{
  CHECK(dna_sequence{ "ACCGTAT" } + dna_sequence{} ==
        dna_sequence{ "ACCGTAT" });
  CHECK(dna_sequence{ "ACCGTAT" } + dna_sequence{ "CCGT" } ==
        dna_sequence{ "ACCGTATCCGT" });
}

TEST_CASE("dna_sequence equality")
{
  CHECK(dna_sequence{} == dna_sequence{});
  CHECK_FALSE(dna_sequence{ "ACCNTAT" } == dna_sequence{});
  CHECK(dna_sequence{ "ACCNTAT" } == dna_sequence{ "ACCNTAT" });
  CHECK_FALSE(dna_sequence{ "ACCNTAT" } == dna_sequence{ "CGTGGTA" });
  CHECK_FALSE(dna_sequence{ "ACCNTAT" } == dna_sequence{ "CC" });
}

} // namespace adapterremoval
