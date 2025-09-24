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
  CHECK(seq.reverse_complement() == "ATANGGT"_dna);
  CHECK(seq.as_string() == "ACCNTAT");
  CHECK(fallbackStringifier(seq) == "dna_sequence{'ACCNTAT'}");
}

TEST_CASE("dna_equence addition", "[dna_sequence]")
{
  CHECK("ACCGTAT"_dna + dna_sequence{} == "ACCGTAT"_dna);
  CHECK("ACCGTAT"_dna + "CCGT"_dna == "ACCGTATCCGT"_dna);
}

TEST_CASE("dna_sequence equality")
{
  CHECK(dna_sequence{} == dna_sequence{});
  CHECK_FALSE("ACCNTAT"_dna == dna_sequence{});
  CHECK("ACCNTAT"_dna == "ACCNTAT"_dna);
  CHECK_FALSE("ACCNTAT"_dna == "CGTGGTA"_dna);
  CHECK_FALSE("ACCNTAT"_dna == "CC"_dna);
}

TEST_CASE("dna_sequence literal")
{
  CHECK(""_dna == dna_sequence{});
  CHECK("ACCNTAT"_dna == "ACCNTAT"_dna);
}

} // namespace adapterremoval
