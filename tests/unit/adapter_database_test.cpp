// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_database.hpp" // declarations
#include "catch.hpp"
#include "commontypes.hpp" // for read_mate
#include "errors.hpp"      // assert_failed
#include "sequence.hpp"    // for dna_sequence
#include "testing.hpp"     // for TEST_CASE, REQUIRE, ...
#include <vector>          // for vector

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// known_adapters

TEST_CASE("simple known adapter requirements")
{
  REQUIRE_NOTHROW((known_adapters{ "src", dna_sequence{} }));
  REQUIRE_NOTHROW((known_adapters{ "src", dna_sequence{ "ACGT" } }));
  REQUIRE_NOTHROW((known_adapters{ "src", dna_sequence{ "ACGT" }, {} }));
  REQUIRE_NOTHROW(
    (known_adapters{ "src", dna_sequence{ "ACGT" }, dna_sequence{ "TATA" } }));

  // Non-empty source is required
  REQUIRE_THROWS_AS((known_adapters{ {}, dna_sequence{ "ACGT" } }),
                    assert_failed);
}

TEST_CASE("complex known adapter requirements")
{
  REQUIRE_NOTHROW((known_adapters{ "src", { "ACGT" } }));
  REQUIRE_NOTHROW((known_adapters{ "src", { "ACGT" }, {} }));
  REQUIRE_NOTHROW((known_adapters{ "src", { "ACGT" }, { "TATA" } }));

  // Non-empty name is required
  REQUIRE_THROWS_AS((known_adapters{ {}, { "ACGT" } }), assert_failed);

  // Non-empty adapter 1 sequences are required
  REQUIRE_THROWS_AS((known_adapters{ "src", {} }), assert_failed);
  REQUIRE_THROWS_AS((known_adapters{ "src", { "" } }), assert_failed);

  // (Optional) adapter 2 sequences must be non-empty
  REQUIRE_THROWS_AS((known_adapters{ "src", { "ACGT" }, { "" } }),
                    assert_failed);
}

TEST_CASE("test known adapters constructor with implicit read 2")
{
  known_adapters adapters{
    "foo",
    { "ACGT", "TTTT" },
  };

  REQUIRE(adapters.name() == "foo");
  REQUIRE(adapters.adapter_1() ==
          std::vector{ dna_sequence{ "ACGT" }, dna_sequence{ "TTTT" } });
  REQUIRE(adapters.adapter_2() ==
          std::vector{ dna_sequence{ "ACGT" }, dna_sequence{ "TTTT" } });
  REQUIRE_FALSE(adapters.user_provided());
}

TEST_CASE("test known adapters constructor with explicitg read 2")
{
  known_adapters adapters{
    "foo",
    { "ACGT", "TTTT" },
    { "TATA", "GTGT" },
  };

  REQUIRE(adapters.name() == "foo");
  REQUIRE(adapters.adapter_1() ==
          std::vector{ dna_sequence{ "ACGT" }, dna_sequence{ "TTTT" } });
  REQUIRE(adapters.adapter_2() ==
          std::vector{ dna_sequence{ "TATA" }, dna_sequence{ "GTGT" } });
  REQUIRE_FALSE(adapters.user_provided());
}

TEST_CASE("test known adapters constructor for user-provided sequences")
{
  known_adapters adapters{
    "foo",
    { "ACGT", "TTTT" },
    { "TATA", "GTGT" },
    true,
  };

  REQUIRE(adapters.name() == "foo");
  REQUIRE(adapters.adapter_1() ==
          std::vector{ dna_sequence{ "ACGT" }, dna_sequence{ "TTTT" } });
  REQUIRE(adapters.adapter_2() ==
          std::vector{ dna_sequence{ "TATA" }, dna_sequence{ "GTGT" } });
  REQUIRE(adapters.user_provided());
}

////////////////////////////////////////////////////////////////////////////////
// identified_adapter

TEST_CASE("identified adapter equality")
{
  const dna_sequence seq_1{ "ACGT" };
  const dna_sequence seq_2{ "TGTA" };

  REQUIRE(identified_adapter{} == identified_adapter{});
  REQUIRE(identified_adapter{ "foo", seq_1, read_mate::_1 } ==
          identified_adapter{ "foo", seq_1, read_mate::_1 });

  REQUIRE_FALSE(identified_adapter{ "foo", seq_1, read_mate::_1 } ==
                identified_adapter{ "bar", seq_1, read_mate::_1 });
  REQUIRE_FALSE(identified_adapter{ "foo", seq_1, read_mate::_1 } ==
                identified_adapter{ "foo", seq_2, read_mate::_1 });
  REQUIRE_FALSE(identified_adapter{ "foo", seq_1, read_mate::_1 } ==
                identified_adapter{ "foo", seq_1, read_mate::_2 });
}

TEST_CASE("identified adapter to string")
{
  const identified_adapter adapter{ "Illumina ScriptSeq",
                                    dna_sequence{ "AGATCGGAAGAGCACACGTCTGAAC" },
                                    read_mate::_1 };

  REQUIRE(Catch::fallbackStringifier(adapter) ==
          "identified_adapter{source='Illumina ScriptSeq', "
          "sequence=dna_sequence{'AGATCGGAAGAGCACACGTCTGAAC'}, mate=read 1}");
}

////////////////////////////////////////////////////////////////////////////////
// adapter_database

TEST_CASE("exact match of empty sequences")
{
  adapter_database db;
  if (GENERATE(true, false)) {
    db.add_known();
  }

  const auto expected = std::pair(identified_adapter{}, identified_adapter{});
  const auto result = db.identify_exact(dna_sequence{ "" }, dna_sequence{ "" });

  REQUIRE(result == expected);
}

TEST_CASE("closest match of empty sequences")
{
  adapter_database db;
  if (GENERATE(true, false)) {
    db.add_known();
  }

  const auto expected = std::pair(identified_adapter{}, identified_adapter{});
  const auto result =
    db.identify_closest(dna_sequence{ "" }, dna_sequence{ "" });

  REQUIRE(result == expected);
}

TEST_CASE("match of SE sequence with exact match")
{
  adapter_database db;
  db.add_known();

  const dna_sequence seq{ "AGATCGGAAGAGCACACGTCTGAAC" };
  const identified_adapter match{ "Illumina ScriptSeq",
                                  dna_sequence{ "AGATCGGAAGAGCACACGTCTGAAC" },
                                  read_mate::_1 };

  SECTION("exact match")
  {
    REQUIRE(db.identify_exact(seq, dna_sequence{}) ==
            std::pair{ match, identified_adapter{} });
    REQUIRE(db.identify_exact(dna_sequence{}, seq) ==
            std::pair{ identified_adapter{}, match });
  }

  SECTION("closest match")
  {
    REQUIRE(db.identify_closest(seq, dna_sequence{}) ==
            std::pair{ match, identified_adapter{} });
    REQUIRE(db.identify_closest(dna_sequence{}, seq) ==
            std::pair{ identified_adapter{}, match });
  }
}

TEST_CASE("prefer related PE adapters")
{
  adapter_database db;
  db.add({
    { { "TTTT" }, { "TTTT" } },
    { { "TTTT" }, { "CCCC" } },
    { { "CCCC" }, { "CCCC" } },
  });

  const identified_adapter expected_1{ "User adapters #2",
                                       dna_sequence{ "TTTT" },
                                       read_mate::_1 };
  const identified_adapter expected_2{ "User adapters #2",
                                       dna_sequence{ "CCCC" },
                                       read_mate::_2 };
  const dna_sequence seq_1{ "TTTT" };
  const dna_sequence seq_2{ "CCCC" };

  REQUIRE(db.identify_exact(seq_1, seq_2) ==
          std::pair{ expected_1, expected_2 });
  REQUIRE(db.identify_closest(seq_1, seq_2) ==
          std::pair{ expected_1, expected_2 });
}

TEST_CASE("prefer related PE adapters with mismatches")
{
  adapter_database db;
  db.add({
    { { "TTTT" }, { "TTTT" } },
    { { "CCCC" }, { "TTTT" } },
    { { "TTTT" }, { "CCCC" } },
    { { "CCCC" }, { "CCCC" } },
  });

  const identified_adapter expected_1{ "User adapters #3",
                                       dna_sequence{ "TTTT" },
                                       read_mate::_1 };
  const identified_adapter expected_2{ "User adapters #3",
                                       dna_sequence{ "CCCC" },
                                       read_mate::_2 };
  const dna_sequence seq_1{ "TTTA" };
  const dna_sequence seq_2{ "CCCA" };

  REQUIRE_THROWS_AS(db.identify_exact(seq_1, seq_2), assert_failed);
  REQUIRE(db.identify_closest(seq_1, seq_2) ==
          std::pair{ expected_1, expected_2 });
}

TEST_CASE("prefer independent PE match if no better match")
{
  adapter_database db;
  db.add({
    { { "TTTT" }, { "TTTT" } },
    { { "CCCC" }, { "CCCC" } },
  });

  const identified_adapter expected_1{ "User adapters #1",
                                       dna_sequence{ "TTTT" },
                                       read_mate::_1 };
  const identified_adapter expected_2{ "User adapters #2",
                                       dna_sequence{ "CCCC" },
                                       read_mate::_2 };
  const std::pair expected{ expected_1, expected_2 };

  const dna_sequence seq_1{ "TTTT" };
  const dna_sequence seq_2{ "CCCC" };

  REQUIRE(db.identify_exact(seq_1, seq_2) == expected);
  REQUIRE(db.identify_closest(seq_1, seq_2) == expected);
}

TEST_CASE("prefer flipped PE match if no better match")
{
  adapter_database db;
  db.add({
    { { "TTTT" }, { "TTTT" } },
    { { "CCCC" }, { "CCCC" } },
    { { "CCCC" }, { "TTTT" } },
  });

  const dna_sequence seq_1{ "TTTT" };
  const dna_sequence seq_2{ "CCCC" };

  const identified_adapter expected_1{ "User adapters #3",
                                       dna_sequence{ "TTTT" },
                                       read_mate::_2 };
  const identified_adapter expected_2{ "User adapters #3",
                                       dna_sequence{ "CCCC" },
                                       read_mate::_1 };
  const std::pair expected{ expected_1, expected_2 };

  REQUIRE(db.identify_exact(seq_1, seq_2) == expected);
  REQUIRE(db.identify_closest(seq_1, seq_2) == expected);
}

TEST_CASE("prefer known SE over user provided")
{
  adapter_database db;
  db.add({ { { "AGATCGGAAGAGCACACGTCT" }, { "TTTT" } } });
  db.add_known();
  db.add({ { { "AGATCGGAAGAGCACACGTCT" }, { "TTTT" } } });

  const dna_sequence seq{ "AGATCGGAAGAGCACACGTCT" };
  const identified_adapter expected_1{ "Illumina TruSeq Ribo Profile",
                                       seq,
                                       read_mate::_1 };
  const identified_adapter expected_2{ "Illumina TruSeq Ribo Profile",
                                       seq,
                                       read_mate::_2 };

  REQUIRE(db.identify_exact(seq, {}) ==
          std::pair{ expected_1, identified_adapter{} });
  REQUIRE(db.identify_exact({}, seq) ==
          std::pair{ identified_adapter{}, expected_2 });

  REQUIRE(db.identify_closest(seq, {}) ==
          std::pair{ expected_1, identified_adapter{} });
  REQUIRE(db.identify_closest({}, seq) ==
          std::pair{ identified_adapter{}, expected_2 });
}

TEST_CASE("prefer known PE over user provided")
{
  adapter_database db;
  db.add({ { { "CTGTCTCTTATACACATCT" }, { "ATGTGTATAAGAGACA" } } });
  db.add_known();
  db.add({ { { "CTGTCTCTTATACACATCT" }, { "ATGTGTATAAGAGACA" } } });

  const dna_sequence seq_1{ "CTGTCTCTTATACACATCT" };
  const dna_sequence seq_2{ "ATGTGTATAAGAGACA" };
  const identified_adapter expected_1{
    "Illumina DNA PCR-Free Prep, Tagmentation",
    seq_1,
    read_mate::_1
  };
  const identified_adapter expected_2{
    "Illumina DNA PCR-Free Prep, Tagmentation",
    seq_2,
    read_mate::_2
  };
  const std::pair expected{ expected_1, expected_2 };

  REQUIRE(db.identify_exact(seq_1, seq_2) == expected);
  REQUIRE(db.identify_closest(seq_1, seq_2) == expected);
}

// This is arbitrary, but prevents identical additions from changing results
TEST_CASE("prefer first match")
{
  adapter_database db;
  db.add({ { { "CTGTCTCTTATACACATCT" }, { "ATGTGTATAAGAGACA" } } });
  db.add({ { { "CTGTCTCTTATACACATCT" }, { "ATGTGTATAAGAGACA" } } });

  const dna_sequence seq_1{ "CTGTCTCTTATACACATCT" };
  const dna_sequence seq_2{ "ATGTGTATAAGAGACA" };
  const identified_adapter expected_1{ "User adapters #1",
                                       seq_1,
                                       read_mate::_1 };
  const identified_adapter expected_2{ "User adapters #1",
                                       seq_2,
                                       read_mate::_2 };

  SECTION("SE")
  {
    REQUIRE(db.identify_exact(seq_1, {}) ==
            std::pair{ expected_1, identified_adapter{} });
    REQUIRE(db.identify_closest(seq_1, {}) ==
            std::pair{ expected_1, identified_adapter{} });

    REQUIRE(db.identify_exact({}, seq_2) ==
            std::pair{ identified_adapter{}, expected_2 });
    REQUIRE(db.identify_closest({}, seq_2) ==
            std::pair{ identified_adapter{}, expected_2 });
  }

  SECTION("PE")
  {
    const std::pair expected{ expected_1, expected_2 };

    REQUIRE(db.identify_exact(seq_1, seq_2) == expected);
    REQUIRE(db.identify_closest(seq_1, seq_2) == expected);
  }
}

TEST_CASE("match partially identified PE adapters")
{
  adapter_database db;
  db.add_known();

  const dna_sequence seq_1{ "AGATCGGAAGAGCACACGTCTGAAC" };
  const dna_sequence seq_2{ "AGATCGGAAGAGCGTCGTGTAGGGA" };

  const identified_adapter expected_1{ "Illumina ScriptSeq",
                                       seq_1,
                                       read_mate::_1 };
  const identified_adapter expected_2{ "Illumina ScriptSeq",
                                       seq_2,
                                       read_mate::_2 };

  REQUIRE(db.identify_exact(seq_1, {}) ==
          std::pair{ expected_1, identified_adapter{} });
  REQUIRE(db.identify_closest(seq_1, {}) ==
          std::pair{ expected_1, identified_adapter{} });

  REQUIRE(db.identify_exact({}, seq_2) ==
          std::pair{ identified_adapter{}, expected_2 });
  REQUIRE(db.identify_closest({}, seq_2) ==
          std::pair{ identified_adapter{}, expected_2 });
}

} // namespace adapterremoval
