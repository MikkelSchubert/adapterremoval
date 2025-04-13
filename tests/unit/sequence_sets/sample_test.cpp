// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "commontypes.hpp"   // for barcode_orientation
#include "errors.hpp"        // for parsing_error
#include "read_group.hpp"    // for read_group
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for read_group
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <sstream>
#include <vector> // for vector

using Contains = Catch::Matchers::StdString::ContainsMatcher;

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// Sample sequences

TEST_CASE("sample_sequences constructor", "[sample_sequences]")
{
  sample_sequences set;
  CHECK(set.barcode_1.empty());
  CHECK(set.barcode_2.empty());
}

TEST_CASE("sample_sequences explicit constructor", "[sample_sequences]")
{
  sample_sequences set{ dna_sequence{ "ACGT" }, dna_sequence{ "TTAA" } };
  CHECK(set.barcode_1 == dna_sequence{ "ACGT" });
  CHECK(set.barcode_2 == dna_sequence{ "TTAA" });
}

TEST_CASE("sample_sequences equality operator", "[sample_sequences]")
{
  sample_sequences set_1{ dna_sequence{ "ACGT" }, dna_sequence{ "TTAA" } };

  SECTION("same")
  {
    CHECK(sample_sequences{ set_1 } == set_1);
  }

  SECTION("has_read_group")
  {
    auto set_2 = set_1;
    set_2.has_read_group = true;
    CHECK_FALSE(set_1 == set_2);
  }

  SECTION("barcode_1")
  {
    auto set_2 = set_1;
    set_2.barcode_1 = dna_sequence("AAAAAAAAAAA");
    CHECK_FALSE(set_1 == set_2);
  }

  SECTION("barcode_2")
  {
    auto set_2 = set_1;
    set_2.barcode_2 = dna_sequence("AAAAAAAAAAA");
    CHECK_FALSE(set_1 == set_2);
  }

  SECTION("adapters")
  {
    auto set_2 = set_1;
    set_2.adapters = adapter_set{ { "AAAA", "TTTT" } };
    CHECK_FALSE(set_1 == set_2);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Sample sequences -- operator<<

TEST_CASE("sample_sequences to string", "[sample_sequences]")
{
  sample_sequences ss{ dna_sequence{ "ACGT" }, dna_sequence{ "TTAA" } };

  ss.has_read_group = true;
  ss.read_group_ = read_group{ "SM:foo" };
  ss.barcode_1 = dna_sequence{ "AGAA" };
  ss.barcode_2 = dna_sequence{ "TCTT" };

  std::ostringstream os;
  os << ss;

  CHECK(os.str() ==
        "sample_sequences{has_read_group=true, read_group=read_group{id='1', "
        "header='@RG\\tID:1\\tSM:foo'}, barcode_1=dna_sequence{'AGAA'}, "
        "barcode_2=dna_sequence{'TCTT'}, adapters=adapter_set{[]}}");
}

////////////////////////////////////////////////////////////////////////////////
// Sample

TEST_CASE("default sample constructor", "[sample]")
{
  sample s;

  CHECK(s == s);
  CHECK(s.name() == "");
  CHECK(s.size() == 1);
  CHECK(s.at(0) == sample_sequences{});
  CHECK(std::vector(s.begin(), s.end()) == std::vector{ sample_sequences{} });
}

TEST_CASE("explicit sample constructor", "[sample]")
{
  sample s{ "foo", dna_sequence{ "TTAC" }, dna_sequence{ "GATG" } };

  sample_sequences ss{ dna_sequence{ "TTAC" }, dna_sequence{ "GATG" } };

  CHECK(s == s);
  CHECK(s.name() == "foo");
  CHECK(s.size() == 1);
  CHECK(s.at(0) == ss);
  CHECK(std::vector(s.begin(), s.end()) == std::vector{ ss });
}

TEST_CASE("barcode2 requires barcode1 in sample constructor ", "[sample]")
{
  dna_sequence s1{};
  dna_sequence s2{ "ACGT" };

  CHECK_THROWS_AS(sample("foo", s1, s2), assert_failed);
}

TEST_CASE("add to sample", "[sample]")
{
  sample s{ "foo", dna_sequence{ "TTAC" }, dna_sequence{ "GATG" } };
  s.add(dna_sequence{ "GTGT" }, dna_sequence{ "GATC" });

  sample_sequences ss_1{ dna_sequence{ "TTAC" }, dna_sequence{ "GATG" } };
  sample_sequences ss_2{ dna_sequence{ "GTGT" }, dna_sequence{ "GATC" } };

  CHECK(s.size() == 2);
  CHECK(std::vector(s.begin(), s.end()) == std::vector{ ss_1, ss_2 });
}

TEST_CASE("set adapters for sample", "[sample]")
{
  sample s{ "foo", dna_sequence{ "TTAC" }, dna_sequence{ "GATG" } };
  s.set_adapters(adapter_set{ { "AAAA", "GGGG" } });

  sample_sequences ss{ dna_sequence{ "TTAC" }, dna_sequence{ "GATG" } };
  ss.adapters = adapter_set{ { "CATCAAAA", "GTAAGGGG" } };

  CHECK(std::vector(s.begin(), s.end()) == std::vector{ ss });
}

TEST_CASE("set read group for sample w/wo barcodes", "[sample]")
{
  SECTION("no barcodes")
  {
    sample s{ "sample", dna_sequence{}, dna_sequence{} };
    s.set_read_group(read_group{ "LB:foo" });
    CHECK(s.at(0).read_group_ == read_group{ "ID:sample\tLB:foo\tSM:sample" });
  }

  SECTION("single barcode")
  {
    sample s{ "sample", dna_sequence{ "ACAT" }, dna_sequence{} };
    s.set_read_group(read_group{ "LB:foo" });
    CHECK(s.at(0).read_group_ ==
          read_group{ "ID:sample\tLB:foo\tSM:sample\tBC:ACAT" });
  }

  SECTION("two barcodes")
  {
    sample s{ "sample", dna_sequence{ "ACAT" }, dna_sequence{ "GTGT" } };
    s.set_read_group(read_group{ "LB:foo" });
    CHECK(s.at(0).read_group_ ==
          read_group{ "ID:sample\tLB:foo\tSM:sample\tBC:ACAT-GTGT" });
  }
}

TEST_CASE("set read group for multiple sample", "[sample]")
{
  sample s{ "sample", dna_sequence{}, dna_sequence{} };
  s.add(dna_sequence{}, dna_sequence{});
  s.set_read_group(read_group{ "LB:foo" });

  CHECK(s.at(0).read_group_ == read_group{ "ID:sample.1\tLB:foo\tSM:sample" });
  CHECK(s.at(1).read_group_ == read_group{ "ID:sample.2\tLB:foo\tSM:sample" });
}

TEST_CASE("sample equality operator", "[sample]")
{
  CHECK(sample{} == sample{ "", dna_sequence{}, dna_sequence{} });
  CHECK_FALSE(sample{} == sample{ "foo", dna_sequence{}, dna_sequence{} });
  CHECK_FALSE(sample{} == sample{ "", dna_sequence{ "ACGT" }, dna_sequence{} });
  CHECK_FALSE(sample{} ==
              sample{ "", dna_sequence{ "ACGT" }, dna_sequence{ "ACGT" } });
}

////////////////////////////////////////////////////////////////////////////////
// Sample  -- operator<<

TEST_CASE("sample to string", "[sample]")
{
  sample ss{ "foo", dna_sequence{ "ACGT" }, dna_sequence{ "TTAA" } };

  std::ostringstream os;
  os << ss;

  CHECK(os.str() ==
        "sample{name='foo', barcodes=[sample_sequences{has_read_group=false, "
        "read_group=read_group{id='1', header='@RG\\tID:1'}, "
        "barcode_1=dna_sequence{'ACGT'}, barcode_2=dna_sequence{'TTAA'}, "
        "adapters=adapter_set{[]}}]}");
}

} // namespace adapterremoval
