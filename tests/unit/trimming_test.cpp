// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "catch.hpp"
#include "errors.hpp"
#include "fastq.hpp"
#include "testing.hpp" // for TEST_CASE, REQUIRE, ...
#include "trimming.hpp"

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// tests for `merged_reads`

TEST_CASE("merged_reads constructor ")
{
  const fastq r1{ "r1", "ACGTT" };
  const fastq r2{ "r2", "TGGGTAG" };

  REQUIRE_NOTHROW(merged_reads(fastq{}, fastq{}, 0));
  REQUIRE_NOTHROW(merged_reads(r1, r2, -1));
  REQUIRE_NOTHROW(merged_reads(r1, r2, 0));
  REQUIRE_NOTHROW(merged_reads(r1, r2, 4));
  REQUIRE_THROWS_AS(merged_reads(r1, r2, 5), assert_failed);

  REQUIRE_THROWS_AS(merged_reads(r2, r1, 0), assert_failed);
}

TEST_CASE("merged_reads trim fully overlapping")
{
  const fastq r1{ "r1", "CTACCCA" };
  const fastq r2{ "r2", "TGGGTAG" };
  merged_reads m{ r1, r2, 0 };

  SECTION("from 5p")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE(m.increment(1, 0));
    REQUIRE(m.increment(2, 0));
    REQUIRE_FALSE(m.increment(0, 0));
  }

  SECTION("from 3p")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE(m.increment(0, 1));
    REQUIRE(m.increment(0, 2));
    REQUIRE_FALSE(m.increment(0, 0));
  }

  SECTION("mixed")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE(m.increment(1, 0));
    REQUIRE(m.increment(0, 2));
    REQUIRE(m.increment(2, 1));
    REQUIRE_FALSE(m.increment(0, 0));
  }
}

TEST_CASE("merged_reads trim partially overlapping, symmetric")
{
  const fastq r1{ "r1", "CTACCCA" };
  const fastq r2{ "r2", "TGGGTAG" };
  merged_reads m{ r1, r2, 3 };

  SECTION("from 5p")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE_FALSE(m.increment(3, 0));
    REQUIRE(m.increment(2, 0));
    REQUIRE(m.increment(2, 0));
    REQUIRE_FALSE(m.increment(0, 0));
  }

  SECTION("from 3p")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE_FALSE(m.increment(0, 3));
    REQUIRE(m.increment(0, 2));
    REQUIRE(m.increment(0, 2));
    REQUIRE_FALSE(m.increment(0, 0));
  }

  SECTION("mixed")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE_FALSE(m.increment(1, 0));
    REQUIRE(m.increment(0, 4));
    REQUIRE_FALSE(m.increment(2, 0));
    REQUIRE(m.increment(2, 0));
    REQUIRE(m.increment(0, 1));
    REQUIRE_FALSE(m.increment(0, 0));
  }
}

TEST_CASE("merged_reads trim partially overlapping, asymmetric")
{
  const fastq r1{ "r1", "CTACC" };
  const fastq r2{ "r2", "TGGGTAG" };
  merged_reads m{ r1, r2, 2 };

  SECTION("from 5p")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE_FALSE(m.increment(2, 0));
    REQUIRE(m.increment(1, 0));
    REQUIRE(m.increment(2, 0));
    REQUIRE_FALSE(m.increment(0, 0));
  }

  SECTION("from 3p")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE_FALSE(m.increment(0, 3));
    REQUIRE_FALSE(m.increment(0, 1));
    REQUIRE(m.increment(0, 2));
    REQUIRE_FALSE(m.increment(0, 0));
  }

  SECTION("mixed")
  {
    REQUIRE_FALSE(m.increment(0, 0));
    REQUIRE_FALSE(m.increment(1, 0));
    REQUIRE_FALSE(m.increment(0, 3));
    REQUIRE(m.increment(2, 0));
    REQUIRE(m.increment(2, 0));
    REQUIRE_FALSE(m.increment(0, 1));
    REQUIRE(m.increment(0, 1));
    REQUIRE_FALSE(m.increment(0, 0));
  }
}

TEST_CASE("merged_reads totals are not checked")
{
  const fastq r1{ "r1", "CTACCCA" };
  const fastq r2{ "r2", "TGGGTAG" };
  merged_reads m{ r1, r2, 3 };

  REQUIRE(m.increment(1024, 0));
  REQUIRE(m.increment(0, 1024));
  REQUIRE(m.increment(1024, 1024));
}

} // namespace adapterremoval
