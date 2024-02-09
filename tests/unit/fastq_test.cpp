/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "catch.hpp"       // for operator""_catch_sr, AssertionHandler
#include "commontypes.hpp" // for fastq_vec
#include "debug.hpp"       // for assert_failed
#include "fastq.hpp"       // for fastq, fastq::ntrimmed, ACGTN, ACGT
#include "fastq_enc.hpp"   // for fastq_error (ptr only), FASTQ_ENCODING_33
#include "linereader.hpp"  // for vec_reader
#include "strutils.hpp"    // for string_vec
#include <limits>          // for numeric_limits
#include <stddef.h>        // for size_t
#include <string>          // for operator==, basic_string, string, operator+
#include <utility>         // for operator==, pair, move
#include <vector>          // for vector, vector<>::const_iterator

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Helper functions

TEST_CASE("ACGT::to_index", "[fastq::*]")
{
  // The exact encoding is unimportant, but it must be unique and 2 bit
  REQUIRE(ACGT::to_index('A') <= 3);
  REQUIRE(ACGT::to_index('C') <= 3);
  REQUIRE(ACGT::to_index('G') <= 3);
  REQUIRE(ACGT::to_index('T') <= 3);

  REQUIRE(ACGT::to_index('A') != ACGT::to_index('C'));
  REQUIRE(ACGT::to_index('A') != ACGT::to_index('G'));
  REQUIRE(ACGT::to_index('A') != ACGT::to_index('T'));
  REQUIRE(ACGT::to_index('C') != ACGT::to_index('G'));
  REQUIRE(ACGT::to_index('C') != ACGT::to_index('T'));
  REQUIRE(ACGT::to_index('G') != ACGT::to_index('T'));
}

TEST_CASE("ACGT::to_value", "[fastq::*]")
{
  REQUIRE(ACGT::to_value(ACGT::to_index('A')) == 'A');
  REQUIRE(ACGT::to_value(ACGT::to_index('C')) == 'C');
  REQUIRE(ACGT::to_value(ACGT::to_index('G')) == 'G');
  REQUIRE(ACGT::to_value(ACGT::to_index('T')) == 'T');
}

TEST_CASE("ACGTN::to_index", "[fastq::*]")
{
  // The exact encoding is unimportant, but it must be unique in the range 0-4
  REQUIRE(ACGTN::to_index('A') <= 4);
  REQUIRE(ACGTN::to_index('C') <= 4);
  REQUIRE(ACGTN::to_index('G') <= 4);
  REQUIRE(ACGTN::to_index('T') <= 4);
  REQUIRE(ACGTN::to_index('N') <= 4);

  REQUIRE(ACGTN::to_index('A') != ACGTN::to_index('C'));
  REQUIRE(ACGTN::to_index('A') != ACGTN::to_index('G'));
  REQUIRE(ACGTN::to_index('A') != ACGTN::to_index('T'));
  REQUIRE(ACGTN::to_index('A') != ACGTN::to_index('N'));
  REQUIRE(ACGTN::to_index('C') != ACGTN::to_index('G'));
  REQUIRE(ACGTN::to_index('C') != ACGTN::to_index('T'));
  REQUIRE(ACGTN::to_index('C') != ACGTN::to_index('N'));
  REQUIRE(ACGTN::to_index('G') != ACGTN::to_index('T'));
  REQUIRE(ACGTN::to_index('G') != ACGTN::to_index('N'));
  REQUIRE(ACGTN::to_index('T') != ACGTN::to_index('N'));
}

TEST_CASE("ACGTN::to_value", "[fastq::*]")
{
  REQUIRE(ACGTN::to_value(ACGTN::to_index('A')) == 'A');
  REQUIRE(ACGTN::to_value(ACGTN::to_index('C')) == 'C');
  REQUIRE(ACGTN::to_value(ACGTN::to_index('G')) == 'G');
  REQUIRE(ACGTN::to_value(ACGTN::to_index('T')) == 'T');
  REQUIRE(ACGTN::to_value(ACGTN::to_index('N')) == 'N');
}

///////////////////////////////////////////////////////////////////////////////
// Default constructor

TEST_CASE("default_constructor", "[fastq::fastq]")
{
  const fastq record;
  REQUIRE(record.header() == "");
  REQUIRE(record.sequence() == "");
  REQUIRE(record.qualities() == "");
}

///////////////////////////////////////////////////////////////////////////////
// Primary constructor

TEST_CASE("constructor_empty_fields", "[fastq::fastq]")
{
  const fastq record("", "", "");
  REQUIRE(record.header() == "");
  REQUIRE(record.sequence() == "");
  REQUIRE(record.qualities() == "");
}

TEST_CASE("constructor_simple_record_phred_33_encoded", "[fastq::fastq]")
{
  const fastq record("record_1", "ACGAGTCA", "!7BF8DGI");
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "ACGAGTCA");
  REQUIRE(record.qualities() == "!7BF8DGI");
}

TEST_CASE("constructor_simple_record_phred_64_encoded", "[fastq::fastq]")
{
  const fastq record("record_2", "ACGAGTCA", "@VaeWcfh", FASTQ_ENCODING_64);
  REQUIRE(record.header() == "record_2");
  REQUIRE(record.sequence() == "ACGAGTCA");
  REQUIRE(record.qualities() == "!7BF8DGI");
}

TEST_CASE("constructor_simple_record_phred_solexa_encoded", "[fastq::fastq]")
{
  const fastq record(
    "record_3", "AAACGAGTCA", ";h>S\\TCDUJ", FASTQ_ENCODING_SOLEXA);
  REQUIRE(record.header() == "record_3");
  REQUIRE(record.sequence() == "AAACGAGTCA");
  REQUIRE(record.qualities() == "\"I#4=5&&6+");
}

TEST_CASE("constructor_simple_record_lowercase_to_uppercase", "[fastq::fastq]")
{
  const fastq record("record_1", "AnGaGtcA", "!7BF8DGI");
  REQUIRE(record.sequence() == "ANGAGTCA");
}

TEST_CASE("constructor_score_boundries_phred_33", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", "!!\"", FASTQ_ENCODING_33));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", " !\"", FASTQ_ENCODING_33),
                    fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "IJN", FASTQ_ENCODING_33));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "IJO", FASTQ_ENCODING_33), fastq_error);
}

TEST_CASE("constructor_score_boundries_phred_64", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", "@@A", FASTQ_ENCODING_64));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "?@A", FASTQ_ENCODING_64), fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "gh~", FASTQ_ENCODING_64));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "gh\x7f", FASTQ_ENCODING_64),
                    fastq_error);
}

TEST_CASE("constructor_score_boundries_phred_64 suggests solexa")
{
  REQUIRE_THROWS_WITH(fastq("Rec", "A", ":", FASTQ_ENCODING_64),
                      !Catch::Contains("Solexa format"));
  REQUIRE_THROWS_WITH(fastq("Rec", "A", ";", FASTQ_ENCODING_64),
                      Catch::Contains("Solexa format"));
}

TEST_CASE("constructor_score_boundries_solexa", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", ";;<", FASTQ_ENCODING_SOLEXA));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", ":;<", FASTQ_ENCODING_SOLEXA),
                    fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "fgh", FASTQ_ENCODING_SOLEXA));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "fgi", FASTQ_ENCODING_SOLEXA),
                    fastq_error);
}

TEST_CASE("constructor_score_boundries_phred_sam", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", "!!\"", FASTQ_ENCODING_SAM));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", " !\"", FASTQ_ENCODING_SAM),
                    fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "gh~", FASTQ_ENCODING_SAM));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "gh\x7f", FASTQ_ENCODING_SAM),
                    fastq_error);
}

TEST_CASE("constructor_score_boundries_ignored", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", "!!\"", FASTQ_ENCODING_SAM));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", " !\"", FASTQ_ENCODING_SAM),
                    fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "gh~", FASTQ_ENCODING_SAM));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "gh\x7f", FASTQ_ENCODING_SAM),
                    fastq_error);
}

TEST_CASE("constructor_field_lengths", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Name", "CAT", "IJJ"));
  // A non-empty sequence must be specified
  REQUIRE_THROWS_AS(fastq("Name", "", "IJJ"), fastq_error);
  // A non-empty quality string must be specified
  REQUIRE_THROWS_AS(fastq("Name", "CAT", ""), fastq_error);
  // And the length of each must be the same
  REQUIRE_THROWS_AS(fastq("Name", "CA", "IJJ"), fastq_error);
  REQUIRE_THROWS_AS(fastq("Name", "CAT", "IJ"), fastq_error);
}

TEST_CASE("constructor_invalid_nucleotides", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Name", "CATT", "IJJI"));
  // Non-alpha characters are not allowed
  REQUIRE_THROWS_AS(fastq("Name", "CAT!", "IJJI"), fastq_error);
  // Numeric charecters are not allowed
  REQUIRE_THROWS_AS(fastq("Name", "CAT7", "IJJI"), fastq_error);
  // But neither are non acgtn/ACGTN allowed
  REQUIRE_THROWS_AS(fastq("Name", "CATS", "IJJI"), fastq_error);
  REQUIRE_THROWS_AS(fastq("Name", "CATs", "IJJI"), fastq_error);
}

///////////////////////////////////////////////////////////////////////////////
// Constructor without qualities

TEST_CASE("constructor_no_qualities", "[fastq::fastq]")
{
  const fastq record("record_1", "ACGT");
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "ACGT");
  REQUIRE(record.qualities() == "!!!!");
}

TEST_CASE("constructor_no_qualities_no_sequence", "[fastq::fastq]")
{
  const fastq record("record_1", "");
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "");
  REQUIRE(record.qualities() == "");
}

///////////////////////////////////////////////////////////////////////////////
// Move constructor

TEST_CASE("move_constructor", "[fastq::fastq]")
{
  fastq record1("record_1", "ACGT", "1234");
  fastq record2(std::move(record1));

  REQUIRE(record1.header() == "");
  REQUIRE(record1.sequence() == "");
  REQUIRE(record1.qualities() == "");

  REQUIRE(record2.header() == "record_1");
  REQUIRE(record2.sequence() == "ACGT");
  REQUIRE(record2.qualities() == "1234");
}

///////////////////////////////////////////////////////////////////////////////
// Move constructor

TEST_CASE("assignment_operator", "[fastq::fastq]")
{
  fastq record1("record_1", "ACGT", "1234");
  fastq record2;
  record2 = record1;

  REQUIRE(record1.header() == "record_1");
  REQUIRE(record1.sequence() == "ACGT");
  REQUIRE(record1.qualities() == "1234");

  REQUIRE(record2.header() == "record_1");
  REQUIRE(record2.sequence() == "ACGT");
  REQUIRE(record2.qualities() == "1234");
}

///////////////////////////////////////////////////////////////////////////////
// misc properties

TEST_CASE("name", "[fastq::fastq]")
{
  REQUIRE(fastq("name", "", "").name() == "name");
  REQUIRE(fastq("name meta", "", "").name() == "name");
  REQUIRE(fastq("name meta more", "", "").name() == "name");
}

TEST_CASE("length", "[fastq::fastq]")
{
  REQUIRE(fastq("record_1", "", "").length() == 0);
  REQUIRE(fastq("record_1", "A", "G").length() == 1);
  REQUIRE(fastq("record_1", "AC", "!B").length() == 2);
  REQUIRE(fastq("record_1", "ACG", "!7B").length() == 3);
}

TEST_CASE("count_ns", "[fastq::fastq]")
{
  REQUIRE(fastq("Rec", "ACGTA", "IJIJI").count_ns() == 0);
  REQUIRE(fastq("Rec", "ANGTA", "IJIJI").count_ns() == 1);
  REQUIRE(fastq("Rec", "ANGTN", "IJIJI").count_ns() == 2);
  REQUIRE(fastq("Rec", "ANGNN", "IJIJI").count_ns() == 3);
  REQUIRE(fastq("Rec", "NNGNN", "IJIJI").count_ns() == 4);
  REQUIRE(fastq("Rec", "NNNNN", "IJIJI").count_ns() == 5);
}

TEST_CASE("complexity", "[fastq::fastq]")
{
  using Catch::WithinAbs;

  REQUIRE_THAT(fastq("c", "").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "A").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "AA").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "AAA").complexity(), WithinAbs(0.0, 1e-6));

  REQUIRE_THAT(fastq("c", "N").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "NA").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "NN").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "NAN").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "NNN").complexity(), WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(fastq("c", "NANT").complexity(), WithinAbs(1.0 / 3.0, 1e-6));

  REQUIRE_THAT(fastq("c", "AT").complexity(), WithinAbs(1.0 / 1.0, 1e-6));
  REQUIRE_THAT(fastq("c", "ACGT").complexity(), WithinAbs(3.0 / 3.0, 1e-6));
  REQUIRE_THAT(fastq("c", "AACGT").complexity(), WithinAbs(3.0 / 4.0, 1e-6));
  REQUIRE_THAT(fastq("c", "ANCGT").complexity(), WithinAbs(3.0 / 4.0, 1e-6));
  REQUIRE_THAT(fastq("c", "AGACCGT").complexity(), WithinAbs(5.0 / 6.0, 1e-6));
  REQUIRE_THAT(fastq("c", "AGACCGGT").complexity(), WithinAbs(5.0 / 7.0, 1e-6));
  REQUIRE_THAT(fastq("c", "AGACCGGTN").complexity(),
               WithinAbs(5.0 / 8.0, 1e-6));
}

///////////////////////////////////////////////////////////////////////////////
// trim_trailing_bases

TEST_CASE("trim_trailing_bases__empty_record", "[fastq::fastq]")
{
  fastq record("Empty", "", "");
  const fastq::ntrimmed expected(0, 0);
  REQUIRE(record.trim_trailing_bases(true, 10) == expected);
  REQUIRE(record == fastq("Empty", "", ""));
}

TEST_CASE("trim_trailing_bases__trim_nothing", "[fastq::fastq]")
{
  const fastq reference("Rec", "NNNNN", "!!!!!");
  const fastq::ntrimmed expected(0, 0);
  fastq record = reference;

  // Trim neither Ns nor low Phred score bases
  REQUIRE(record.trim_trailing_bases(false, -1) == expected);
  REQUIRE(record == reference);
}

TEST_CASE("trim_trailing_bases__trim_ns", "[fastq::fastq]")
{
  fastq record("Rec", "NNANT", "23456");
  const fastq expected_record("Rec", "ANT", "456");
  const fastq::ntrimmed expected_ntrim(2, 0);

  REQUIRE(record.trim_trailing_bases(true, -1) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_trailing_bases__trim_trailing_bases", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "TN", "%$");
  const fastq::ntrimmed expected_ntrim(0, 3);
  fastq record("Rec", "TNANT", "%$#!\"");

  REQUIRE(record.trim_trailing_bases(false, 2) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_trailing_bases__trim_mixed", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "TAG", "$12");
  const fastq::ntrimmed expected_ntrim(3, 2);
  fastq record("Rec", "NTNTAGNT", "1!#$12#\"");

  REQUIRE(record.trim_trailing_bases(true, 2) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_trailing_bases__trim_mixed__no_low_quality_bases",
          "[fastq::fastq]")
{
  const fastq expected_record("Rec", "ACTTAG", "12I$12");
  const fastq::ntrimmed expected_ntrim(0, 0);
  fastq record = expected_record;

  REQUIRE(record.trim_trailing_bases(true, 2) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_trailing_bases__trim_everything", "[fastq::fastq]")
{
  fastq record("Rec", "TAG", "!!!");
  const fastq expected_record = fastq("Rec", "", "");
  const fastq::ntrimmed expected_ntrim(0, 3);
  REQUIRE(record.trim_trailing_bases(true, 2) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_trailing_bases__trim_3p__trim_ns", "[fastq::fastq]")
{
  fastq record("Rec", "NNATN", "23456");
  const fastq expected_record("Rec", "NNAT", "2345");
  const fastq::ntrimmed expected_ntrim(0, 1);

  REQUIRE(record.trim_trailing_bases(true, -1, true) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_trailing_bases__trim_3p__trim_trailing_bases", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "TN", "!$");
  const fastq::ntrimmed expected_ntrim(0, 3);
  fastq record("Rec", "TNANT", "!$#!\"");

  REQUIRE(record.trim_trailing_bases(false, 2, true) == expected_ntrim);
  REQUIRE(record == expected_record);
}

///////////////////////////////////////////////////////////////////////////////
// trim_windowed_bases

TEST_CASE("Window trimming with invalid parameters", "[fastq::windows]")
{
  const std::vector<double> values = {
    -1.0, std::numeric_limits<double>::quiet_NaN()
  };
  for (const auto& value : values) {
    fastq record("Rec", "TAGTGACAT", "111111111");
    REQUIRE_THROWS_AS(record.trim_windowed_bases(false, -1, value),
                      assert_failed);
  }
}

TEST_CASE("Window trimming empty reads", "[fastq::windows]")
{
  const std::vector<double> values = { 1, 0.1, 3 };
  for (const auto& value : values) {
    fastq record("Empty", "", "");
    const fastq::ntrimmed expected(0, 0);
    REQUIRE(record.trim_windowed_bases(true, 10, value) == expected);
    REQUIRE(record == fastq("Empty", "", ""));
  }
}

TEST_CASE("Window trimming entire read", "[fastq::windows]")
{
  const std::vector<double> values = { 1, 0.2, 4, 10 };
  for (const auto& value : values) {
    fastq record("Rec", "TAGTGACAT", "111111111");
    const fastq expected_record = fastq("Rec", "", "");
    const fastq::ntrimmed expected_ntrim(9, 0);
    REQUIRE(record.trim_windowed_bases(false, '2' - '!', value) ==
            expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming nothing", "[fastq::windows]")
{
  const std::vector<double> values = { 0, 1, 0.2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  for (const auto& value : values) {
    fastq record("Rec", "TAGTGACAT", "111111111");
    const fastq expected_record = record;
    const fastq::ntrimmed expected_ntrim(0, 0);
    REQUIRE(record.trim_windowed_bases(false, -1, value) == expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming of Ns", "[fastq::windows]")
{
  const std::vector<double> values = { 1, 0.2 };
  for (const auto& value : values) {
    fastq record("Rec", "NNATNT", "234567");
    const fastq expected_record("Rec", "ATNT", "4567");
    const fastq::ntrimmed expected_ntrim(2, 0);

    REQUIRE(record.trim_windowed_bases(true, -1, value) == expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming of Ns, final window contains Ns", "[fastq::windows]")
{
  const std::vector<double> values = { 2, 3, 4 };
  for (const auto& value : values) {
    fastq record("Rec", "NNATNT", "234567");
    const fastq expected_record("Rec", "AT", "45");
    const fastq::ntrimmed expected_ntrim(2, 2);

    REQUIRE(record.trim_windowed_bases(true, -1, value) == expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming of Ns, no valid window found", "[fastq::fastq]")
{
  fastq record("Rec", "NNATNT", "234567");
  const fastq expected_record("Rec", "", "");
  const fastq::ntrimmed expected_ntrim(6, 0);

  REQUIRE(record.trim_windowed_bases(true, -1, 5) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("Window trimming 5p only, 1bp", "[fastq::windows]")
{
  const std::vector<double> values = { 1, 0.1 };
  for (const auto& value : values) {
    fastq record("Rec", "TAACGATCCG", "0123456789");
    const fastq expected_record("Rec", "CGATCCG", "3456789");
    const fastq::ntrimmed expected_ntrim(3, 0);

    REQUIRE(record.trim_windowed_bases(false, '2' - '!', value) ==
            expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming 5p only, 2bp", "[fastq::windows]")
{
  const std::vector<double> values = { 2, 0.2 };
  for (const auto& value : values) {
    fastq record("Rec", "TAACGATCCG", "0123456789");
    const fastq expected_record("Rec", "CGATCCG", "3456789");
    const fastq::ntrimmed expected_ntrim(3, 0);

    REQUIRE(record.trim_windowed_bases(false, '2' - '!', value) ==
            expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming 5p only, inclusive lower bound", "[fastq::windows]")
{
  const std::vector<double> values = { 2, 3 };
  for (const auto& value : values) {
    fastq record("Rec", "TAACGATCCG", "0123126789");
    const fastq expected_record("Rec", "TCCG", "6789");
    const fastq::ntrimmed expected_ntrim(6, 0);

    REQUIRE(record.trim_windowed_bases(false, '2' - '!', value) ==
            expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming 3p only, inclusive lower bound", "[fastq::windows]")
{
  const std::vector<double> values = { 3 };
  for (const auto& value : values) {
    fastq record("Rec", "TAACGATCCG", "9876312333");
    const fastq expected_record("Rec", "TAACG", "98763");
    const fastq::ntrimmed expected_ntrim(0, 5);

    REQUIRE(record.trim_windowed_bases(false, '2' - '!', value) ==
            expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming big and small windows #1", "[fastq::windows]")
{
  const std::vector<double> values = { 0, 0.01, 20 };
  for (const auto& value : values) {
    fastq record("Rec", "TAACGATC", "23456789");
    const fastq expected_record("Rec", "", "");
    const fastq::ntrimmed expected_ntrim(8, 0);

    CHECK(record.trim_windowed_bases(false, '2' - '!', value) ==
          expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming big and small windows #2", "[fastq::windows]")
{
  const std::vector<double> values = { 0, 0.01, 20 };
  for (const auto& value : values) {
    fastq record("Rec", "TAACGATC", "23456780");
    const fastq expected_record("Rec", "TAACGAT", "2345678");
    const fastq::ntrimmed expected_ntrim(0, 1);

    REQUIRE(record.trim_windowed_bases(false, '1' - '!', value) ==
            expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("Window trimming last trailing window", "[fastq::windows]")
{
  const std::vector<double> values = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  for (const auto& value : values) {
    fastq record("Rec", "TAACGATCC", "234567811");
    const fastq expected_record("Rec", "TAACGAT", "2345678");
    const fastq::ntrimmed expected_ntrim(0, 2);

    REQUIRE(record.trim_windowed_bases(false, '1' - '!', value) ==
            expected_ntrim);
    REQUIRE(record == expected_record);
  }
}

TEST_CASE("trim_windowed_bases__trim_window", "[fastq::fastq]")
{
  // Should trim starting at the window of low quality bases in the middle
  // even with high qual bases at the end.
  fastq record("Rec", "NNAAAAAAAAATNNNNNNNA", "##EEEEEEEEEE#######E");
  const fastq expected_record = fastq("Rec", "AAAAAAAAAT", "EEEEEEEEEE");
  const fastq::ntrimmed expected_ntrim(2, 8);
  REQUIRE(record.trim_windowed_bases(true, 10, 5) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_windowed_bases__trim_3p", "[fastq::fastq]")
{
  fastq record("Rec", "NNAAAAAAAAATNNNNNNN", "##EEEEEEEEEE#######");
  const fastq expected_record = fastq("Rec", "NNAAAAAAAAAT", "##EEEEEEEEEE");
  const fastq::ntrimmed expected_ntrim(0, 7);
  REQUIRE(record.trim_windowed_bases(false, 10, 5, true) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_windowed_bases__trim_3p__trim_ns", "[fastq::fastq]")
{
  fastq record("Rec", "NNAAAAAAAAATNNNNNNN", "##EEEEEEEEEE#######");
  const fastq expected_record = fastq("Rec", "NNAAAAAAAAAT", "##EEEEEEEEEE");
  const fastq::ntrimmed expected_ntrim(0, 7);
  REQUIRE(record.trim_windowed_bases(true, 10, 5, true) == expected_ntrim);
  REQUIRE(record == expected_record);
}

TEST_CASE("trim_windowed_bases__trim_3p__reversed", "[fastq::fastq]")
{
  fastq record("Rec", "NNNNNNNTAAAAAAAAANN", "#######EEEEEEEEEE##");
  const fastq expected_record =
    fastq("Rec", "NNNNNNNTAAAAAAAAA", "#######EEEEEEEEEE");
  const fastq::ntrimmed expected_ntrim(0, 2);
  REQUIRE(record.trim_windowed_bases(true, 10, 5, true) == expected_ntrim);
  REQUIRE(record == expected_record);
}

///////////////////////////////////////////////////////////////////////////////
// mott_trimming

TEST_CASE("Mott trimming empty sequence yields empty sequence")
{
  const fastq expected_record("Rec", "", "");
  fastq record = expected_record;

  REQUIRE(record.mott_trimming(0.05) == fastq::ntrimmed());
  REQUIRE(record == expected_record);
}

TEST_CASE("Mott doesn't trim high quality bases")
{
  const fastq expected("Rec", "GGAACGTG", "JJJJJJJJ");
  auto record = expected;

  REQUIRE(record.mott_trimming(0.05) == fastq::ntrimmed());
  REQUIRE(record == expected);
}

TEST_CASE("Mott trimming treats Ns as low quality")
{
  fastq record("Rec", "NNNNNNNN", "JJJJJJJJ");

  REQUIRE(record.mott_trimming(0.05) == fastq::ntrimmed(0, 8));
  REQUIRE(record == fastq("Rec", ""));
}

TEST_CASE("Mott trimming both ends by default")
{
  fastq record("Rec", "GGAACGTGTACTAT", "####JJJJJ#####");

  REQUIRE(record.mott_trimming(0.05) == fastq::ntrimmed(4, 5));
  REQUIRE(record == fastq("Rec", "CGTGT", "JJJJJ"));
}

TEST_CASE("Mott trimming trims only 3' if preserve3p")
{
  fastq record("Rec", "GGAACGTGTACTAT", "####JJJJJ#####");

  REQUIRE(record.mott_trimming(0.05, true) == fastq::ntrimmed(0, 5));
  REQUIRE(record == fastq("Rec", "GGAACGTGT", "####JJJJJ"));
}

TEST_CASE("Mott trimming shorter segments if first")
{
  fastq record("Rec", "GGAACGTGTACTATGGAACGTG", "####JJJJ#####JJJJJ####");

  REQUIRE(record.mott_trimming(0.05) == fastq::ntrimmed(13, 4));
  REQUIRE(record == fastq("Rec", "TGGAA", "JJJJJ"));
}

TEST_CASE("Mott trimming shorter segments if last")
{
  fastq record("Rec", "GGAACGTGTACTATGGAACGTG", "####JJJJJ####JJJJ#####");

  REQUIRE(record.mott_trimming(0.05) == fastq::ntrimmed(4, 13));
  REQUIRE(record == fastq("Rec", "CGTGT", "JJJJJ"));
}

///////////////////////////////////////////////////////////////////////////////
// poly_x_trimming

std::string
poly_x_trimming(const std::string& nucleotides,
                const size_t min_length,
                const std::string sequence)
{
  fastq record("Test", sequence);
  const auto result = record.poly_x_trimming(nucleotides, min_length);
  REQUIRE(record.length() + result.second == sequence.length());

  return record.sequence();
}

TEST_CASE("trim_poly_x on empty sequence and with empty X")
{
  REQUIRE(poly_x_trimming("", 10, "") == "");
}

TEST_CASE("trim_poly_x on empty sequence and with X")
{
  REQUIRE(poly_x_trimming("ACGT", 10, "") == "");
}

TEST_CASE("trim_poly_x on sequence with empty X")
{
  const std::string seq = "TTTTTTTTTTTTTTTTTTTT";
  REQUIRE(poly_x_trimming("", 10, seq) == seq);
}

TEST_CASE("trim_poly_x no trimming of too short sequences")
{
  REQUIRE(poly_x_trimming("T", 5, "") == "");
  REQUIRE(poly_x_trimming("T", 5, "T") == "T");
  REQUIRE(poly_x_trimming("T", 5, "TT") == "TT");
  REQUIRE(poly_x_trimming("T", 5, "TTT") == "TTT");
  REQUIRE(poly_x_trimming("T", 5, "TTTT") == "TTTT");
  REQUIRE(poly_x_trimming("T", 5, "TTTTT") == "");
}

TEST_CASE("trim_poly_x trimming until, but not including last mismtches")
{
  REQUIRE(poly_x_trimming("T", 10, "TATTTTATTTTTTTT") == "TA");
  REQUIRE(poly_x_trimming("T", 10, "TAATTTTTTTTTTTT") == "TAA");
}

TEST_CASE("trim_poly_x trims only exact matches")
{
  const std::string seq = "TATTTTATTTTTTTT";
  REQUIRE(poly_x_trimming("A", 10, seq) == seq);
  REQUIRE(poly_x_trimming("C", 10, seq) == seq);
  REQUIRE(poly_x_trimming("G", 10, seq) == seq);
  REQUIRE(poly_x_trimming("N", 10, seq) == seq);
  REQUIRE(poly_x_trimming("t", 10, seq) == seq);
  REQUIRE(poly_x_trimming("T", 10, seq) == "TA");
}

TEST_CASE("trim_poly_x accept one mismatch per 8 bases")
{
  REQUIRE(poly_x_trimming("A", 10, "AAAAAAAAAA") == "");
  REQUIRE(poly_x_trimming("A", 10, "AAAAAATAAA") == "");
  REQUIRE(poly_x_trimming("A", 10, "TTAAAAAAAAA") == "TT");
  REQUIRE(poly_x_trimming("A", 10, "AATAAATAAA") == "AATAAATAAA");
  REQUIRE(poly_x_trimming("A", 10, "AAATAAAAAAAAATAAA") == "AAAT");
  REQUIRE(poly_x_trimming("A", 10, "AATAAAAAAAAAATAAA") == "AAT");
  REQUIRE(poly_x_trimming("A", 10, "ATAAAAAAAAAAATAAA") == "");
}

TEST_CASE("trim_poly_x accept early mismatches if in window")
{
  REQUIRE(poly_x_trimming("A", 10, "AAAAAAAAAAAAATAAA") == "");
  REQUIRE(poly_x_trimming("A", 10, "AAAAAAAAAAAAATAAA") == "");
}

TEST_CASE("trim_poly_x stops trims best candidate")
{
  REQUIRE(poly_x_trimming("ACGT", 10, "AAAAAAAATT") == "AAAAAAAATT");
  REQUIRE(poly_x_trimming("ACGT", 10, "AAAAAAAAAT") == "");

  REQUIRE(poly_x_trimming("ACGT", 10, "GCCCCCTCCC") == "GCCCCCTCCC");
  REQUIRE(poly_x_trimming("ACGT", 10, "CCCCCCTCCC") == "");

  REQUIRE(poly_x_trimming("ACGT", 10, "GGGGGGGGAA") == "GGGGGGGGAA");
  REQUIRE(poly_x_trimming("ACGT", 10, "AGGGGGGGGG") == "A");

  REQUIRE(poly_x_trimming("ACGT", 10, "AATTTTTTTT") == "AATTTTTTTT");
  REQUIRE(poly_x_trimming("ACGT", 10, "TATTTTTTTT") == "");
}

TEST_CASE("trim_poly_x doesn't count Ns in alignnment")
{
  REQUIRE(poly_x_trimming("C", 5, "CCCCCC") == "");
  REQUIRE(poly_x_trimming("C", 5, "CTCCCC") == "CTCCCC");
  REQUIRE(poly_x_trimming("C", 5, "CNCCCC") == "");
  REQUIRE(poly_x_trimming("C", 5, "CNCCNC") == "CNCCNC");

  REQUIRE(poly_x_trimming("C", 10, "CCCCNCCCCC") == "CCCCNCCCCC");
  REQUIRE(poly_x_trimming("C", 10, "CCCCCNCCCCC") == "");
  REQUIRE(poly_x_trimming("C", 10, "TCCCCNCCCCC") == "T");

  REQUIRE(poly_x_trimming("C", 10, "CCCCTCCCCCCCCCCCTCCC") == "");
  REQUIRE(poly_x_trimming("C", 10, "CCCCTCCCCCNCCCCCTCCC") == "CCCCT");
}

TEST_CASE("trim_poly_x only trims trailing Ns after matches")
{
  REQUIRE(poly_x_trimming("C", 5, "NCCCNCCCNCC") == "");
  REQUIRE(poly_x_trimming("C", 5, "NNCCCCCCCC") == "");
  REQUIRE(poly_x_trimming("C", 5, "NNTCCCCCCC") == "NNT");
}

TEST_CASE("trim_poly_x returns nucleotide trimmed or N")
{
  fastq record("Test", "CCCCC");

  using pair = std::pair<char, size_t>;

  REQUIRE(fastq(record).poly_x_trimming("A", 4) == pair{ 'N', 0 });
  REQUIRE(fastq(record).poly_x_trimming("C", 4) == pair{ 'C', 5 });
  REQUIRE(fastq(record).poly_x_trimming("CAGT", 4) == pair{ 'C', 5 });
  REQUIRE(fastq(record).poly_x_trimming("AGTC", 4) == pair{ 'C', 5 });
  REQUIRE(fastq(record).poly_x_trimming("ATG", 4) == pair{ 'N', 0 });
}

///////////////////////////////////////////////////////////////////////////////
// Truncate

TEST_CASE("truncate_empty", "[fastq::fastq]")
{
  fastq record("Empty", "", "");
  record.truncate(0, 10);
  REQUIRE(record == fastq("Empty", "", ""));
}

TEST_CASE("truncate_zero_bases", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "ACTTAG", "12I$12");
  fastq current_record = expected_record;
  current_record.truncate();
  REQUIRE(current_record == expected_record);
}

TEST_CASE("truncate_all_bases", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "", "");
  fastq current_record("Rec", "ACTTAG", "12I$12");
  current_record.truncate(1, 0);
  REQUIRE(current_record == expected_record);
}

TEST_CASE("truncate_5p", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "TTAG", "I$12");
  fastq current_record("Rec", "ACTTAG", "12I$12");
  current_record.truncate(2);
  REQUIRE(current_record == expected_record);
}

TEST_CASE("truncate_3p", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "ACT", "12I");
  fastq current_record("Rec", "ACTTAG", "12I$12");
  current_record.truncate(0, 3);
  REQUIRE(current_record == expected_record);
}

TEST_CASE("truncate_middle", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "TTA", "I$1");
  fastq current_record("Rec", "ACTTAG", "12I$12");
  current_record.truncate(2, 3);
  REQUIRE(current_record == expected_record);
}

TEST_CASE("truncate_len_higher_than_n_bases", "[fastq::fastq]")
{
  const fastq expected_record("Rec", "TTAG", "I$12");
  fastq current_record("Rec", "ACTTAG", "12I$12");
  current_record.truncate(2, 1024);
  REQUIRE(current_record == expected_record);
}

TEST_CASE("truncate_pos_after_last_base", "[fastq::fastq]")
{
  // Same behavior as string::substr
  fastq current_record("Rec", "ACTTAG", "12I$12");
  REQUIRE_NOTHROW(current_record.truncate(6));
  REQUIRE_THROWS_AS(current_record.truncate(7), assert_failed);
}

///////////////////////////////////////////////////////////////////////////////
// Reverse complement

TEST_CASE("reverse_complement__empty", "[fastq::fastq]")
{
  const fastq expected = fastq("Empty", "", "");
  fastq result = fastq("Empty", "", "");
  result.reverse_complement();
  REQUIRE(result == expected);
}

TEST_CASE("reverse_complement", "[fastq::fastq]")
{
  const fastq expected = fastq("Rec", "TACAGANGTN", "0123456789");
  fastq result = fastq("Rec", "NACNTCTGTA", "9876543210");
  result.reverse_complement();
  REQUIRE(result == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Adding prefixes to the header

TEST_CASE("add_prefix_to_name", "[fastq::fastq]")
{
  const fastq expected("not_my_header", "ACGTA", "12345");
  fastq record("my_header", "ACGTA", "12345");
  record.add_prefix_to_name("not_");
  REQUIRE(record == expected);
}

TEST_CASE("add_prefix_to_header__empty_prefix", "[fastq::fastq]")
{
  const fastq expected("my_header", "ACGTA", "12345");
  fastq record = expected;
  record.add_prefix_to_name("");
  REQUIRE(record == expected);
}

TEST_CASE("add_prefix_to_header__header", "[fastq::fastq]")
{
  const fastq expected("new_header", "ACGTA", "12345");
  fastq record("", "ACGTA", "12345");
  record.add_prefix_to_name("new_header");
  REQUIRE(record == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Reading from stream

TEST_CASE("simple_fastq_record_1", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("ACGAGTCA");
  lines.push_back("+");
  lines.push_back("!7BF8DGI");
  vec_reader reader(lines);

  fastq record;
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "ACGAGTCA");
  REQUIRE(record.qualities() == "!7BF8DGI");
}

TEST_CASE("simple_fastq_record_2", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("ACGAGTCA");
  lines.push_back("+");
  lines.push_back("!7BF8DGI");
  lines.push_back("@record_2");
  lines.push_back("GTCAGGAT");
  lines.push_back("+");
  lines.push_back("D7BIG!F8");
  vec_reader reader(lines);

  fastq record;
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "ACGAGTCA");
  REQUIRE(record.qualities() == "!7BF8DGI");
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_2");
  REQUIRE(record.sequence() == "GTCAGGAT");
  REQUIRE(record.qualities() == "D7BIG!F8");
  CHECK(!record.read(reader, FASTQ_ENCODING_33));
}

TEST_CASE("simple_fastq_record__with_extra_header_1", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1 Extra header here");
  lines.push_back("ACGAGTCA");
  lines.push_back("+");
  lines.push_back("!7BF8DGI");
  vec_reader reader(lines);

  fastq record;
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_1 Extra header here");
  REQUIRE(record.sequence() == "ACGAGTCA");
  REQUIRE(record.qualities() == "!7BF8DGI");
  CHECK(!record.read(reader, FASTQ_ENCODING_33));
}

TEST_CASE("simple_fastq_record__with_extra_header_2", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("ACGAGTCA");
  lines.push_back("+Extra header here");
  lines.push_back("!7BF8DGI");
  vec_reader reader(lines);

  fastq record;
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "ACGAGTCA");
  REQUIRE(record.qualities() == "!7BF8DGI");
  CHECK(!record.read(reader, FASTQ_ENCODING_33));
}

TEST_CASE("simple_fastq_record__no_header", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@");
  lines.push_back("ACGAGTCA");
  lines.push_back("+");
  lines.push_back("!7BF8DGI");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}

TEST_CASE("simple_fastq_record__no_sequence", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("");
  lines.push_back("+");
  lines.push_back("!7BF8DGI");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}

TEST_CASE("simple_fastq_record__no_qualities", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("ACGAGTCA");
  lines.push_back("+");
  lines.push_back("");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}

TEST_CASE("simple_fastq_record__no_qualities_or_sequence", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("");
  lines.push_back("+");
  lines.push_back("");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}

TEST_CASE("simple_fastq_record__mismatching_seq_qual_length", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("ACGAGTCA");
  lines.push_back("+");
  lines.push_back("!!!!!!!");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}

TEST_CASE("eof_when_starting_to_read_record", "[fastq::fastq]")
{
  string_vec lines;
  vec_reader reader(lines);

  fastq record;
  CHECK(!record.read(reader));
}

TEST_CASE("eof_after_header", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader), fastq_error);
}

TEST_CASE("eof_after_sequence_1", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record");
  lines.push_back("ACGTA");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader), fastq_error);
}

TEST_CASE("eof_after_sequence_2", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record");
  lines.push_back("ACGTA");
  lines.push_back("");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader), fastq_error);
}

TEST_CASE("eof_after_sep_1", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record");
  lines.push_back("ACGTA");
  lines.push_back("+");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader), fastq_error);
}

TEST_CASE("eof_after_sep_2", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record");
  lines.push_back("ACGTA");
  lines.push_back("+");
  lines.push_back("");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_THROWS_AS(record.read(reader), fastq_error);
}

TEST_CASE("eof_after_qualities_following_previous_read_1", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("ACGTA");
  lines.push_back("+");
  lines.push_back("!!!!!");
  lines.push_back("@record_2");
  lines.push_back("ACGTA");
  lines.push_back("+");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_NOTHROW(record.read(reader));
  REQUIRE_THROWS_AS(record.read(reader), fastq_error);
}

TEST_CASE("eof_after_qualities_following_previous_read_2", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("@record_1");
  lines.push_back("ACGTA");
  lines.push_back("+");
  lines.push_back("!!!!!");
  lines.push_back("@record_2");
  lines.push_back("ACGTA");
  lines.push_back("+");
  lines.push_back("");
  vec_reader reader(lines);

  fastq record;
  REQUIRE_NOTHROW(record.read(reader));
  REQUIRE_THROWS_AS(record.read(reader), fastq_error);
}

TEST_CASE("ignores_trailing_newlines", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("");
  lines.push_back("@record_1");
  lines.push_back("ACGTA");
  lines.push_back("+");
  lines.push_back("!!!!!");
  lines.push_back("");
  lines.push_back("");
  vec_reader reader(lines);

  fastq record;
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "ACGTA");
  REQUIRE(record.qualities() == "!!!!!");
  REQUIRE_NOTHROW(!record.read(reader));
}

TEST_CASE("ignores_newlines_between_records", "[fastq::fastq]")
{
  string_vec lines;
  lines.push_back("");
  lines.push_back("@record_1");
  lines.push_back("ACGAGTCA");
  lines.push_back("+");
  lines.push_back("!7BF8DGI");
  lines.push_back("");
  lines.push_back("");
  lines.push_back("@record_2");
  lines.push_back("GTCAGGAT");
  lines.push_back("+");
  lines.push_back("D7BIG!F8");
  lines.push_back("");
  vec_reader reader(lines);

  fastq record;
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_1");
  REQUIRE(record.sequence() == "ACGAGTCA");
  REQUIRE(record.qualities() == "!7BF8DGI");
  CHECK(record.read(reader, FASTQ_ENCODING_33));
  REQUIRE(record.header() == "record_2");
  REQUIRE(record.sequence() == "GTCAGGAT");
  REQUIRE(record.qualities() == "D7BIG!F8");
  CHECK(!record.read(reader, FASTQ_ENCODING_33));
}

///////////////////////////////////////////////////////////////////////////////
// Writing to stream

TEST_CASE("Writing_to_stream_phred_33", "[fastq::fastq]")
{
  const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
  std::string str;
  record.into_string(str);

  REQUIRE(str == "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n");
}

TEST_CASE("Writing_to_stream_phred_33_explicit", "[fastq::fastq]")
{
  const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
  std::string str;
  record.into_string(str);

  REQUIRE(str == "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n");
}

///////////////////////////////////////////////////////////////////////////////
// Guessing the mate separator

TEST_CASE("guess_mate_separator__empty_lists", "[fastq::guess_mate_separator]")
{
  REQUIRE(fastq::guess_mate_separator(fastq_vec(), fastq_vec()) == '/');
}

TEST_CASE("guess_mate_separator asserts on mismatching input lengths")
{
  const fastq_vec reads_1{ fastq("foo", "") };
  const fastq_vec reads_2;

  REQUIRE_THROWS_AS(fastq::guess_mate_separator(reads_1, reads_2),
                    assert_failed);
  REQUIRE_THROWS_AS(fastq::guess_mate_separator(reads_2, reads_1),
                    assert_failed);
}

TEST_CASE("guess_mate_separator should default to / if there are no separators")
{
  const fastq_vec reads{ fastq("foo", "") };

  REQUIRE(fastq::guess_mate_separator(reads, reads) == '/');
}

TEST_CASE("guess_mate_separator should return valid separator")
{
  const std::string name = "foo";
  const char sep = GENERATE('/', '.', ':');
  const fastq_vec reads_1{ fastq(name + sep + '1', "A") };
  const fastq_vec reads_2{ fastq(name + sep + '2', "A") };

  REQUIRE(fastq::guess_mate_separator(reads_1, reads_2) == sep);
}

TEST_CASE("guess_mate_separator fails on mismatching separators")
{
  const fastq_vec reads_1{ fastq("foo/1", "") };
  const fastq_vec reads_2{ fastq("foo.2", "") };

  REQUIRE(fastq::guess_mate_separator(reads_1, reads_2) == 0);
}

TEST_CASE("guess_mate_separator fails on mismatching mates")
{
  const fastq_vec reads_1{ fastq("foo/1", "A") };
  const fastq_vec reads_2{ fastq("foo/1", "A") };

  REQUIRE(fastq::guess_mate_separator(reads_1, reads_2) == 0);
}

TEST_CASE("guess_mate_separator returns separator for partial information")
{
  const fastq_vec reads_1{ fastq("foo/1", "") };
  const fastq_vec reads_2{ fastq("foo", "") };

  REQUIRE(fastq::guess_mate_separator(reads_1, reads_2) == '/');
}

TEST_CASE("guess_mate_separator fails on inconsistent data")
{
  const fastq_vec reads_1{ fastq("foo/1", ""), fastq("foo.1", "") };
  const fastq_vec reads_2{ fastq("foo/2", ""), fastq("foo.2", "") };

  REQUIRE(fastq::guess_mate_separator(reads_1, reads_2) == 0);
}

TEST_CASE("guess_mate_separator fails on malformed data")
{
  const fastq_vec reads_1{ fastq("foo/1", ""), fastq("foo/1", "") };
  const fastq_vec reads_2{ fastq("foo/2", ""), fastq("foox", "") };

  REQUIRE(fastq::guess_mate_separator(reads_1, reads_2) == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Normalizing pairs

TEST_CASE("normalize_paired_reads__throws_if_order_or_number_is_wrong",
          "[fastq::fastq]")
{
  const fastq ref_mate0 = fastq("Mate/0", "ACGT", "!!#$");
  const fastq ref_mate1 = fastq("Mate/1", "ACGT", "!!#$");
  const fastq ref_mate2 = fastq("Mate/2", "ACGT", "!!#$");
  const fastq ref_mate3 = fastq("Mate/3", "ACGT", "!!#$");
  const fastq ref_matea = fastq("Mate/A", "ACGT", "!!#$");
  const fastq ref_mateb = fastq("Mate/B", "ACGT", "!!#$");

  {
    fastq mate0 = ref_mate0;
    fastq mate1 = ref_mate1;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate0, mate1), fastq_error);
  }

  {
    fastq mate0 = ref_mate0;
    fastq mate1 = ref_mate1;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate0), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    fastq::normalize_paired_reads(mate1, mate2);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate2, mate1), fastq_error);
  }

  {
    fastq mate2 = ref_mate2;
    fastq mate3 = ref_mate3;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate2, mate3), fastq_error);
  }

  {
    fastq mate2 = ref_mate2;
    fastq mate3 = ref_mate3;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate3, mate2), fastq_error);
  }

  {
    fastq matea = ref_matea;
    fastq mateb = ref_mateb;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(matea, mateb), fastq_error);
  }

  {
    fastq matea = ref_matea;
    fastq mateb = ref_mateb;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mateb, matea), fastq_error);
  }
}

TEST_CASE("normalize_paired_reads__allows_other_separators", "[fastq::fastq]")
{
  const fastq ref_mate1 = fastq("Mate:1", "ACGT", "!!#$");
  const fastq ref_mate2 = fastq("Mate:2", "GCTAA", "$!@#$");

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    fastq::normalize_paired_reads(mate1, mate2, ':');
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate2, mate1), fastq_error);
  }
}

TEST_CASE("normalize_paired_reads__mate_separator_is_updated", "[fastq::fastq]")
{
  const fastq ref_mate_1 = fastq("Mate/1", "ACGT", "!!#$");
  const fastq ref_mate_2 = fastq("Mate/2", "GCTAA", "$!@#$");

  fastq mate1 = fastq("Mate:1", "ACGT", "!!#$");
  fastq mate2 = fastq("Mate:2", "GCTAA", "$!@#$");
  fastq::normalize_paired_reads(mate1, mate2, ':');

  REQUIRE(ref_mate_1 == mate1);
  REQUIRE(ref_mate_2 == mate2);
}

TEST_CASE("normalize_paired_reads__throws_if_mate_is_empty", "[fastq::fastq]")
{
  const fastq ref_mate1 = fastq("Mate", "", "");
  const fastq ref_mate2 = fastq("Mate", "ACGT", "!!#$");
  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate2), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate2, mate1), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate1), fastq_error);
  }
}

TEST_CASE("normalize_paired_reads__throws_if_only_mate_1_is_numbered",
          "[fastq::fastq]")
{
  const fastq ref_mate2 = fastq("Mate/1", "GCTAA", "$!@#$");
  const fastq ref_mate1 = fastq("Mate", "ACGT", "!!#$");

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate2), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;

    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate2, mate1), fastq_error);
  }
}

TEST_CASE("normalize_paired_reads__throws_if_only_mate_2_is_numbered",
          "[fastq::fastq]")
{
  const fastq ref_mate1 = fastq("Mate", "GCTAA", "$!@#$");
  const fastq ref_mate2 = fastq("Mate/2", "ACGT", "!!#$");

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate2), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate2, mate1), fastq_error);
  }
}

TEST_CASE("normalize_paired_reads__throws_if_mate_is_misnumbered",
          "[fastq::fastq]")
{
  fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
  fastq mate2 = fastq("Mate/3", "ACGT", "!!#$");
  REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate2), fastq_error);
}

TEST_CASE("normalize_paired_reads__throws_if_same_mate_numbers",
          "[fastq::fastq]")
{
  fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
  fastq mate2 = fastq("Mate/1", "ACGT", "!!#$");
  REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate2), fastq_error);
}

TEST_CASE("normalize_paired_reads__throws_if_name_differs", "[fastq::fastq]")
{
  fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
  fastq mate2 = fastq("WrongName/2", "ACGT", "!!#$");
  REQUIRE_THROWS_AS(fastq::normalize_paired_reads(mate1, mate2), fastq_error);
}

TEST_CASE("normalize_paired_reads doesn't modify reads without mate numbers")
{
  const auto name = GENERATE("Name", "Name and meta");
  fastq mate1 = fastq(name, "GCTAA", "$!@#$");
  fastq mate2 = fastq(name, "ACGTA", "@$!$#");

  fastq::normalize_paired_reads(mate1, mate2);

  REQUIRE(mate1.header() == name);
  REQUIRE(mate2.header() == name);
}

} // namespace adapterremoval
