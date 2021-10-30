/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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

#include "debug.hpp"
#include "fastq.hpp"
#include "linereader.hpp"
#include "testing.hpp"

class vec_reader : public line_reader_base
{
public:
  vec_reader(const string_vec& lines)
    : m_lines(lines)
    , m_it(m_lines.begin())
  {}

  bool getline(std::string& dst)
  {
    if (m_it == m_lines.end()) {
      return false;
    }

    dst = *m_it++;
    return true;
  }

private:
  string_vec m_lines;
  string_vec::const_iterator m_it;
};

///////////////////////////////////////////////////////////////////////////////
// Helper functions

TEST_CASE("ACGT_TO_IDX", "[fastq::*]")
{
  // The exact encoding is unimportant, but it must be unique and 2 bit
  REQUIRE(ACGT_TO_IDX('A') <= 3);
  REQUIRE(ACGT_TO_IDX('C') <= 3);
  REQUIRE(ACGT_TO_IDX('G') <= 3);
  REQUIRE(ACGT_TO_IDX('T') <= 3);

  REQUIRE(ACGT_TO_IDX('A') != ACGT_TO_IDX('C'));
  REQUIRE(ACGT_TO_IDX('A') != ACGT_TO_IDX('G'));
  REQUIRE(ACGT_TO_IDX('A') != ACGT_TO_IDX('T'));
  REQUIRE(ACGT_TO_IDX('C') != ACGT_TO_IDX('G'));
  REQUIRE(ACGT_TO_IDX('C') != ACGT_TO_IDX('T'));
  REQUIRE(ACGT_TO_IDX('G') != ACGT_TO_IDX('T'));
}

TEST_CASE("IDX_TO_ACGT", "[fastq::*]")
{
  REQUIRE(IDX_TO_ACGT(ACGT_TO_IDX('A')) == 'A');
  REQUIRE(IDX_TO_ACGT(ACGT_TO_IDX('C')) == 'C');
  REQUIRE(IDX_TO_ACGT(ACGT_TO_IDX('G')) == 'G');
  REQUIRE(IDX_TO_ACGT(ACGT_TO_IDX('T')) == 'T');
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

TEST_CASE("constructor_simple_record_dots_to_n", "[fastq::fastq]")
{
  const fastq record("record_1", "AC.AG.C.", "!7BF8DGI");
  REQUIRE(record.sequence() == "ACNAGNCN");
}

TEST_CASE("constructor_score_boundries_phred_33", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", "!!\"", FASTQ_ENCODING_33));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", " !\"", FASTQ_ENCODING_33),
                    fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "IJJ", FASTQ_ENCODING_33));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "IJK", FASTQ_ENCODING_33), fastq_error);
}

TEST_CASE("constructor_score_boundries_phred_64", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", "@@A", FASTQ_ENCODING_64));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "?@A", FASTQ_ENCODING_64), fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "ghi", FASTQ_ENCODING_64));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "ghj", FASTQ_ENCODING_64), fastq_error);
}

TEST_CASE("constructor_score_boundries_solexa", "[fastq::fastq]")
{
  REQUIRE_NOTHROW(fastq("Rec", "CAT", ";;<", FASTQ_ENCODING_SOLEXA));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", ":;<", FASTQ_ENCODING_SOLEXA),
                    fastq_error);

  REQUIRE_NOTHROW(fastq("Rec", "CAT", "ghi", FASTQ_ENCODING_SOLEXA));
  REQUIRE_THROWS_AS(fastq("Rec", "CAT", "ghj", FASTQ_ENCODING_SOLEXA),
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

TEST_CASE("add_prefix_to_header", "[fastq::fastq]")
{
  const fastq expected("not_my_header", "ACGTA", "12345");
  fastq record("my_header", "ACGTA", "12345");
  record.add_prefix_to_header("not_");
  REQUIRE(record == expected);
}

TEST_CASE("add_prefix_to_header__empty_prefix", "[fastq::fastq]")
{
  const fastq expected("my_header", "ACGTA", "12345");
  fastq record = expected;
  record.add_prefix_to_header("");
  REQUIRE(record == expected);
}

TEST_CASE("add_prefix_to_header__header", "[fastq::fastq]")
{
  const fastq expected("new_header", "ACGTA", "12345");
  fastq record("", "ACGTA", "12345");
  record.add_prefix_to_header("new_header");
  REQUIRE(record == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Adding postfixes to the header

TEST_CASE("add_postfix_to_header", "[fastq::fastq]")
{
  const fastq expected("my_header new postfix", "ACGTA", "12345");
  fastq record("my_header", "ACGTA", "12345");
  record.add_postfix_to_header(" new postfix");
  REQUIRE(record == expected);
}

TEST_CASE("add_postfix_to_header__empty_prefix", "[fastq::fastq]")
{
  const fastq expected("my_header", "ACGTA", "12345");
  fastq record = expected;
  record.add_postfix_to_header("");
  REQUIRE(record == expected);
}

TEST_CASE("add_postfix_to_header__header", "[fastq::fastq]")
{
  const fastq expected("new_header", "ACGTA", "12345");
  fastq record("", "ACGTA", "12345");
  record.add_postfix_to_header("new_header");
  REQUIRE(record == expected);
}

///////////////////////////////////////////////////////////////////////////////
// Discarding read, setting seq to N and qual to '!'

TEST_CASE("discard_read", "[fastq::fastq]")
{
  const fastq expected("my_header", "N", "!");
  fastq record("my_header", "ACGTA", "12345");
  record.discard();
  REQUIRE(record == expected);
}

TEST_CASE("discard_discarded_read", "[fastq::fastq]")
{
  const fastq expected("my_header", "N", "!");
  fastq record("my_header", "N", "!");
  record.discard();
  REQUIRE(record == expected);
}

TEST_CASE("discard_empty_read", "[fastq::fastq]")
{
  const fastq expected("my_header", "N", "!");
  fastq record("my_header", "", "");
  record.discard();
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
  lines.push_back("@");
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
  lines.push_back("@");
  lines.push_back("");
  lines.push_back("+");
  lines.push_back("");
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

TEST_CASE("Writing_to_stream_phred_64_explicit", "[fastq::fastq]")
{
  const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
  std::string str;
  record.into_string(str, FASTQ_ENCODING_64);

  REQUIRE(str == "@record_1\nACGTACGATA\n+\n@CBCIUWbfi\n");
}

///////////////////////////////////////////////////////////////////////////////
// Validating pairs

TEST_CASE("validate_paired_reads__throws_if_order_or_number_is_wrong",
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
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate0, mate1), fastq_error);
  }

  {
    fastq mate0 = ref_mate0;
    fastq mate1 = ref_mate1;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate0), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    fastq::validate_paired_reads(mate1, mate2);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate2, mate1), fastq_error);
  }

  {
    fastq mate2 = ref_mate2;
    fastq mate3 = ref_mate3;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate2, mate3), fastq_error);
  }

  {
    fastq mate2 = ref_mate2;
    fastq mate3 = ref_mate3;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate3, mate2), fastq_error);
  }

  {
    fastq matea = ref_matea;
    fastq mateb = ref_mateb;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(matea, mateb), fastq_error);
  }

  {
    fastq matea = ref_matea;
    fastq mateb = ref_mateb;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mateb, matea), fastq_error);
  }
}

TEST_CASE("validate_paired_reads__allows_other_separators", "[fastq::fastq]")
{
  const fastq ref_mate1 = fastq("Mate:1", "ACGT", "!!#$");
  const fastq ref_mate2 = fastq("Mate:2", "GCTAA", "$!@#$");

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    fastq::validate_paired_reads(mate1, mate2, ':');
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate2, mate1), fastq_error);
  }
}

TEST_CASE("validate_paired_reads__mate_separator_is_updated", "[fastq::fastq]")
{
  const fastq ref_mate_1 = fastq("Mate/1", "ACGT", "!!#$");
  const fastq ref_mate_2 = fastq("Mate/2", "GCTAA", "$!@#$");

  fastq mate1 = fastq("Mate:1", "ACGT", "!!#$");
  fastq mate2 = fastq("Mate:2", "GCTAA", "$!@#$");
  fastq::validate_paired_reads(mate1, mate2, ':');

  REQUIRE(ref_mate_1 == mate1);
  REQUIRE(ref_mate_2 == mate2);
}

TEST_CASE("validate_paired_reads__throws_if_mate_is_empty", "[fastq::fastq]")
{
  const fastq ref_mate1 = fastq("Mate", "", "");
  const fastq ref_mate2 = fastq("Mate", "ACGT", "!!#$");
  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate2), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate2, mate1), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate1), fastq_error);
  }
}

TEST_CASE("validate_paired_reads__throws_if_only_mate_1_is_numbered",
          "[fastq::fastq]")
{
  const fastq ref_mate2 = fastq("Mate/1", "GCTAA", "$!@#$");
  const fastq ref_mate1 = fastq("Mate", "ACGT", "!!#$");

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate2), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;

    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate2, mate1), fastq_error);
  }
}

TEST_CASE("validate_paired_reads__throws_if_only_mate_2_is_numbered",
          "[fastq::fastq]")
{
  const fastq ref_mate1 = fastq("Mate", "GCTAA", "$!@#$");
  const fastq ref_mate2 = fastq("Mate/2", "ACGT", "!!#$");

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate2), fastq_error);
  }

  {
    fastq mate1 = ref_mate1;
    fastq mate2 = ref_mate2;
    REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate2, mate1), fastq_error);
  }
}

TEST_CASE("validate_paired_reads__throws_if_mate_is_misnumbered",
          "[fastq::fastq]")
{
  fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
  fastq mate2 = fastq("Mate/3", "ACGT", "!!#$");
  REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}

TEST_CASE("validate_paired_reads__throws_if_same_mate_numbers",
          "[fastq::fastq]")
{
  fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
  fastq mate2 = fastq("Mate/1", "ACGT", "!!#$");
  REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}

TEST_CASE("validate_paired_reads__throws_if_name_differs", "[fastq::fastq]")
{
  fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
  fastq mate2 = fastq("WrongName/2", "ACGT", "!!#$");
  REQUIRE_THROWS_AS(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}
