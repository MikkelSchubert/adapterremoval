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
#include <gtest/gtest.h>

#include "testing.hpp"
#include "debug.hpp"
#include "fastq.hpp"
#include "linereader.hpp"

namespace ar
{

class vec_reader : public line_reader_base
{
public:
    vec_reader(const string_vec& lines)
        : m_lines(lines)
        , m_it(m_lines.begin())
    {
    }

    bool getline(std::string& dst) {
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
// Default constructor

TEST(fastq, default_constructor)
{
    const fastq record;
    ASSERT_EQ("", record.header());
    ASSERT_EQ("", record.sequence());
    ASSERT_EQ("", record.qualities());
}


///////////////////////////////////////////////////////////////////////////////
// Primary constructor

TEST(fastq, constructor_empty_fields)
{
    const fastq record("", "", "");
    ASSERT_EQ("", record.header());
    ASSERT_EQ("", record.sequence());
    ASSERT_EQ("", record.qualities());
}


TEST(fastq, constructor_simple_record_phred_33_encoded)
{
    const fastq record("record_1", "ACGAGTCA", "!7BF8DGI");
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
}


TEST(fastq, constructor_simple_record_phred_64_encoded)
{
    const fastq record("record_2", "ACGAGTCA", "@VaeWcfh", FASTQ_ENCODING_64);
    ASSERT_EQ("record_2", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
}


TEST(fastq, constructor_simple_record_phred_solexa_encoded)
{
    const fastq record("record_3", "AAACGAGTCA", ";h>S\\TCDUJ", FASTQ_ENCODING_SOLEXA);
    ASSERT_EQ("record_3", record.header());
    ASSERT_EQ("AAACGAGTCA", record.sequence());
    ASSERT_EQ("\"I#4=5&&6+", record.qualities());
}


TEST(fastq, constructor_simple_record_lowercase_to_uppercase)
{
    const fastq record("record_1", "AnGaGtcA", "!7BF8DGI");
    ASSERT_EQ("ANGAGTCA", record.sequence());
}


TEST(fastq, constructor_simple_record_dots_to_n)
{
    const fastq record("record_1", "AC.AG.C.", "!7BF8DGI");
    ASSERT_EQ("ACNAGNCN", record.sequence());
}


TEST(fastq, constructor_score_boundries_phred_33)
{
    ASSERT_NO_THROW(fastq("Rec", "CAT", "!!\"", FASTQ_ENCODING_33));
    ASSERT_THROW(fastq("Rec", "CAT", " !\"", FASTQ_ENCODING_33), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "IJJ", FASTQ_ENCODING_33));
    ASSERT_THROW(fastq("Rec", "CAT", "IJK", FASTQ_ENCODING_33), fastq_error);
}


TEST(fastq, constructor_score_boundries_phred_64)
{
    ASSERT_NO_THROW(fastq("Rec", "CAT", "@@A", FASTQ_ENCODING_64));
    ASSERT_THROW(fastq("Rec", "CAT", "?@A", FASTQ_ENCODING_64), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "ghi", FASTQ_ENCODING_64));
    ASSERT_THROW(fastq("Rec", "CAT", "ghj", FASTQ_ENCODING_64), fastq_error);
}


TEST(fastq, constructor_score_boundries_solexa)
{
    ASSERT_NO_THROW(fastq("Rec", "CAT", ";;<", FASTQ_ENCODING_SOLEXA));
    ASSERT_THROW(fastq("Rec", "CAT", ":;<", FASTQ_ENCODING_SOLEXA), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "ghi", FASTQ_ENCODING_SOLEXA));
    ASSERT_THROW(fastq("Rec", "CAT", "ghj", FASTQ_ENCODING_SOLEXA), fastq_error);
}

TEST(fastq, constructor_score_boundries_ignored)
{
    ASSERT_NO_THROW(fastq("Rec", "CAT", "!!\"", FASTQ_ENCODING_SAM));
    ASSERT_THROW(fastq("Rec", "CAT", " !\"", FASTQ_ENCODING_SAM), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "gh~", FASTQ_ENCODING_SAM));
    ASSERT_THROW(fastq("Rec", "CAT", "gh\x7f", FASTQ_ENCODING_SAM), fastq_error);
}


TEST(fastq, constructor_field_lengths)
{
    ASSERT_NO_THROW(fastq("Name", "CAT", "IJJ"));
    // A non-empty sequence must be specified
    ASSERT_THROW(fastq("Name", "", "IJJ"), fastq_error);
    // A non-empty quality string must be specified
    ASSERT_THROW(fastq("Name", "CAT", ""), fastq_error);
    // And the length of each must be the same
    ASSERT_THROW(fastq("Name", "CA", "IJJ"), fastq_error);
    ASSERT_THROW(fastq("Name", "CAT", "IJ"), fastq_error);
}


TEST(fastq, constructor_invalid_nucleotides)
{
    ASSERT_NO_THROW(fastq("Name", "CATT", "IJJI"));
    // Non-alpha characters are not allowed
    ASSERT_THROW(fastq("Name", "CAT!", "IJJI"), fastq_error);
    // Numeric charecters are not allowed
    ASSERT_THROW(fastq("Name", "CAT7", "IJJI"), fastq_error);
    // But neither are non acgtn/ACGTN allowed
    ASSERT_THROW(fastq("Name", "CATS", "IJJI"), fastq_error);
    ASSERT_THROW(fastq("Name", "CATs", "IJJI"), fastq_error);
}


///////////////////////////////////////////////////////////////////////////////
// Constructor without qualities

TEST(fastq, constructor_no_qualities)
{
    const fastq record("record_1", "ACGT");
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGT", record.sequence());
    ASSERT_EQ("!!!!", record.qualities());
}


TEST(fastq, constructor_no_qualities_no_sequence)
{
    const fastq record("record_1", "");
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("", record.sequence());
    ASSERT_EQ("", record.qualities());
}


///////////////////////////////////////////////////////////////////////////////
// misc properties

TEST(fastq, name)
{
    ASSERT_EQ("name", fastq("name", "", "").name());
    ASSERT_EQ("name", fastq("name meta", "", "").name());
    ASSERT_EQ("name", fastq("name meta more", "", "").name());
}


TEST(fastq, length)
{
    ASSERT_EQ(0, fastq("record_1", "", "").length());
    ASSERT_EQ(1, fastq("record_1", "A", "G").length());
    ASSERT_EQ(2, fastq("record_1", "AC", "!B").length());
    ASSERT_EQ(3, fastq("record_1", "ACG", "!7B").length());
}


TEST(fastq, count_ns)
{
    ASSERT_EQ(0, fastq("Rec", "ACGTA", "IJIJI").count_ns());
    ASSERT_EQ(1, fastq("Rec", "ANGTA", "IJIJI").count_ns());
    ASSERT_EQ(2, fastq("Rec", "ANGTN", "IJIJI").count_ns());
    ASSERT_EQ(3, fastq("Rec", "ANGNN", "IJIJI").count_ns());
    ASSERT_EQ(4, fastq("Rec", "NNGNN", "IJIJI").count_ns());
    ASSERT_EQ(5, fastq("Rec", "NNNNN", "IJIJI").count_ns());
}


///////////////////////////////////////////////////////////////////////////////
// trim_trailing_bases

TEST(fastq, trim_trailing_bases__empty_record)
{
    fastq record("Empty", "", "");
    const fastq::ntrimmed expected(0, 0);
    ASSERT_EQ(expected, record.trim_trailing_bases(true, 10));
    ASSERT_EQ(fastq("Empty", "", ""), record);
}


TEST(fastq, trim_trailing_bases__trim_nothing)
{
    const fastq reference("Rec", "NNNNN", "!!!!!");
    const fastq::ntrimmed expected(0, 0);
    fastq record = reference;

    // Trim neither Ns nor low Phred score bases
    ASSERT_EQ(expected, record.trim_trailing_bases(false, -1));
    ASSERT_EQ(reference, record);
}


TEST(fastq, trim_trailing_bases__trim_ns)
{
    fastq record("Rec", "NNANT", "23456");
    const fastq expected_record("Rec", "ANT", "456");
    const fastq::ntrimmed expected_ntrim(2, 0);

    ASSERT_EQ(expected_ntrim, record.trim_trailing_bases(true, -1));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_trailing_bases__trim_trailing_bases)
{
    const fastq expected_record("Rec", "TN", "%$");
    const fastq::ntrimmed expected_ntrim(0, 3);
    fastq record("Rec", "TNANT", "%$#!\"");

    ASSERT_EQ(expected_ntrim, record.trim_trailing_bases(false, 2));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_trailing_bases__trim_mixed)
{
    const fastq expected_record("Rec", "TAG", "$12");
    const fastq::ntrimmed expected_ntrim(3, 2);
    fastq record("Rec", "NTNTAGNT", "1!#$12#\"");

    ASSERT_EQ(expected_ntrim, record.trim_trailing_bases(true, 2));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_trailing_bases__trim_mixed__no_low_quality_bases)
{
    const fastq expected_record("Rec", "ACTTAG", "12I$12");
    const fastq::ntrimmed expected_ntrim(0, 0);
    fastq record = expected_record;

    ASSERT_EQ(expected_ntrim, record.trim_trailing_bases(true, 2));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_trailing_bases__trim_everything)
{
    fastq record("Rec", "TAG", "!!!");
    const fastq expected_record = fastq("Rec", "", "");
    const fastq::ntrimmed expected_ntrim(0, 3);
    ASSERT_EQ(expected_ntrim, record.trim_trailing_bases(true, 2));
    ASSERT_EQ(expected_record, record);
}


///////////////////////////////////////////////////////////////////////////////
// trim_windowed_bases


#define PARAMETERIZED_TEST(name, values) \
    class name : public ::testing::TestWithParam<double> {}; \
    INSTANTIATE_TEST_CASE_P(fastq_windowed_trimming, name, values); \
    TEST_P(name, test)


// Test for invalid parameters
PARAMETERIZED_TEST(invalid_parameters, ::testing::Values(-1.0, std::numeric_limits<double>::quiet_NaN()))
{
    fastq record("Rec", "TAGTGACAT", "111111111");
    ASSERT_THROW(record.trim_windowed_bases(false, -1, GetParam()), assert_failed);
}


// Test for trimming empty reads
PARAMETERIZED_TEST(empty_reads, ::testing::Values(1, 0.1, 3))
{
    fastq record("Empty", "", "");
    const fastq::ntrimmed expected(0, 0);
    ASSERT_EQ(expected, record.trim_windowed_bases(true, 10, GetParam()));
    ASSERT_EQ(fastq("Empty", "", ""), record);
}


// Test for when entire read is trimmed
PARAMETERIZED_TEST(trim_everything, ::testing::Values(1, 0.2, 4, 10))
{
    fastq record("Rec", "TAGTGACAT", "111111111");
    const fastq expected_record = fastq("Rec", "", "");
    const fastq::ntrimmed expected_ntrim(9, 0);
    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '2' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}


// Test for when nothing is trimmed
PARAMETERIZED_TEST(trim_nothing, ::testing::Values(0, 1, 0.2, 3, 4, 5, 6, 7, 8, 9, 10, 11))
{
    fastq record("Rec", "TAGTGACAT", "111111111");
    const fastq expected_record = record;
    const fastq::ntrimmed expected_ntrim(0, 0);
    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, -1, GetParam()));
    ASSERT_EQ(expected_record, record);
}


// Test for trimming of Ns
PARAMETERIZED_TEST(trim_ns_1, ::testing::Values(1, 0.2))
{
    fastq record("Rec", "NNATNT", "234567");
    const fastq expected_record("Rec", "ATNT", "4567");
    const fastq::ntrimmed expected_ntrim(2, 0);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(true, -1, GetParam()));
    ASSERT_EQ(expected_record, record);
}

// - trimming of Ns - the final window contains Ns and is therefore truncated
PARAMETERIZED_TEST(trim_ns_2, ::testing::Values(2, 3, 4))
{
    fastq record("Rec", "NNATNT", "234567");
    const fastq expected_record("Rec", "AT", "45");
    const fastq::ntrimmed expected_ntrim(2, 2);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(true, -1, GetParam()));
    ASSERT_EQ(expected_record, record);
}

// - trimming of Ns - No valid window is found
TEST(fastq, trim_windowed_bases__trim_ns_5)
{
    fastq record("Rec", "NNATNT", "234567");
    const fastq expected_record("Rec", "", "");
    const fastq::ntrimmed expected_ntrim(6, 0);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(true, -1, 5));
    ASSERT_EQ(expected_record, record);
}


// Trimming of 5' only
PARAMETERIZED_TEST(trim_5p_1bp, ::testing::Values(1, 0.1))
{
    fastq record("Rec", "TAACGATCCG", "0123456789");
    const fastq expected_record("Rec", "CGATCCG", "3456789");
    const fastq::ntrimmed expected_ntrim(3, 0);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '2' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}


// Trimming of 5' only
PARAMETERIZED_TEST(trim_5p_2bp, ::testing::Values(2, 0.2))
{
    fastq record("Rec", "TAACGATCCG", "0123456789");
    const fastq expected_record("Rec", "CGATCCG", "3456789");
    const fastq::ntrimmed expected_ntrim(3, 0);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '2' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}



// Trimming of 5' only - lowquality is inclusive
PARAMETERIZED_TEST(trim_5p_inclusive_low_quality, ::testing::Values(2, 3))
{
    fastq record("Rec", "TAACGATCCG", "0123126789");
    const fastq expected_record("Rec", "TCCG", "6789");
    const fastq::ntrimmed expected_ntrim(6, 0);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '2' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}


// Trimming of 3' only - lowquality is inclusive
PARAMETERIZED_TEST(trim_3p_inclusive_low_quality, ::testing::Values(3))
{
    fastq record("Rec", "TAACGATCCG", "9876312333");
    const fastq expected_record("Rec", "TAACG", "98763");
    const fastq::ntrimmed expected_ntrim(0, 5);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '2' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}


PARAMETERIZED_TEST(tiny_and_huge_window_sizes_1, ::testing::Values(0, 0.01, 20))
{
    fastq record("Rec", "TAACGATC", "23456789");
    const fastq expected_record("Rec", "", "");
    const fastq::ntrimmed expected_ntrim(8, 0);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '2' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}


PARAMETERIZED_TEST(tiny_and_huge_window_sizes_2, ::testing::Values(0, 0.01, 20))
{
    fastq record("Rec", "TAACGATC", "23456780");
    const fastq expected_record("Rec", "TAACGAT", "2345678");
    const fastq::ntrimmed expected_ntrim(0, 1);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '1' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}


PARAMETERIZED_TEST(last_trailing_window, ::testing::Values(1, 2, 3, 4, 5, 6, 7, 8, 9))
{
    fastq record("Rec", "TAACGATCC", "234567811");
    const fastq expected_record("Rec", "TAACGAT", "2345678");
    const fastq::ntrimmed expected_ntrim(0, 2);

    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(false, '1' - '!',
                                                         GetParam()));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_windowed_bases__trim_window)
{
    // Should trim starting at the window of low quality bases in the middle
    // even with high qual bases at the end.
    fastq record("Rec", "NNAAAAAAAAATNNNNNNNA", "##EEEEEEEEEE#######E");
    const fastq expected_record = fastq("Rec", "AAAAAAAAAT", "EEEEEEEEEE");
    const fastq::ntrimmed expected_ntrim(2, 8);
    ASSERT_EQ(expected_ntrim, record.trim_windowed_bases(true, 10, 5));
    ASSERT_EQ(expected_record, record);
}


///////////////////////////////////////////////////////////////////////////////
// Truncate

TEST(fastq, truncate_empty)
{
    fastq record("Empty", "", "");
    record.truncate(0, 10);
    ASSERT_EQ(fastq("Empty", "", ""), record);
}


TEST(fastq, truncate_zero_bases)
{
    const fastq expected_record("Rec", "ACTTAG", "12I$12");
    fastq current_record = expected_record;
    current_record.truncate();
    ASSERT_EQ(expected_record, current_record);
}


TEST(fastq, truncate_all_bases)
{
    const fastq expected_record("Rec", "", "");
    fastq current_record("Rec", "ACTTAG", "12I$12");
    current_record.truncate(1, 0);
    ASSERT_EQ(expected_record, current_record);
}


TEST(fastq, truncate_5p)
{
    const fastq expected_record("Rec", "TTAG", "I$12");
    fastq current_record("Rec", "ACTTAG", "12I$12");
    current_record.truncate(2);
    ASSERT_EQ(expected_record, current_record);
}


TEST(fastq, truncate_3p)
{
    const fastq expected_record("Rec", "ACT", "12I");
    fastq current_record("Rec", "ACTTAG", "12I$12");
    current_record.truncate(0, 3);
    ASSERT_EQ(expected_record, current_record);
}


TEST(fastq, truncate_middle)
{
    const fastq expected_record("Rec", "TTA", "I$1");
    fastq current_record("Rec", "ACTTAG", "12I$12");
    current_record.truncate(2, 3);
    ASSERT_EQ(expected_record, current_record);
}


TEST(fastq, truncate_len_higher_than_n_bases)
{
    const fastq expected_record("Rec", "TTAG", "I$12");
    fastq current_record("Rec", "ACTTAG", "12I$12");
    current_record.truncate(2, 1024);
    ASSERT_EQ(expected_record, current_record);
}


TEST(fastq, truncate_pos_after_last_base)
{
    // Same behavior as string::substr
    fastq current_record("Rec", "ACTTAG", "12I$12");
    ASSERT_NO_THROW(current_record.truncate(6));
    ASSERT_THROW(current_record.truncate(7), std::out_of_range);
}


///////////////////////////////////////////////////////////////////////////////
// Reverse complement

TEST(fastq, reverse_complement__empty)
{
    const fastq expected = fastq("Empty", "", "");
    fastq result = fastq("Empty", "", "");
    result.reverse_complement();
    ASSERT_EQ(expected, result);
}


TEST(fastq, reverse_complement)
{
    const fastq expected = fastq("Rec", "TACAGANGTN", "0123456789");
    fastq result = fastq("Rec", "NACNTCTGTA", "9876543210");
    result.reverse_complement();
    ASSERT_EQ(expected, result);
}


///////////////////////////////////////////////////////////////////////////////
// Adding prefixes to the header

TEST(fastq, add_prefix_to_header)
{
    const fastq expected("not_my_header", "ACGTA", "12345");
    fastq record("my_header", "ACGTA", "12345");
    record.add_prefix_to_header("not_");
    ASSERT_EQ(expected, record);
}


TEST(fastq, add_prefix_to_header__empty_prefix)
{
    const fastq expected("my_header", "ACGTA", "12345");
    fastq record = expected;
    record.add_prefix_to_header("");
    ASSERT_EQ(expected, record);
}


TEST(fastq, add_prefix_to_header__header)
{
    const fastq expected("new_header", "ACGTA", "12345");
    fastq record("", "ACGTA", "12345");
    record.add_prefix_to_header("new_header");
    ASSERT_EQ(expected, record);
}


///////////////////////////////////////////////////////////////////////////////
// Adding postfixes to the header

TEST(fastq, add_postfix_to_header)
{
    const fastq expected("my_header new postfix", "ACGTA", "12345");
    fastq record("my_header", "ACGTA", "12345");
    record.add_postfix_to_header(" new postfix");
    ASSERT_EQ(expected, record);
}


TEST(fastq, add_postfix_to_header__empty_prefix)
{
    const fastq expected("my_header", "ACGTA", "12345");
    fastq record = expected;
    record.add_postfix_to_header("");
    ASSERT_EQ(expected, record);
}


TEST(fastq, add_postfix_to_header__header)
{
    const fastq expected("new_header", "ACGTA", "12345");
    fastq record("", "ACGTA", "12345");
    record.add_postfix_to_header("new_header");
    ASSERT_EQ(expected, record);
}


///////////////////////////////////////////////////////////////////////////////
// Discarding read, setting seq to N and qual to '!'

TEST(fastq, discard_read)
{
    const fastq expected("my_header", "N", "!");
    fastq record("my_header", "ACGTA", "12345");
    record.discard();
    ASSERT_EQ(expected, record);
}


TEST(fastq, discard_discarded_read)
{
    const fastq expected("my_header", "N", "!");
    fastq record("my_header", "N", "!");
    record.discard();
    ASSERT_EQ(expected, record);
}


TEST(fastq, discard_empty_read)
{
    const fastq expected("my_header", "N", "!");
    fastq record("my_header", "", "");
    record.discard();
    ASSERT_EQ(expected, record);
}


///////////////////////////////////////////////////////////////////////////////
// Reading from stream

TEST(fastq, simple_fastq_record_1)
{
    string_vec lines;
    lines.push_back("@record_1");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    vec_reader reader(lines);

    fastq record;
    ASSERT_TRUE(record.read(reader, FASTQ_ENCODING_33));
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
}

TEST(fastq, simple_fastq_record_2)
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
    ASSERT_TRUE(record.read(reader, FASTQ_ENCODING_33));
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
    ASSERT_TRUE(record.read(reader, FASTQ_ENCODING_33));
    ASSERT_EQ("record_2", record.header());
    ASSERT_EQ("GTCAGGAT", record.sequence());
    ASSERT_EQ("D7BIG!F8", record.qualities());
    ASSERT_FALSE(record.read(reader, FASTQ_ENCODING_33));
}


TEST(fastq, simple_fastq_record__with_extra_header_1)
{
    string_vec lines;
    lines.push_back("@record_1 Extra header here");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    vec_reader reader(lines);

    fastq record;
    ASSERT_TRUE(record.read(reader, FASTQ_ENCODING_33));
    ASSERT_EQ("record_1 Extra header here", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
    ASSERT_FALSE(record.read(reader, FASTQ_ENCODING_33));
}


TEST(fastq, simple_fastq_record__with_extra_header_2)
{
    string_vec lines;
    lines.push_back("@record_1");
    lines.push_back("ACGAGTCA");
    lines.push_back("+Extra header here");
    lines.push_back("!7BF8DGI");
    vec_reader reader(lines);

    fastq record;
    ASSERT_TRUE(record.read(reader, FASTQ_ENCODING_33));
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
    ASSERT_FALSE(record.read(reader, FASTQ_ENCODING_33));
}


TEST(fastq, simple_fastq_record__no_header)
{
    string_vec lines;
    lines.push_back("@");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}


TEST(fastq, simple_fastq_record__no_sequence)
{
    string_vec lines;
    lines.push_back("@record_1");
    lines.push_back("");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}


TEST(fastq, simple_fastq_record__no_qualities)
{
    string_vec lines;
    lines.push_back("@");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}


TEST(fastq, simple_fastq_record__no_qualities_or_sequence)
{
    string_vec lines;
    lines.push_back("@");
    lines.push_back("");
    lines.push_back("+");
    lines.push_back("");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader, FASTQ_ENCODING_33), fastq_error);
}


TEST(fastq, eof_when_starting_to_read_record)
{
    string_vec lines;
    vec_reader reader(lines);

    fastq record;
    ASSERT_FALSE(record.read(reader));
}


TEST(fastq, eof_after_header)
{
    string_vec lines;
    lines.push_back("@record");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader), fastq_error);
}


TEST(fastq, eof_after_sequence_1)
{
    string_vec lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader), fastq_error);
}


TEST(fastq, eof_after_sequence_2)
{
    string_vec lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    lines.push_back("");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader), fastq_error);
}


TEST(fastq, eof_after_sep_1)
{
    string_vec lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    lines.push_back("+");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader), fastq_error);
}


TEST(fastq, eof_after_sep_2)
{
    string_vec lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    lines.push_back("+");
    lines.push_back("");
    vec_reader reader(lines);

    fastq record;
    ASSERT_THROW(record.read(reader), fastq_error);
}


TEST(fastq, eof_after_qualities_following_previous_read_1)
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
    ASSERT_NO_THROW(record.read(reader));
    ASSERT_THROW(record.read(reader), fastq_error);
}


TEST(fastq, eof_after_qualities_following_previous_read_2)
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
    ASSERT_NO_THROW(record.read(reader));
    ASSERT_THROW(record.read(reader), fastq_error);
}

///////////////////////////////////////////////////////////////////////////////
// Writing to stream

TEST(fastq, Writing_to_stream_phred_33)
{
    const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
    ASSERT_EQ("@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n", record.to_str());
}

TEST(fastq, Writing_to_stream_phred_33_explicit)
{
    const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
    ASSERT_EQ("@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n", record.to_str());
}

TEST(fastq, Writing_to_stream_phred_64_explicit)
{
    const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
    ASSERT_EQ("@record_1\nACGTACGATA\n+\n@CBCIUWbfi\n", record.to_str(FASTQ_ENCODING_64));
}


///////////////////////////////////////////////////////////////////////////////
// Validating pairs

TEST(fastq, validate_paired_reads__throws_if_order_or_number_is_wrong)
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
        ASSERT_THROW(fastq::validate_paired_reads(mate0, mate1), fastq_error);
    }

    {
        fastq mate0 = ref_mate0;
        fastq mate1 = ref_mate1;
        ASSERT_THROW(fastq::validate_paired_reads(mate1, mate0), fastq_error);
    }

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        fastq::validate_paired_reads(mate1, mate2);
    }

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
    }

    {
        fastq mate2 = ref_mate2;
        fastq mate3 = ref_mate3;
        ASSERT_THROW(fastq::validate_paired_reads(mate2, mate3), fastq_error);
    }

    {
        fastq mate2 = ref_mate2;
        fastq mate3 = ref_mate3;
        ASSERT_THROW(fastq::validate_paired_reads(mate3, mate2), fastq_error);
    }

    {
        fastq matea = ref_matea;
        fastq mateb = ref_mateb;
        ASSERT_THROW(fastq::validate_paired_reads(matea, mateb), fastq_error);
    }

    {
        fastq matea = ref_matea;
        fastq mateb = ref_mateb;
        ASSERT_THROW(fastq::validate_paired_reads(mateb, matea), fastq_error);
    }
}


TEST(fastq, validate_paired_reads__allows_other_separators)
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
        ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
    }
}


TEST(fastq, validate_paired_reads__mate_separator_is_updated)
{
    const fastq ref_mate_1 = fastq("Mate/1", "ACGT", "!!#$");
    const fastq ref_mate_2 = fastq("Mate/2", "GCTAA", "$!@#$");

    fastq mate1 = fastq("Mate:1", "ACGT", "!!#$");
    fastq mate2 = fastq("Mate:2", "GCTAA", "$!@#$");
    fastq::validate_paired_reads(mate1, mate2, ':');

    ASSERT_EQ(mate1, ref_mate_1);
    ASSERT_EQ(mate2, ref_mate_2);
}


TEST(fastq, validate_paired_reads__throws_if_mate_is_empty)
{
    const fastq ref_mate1 = fastq("Mate", "", "");
    const fastq ref_mate2 = fastq("Mate", "ACGT", "!!#$");
    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
    }

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
    }

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        ASSERT_THROW(fastq::validate_paired_reads(mate1, mate1), fastq_error);
    }
}


TEST(fastq, validate_paired_reads__throws_if_only_mate_1_is_numbered)
{
    const fastq ref_mate2 = fastq("Mate/1", "GCTAA", "$!@#$");
    const fastq ref_mate1 = fastq("Mate", "ACGT", "!!#$");

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
    }

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;

        ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
    }
}


TEST(fastq, validate_paired_reads__throws_if_only_mate_2_is_numbered)
{
    const fastq ref_mate1 = fastq("Mate", "GCTAA", "$!@#$");
    const fastq ref_mate2 = fastq("Mate/2", "ACGT", "!!#$");

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
    }

    {
        fastq mate1 = ref_mate1;
        fastq mate2 = ref_mate2;
        ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
    }
}


TEST(fastq, validate_paired_reads__throws_if_mate_is_misnumbered)
{
   fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
   fastq mate2 = fastq("Mate/3", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_same_mate_numbers)
{
   fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
   fastq mate2 = fastq("Mate/1", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_name_differs)
{
   fastq mate1 = fastq("Mate/1", "GCTAA", "$!@#$");
   fastq mate2 = fastq("WrongName/2", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}

} // namespace ar
