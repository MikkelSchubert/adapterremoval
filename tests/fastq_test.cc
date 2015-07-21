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

#include "fastq.h"



inline std::istream& operator>>(std::istream& stream, fastq& record)
{
    record.read(stream);
    return stream;
}


inline std::ostream& operator<<(std::ostream& stream, const fastq& record)
{
    record.write(stream);
    return stream;
}


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
    ASSERT_EQ(std::string("!7BF8DGI", 8), record.qualities());
}


TEST(fastq, constructor_simple_record_phred_64_encoded)
{
    const fastq record("record_2", "ACGAGTCA", "@VaeWcfh", fastq::phred_64);
    ASSERT_EQ("record_2", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ(std::string("!7BF8DGI", 8), record.qualities());
}


TEST(fastq, constructor_simple_record_phred_solexa_encoded)
{
    const fastq record("record_3", "AAACGAGTCA", ";h>S\\TCDUJ", fastq::solexa);
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
    ASSERT_NO_THROW(fastq("Rec", "CAT", "!!\"", fastq::phred_33));
    ASSERT_THROW(fastq("Rec", "CAT", " !\"", fastq::phred_33), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "IJJ", fastq::phred_33));
    ASSERT_THROW(fastq("Rec", "CAT", "IJK", fastq::phred_33), fastq_error);
}


TEST(fastq, constructor_score_boundries_phred_64)
{
    ASSERT_NO_THROW(fastq("Rec", "CAT", "@@A", fastq::phred_64));
    ASSERT_THROW(fastq("Rec", "CAT", "?@A", fastq::phred_64), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "ghf", fastq::phred_64));
    ASSERT_THROW(fastq("Rec", "CAT", "ghi", fastq::phred_64), fastq_error);
}


TEST(fastq, constructor_score_boundries_solexa)
{
    ASSERT_NO_THROW(fastq("Rec", "CAT", ";;<", fastq::solexa));
    ASSERT_THROW(fastq("Rec", "CAT", ":;<", fastq::solexa), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "ghf", fastq::solexa));
    ASSERT_THROW(fastq("Rec", "CAT", "ghi", fastq::solexa), fastq_error);
}

TEST(fastq, constructor_score_boundries_ignored)
{
    ASSERT_NO_THROW(fastq("Rec", "CAT", "!!\"", fastq::ignored));
    ASSERT_THROW(fastq("Rec", "CAT", " !\"", fastq::ignored), fastq_error);

    ASSERT_NO_THROW(fastq("Rec", "CAT", "ghf", fastq::ignored));
    ASSERT_THROW(fastq("Rec", "CAT", "ghi", fastq::ignored), fastq_error);
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
// misc properties

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
// trim_low_quality_bases

TEST(fastq, trim_low_quality_bases__empty_record)
{
    fastq record("Empty", "", "");
    const fastq::ntrimmed expected(0, 0);
    ASSERT_EQ(expected, record.trim_low_quality_bases(true, 10));
    ASSERT_EQ(fastq("Empty", "", ""), record);
}


TEST(fastq, trim_low_quality_bases__trim_nothing)
{
    const fastq reference("Rec", "NNNNN", "!!!!!");
    const fastq::ntrimmed expected(0, 0);
    fastq record = reference;
    // Trim neither Ns nor low Phred score bases

    ASSERT_EQ(expected, record.trim_low_quality_bases(false, -1));
    ASSERT_EQ(reference, record);
}


TEST(fastq, trim_low_quality_bases__trim_ns)
{
    const fastq expected_record("Rec", "ANT", "456");
    const fastq::ntrimmed expected_ntrim(2, 0);
    fastq record("Rec", "NNANT", "23456");
    // Trim neither Ns nor low Phred score bases
    ASSERT_EQ(expected_ntrim, record.trim_low_quality_bases(true, -1));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_low_quality_bases__trim_low_quality_bases)
{
    const fastq expected_record("Rec", "TN", "%$");
    const fastq::ntrimmed expected_ntrim(0, 3);
    fastq record("Rec", "TNANT", "%$#!\"");
    // Trim neither Ns nor low Phred score bases
    ASSERT_EQ(expected_ntrim, record.trim_low_quality_bases(false, 2));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_low_quality_bases__trim_mixed)
{
    const fastq expected_record("Rec", "TAG", "$12");
    const fastq::ntrimmed expected_ntrim(3, 2);
    fastq record("Rec", "NTNTAGNT", "1!#$12#\"");
    // Trim neither Ns nor low Phred score bases
    ASSERT_EQ(expected_ntrim, record.trim_low_quality_bases(true, 2));
    ASSERT_EQ(expected_record, record);
}


TEST(fastq, trim_low_quality_bases__trim_mixed__no_low_quality_bases)
{
    const fastq expected_record("Rec", "ACTTAG", "12I$12");
    const fastq::ntrimmed expected_ntrim(0, 0);
    fastq record = expected_record;
    // Trim neither Ns nor low Phred score bases
    ASSERT_EQ(expected_ntrim, record.trim_low_quality_bases(true, 2));
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
// Reading from stream

TEST(fastq, simple_fastq_record_1)
{
    string_list lines;
    lines.push_back("@record_1");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_TRUE(record.read(it, lines.end(), fastq::phred_33));
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
    ASSERT_EQ(it, lines.end());
}

TEST(fastq, simple_fastq_record_2)
{
    string_list lines;
    lines.push_back("@record_1");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    lines.push_back("@record_2");
    lines.push_back("GTCAGGAT");
    lines.push_back("+");
    lines.push_back("D7BIG!F8");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_TRUE(record.read(it, lines.end(), fastq::phred_33));
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
    ASSERT_TRUE(record.read(it, lines.end(), fastq::phred_33));
    ASSERT_EQ("record_2", record.header());
    ASSERT_EQ("GTCAGGAT", record.sequence());
    ASSERT_EQ("D7BIG!F8", record.qualities());
    ASSERT_EQ(it, lines.end());
    ASSERT_FALSE(record.read(it, lines.end(), fastq::phred_33));
}


TEST(fastq, simple_fastq_record__with_extra_header_1)
{
    string_list lines;
    lines.push_back("@record_1 Extra header here");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_TRUE(record.read(it, lines.end(), fastq::phred_33));
    ASSERT_EQ("record_1 Extra header here", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
    ASSERT_EQ(it, lines.end());
    ASSERT_FALSE(record.read(it, lines.end(), fastq::phred_33));
}


TEST(fastq, simple_fastq_record__with_extra_header_2)
{
    string_list lines;
    lines.push_back("@record_1");
    lines.push_back("ACGAGTCA");
    lines.push_back("+Extra header here");
    lines.push_back("!7BF8DGI");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_TRUE(record.read(it, lines.end(), fastq::phred_33));
    ASSERT_EQ("record_1", record.header());
    ASSERT_EQ("ACGAGTCA", record.sequence());
    ASSERT_EQ("!7BF8DGI", record.qualities());
    ASSERT_EQ(it, lines.end());
    ASSERT_FALSE(record.read(it, lines.end(), fastq::phred_33));
}


TEST(fastq, simple_fastq_record__no_header)
{
    string_list lines;
    lines.push_back("@");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end(), fastq::phred_33), fastq_error);
}


TEST(fastq, simple_fastq_record__no_sequence)
{
    string_list lines;
    lines.push_back("@record_1");
    lines.push_back("");
    lines.push_back("+");
    lines.push_back("!7BF8DGI");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end(), fastq::phred_33), fastq_error);
}


TEST(fastq, simple_fastq_record__no_qualities)
{
    string_list lines;
    lines.push_back("@");
    lines.push_back("ACGAGTCA");
    lines.push_back("+");
    lines.push_back("");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end(), fastq::phred_33), fastq_error);
}


TEST(fastq, simple_fastq_record__no_qualities_or_sequence)
{
    string_list lines;
    lines.push_back("@");
    lines.push_back("");
    lines.push_back("+");
    lines.push_back("");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end(), fastq::phred_33), fastq_error);
}


TEST(fastq, eof_when_starting_to_read_record)
{
    string_list lines;
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_FALSE(record.read(it, lines.end()));
    ASSERT_EQ(it, lines.end());
}


TEST(fastq, eof_after_header)
{
    string_list lines;
    lines.push_back("@record");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end()), fastq_error);
}


TEST(fastq, eof_after_sequence_1)
{
    string_list lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end()), fastq_error);
}


TEST(fastq, eof_after_sequence_2)
{
    string_list lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    lines.push_back("");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end()), fastq_error);
}


TEST(fastq, eof_after_sep_1)
{
    string_list lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    lines.push_back("+");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end()), fastq_error);
}


TEST(fastq, eof_after_sep_2)
{
    string_list lines;
    lines.push_back("@record");
    lines.push_back("ACGTA");
    lines.push_back("+");
    lines.push_back("");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_THROW(record.read(it, lines.end()), fastq_error);
}


TEST(fastq, eof_after_qualities_following_previous_read_1)
{
    string_list lines;
    lines.push_back("@record_1");
    lines.push_back("ACGTA");
    lines.push_back("+");
    lines.push_back("!!!!!");
    lines.push_back("@record_2");
    lines.push_back("ACGTA");
    lines.push_back("+");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_NO_THROW(record.read(it, lines.end()));
    ASSERT_THROW(record.read(it, lines.end()), fastq_error);
}


TEST(fastq, eof_after_qualities_following_previous_read_2)
{
    string_list lines;
    lines.push_back("@record_1");
    lines.push_back("ACGTA");
    lines.push_back("+");
    lines.push_back("!!!!!");
    lines.push_back("@record_2");
    lines.push_back("ACGTA");
    lines.push_back("+");
    lines.push_back("");
    string_list_citer it = lines.begin();

    fastq record;
    ASSERT_NO_THROW(record.read(it, lines.end()));
    ASSERT_THROW(record.read(it, lines.end()), fastq_error);
}



///////////////////////////////////////////////////////////////////////////////
// Writing to stream

TEST(fastq, Writing_to_stream_phred_33)
{
    std::stringstream outstream;
    const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
    outstream << record;
    ASSERT_EQ("@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n", outstream.str());
}

TEST(fastq, Writing_to_stream_phred_33_explicit)
{
    std::stringstream outstream;
    const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
    record.write(outstream);
    ASSERT_EQ("@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n", outstream.str());
}

TEST(fastq, Writing_to_stream_phred_64_explicit)
{
    std::stringstream outstream;
    const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
    record.write(outstream, fastq::phred_64);
    ASSERT_EQ("@record_1\nACGTACGATA\n+\n@CBCIUWbfh\n", outstream.str());
}


///////////////////////////////////////////////////////////////////////////////
// Clean sequence

TEST(fastq, clean_empty)
{
    std::string sequence;
    fastq::clean_sequence(sequence);
    ASSERT_EQ("", sequence);
}


TEST(fastq, clean_lowercase)
{
    std::string sequence = "acGtAcngN";
    fastq::clean_sequence(sequence);
    ASSERT_EQ("ACGTACNGN", sequence);
}


TEST(fastq, clean_dots)
{
    std::string sequence = "ACGTAC.G.";
    fastq::clean_sequence(sequence);
    ASSERT_EQ("ACGTACNGN", sequence);
}


TEST(fastq, reject_non_nucleotides_1)
{
    std::string sequence = "AsTACNGN";
    ASSERT_THROW(fastq::clean_sequence(sequence), fastq_error);
}


TEST(fastq, reject_non_nucleotides_2)
{
    std::string sequence = "ACGTAC1GN";
    ASSERT_THROW(fastq::clean_sequence(sequence), fastq_error);
}


///////////////////////////////////////////////////////////////////////////////
// Validating pairs

TEST(fastq, validate_paired_reads__throws_if_order_is_wrong)
{
    const fastq mate1 = fastq("Mate/1", "ACGT", "!!#$");
    const fastq mate2 = fastq("Mate/2", "GCTAA", "$!@#$");
    fastq::validate_paired_reads(mate1, mate2);
    ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_mate_is_empty)
{
    const fastq mate1 = fastq("Mate", "", "");
    const fastq mate2 = fastq("Mate", "ACGT", "!!#$");
    ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
    ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
    ASSERT_THROW(fastq::validate_paired_reads(mate1, mate1), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_only_mate_1_is_numbered)
{
   const fastq mate2 = fastq("Mate/1", "GCTAA", "$!@#$");
   const fastq mate1 = fastq("Mate", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
   ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_only_mate_2_is_numbered)
{
   const fastq mate2 = fastq("Mate", "GCTAA", "$!@#$");
   const fastq mate1 = fastq("Mate/2", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
   ASSERT_THROW(fastq::validate_paired_reads(mate2, mate1), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_mate_is_misnumbered)
{
   const fastq mate2 = fastq("Mate/1", "GCTAA", "$!@#$");
   const fastq mate1 = fastq("Mate/3", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_same_mate_numbers)
{
   const fastq mate2 = fastq("Mate/1", "GCTAA", "$!@#$");
   const fastq mate1 = fastq("Mate/1", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}


TEST(fastq, validate_paired_reads__throws_if_name_differs)
{
   const fastq mate2 = fastq("Mate/1", "GCTAA", "$!@#$");
   const fastq mate1 = fastq("WrongName/2", "ACGT", "!!#$");
   ASSERT_THROW(fastq::validate_paired_reads(mate1, mate2), fastq_error);
}
