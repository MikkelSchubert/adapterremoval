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
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "fastq.h"


///////////////////////////////////////////////////////////////////////////////
// Utility functions

std::string calc_solexa_to_phred()
{
    const char scores[] = {
    1,  1,  2,  2,  3, // -5 .. -1
    3,  4,  4,  5,  5, //  0 ..  4
    6,  7,  8,  9, 10, //  5 ..  9
   10, 11, 12, 13, 14, // 10 .. 14
   15, 16, 17, 18, 19, // 15 .. 19
   20, 21, 22, 23, 24, // 20 .. 24
   25, 26, 27, 28, 29, // 25 .. 29
   30, 31, 32, 33, 34, // 30 .. 34
   35, 36, 37, 38, 39, // 35 .. 39
   40, 41, '\0'};      // 40 .. 41

   return std::string(scores);
}


//! The corresponding Phred score for each solexa score
const std::string g_solexa_to_phred = calc_solexa_to_phred();


char process_scores_phred_33(char quality)
{
    if (quality < '!') { // Neither solexa nor Phred+33
        throw fastq_error("ASCII value of Phred score is less than 33 = '!'; "
                          "input is corrupt or not FASTQ format!");
    } else if (quality > 'J') {
        throw fastq_error("Phred+33 score is greater than 41 (ASCII = 'J'); "
                          "Are these FASTQ reads actually in Phred+64 format?\n"
                          "If so, use the option \"--qualitybase 64\"");
    }

    return quality;
}


char process_scores_phred_64(char quality)
{
    // The scores from '@' (0) to 'h' (40) are expected for Phred+64 reads
    // As an older format, scores outside of this range is not expected
    if (quality < '@') {
        if (quality < ';') {
            throw fastq_error("Phred+64 score is less than 0 (ASCII = '@');\n"
                              "Are these FASTQ reads actually in Phred+33 format?\n"
                              "If so, use the option \"--qualitybase 33\"");
        } else {
            // Value not less than -5, which is the lowest Solexa score
            throw fastq_error("Phred+64 score is less than 0 (ASCII = '@');\n"
                              "Are these FASTQ reads actually in Solexa format?\n"
                              "If so, use the option \"--qualitybase solexa\"");
        }
    } else if (quality > 'h') {
        throw fastq_error("Phred+64 score is greater than 40 (ASCII = 'h'); "
                          "input is corrupt or not FASTQ format!");
    }

    return (quality - '@') + PHRED_OFFSET_33;
}


char process_scores_solexa(char quality)
{
    // The scores from ';' (-5) to 'h' (40) are expected for Solexa reads
    // As an older format, scores outside of this range is not expected
    if (quality < ';') { // TODO: -5, maybe Phred+33
        throw fastq_error("Solexa score is less than -5 (ASCII = ';'); "
                          "Is this actually Phred+33 data?");
    } else if (quality > 'h') {
        throw fastq_error("Solexa score is greater than 40 (ASCII = 'h'); "
                          "This is not a valid Solexa score for FASTQ reads.");
    }

    return g_solexa_to_phred.at(quality - ';') + PHRED_OFFSET_33;
}


char process_scores_ignored(char quality)
{
    if (quality < '!') { // Neither solexa nor Phred+33
        throw fastq_error("ASCII value of Phred score is less than 33 = '!'; "
                          "input is corrupt or not FASTQ format!");
    } else if (quality > 'h') {
        throw fastq_error("Quality score is greater than 'h' (ASCII); "
                          "This is not a valid quality score for FASTQ reads.");
    }

    return PHRED_OFFSET_33;
}



inline void convert_qualities(std::string& qualities, const fastq::quality_format base)
{
    char (*func)(char) = NULL;
    switch (base) {
        case fastq::phred_33:
            func = &process_scores_phred_33;
            break;

        case fastq::phred_64:
            func = &process_scores_phred_64;
            break;

        case fastq::solexa:
            func = &process_scores_solexa;
            break;

        case fastq::ignored:
            func = &process_scores_ignored;
            break;

        default:
            throw std::logic_error("Unhandled FASTQ quality score format in 'convert_qualities'");
    }

    std::transform(qualities.begin(),
                   qualities.end(),
                   qualities.begin(),
                   func);
}


struct mate_info
{
    mate_info()
      : name()
      , mate(unknown)
    {}

    std::string name;
    enum { unknown, mate1, mate2 } mate;
};


inline mate_info get_mate_information(const fastq& read)
{
    mate_info info;
    const std::string& header = read.header();

    size_t pos = header.find_first_of(' ');
    if (pos == std::string::npos) {
        pos = header.length();
    }

    if (pos >= 2) {
        const std::string substr = header.substr(0, pos);
        if (substr.substr(pos - 2) == "/1") {
            info.mate = mate_info::mate1;
            pos -= 2;
        } else if (substr.substr(pos - 2) == "/2") {
            info.mate = mate_info::mate2;
            pos -= 2;
        }
    }

    info.name = header.substr(0, pos);
    return info;
}


///////////////////////////////////////////////////////////////////////////////
// fastq_error

fastq_error::fastq_error(const std::string& message)
    : std::exception()
    , m_message(message)
{
}


fastq_error::~fastq_error() throw()
{
}


const char* fastq_error::what() const throw()
{
    return m_message.c_str();
}


///////////////////////////////////////////////////////////////////////////////
// fastq

fastq::fastq()
    : m_header()
	, m_sequence()
	, m_qualities()
{
}

fastq::fastq(const std::string& header,
             const std::string& sequence,
             const std::string& qualities,
             quality_format encoding)
    : m_header(header)
    , m_sequence(sequence)
    , m_qualities(qualities)
{
    process_record(encoding);
}


bool fastq::operator==(const fastq& other) const
{
    return (m_header == other.m_header)
        && (m_sequence == other.m_sequence)
        && (m_qualities == other.m_qualities);
}


size_t fastq::length() const
{
    return m_sequence.length();
}


size_t fastq::count_ns() const
{
    return static_cast<size_t>(std::count(m_sequence.begin(), m_sequence.end(), 'N'));
}


fastq::ntrimmed fastq::trim_low_quality_bases(bool trim_ns, char low_quality)
{
    low_quality += PHRED_OFFSET_33;
    if (m_sequence.empty()) {
        return ntrimmed();
    }

	size_t lower_bound = m_sequence.length();
	for (size_t i = 0; i < m_sequence.length(); ++i) {
		if ((!trim_ns || m_sequence.at(i) != 'N') && (m_qualities.at(i) > low_quality)) {
			lower_bound = i;
			break;
		}
	}

	size_t upper_bound = lower_bound;
	for (size_t i = m_sequence.length() - 1; i > lower_bound; --i) {
		if ((!trim_ns || m_sequence.at(i) != 'N') && (m_qualities.at(i) > low_quality)) {
			upper_bound = i;
			break;
		}
	}

	const size_t retained = upper_bound - lower_bound + 1;
    const ntrimmed summary(lower_bound, m_sequence.length() - upper_bound - 1);

    if (summary.first || summary.second) {
    	m_sequence = m_sequence.substr(lower_bound, retained);
    	m_qualities = m_qualities.substr(lower_bound, retained);
    }

    return summary;
}


void fastq::truncate(size_t pos, size_t len)
{
    if (pos || len < length()) {
        m_sequence = m_sequence.substr(pos, len);
        m_qualities = m_qualities.substr(pos, len);
    }
}


void fastq::reverse_complement()
{
    std::reverse(m_sequence.begin(), m_sequence.end());
    std::reverse(m_qualities.begin(), m_qualities.end());

    // Lookup table for complementary bases based only on the last 4 bits
    static const char complements[] = "-T-GA--C------N-";
    for(std::string::iterator it = m_sequence.begin(); it != m_sequence.end(); ++it) {
        *it = complements[*it & 0xf];
    }
}


void fastq::add_prefix_to_header(const std::string& prefix)
{
    m_header.insert(0, prefix);
}


bool fastq::read(string_vec_citer& it, const string_vec_citer& end, quality_format encoding)
{
    if (it == end) {
        return false;
    } else {
        const std::string& header_line = *it;
        if (header_line.empty() || header_line.at(0) != '@') {
            throw fastq_error("FASTQ header did not start with '@'  ");
        }

        m_header = it->substr(1);
        if (m_header.empty()) {
            throw fastq_error("FASTQ header is empty");
        }

        ++it;
    }

    if (it == end) {
        throw fastq_error("partial FASTQ record; cut off after header");
    } else {
        m_sequence = *it++;
        if (m_sequence.empty()) {
            throw fastq_error("sequence is empty");
        }
    }

    if (it == end) {
        throw fastq_error("partial FASTQ record; cut off after sequence");
    } else {
        const std::string& separator = *it++;
        if (separator.empty() || separator.at(0) != '+') {
            throw fastq_error("FASTQ record lacks seperator character (+)");
        }
    }

    if (it == end) {
        throw fastq_error("partial FASTQ record; cut off after separator");
    } else {
        m_qualities = *it++;
        if (m_sequence.empty()) {
            throw fastq_error("sequence is empty");
        }
    }

    process_record(encoding);
    return true;
}


std::string fastq::to_str(quality_format encoding) const
{
    std::string result;
    // Size of header, sequence, qualities, 4 new-lines, @ and +
    result.reserve(m_header.size() + m_sequence.size() * 2 + 6);

    result.push_back('@');
    result += m_header;
    result.push_back('\n');
    result += m_sequence;
    result += "\n+\n";

    if (encoding == phred_33) {
        result += m_qualities;
    } else if (encoding == phred_64) {
        for (std::string::const_iterator it = m_qualities.begin(); it != m_qualities.end(); ++it) {
            // Phred+64 is limite to scores in the range 64 + (0 .. 40)
            result.push_back(std::min<char>(64 + 40, *it + 31));
        }
    } else {
        throw std::invalid_argument("writing solexa / ignored scores not supported");
    }

    result.push_back('\n');

    return result;
}



///////////////////////////////////////////////////////////////////////////////
// Public helper functions

void fastq::clean_sequence(std::string& sequence)
{
    for (std::string::iterator it = sequence.begin(); it != sequence.end(); ++it) {
        switch (*it) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                break;

            case 'a':
            case 'c':
            case 'g':
            case 't':
            case 'n':
                *it += 'A' - 'a';
                break;

            case '.':
                *it = 'N';
                break;

            default:
                throw fastq_error("invalid character in FASTQ sequence");
        }
    }
}


char fastq::p_to_phred_33(double p)
{
    const int raw_score = static_cast<int>(-10.0 * std::log10(p));
    const char phred_score = static_cast<char>(std::min<int>(MAX_PHRED_SCORE, raw_score));
    return phred_score + PHRED_OFFSET_33;
}


void fastq::validate_paired_reads(const fastq& mate1, const fastq& mate2)
{
    if (mate1.length() == 0 || mate2.length() == 0) {
        throw fastq_error("Pair contains empty reads");
    }

    const mate_info info1 = get_mate_information(mate1);
    const mate_info info2 = get_mate_information(mate2);

    if (info1.name != info2.name) {
        std::string message = "Pair contains reads with mismatching names: '";
        message += info1.name;
        message += "' and '";
        message += info2.name;
        message += "'";

        throw fastq_error(message);
    }

    if (info1.mate != mate_info::unknown || info2.mate != mate_info::unknown) {
        if (info1.mate != mate_info::mate1 || info2.mate != mate_info::mate2) {
            throw fastq_error("Inconsistent mate numbering");
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// Private helper functions

void fastq::process_record(quality_format encoding)
{
    if (m_qualities.length() != m_sequence.length()) {
        throw fastq_error("invalid FASTQ record; sequence/quality length does not match");
    }

    clean_sequence(m_sequence);
    convert_qualities(m_qualities, encoding);
}
