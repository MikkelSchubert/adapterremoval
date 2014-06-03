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
#include <istream>
#include <stdexcept>

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


inline void convert_qualities(std::string& qualities, quality_format base)
{
    for (std::string::iterator it = qualities.begin(); it != qualities.end(); ++it) {
        int quality = *it;
        if (base == solexa) {
            // Only his range of values are mapped to Phred scores (see above)
            quality -= PHRED_OFFSET_64;
            if (quality < MIN_SOLEXA_SCORE || quality > MAX_SOLEXA_SCORE) {
                throw fastq_error("invalid solexa quality score");
            }

            quality = g_solexa_to_phred.at(quality - MIN_SOLEXA_SCORE);
        } else if (base == phred_64) {
            quality -= PHRED_OFFSET_64;
        } else {
            quality -= PHRED_OFFSET_33;
        }

        if (quality < MIN_PHRED_SCORE) {
            throw fastq_error("phred score less than 0; are these solexa scores?");
        } else if (quality > MAX_PHRED_SCORE) {
            throw fastq_error("phred score greater than the maximum allowed score");
        }

        *it = static_cast<char>(PHRED_OFFSET_33 + quality);
    }
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


bool fastq::read(std::istream& instream, quality_format encoding)
{
    char at_header = -1;
    instream >> at_header;
    if (instream.eof()) {
        return false;
    } else if (instream.fail()) {
        throw fastq_error("IO error while reading FASTQ header");
    } else if (at_header != '@') {
        throw fastq_error("FASTQ header did not start with '@'  ");
    }

    std::getline(instream, m_header);
    if (instream.eof()) {
        // partial header read before EOF
        throw fastq_error("EOF while reading FASTQ header");
    } else if (instream.fail()) {
        throw fastq_error("IO error while reading FASTQ header");
    } else if (m_header.empty()) {
        throw fastq_error("FASTQ header is empty");
    }

    std::getline(instream, m_sequence);
    if (instream.eof()) {
        throw fastq_error("partial FASTQ record; eof encountered instead of sequence");
    } else if (instream.fail()) {
        throw fastq_error("IO error reading FASTQ sequence");
    } else if (m_sequence.empty()) {
        throw fastq_error("sequence is empty");
    }

    std::string separator;
    getline(instream, separator);
    if (instream.eof()) {
        throw fastq_error("partial FASTQ record; eof encountered instead of separator");
    } else if (instream.fail()) {
        throw fastq_error("IO error reading FASTQ separator");
    } else if (separator.empty() || separator.at(0) != '+') {
        throw fastq_error("partial FASTQ record; expected separator (+) not found");
    }

    std::getline(instream, m_qualities);
    if (instream.eof() && m_qualities.empty()) {
        throw fastq_error("partial FASTQ record; eof encountered instead of qualities");
    } else if (instream.fail()) {
        throw fastq_error("IO error reading FASTQ qualities");
    } else if (m_qualities.empty()) {
        throw fastq_error("sequence is empty");
    }

    process_record(encoding);
    return true;
}


bool fastq::write(std::ostream& outstream, quality_format encoding) const
{
    std::string qualities = m_qualities;
    if (encoding == phred_64) {
        for (std::string::iterator it = qualities.begin(); it != qualities.end(); ++it) {
            // Phred+64 is limite to scores in the range 64 + (0 .. 40)
            *it = std::min<char>(64 + 40, *it + 31);
        }
    } else if (encoding == solexa) {
        throw std::invalid_argument("writing solexa scores not supported");
    }

    outstream
        << "@" << m_header << "\n"
        << m_sequence << "\n"
        << "+\n"
        << qualities << "\n";

    return !outstream.fail();
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
