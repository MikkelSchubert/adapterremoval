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
#include <sstream>

#include "fastq.h"
#include "linereader.h"

namespace ar
{

struct mate_info
{
    mate_info()
      : name()
      , mate(unknown)
    {}

    std::string desc() const
    {
        switch (mate) {
            case unknown: return "unknown";
            case mate1: return "mate 1";
            case mate2: return "mate 2";
            default:
                throw std::invalid_argument("Invalid mate in mate_info::desc");
        }
    }

    std::string name;
    enum { unknown, mate1, mate2 } mate;
};


inline mate_info get_and_fix_mate_info(fastq& read, char mate_separator)
{
    mate_info info;
    std::string& header = read.m_header;

    size_t pos = header.find_first_of(' ');
    if (pos == std::string::npos) {
        pos = header.length();
    }

    if (pos >= 2 && header.at(pos - 2) == mate_separator) {
        const char digit = header.at(pos - 1);

        if (digit == '1') {
            header[pos - 2] = MATE_SEPARATOR;
            info.mate = mate_info::mate1;
            pos -= 2;
        } else if (digit == '2') {
            header[pos - 2] = MATE_SEPARATOR;
            info.mate = mate_info::mate2;
            pos -= 2;
        }
    }

    info.name = header.substr(0, pos);
    return info;
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
             const fastq_encoding& encoding)
    : m_header(header)
    , m_sequence(sequence)
    , m_qualities(qualities)
{
    process_record(encoding);
}


fastq::fastq(const std::string& header,
             const std::string& sequence)
    : m_header(header)
    , m_sequence(sequence)
    , m_qualities(std::string(sequence.length(), '!'))
{
    process_record(FASTQ_ENCODING_33);
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

    size_t right_exclusive = 0;
    for (size_t i = m_sequence.length(); i; --i) {
        if ((!trim_ns || m_sequence.at(i - 1) != 'N') && (m_qualities.at(i - 1) > low_quality)) {
            right_exclusive = i;
            break;
        }
    }

    size_t left_inclusive = 0;
    for (size_t i = 0; i < right_exclusive; ++i) {
        if ((!trim_ns || m_sequence.at(i) != 'N') && (m_qualities.at(i) > low_quality)) {
            left_inclusive = i;
            break;
        }
    }

    const ntrimmed summary(left_inclusive, m_sequence.length() - right_exclusive);

    if (summary.first || summary.second) {
        const size_t retained = right_exclusive - left_inclusive;
        m_sequence = m_sequence.substr(left_inclusive, retained);
        m_qualities = m_qualities.substr(left_inclusive, retained);
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


void fastq::add_postfix_to_header(const std::string& postfix)
{
    m_header += postfix;
}


void fastq::discard()
{
    m_sequence = 'N';
    m_qualities = '!';
}



bool fastq::read(line_reader_base& reader, const fastq_encoding& encoding)
{
    std::string line;
    if (!reader.getline(line)) {
        // End of file; terminate gracefully
        return false;
    }

    m_header = line.substr(1);
    if (m_header.empty() || line.at(0) != '@') {
        throw fastq_error("Malformed or empty FASTQ header");
    }

    if (!reader.getline(m_sequence)) {
        throw fastq_error("partial FASTQ record; cut off after header");
    } else if (m_sequence.empty()) {
        throw fastq_error("sequence is empty");
    }

    if (!reader.getline(line)) {
        throw fastq_error("partial FASTQ record; cut off after sequence");
    } else if (line.empty() || line.at(0) != '+') {
        throw fastq_error("FASTQ record lacks separator character (+)");
    }

    if (!reader.getline(m_qualities)) {
        throw fastq_error("partial FASTQ record; cut off after separator");
    } else if (m_sequence.empty()) {
        throw fastq_error("sequence is empty");
    }

    process_record(encoding);
    return true;
}


std::string fastq::to_str(const fastq_encoding& encoding) const
{
    std::string result;
    // Size of header, sequence, qualities, 4 new-lines, '@' and '+'
    result.reserve(m_header.size() + m_sequence.size() * 2 + 6);

    result.push_back('@');
    result.append(m_header);
    result.push_back('\n');
    result.append(m_sequence);
    result.append("\n+\n", 3);
    result.append(m_qualities);
    result.push_back('\n');

    // Encode quality-scores in place
    size_t quality_start = m_header.size() + m_sequence.size() + 5;
    size_t quality_end = quality_start + m_sequence.size();
    encoding.encode_string(result.begin() + quality_start,
                           result.begin() + quality_end);

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
                throw fastq_error("invalid character in FASTQ sequence; "
                                  "only A, C, G, T and N are expected!");
        }
    }
}


char fastq::p_to_phred_33(double p)
{
    // Lowest possible error rate representable is '~' (~5e-10)
    const auto min_p = std::max(5e-10, p);
    const auto raw_score = static_cast<int>(-10.0 * std::log10(min_p));
    return std::min<int>('~', raw_score + PHRED_OFFSET_33);
}


void fastq::validate_paired_reads(fastq& mate1, fastq& mate2,
                                  char mate_separator)
{
    if (mate1.length() == 0 || mate2.length() == 0) {
        throw fastq_error("Pair contains empty reads");
    }

    const mate_info info1 = get_and_fix_mate_info(mate1, mate_separator);
    const mate_info info2 = get_and_fix_mate_info(mate2, mate_separator);

    if (info1.name != info2.name) {
        std::stringstream error;
        error << "Pair contains reads with mismatching names:\n"
              << " - '" << info1.name << "'\n"
              << " - '" << info2.name << "'";

        if (info1.mate == mate_info::unknown || info2.mate == mate_info::unknown) {
            error << "\n\nNote that AdapterRemoval by determines the mate "
                     "numbers as the digit found at the end of the read name, "
                     "if this is preceeded by the character '"
                  << mate_separator
                  << "'; if these data makes use of a different character to "
                     "separate the mate number from the read name, then you "
                     "will need to set the --mate-separator command-line "
                     "option to the appropriate character.";
        }

        throw fastq_error(error.str());
    }

    if (info1.mate != mate_info::unknown || info2.mate != mate_info::unknown) {
        if (info1.mate != mate_info::mate1 || info2.mate != mate_info::mate2) {
            std::stringstream error;
            error << "Inconsistent mate numbering; please verify data:\n"
                  << "\nRead 1 identified as " << info1.desc() << ": " << mate1.name()
                  << "\nRead 2 identified as " << info2.desc() << ": " << mate2.name();

            throw fastq_error(error.str());
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// Private helper functions

void fastq::process_record(const fastq_encoding& encoding)
{
    if (m_qualities.length() != m_sequence.length()) {
        throw fastq_error("invalid FASTQ record; sequence/quality length does not match");
    }

    clean_sequence(m_sequence);
    encoding.decode_string(m_qualities.begin(), m_qualities.end());
}

} // namespace ar
