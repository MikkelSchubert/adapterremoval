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
#include <sstream>
#include <iostream>

#include "fastq_enc.h"

namespace ar
{

///////////////////////////////////////////////////////////////////////////////
// fastq_error

fastq_error::fastq_error(const std::string& message)
    : std::exception()
    , m_message(message)
{
}


fastq_error::fastq_error(const fastq_error& error)
    : std::exception()
    , m_message(error.m_message)
{
}


fastq_error::~fastq_error() noexcept
{
}


const char* fastq_error::what() const noexcept
{
    return m_message.c_str();
}


///////////////////////////////////////////////////////////////////////////////
// Pre-calculation of Solexa <-> Phred conversions

std::string calc_solexa_to_phred()
{
    std::string scores;
    scores.resize(MAX_PHRED_SCORE - MIN_SOLEXA_SCORE + 1);

    for (int i = MIN_SOLEXA_SCORE; i <= MAX_PHRED_SCORE; ++i) {
        const double score = round(10.0 * log10(1.0 + pow(10, (i / 10.0))));
        const int transformed = std::max<int>(MIN_PHRED_SCORE, std::min<int>(MAX_PHRED_SCORE, score));
        scores.at(i - MIN_SOLEXA_SCORE) = transformed;
    }

    return scores;
}


std::string calc_phred_to_solexa()
{
    std::string scores;
    scores.resize(MAX_PHRED_SCORE - MIN_PHRED_SCORE + 1);

    for (int i = MIN_PHRED_SCORE; i <= MAX_PHRED_SCORE; ++i) {
        const auto min_i = std::max(1, i);
        const auto score = round(10.0 * log10(pow(10.0, min_i / 10.0) - 1.0));
        const auto transformed = std::max<int>(MIN_SOLEXA_SCORE, std::min<int>(MAX_PHRED_SCORE, score));

        scores.at(i) = transformed;
    }

    return scores;
}


const std::string g_solexa_to_phred = calc_solexa_to_phred();
const std::string g_phred_to_solexa = calc_phred_to_solexa();


///////////////////////////////////////////////////////////////////////////////

void invalid_phred(const char offset, const char max_score, const char raw)
{
    if (raw < offset) {
        if (offset == 33) {
            throw fastq_error("ASCII value of quality score is less than 33 "
                              "(ASCII < '!'); input is corrupt or not in "
                              "FASTQ format!");
        } else if (offset == 64) {
            if (raw < ';') {
                throw fastq_error("Phred+64 encoded quality score is less than 0 "
                                  "(ASCII < '@'); Are these FASTQ reads actually in "
                                  "Phred+33 format? If so, use the command-line "
                                  "option \"--qualitybase 33\"\n\n"

                                  "See README for more information.");
            } else if (raw < '@') {
                // Value not less than -5, which is the lowest Solexa score
                throw fastq_error("Phred+64 encoded quality score is less than 0 "
                                  "(ASCII < '@'); Are these FASTQ reads actually in "
                                  "Phred+33 or Solexa format? If so, use the "
                                  "command-line option \"--qualitybase 33\" or "
                                  "\"--qualitybase solexa\"\n\n"

                                  "See README for more information.");
            }
        } else {
            throw std::logic_error("Unexpected offset in fastq_encoding::decode");
        }
    } else if (raw > max_score) {
        if (raw > '~') {
            throw fastq_error("ASCII value of quality score is greater than "
                              "126 (ASCII > '~'); input is corrupt or not in "
                              "FASTQ format!");
        } else if (offset == 33) {
            std::stringstream ss;

            ss << "Phred+33 encoded quality score is greater than the "
               << "expected maximum (" << max_score << " = "
               << static_cast<char>(offset + max_score) << "). Please "
               << "verify the format of these files.\n\n"

               << "If the quality scores are actually Phred+64 encoded, then "
               << "use the '--qualitybase 64' command-line option.\n\n"

               << "If the quality scores are Phred+33 encoded, but includes "
               << "scores in a greater range than expected, then use the "
               << "'--maxquality' option. Note that this option effects both "
               << "reading and writing of FASTQ files.\n\n"

               << "See README for more information.";

            throw fastq_error(ss.str());
        } else if (offset == 64) {
            std::stringstream ss;

            ss << "Phred+64 encoded quality score is greater than the "
               << "expected maximum (" << max_score << " = "
               << static_cast<char>(offset + max_score) << "). Please "
               << "verify the format of these files.\n\n"

               << "If the quality scores are Phred+64 encoded, but includes "
               << "scores in a greater range than expected, then use the "
               << "'--maxquality' command-line option. Note that this option "
               << "effects both reading and writing of FASTQ files.\n\n"

               << "See README for more information.";

            throw fastq_error(ss.str());
        } else {
            throw std::logic_error("Unexpected offset in fastq_encoding::decode");
        }
    }
}


void invalid_solexa(const char offset, const char max_score, const char raw)
{
    if (raw < ';') {
        if (raw < '!') {
            throw fastq_error("ASCII value of quality score is less than 33 "
                              "(ASCII < '!'); input is corrupt or not in "
                              "FASTQ format!");
        } else {
            throw fastq_error("Solexa score is less than -5 (ASCII = ';'); "
                              "Is this actually Phred+33 data? If so, use "
                              "the '--qualitybase 33' command-line option.\n\n"
                              "See the README for more information.");
        }
    } else if (raw > offset + max_score) {
        if (raw > '~') {
            throw fastq_error("ASCII value of quality score is greater than "
                              "126 (ASCII > '~'); input is corrupt or not in "
                              "FASTQ format!");
        } else {
            std::stringstream ss;

            ss << "Solaxa encoded quality score is greater than the "
               << "expected maximum (" << max_score << " = "
               << static_cast<char>(offset + max_score) << "). Please "
               << "verify the format of these files.\n\n"

               << "If the quality scores are Solexa encoded, but includes "
               << "scores in a greater range than expected, then use the "
               << "'--maxquality' command-line option. Note that this option "
               << "effects both reading and writing of FASTQ files.\n\n"

               << "See README for more information.";

            throw fastq_error(ss.str());
        }
    }
}


///////////////////////////////////////////////////////////////////////////////

fastq_encoding::fastq_encoding(char offset, char max_score)
  : m_offset(offset)
  , m_max_score(std::min<size_t>('~' - offset, max_score))
{
    if (offset != 33 && offset != 64) {
        throw std::invalid_argument("Phred offset must be 33 or 64");
    } else if (max_score < 0) {
        throw std::invalid_argument("Max ASCII encoded Phred score less than 0");
    } else if (max_score > '~' - PHRED_OFFSET_33) {
        throw std::invalid_argument("ASCII value cutoff for quality scores "
                                    "lies after printable characters");
    }
}


fastq_encoding::~fastq_encoding()
{
}


void fastq_encoding::encode_string(std::string::iterator it,
                                   const std::string::iterator& end) const
{
    const char ascii_max = m_offset + m_max_score;
    const char offset = m_offset - '!';

    for (; it != end; ++it) {
        *it = std::min<int>(ascii_max, *it + offset);
    }
}


void fastq_encoding::decode_string(std::string::iterator it,
                                   const std::string::iterator& end) const
{
    const char max_score = m_offset + m_max_score;
    for (; it != end; ++it) {
        const char raw = *it;

        if (raw < m_offset || raw > max_score) {
            invalid_phred(m_offset, m_max_score, raw);
        }

        *it = raw - m_offset + PHRED_OFFSET_33;
    }
}


std::string fastq_encoding::name() const
{
    if (m_offset == 33) {
        return "Phred+33";
    } else if (m_offset == 64) {
        return "Phred+64";
    } else {
        throw std::logic_error("Unexpected offset in fastq_encoding::name");
    }
}


size_t fastq_encoding::max_score() const
{
    return m_max_score;
}


fastq_encoding_solexa::fastq_encoding_solexa(unsigned max_score)
  : fastq_encoding(PHRED_OFFSET_64, max_score)
{
}


void fastq_encoding_solexa::encode_string(std::string::iterator it,
                                   const std::string::iterator& end) const
{
    const char ascii_max = m_offset + m_max_score;

    for (; it != end; ++it) {
        *it = std::min<int>(ascii_max, g_phred_to_solexa.at(*it - '!') + '@');
    }
}


void fastq_encoding_solexa::decode_string(std::string::iterator it,
                                   const std::string::iterator& end) const
{
    const char max_score = m_offset + m_max_score;
    for (; it != end; ++it) {
        const char raw = *it;

        if (raw < ';' || raw > max_score) {
            invalid_phred(m_offset, m_max_score, raw);
        }

        *it = g_solexa_to_phred.at(raw - ';') + '!';
    }
}


std::string fastq_encoding_solexa::name() const
{
    return "Solexa";
}

} // namespace ar
