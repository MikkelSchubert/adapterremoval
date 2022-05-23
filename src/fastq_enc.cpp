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
#include <algorithm> // for min, max
#include <cmath>     // for log10, pow, round
#include <iostream>  // for operator<<, basic_ostream, basic_ostream::opera...
#include <sstream>   // for stringstream
#include <stdexcept> // for invalid_argument

#include "debug.hpp" // for AR_FAIL, AR_REQUIRE
#include "fastq_enc.hpp"

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

fastq_error::~fastq_error() {}

const char*
fastq_error::what() const noexcept
{
  return m_message.c_str();
}

///////////////////////////////////////////////////////////////////////////////
// Pre-calculation of Solexa <-> Phred conversions

std::string
calc_solexa_to_phred()
{
  std::string scores;
  scores.resize(MAX_PHRED_SCORE - MIN_SOLEXA_SCORE + 1);

  for (int i = MIN_SOLEXA_SCORE; i <= MAX_PHRED_SCORE; ++i) {
    const double score = round(10.0 * log10(1.0 + pow(10, (i / 10.0))));
    const int transformed =
      std::max<int>(MIN_PHRED_SCORE, std::min<int>(MAX_PHRED_SCORE, score));
    scores.at(i - MIN_SOLEXA_SCORE) = transformed;
  }

  return scores;
}

const std::string g_solexa_to_phred = calc_solexa_to_phred();

///////////////////////////////////////////////////////////////////////////////

void
validate_phred_format(char raw_score)
{
  if (raw_score < '!' || raw_score > '~') {
    std::stringstream ss;

    ss << "Found FASTQ quality score outside of the range of valid ASCII "
       << "encoded values (min = '!', max = '~'). Input file is either corrupt "
          "not in FASTQ format!";

    throw fastq_error(ss.str());
  }
}

[[noreturn]] void
invalid_phred_33(const int max_score, const char raw_score)
{
  validate_phred_format(raw_score);

  const int score = raw_score - PHRED_OFFSET_33;
  const char max_score_ascii = max_score + PHRED_OFFSET_33;
  AR_REQUIRE(score > max_score, "invalid_phred called on valid PHRED score");

  std::stringstream ss;
  ss << "Found Phred+33 encoded quality score of " << score << " ('"
     << raw_score << "'), which is greater than the expected maximum score of "
     << max_score << " ('" << max_score_ascii << "'). Please verify the format "
     << "of these files.\n\n"

     << "If the quality scores are actually Phred+64 encoded, which would mean "
     << "that the Phred quality score is " << raw_score - PHRED_OFFSET_64
     << ", then use the '--qualitybase 64' command-line option.\n\n"

     << "If the quality scores are Phred+33 encoded, but have higher than "
        "expected scores, then increase the '--qualitymax' value.\n\n"

     << "See the documentation for more information.";

  throw fastq_error(ss.str());
}

[[noreturn]] void
invalid_phred_64(const int max_score, const char raw_score)
{
  validate_phred_format(raw_score);

  const int score = raw_score - PHRED_OFFSET_64;
  const char max_score_ascii = max_score + PHRED_OFFSET_64;
  std::stringstream ss;

  if (score < MIN_SOLEXA_SCORE) {
    ss << "Found Phred+64 encoded quality score of " << score << " ('"
       << raw_score << "'), which is less than the expected minimum score "
       << "of " << MIN_PHRED_SCORE << " ('@'). Please verify the format of "
       << "these files.\n\n"

       << "If the quality scores are actually Phred+33 encoded, which would "
       << "mean that the encoded Phred quality score is "
       << raw_score - PHRED_OFFSET_33 << ", then use the '--qualitybase 33' "
       << "command-line option.\n\n"

       << "See the documentation for more information.";
  } else if (score < MIN_PHRED_SCORE) {
    ss << "Found Phred+64 encoded quality score of " << score << " ('"
       << raw_score << "'), which is less than the expected minimum score "
       << "of " << MIN_PHRED_SCORE << " ('@'). Please verify the format of "
       << "these files.\n\n"

       << "If the quality scores are actually Phred+33 encoded, which would "
       << "mean that the encoded Phred quality score is "
       << raw_score - PHRED_OFFSET_33 << ", then use the '--qualitybase 33' "
       << "command-line option.\n\n"

       << "The quality score could also be the older Solexa format, which has "
       << "a minimum score of -5, but data of this type is rare. If it is "
       << "actually Solexa encoded FASTQ data, then use the '--qualitybase "
       << "solexa' command-line option.\n\n"

       << "See the documentation for more information.";
  } else if (score > max_score) {
    ss << "Found Phred+64 encoded quality score of " << score << " ('"
       << raw_score << "'), which is greater than the expected "
       << "maximum score of " << max_score << " ('" << max_score_ascii << "'). "
       << "Please verify the format of these files.\n\n"

       << "If the quality scores are Phred+64 encoded, but have higher than "
       << "expected scores, then increase the '--qualitymax' value.\n\n"

       << "See the documentation for more information.";
  } else {
    AR_FAIL("invalid_phred called on valid PHRED score");
  }

  throw fastq_error(ss.str());
}

[[noreturn]] void
invalid_solexa(const int max_score, const char raw_score)
{
  validate_phred_format(raw_score);

  const int score = raw_score - PHRED_OFFSET_64;
  const char max_score_ascii = max_score + PHRED_OFFSET_64;
  std::stringstream ss;

  if (score < MIN_SOLEXA_SCORE) {
    ss << "Found Solexa encoded quality score of " << score << " ('"
       << raw_score << "'), which is less than the expected minimum score "
       << "of " << MIN_SOLEXA_SCORE << " (';'). Please verify the format of "
       << "these files.\n\n"

       << "If the quality scores are actually Phred+33 encoded, which would "
       << "mean that the encoded Phred quality score is "
       << raw_score - PHRED_OFFSET_33 << ", then use the '--qualitybase 33' "
       << "command-line option.\n\n"

       << "See the documentation for more information.";
  } else if (score > max_score) {
    ss << "Found Solexa encoded quality score of " << score << " ('"
       << raw_score << "'), which is greater than the expected "
       << "maximum score of " << max_score << " ('" << max_score_ascii << "'). "
       << "Please verify the format of these files.\n\n"

       << "If the quality scores are Solexa encoded, but have higher than "
       << "expected scores, then increase the '--qualitymax' value.\n\n"

       << "See the documentation for more information.";
  } else {
    AR_FAIL("invalid_phred called on valid PHRED score");
  }

  throw fastq_error(ss.str());
}

///////////////////////////////////////////////////////////////////////////////

fastq_encoding::fastq_encoding(quality_encoding encoding, char max_score)
  : m_encoding(encoding)
  , m_offset(PHRED_OFFSET_33)
  , m_max_score()
{
  switch (encoding) {
    case quality_encoding::phred_33:
    case quality_encoding::phred_64:
      m_offset = static_cast<char>(encoding);
      break;

    case quality_encoding::solexa:
      m_offset = 64;
      break;

    default:
      AR_FAIL("unknown encoding");
  }

  m_max_score = std::min<char>(max_score, MAX_PHRED_SCORE - (m_offset - '!'));

  AR_REQUIRE(max_score >= MIN_PHRED_SCORE && max_score <= MAX_PHRED_SCORE);
}

void
fastq_encoding::encode(const std::string& qualities, std::string& dst)
{
  dst.append(qualities);
}

void
fastq_encoding::decode(std::string& qualities) const
{
  const char max_score = m_offset + m_max_score;

  if (m_encoding == quality_encoding::solexa) {
    for (auto& quality : qualities) {
      if (quality < ';' || quality > max_score) {
        invalid_solexa(m_max_score, quality);
      }

      quality = g_solexa_to_phred.at(quality - ';') + PHRED_OFFSET_33;
    }
  } else {
    for (auto& quality : qualities) {
      if (quality < m_offset || quality > max_score) {
        if (m_offset == 33) {
          invalid_phred_33(m_max_score, quality);
        } else {
          invalid_phred_64(m_max_score, quality);
        }
      }

      quality -= m_offset - PHRED_OFFSET_33;
    }
  }
}
