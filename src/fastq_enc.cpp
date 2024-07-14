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
#include "fastq_enc.hpp" // header
#include "debug.hpp"     // for AR_FAIL, AR_REQUIRE
#include "errors.hpp"    // for fastq_error
#include "strutils.hpp"  // for shell_escape
#include <algorithm>     // for min, max
#include <array>         // for array
#include <cmath>         // for log10, pow, round
#include <sstream>       // for ostringstream
#include <stdexcept>     // for invalid_argument

namespace adapterremoval {

//! Offset used by Phred+33 and SAM encodings
const int PHRED_33_OFFSET_MIN = '!';
//! The maximum ASCII value allowed for Phred+33 encoded quality scores. This is
//! a generous maximum, taking into account that instruments are assigning
//! higher and higher quality scores for Phred+33 data.
const int PHRED_33_OFFSET_MAX = 'N';
//! Minimum Phred score allowed
const int PHRED_33_SCORE_MIN = 0;
//! Maximum Phred score allowed
const int PHRED_33_SCORE_MAX = PHRED_33_OFFSET_MAX - PHRED_33_OFFSET_MIN;

//! Offset used by Phred+64 encodings
const int PHRED_64_OFFSET_MIN = '@';
//! The maximum ASCII value allowed for encoded Phred scores
const int PHRED_64_OFFSET_MAX = '~';
//! Minimum Phred+64 score allowed
const int PHRED_64_SCORE_MIN = 0;
//! Maximum Phred+64 score allowed
const int PHRED_64_SCORE_MAX = PHRED_64_OFFSET_MAX - PHRED_64_OFFSET_MIN;

//! Offset used by Solexa encoding quality scores
const int SOLEXA_OFFSET_MIN = '@';
//! The maximum ASCII value allowed for encoded Solexa scores
const int SOLEXA_OFFSET_MAX = 'h';
//! Minimum Phred encoded score allowed
const int SOLEXA_SCORE_MIN = -5;
//! Maximum Phred encoded score allowed
const int SOLEXA_SCORE_MAX = SOLEXA_OFFSET_MAX - SOLEXA_OFFSET_MIN;

///////////////////////////////////////////////////////////////////////////////
// Pre-calculation of Solexa <-> Phred conversions

std::array<char, 256>
calc_solexa_to_phred33()
{
  std::array<char, 256> scores = { 0 };

  for (int i = SOLEXA_SCORE_MIN; i <= SOLEXA_SCORE_MAX; ++i) {
    const double score = round(10.0 * log10(1.0 + pow(10, (i / 10.0))));
    const int transformed =
      std::max<int>(PHRED_SCORE_MIN, std::min<int>(PHRED_SCORE_MAX, score));
    scores.at(i + SOLEXA_OFFSET_MIN) = PHRED_OFFSET_MIN + transformed;
  }

  return scores;
}

const auto g_solexa_to_phred33 = calc_solexa_to_phred33();

///////////////////////////////////////////////////////////////////////////////

/** Returns the quality score in a form that is human readable */
std::string
escape_raw_score(char raw)
{
  return shell_escape(std::string(1, raw));
}

[[noreturn]] void
throw_invalid_phred_33(const char raw)
{
  const int score = raw - PHRED_33_OFFSET_MIN;
  const int alt_score = raw - PHRED_64_OFFSET_MIN;
  const int max_score = PHRED_33_SCORE_MAX;
  const auto esc = escape_raw_score(raw);
  const auto esc_max = escape_raw_score(PHRED_33_OFFSET_MAX);

  std::ostringstream ss;
  ss << "Found Phred+33 encoded quality score of " << score << " (encoded as "
     << esc << "), which is greater than the expected maximum score of "
     << max_score << " (encoded as " << esc_max << "). This suggests that "
     << "input may be Phred+64 encoded.\n\n"

     << "If the quality scores are actually Phred+64 encoded, which would "
     << "mean that the Phred quality score is " << alt_score
     << ", then use the '--quality-format 64' command-line option.\n\n"

     << "If the quality scores are Phred+33 encoded, but have higher than "
        "expected scores, then use '--quality-format sam' to permit the full "
        "range of quality scores.\n\n"

     << "See the documentation for more information.";

  throw fastq_error(ss.str());
}

[[noreturn]] void
throw_invalid_phred_64(const char raw)
{
  AR_REQUIRE(raw < SOLEXA_OFFSET_MIN,
             "invalid_phred called on valid PHRED score");

  const int score = raw - PHRED_64_OFFSET_MIN;
  const int alt_score = raw - PHRED_33_OFFSET_MIN;
  const int min_score = PHRED_64_SCORE_MIN;
  const auto esc = escape_raw_score(raw);
  const auto esc_min = escape_raw_score(PHRED_64_OFFSET_MIN);

  std::ostringstream ss;
  ss << "Found Phred+64 encoded quality score of " << score << " (encoded as "
     << esc << "), which is less than the expected minimum score of "
     << min_score << " (encoded as " << esc_min << "). This suggests that "
     << "input may be Phred+33 encoded.\n\n"

     << "If the quality scores are actually Phred+33 encoded, which would "
     << "mean that the Phred quality score is " << alt_score
     << ", then omit the '--quality-format' command-line option or use "
     << "'--quality-format 33' to explicitly set the format.\n\n";

  if (score >= SOLEXA_SCORE_MIN) {
    ss << "The quality score could also be the older Solexa format, which has "
       << "a minimum score of -5, but data of this type is rare. If it is "
       << "actually Solexa encoded FASTQ data, then use the '--quality-format "
       << "solexa' command-line option.\n\n";
  }

  ss << "See the documentation for more information.";

  throw fastq_error(ss.str());
}

[[noreturn]] void
throw_invalid_solexa(const char raw)
{
  const int score = raw - SOLEXA_OFFSET_MIN;
  const int alt_score = raw - PHRED_33_OFFSET_MIN;
  const int min_score = SOLEXA_SCORE_MIN;
  const auto esc = escape_raw_score(raw);
  const auto esc_min = escape_raw_score(SOLEXA_OFFSET_MIN + SOLEXA_SCORE_MIN);
  const auto esc_max = escape_raw_score(SOLEXA_OFFSET_MAX);

  std::ostringstream ss;
  if (score < SOLEXA_SCORE_MIN) {
    ss << "Found Solexa encoded quality score of " << score << " (encoded as "
       << esc << "), which is less than the expected minimum score of "
       << min_score << " (encoded as " << esc_min << "). This suggests that "
       << "input may be Phred+33 encoded.\n\n"

       << "If the quality scores are actually Phred+33 encoded, which would "
       << "mean that the Phred quality score is " << alt_score
       << ", then omit the '--quality-format' command-line option or use "
       << "'--quality-format 33' to explicitly set the format.\n\n"

       << "See the documentation for more information.";
  } else if (raw > SOLEXA_OFFSET_MAX) {
    ss << "Found Solexa encoded quality score of " << score << " (encoded as "
       << esc << "), which is greater than the expected maximum score of "
       << SOLEXA_SCORE_MAX << " (encoded as " << esc_max << ").\n\n"

       << "See the documentation for more information.";
  } else {
    AR_FAIL("invalid_phred called on valid PHRED score");
  }

  throw fastq_error(ss.str());
}

[[noreturn]] void
throw_invalid_score(const quality_encoding encoding, const char raw_score)
{
  if (raw_score < PHRED_OFFSET_MIN || raw_score > PHRED_OFFSET_MAX) {
    std::ostringstream ss;

    ss << "Found raw FASTQ quality score of " << static_cast<int>(raw_score)
       << " (" << escape_raw_score(raw_score) << "). This is outside of the "
       << "range of valid ASCII encoded quality scores (" << PHRED_OFFSET_MIN
       << " to " << PHRED_OFFSET_MAX << "), meaning that the input file is "
       << "either corrupt or not in FASTQ format!";

    throw fastq_error(ss.str());
  }

  switch (encoding) {
    case quality_encoding::phred_33:
      throw_invalid_phred_33(raw_score);

    case quality_encoding::phred_64:
      throw_invalid_phred_64(raw_score);

    case quality_encoding::solexa:
      throw_invalid_solexa(raw_score);

    default:
      AR_FAIL("This case should have been handled by initial check");
  }
}

///////////////////////////////////////////////////////////////////////////////

fastq_encoding::fastq_encoding(quality_encoding encoding) noexcept
  : m_encoding(encoding)
  , m_offset_min()
  , m_offset_max()
{
  switch (encoding) {
    case quality_encoding::phred_33:
      m_offset_min = PHRED_33_OFFSET_MIN;
      m_offset_max = PHRED_33_OFFSET_MAX;
      break;

    case quality_encoding::phred_64:
      m_offset_min = PHRED_64_OFFSET_MIN;
      m_offset_max = PHRED_64_OFFSET_MAX;
      break;

    case quality_encoding::solexa:
      m_offset_min = SOLEXA_OFFSET_MIN;
      m_offset_max = SOLEXA_OFFSET_MAX;
      break;

    case quality_encoding::sam:
      m_offset_min = PHRED_OFFSET_MIN;
      m_offset_max = PHRED_OFFSET_MAX;
      break;

    default:
      AR_FAIL("unknown encoding");
  }
}

void
fastq_encoding::decode(std::string& qualities) const
{
  if (m_encoding == quality_encoding::solexa) {
    for (auto& quality : qualities) {
      const char current = quality;
      // TODO: Handle negative values, e.g. è
      if (!(quality = g_solexa_to_phred33.at(quality))) {
        throw_invalid_score(m_encoding, current);
      }
    }
  } else {
    for (auto& quality : qualities) {
      if (quality < m_offset_min || quality > m_offset_max) {
        throw_invalid_score(m_encoding, quality);
      }

      // Convert Phred+64 to Phred+33 if needed
      quality -= m_offset_min - PHRED_33_OFFSET_MIN;
    }
  }
}

} // namespace adapterremoval
