// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "fastq_enc.hpp" // header
#include "debug.hpp"     // for AR_FAIL, AR_REQUIRE
#include "errors.hpp"    // for fastq_error
#include "strutils.hpp"  // for shell_escape
#include <algorithm>     // for min, max
#include <array>         // for array
#include <cmath>         // for log10, pow, round
#include <sstream>       // for ostringstream
#include <string>        // for string

namespace adapterremoval {

namespace {

//! Offset used by Phred+33 and SAM encodings
constexpr int PHRED_33_OFFSET_MIN = '!';
//! The maximum ASCII value allowed for Phred+33 encoded quality scores. This is
//! a generous maximum, taking into account that instruments are assigning
//! higher and higher quality scores for Phred+33 data.
constexpr int PHRED_33_OFFSET_MAX = 'N';
//! Minimum Phred score allowed
// constexpr int PHRED_33_SCORE_MIN = 0;
//! Maximum Phred score allowed
constexpr int PHRED_33_SCORE_MAX = PHRED_33_OFFSET_MAX - PHRED_33_OFFSET_MIN;

//! Offset used by Phred+64 encodings
constexpr int PHRED_64_OFFSET_MIN = '@';
//! The maximum ASCII value allowed for encoded Phred scores
constexpr int PHRED_64_OFFSET_MAX = '~';
//! Minimum Phred+64 score allowed
constexpr int PHRED_64_SCORE_MIN = 0;
//! Maximum Phred+64 score allowed
// constexpr int PHRED_64_SCORE_MAX = PHRED_64_OFFSET_MAX - PHRED_64_OFFSET_MIN;

//! Offset used by Solexa encoding quality scores
constexpr int SOLEXA_OFFSET_MIN = '@';
//! The maximum ASCII value allowed for encoded Solexa scores
constexpr int SOLEXA_OFFSET_MAX = 'h';
//! Minimum Phred encoded score allowed
constexpr int SOLEXA_SCORE_MIN = -5;
//! Maximum Phred encoded score allowed
constexpr int SOLEXA_SCORE_MAX = SOLEXA_OFFSET_MAX - SOLEXA_OFFSET_MIN;

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
    const auto idx = static_cast<unsigned char>(i + SOLEXA_OFFSET_MIN);

    scores.at(idx) = PHRED_OFFSET_MIN + transformed;
  }

  return scores;
}

const auto g_solexa_to_phred33 = calc_solexa_to_phred33();

///////////////////////////////////////////////////////////////////////////////

/**
 * Uppercase letters in the range a-z, but mangles other characters. Resulting
 * values outside of A-Z should be compared against/reported, since they may not
 * reflect the original value.
 */
constexpr char
mangle_to_upper(const char c)
{
  return c & 0xDF;
}

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
  AR_REQUIRE(raw < PHRED_64_OFFSET_MIN,
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

    case quality_encoding::sam:
      AR_FAIL("This case should have been handled by initial check");

    default:
      AR_FAIL("Invalid quality encoding");
  }
}

[[noreturn]] void
throw_invalid_base(char c)
{
  std::ostringstream stream;

  switch (mangle_to_upper(c)) {
    case 'B': // C / G / T
    case 'D': // A / G / T
    case 'H': // A / C / T
    case 'K': // G / T
    case 'M': // A / C
    case 'R': // A / G
    case 'S': // C / G
    case 'V': // A / C / G
    case 'W': // A / T
    case 'Y': // C / T
      stream << "found degenerate base '" << c << "' in FASTQ sequence, but "
             << "only bases A, C, G, T and N are supported. Use the option "
                "--mask-degenerate-bases to convert degenerate bases to N";
      break;

    case 'U': // Uracils
      stream << "found uracil (U) in FASTQ sequence, but only bases A, C, G, "
                "T and N are supported. Use the option --convert-uracils to "
                "convert uracils (U) to thymine (T)";
      break;

    default:
      stream << "invalid character " << log_escape(std::string(1, c))
             << " found in FASTQ sequence";
      break;
  }

  throw fastq_error(stream.str());
}

// Build once at startup / compile time
constexpr std::array<char, 256>
build_lookup_table(const bool degenerate, const bool uracil)
{
  std::array<char, 256> lut{}; // zero-initialize: 0 = invalid

  const std::array plain_nucleotides = { 'A', 'C', 'G', 'T', 'N' };
  for (const auto c : plain_nucleotides) {
    lut[static_cast<unsigned char>(c)] = c;
    lut[static_cast<unsigned char>(c) ^ 0x20] = c; // from lowercase
  }

  lut[static_cast<unsigned char>('.')] = 'N'; // legacy data

  if (uracil) {
    lut[static_cast<unsigned char>('u')] = 'T';
    lut[static_cast<unsigned char>('U')] = 'T';
  }

  if (degenerate) {
    // https://en.wikipedia.org/wiki/Nucleic_acid_notation
    const std::array degenerate_nucleotides = { 'B', 'D', 'H', 'K', 'M',
                                                'R', 'S', 'V', 'W', 'Y' };
    for (const auto c : degenerate_nucleotides) {
      lut[static_cast<unsigned char>(c)] = 'N';
      lut[static_cast<unsigned char>(c) ^ 0x20] = 'N'; // from lowercase
    }
  }

  return lut;
}

constexpr auto NUCLEOTIDES_LUT_PLAIN = build_lookup_table(false, false);
constexpr auto NUCLEOTIDES_LUT_URACILS = build_lookup_table(false, true);
constexpr auto NUCLEOTIDES_LUT_DEGENERATE = build_lookup_table(true, false);
constexpr auto NUCLEOTIDES_LUT_BOTH = build_lookup_table(true, true);

} // namespace

///////////////////////////////////////////////////////////////////////////////

fastq_encoding::fastq_encoding(quality_encoding encoding,
                               degenerate_encoding degenerate,
                               uracil_encoding uracils) noexcept
  : m_encoding(encoding)

{
  if (degenerate == degenerate_encoding::reject &&
      uracils == uracil_encoding::reject) {
    m_nucleotides_lut = NUCLEOTIDES_LUT_PLAIN.data();
  } else if (degenerate == degenerate_encoding::reject &&
             uracils == uracil_encoding::convert) {
    m_nucleotides_lut = NUCLEOTIDES_LUT_URACILS.data();
  } else if (degenerate == degenerate_encoding::mask &&
             uracils == uracil_encoding::reject) {
    m_nucleotides_lut = NUCLEOTIDES_LUT_DEGENERATE.data();
  } else if (degenerate == degenerate_encoding::mask &&
             uracils == uracil_encoding::convert) {
    m_nucleotides_lut = NUCLEOTIDES_LUT_BOTH.data();
  } else {
    AR_FAIL("invalid degenerate or uracil encoding value");
  }

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
      AR_FAIL("unknown quality encoding");
  }
}

void
fastq_encoding::process_nucleotides(std::string& sequence) const
{
  for (char& nuc : sequence) {
    const char original = nuc;
    nuc = m_nucleotides_lut[static_cast<unsigned char>(nuc)];
    if (AR_UNLIKELY(nuc == '\0')) {
      throw_invalid_base(original);
    }
  }
}

void
fastq_encoding::process_qualities(std::string& qualities) const
{
  if (m_encoding == quality_encoding::solexa) {
    for (auto& quality : qualities) {
      const char current = quality;
      const auto idx = static_cast<unsigned char>(quality);

      quality = g_solexa_to_phred33.at(idx);
      if (quality == '\0') {
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

double
fastq_encoding::phred_to_p(double phred)
{
  AR_REQUIRE(phred >= 0.0);
  return std::pow(10.0, phred / -10.0);
}

double
fastq_encoding::p_to_phred(double p)
{
  AR_REQUIRE(p >= 0.0 && p <= 1.0);
  // max(0.0, ...) to avoid returning -0.0 for p = 1
  return std::max(0.0, -10.0 * std::log10(p));
}

char
fastq_encoding::p_to_phred_33(double p)
{
  AR_REQUIRE(p >= 0.0 && p <= 1.0);

  // Lowest error rate that can be represented is 93 (~5e-10), encoded as '~'
  const auto raw_score = p_to_phred(std::max(5e-10, p));
  return std::min<int>(PHRED_OFFSET_MAX,
                       PHRED_OFFSET_MIN +
                         static_cast<int>(std::round(raw_score)));
}

} // namespace adapterremoval
