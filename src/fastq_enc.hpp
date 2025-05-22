// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <array>   // for array
#include <cstddef> // for size_t
#include <string>  // for string

namespace adapterremoval {

//! Offset used by Phred scores in SAM files
const int PHRED_OFFSET_MIN = '!';
//! The maximum ASCII value allowed for encoded Phred scores in SAM files
const int PHRED_OFFSET_MAX = '~';
//! Minimum Phred score allowed for SAM files; encodes to '!'
const int PHRED_SCORE_MIN = 0;
//! Maximum Phred score allowed for SAM files; encodes to '~'
const int PHRED_SCORE_MAX = PHRED_OFFSET_MAX - PHRED_OFFSET_MIN;

//! Default character used to separate mate number
const char MATE_SEPARATOR = '/';

enum class quality_encoding
{
  phred_33,
  phred_64,
  solexa,
  sam,
};

//! How to handle degenerate bases (IUPAC notation)
enum class degenerate_encoding
{
  //! Replace any degenerate bases with 'N'
  mask,
  //! Report an error if any degenerate bases are observed
  reject,
};

//! How to handle uracil ('U')
enum class uracil_encoding
{
  //! Replace any degenerate bases with 'N'
  convert,
  //! Report an error if a uracil is observed
  reject,
};

class fastq_encoding
{
public:
  /**
   * Create FASTQ encoding with a given offset (33 or 64), allowing for
   * quality-scores up to a given value (0 - N). Input with higher scores
   * is rejected, and output is truncated to this score.
   */
  explicit fastq_encoding(
    quality_encoding encoding,
    degenerate_encoding degenerate = degenerate_encoding::reject,
    uracil_encoding uracils = uracil_encoding::reject) noexcept;

  /** Validates/normalizes a string of nucleotides in-place. */
  void process_nucleotides(std::string& sequence) const;

  /** Validates/converts a string of ASCII encoded PHRED values in-place */
  void process_qualities(std::string& qualities) const;

  /** Converts a Phred score (>= 0) to an error probability */
  static double phred_to_p(double phred);

  /** Converts a error probability (>= 0) to a Phred score */
  static double p_to_phred(double p);

  /** Converts a error probability (>= 0) to a Phred+33 encoded score */
  static char p_to_phred_33(double p);

private:
  //! Mask or reject degenerate bases
  bool m_mask_degenerate;
  //! Convert or reject uracil
  bool m_convert_uracil;
  //! Quality score encoding expected when decoding data
  quality_encoding m_encoding;
  //! Offset of the lowest ASCII value used by the given encoding
  char m_offset_min;
  //! Offset of the maximum ASCII value used by the given encoding
  char m_offset_max;
};

static const fastq_encoding FASTQ_ENCODING_33{ quality_encoding::phred_33 };
static const fastq_encoding FASTQ_ENCODING_64{ quality_encoding::phred_64 };
static const fastq_encoding FASTQ_ENCODING_SAM{ quality_encoding::sam };
static const fastq_encoding FASTQ_ENCODING_SOLEXA{ quality_encoding::solexa };

///////////////////////////////////////////////////////////////////////////////

struct ACGT
{
  using value_type = char;

  //! The number of possible index values
  static constexpr size_t indices = 4;
  //! Nucleotides supported by hashing function
  static constexpr std::array<value_type, 4> values{ 'A', 'C', 'G', 'T' };

  /**
   * Simple hashing function for nucleotides 'A', 'C', 'G', 'T', returning
   * numbers in the range 0-3. Passing characters other than "ACGT" (uppercase
   * only) will result in hash collisions.
   */
  static constexpr auto to_index(value_type nt) { return (nt >> 1) & 0x3; }

  /**
   * Inverse of to_index. Only values in the range 0 to 3 are allowed.
   */
  static constexpr value_type to_value(size_t idx) { return "ACTG"[idx]; }
};

struct ACGTN
{
  using value_type = char;

  //! The number of possible index values
  static constexpr size_t indices = 8;
  //! Nucleotides supported by hashing function
  static constexpr std::array<value_type, 5> values{ 'A', 'C', 'G', 'T', 'N' };

  /**
   * Simple hashing function for nucleotides 'A', 'C', 'G', 'T', 'N', returning
   * numbers in the range 0-4. Passing characters other than "ACGTN" (uppercase
   * only) will result in hash collisions.
   */
  static constexpr auto to_index(value_type nt) { return nt & 0x7; }

  /**
   * Inverse of to_index. Only values in the range 0 to 4 are allowed.
   */
  static constexpr value_type to_value(size_t idx) { return "-A-CT-NG"[idx]; }
};

} // namespace adapterremoval
