/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#pragma once

#include <exception> // for exception
#include <stddef.h>  // for size_t
#include <string>    // for string

namespace adapterremoval {

//! Offset used by Phred+33 and SAM encodings
const int PHRED_OFFSET_33 = '!';
//! Offset used by Phred+64 and Solexa encodings
const int PHRED_OFFSET_64 = '@';

//! Minimum Phred score allowed; encodes to '!'
const int MIN_PHRED_SCORE = 0;
//! Maximum Phred score allowed by default, to ensure backwards compatibility
//! with AdapterRemoval v1.x.
const int MAX_PHRED_SCORE_DEFAULT = 41;
//! Maximum Phred score allowed, as this encodes to the last printable
//! character '~', when using an offset of 33.
const int MAX_PHRED_SCORE = '~' - '!';

//! Minimum Solexa score allowed; encodes to ';' with an offset of 64
const int MIN_SOLEXA_SCORE = -5;
//! Maximum Solexa score allowed; encodes to 'h' with an offset of 64
const int MAX_SOLEXA_SCORE = 40;

//! Default character used to separate mate number
const char MATE_SEPARATOR = '/';

enum class quality_encoding
{
  solexa = -1,
  phred_33 = 33,
  phred_64 = 64,
};

/** Exception raised for FASTQ parsing and validation errors. */
class fastq_error : public std::exception
{
public:
  fastq_error(const std::string& message);
  fastq_error(const fastq_error& error);

  virtual ~fastq_error() override;

  /** Returns error message; string is owned by exception. */
  virtual const char* what() const noexcept override;

private:
  //! Error message associated with exception.
  std::string m_message;
};

class fastq_encoding
{
public:
  /**
   * Create FASTQ encoding with a given offset (33 or 64), allowing for
   * quality-scores up to a given value (0 - N). Input with higher scores
   * is rejected, and output is truncated to this score.
   */
  fastq_encoding(quality_encoding encoding, char max_score);

  /** Appends Phred+33 encoded qualities to dst. */
  static void encode(const std::string& qualities, std::string& dst);
  /** Decodes a string of ASCII values in-place. */
  void decode(std::string& qualities) const;

protected:
  //! Quality score encoding expected when decoding data
  quality_encoding m_encoding;
  //! Offset used by the given encoding
  char m_offset;
  //! Maximum allowed score; used for checking input / truncating output
  char m_max_score;
};

static const fastq_encoding FASTQ_ENCODING_33(quality_encoding::phred_33,
                                              MAX_PHRED_SCORE_DEFAULT);
static const fastq_encoding FASTQ_ENCODING_64(quality_encoding::phred_64,
                                              MAX_PHRED_SCORE_DEFAULT);
static const fastq_encoding FASTQ_ENCODING_SAM(quality_encoding::phred_33,
                                               MAX_PHRED_SCORE);
static const fastq_encoding FASTQ_ENCODING_SOLEXA(quality_encoding::solexa,
                                                  MAX_SOLEXA_SCORE);

} // namespace adapterremoval
