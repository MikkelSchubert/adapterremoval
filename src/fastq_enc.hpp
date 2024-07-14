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
#pragma once

#include <stdexcept> // for runtime_error
#include <string>    // for string

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

class fastq_encoding
{
public:
  /**
   * Create FASTQ encoding with a given offset (33 or 64), allowing for
   * quality-scores up to a given value (0 - N). Input with higher scores
   * is rejected, and output is truncated to this score.
   */
  explicit fastq_encoding(quality_encoding encoding) noexcept;

  /** Decodes a string of ASCII values in-place. */
  void decode(std::string& qualities) const;

private:
  //! Quality score encoding expected when decoding data
  quality_encoding m_encoding;
  //! Offset of the lowest ASCII value used by the given encoding
  char m_offset_min;
  //! Offset of the maximum ASCII value used by the given encoding
  char m_offset_max;
};

static const fastq_encoding FASTQ_ENCODING_33(quality_encoding::phred_33);
static const fastq_encoding FASTQ_ENCODING_64(quality_encoding::phred_64);
static const fastq_encoding FASTQ_ENCODING_SAM(quality_encoding::sam);
static const fastq_encoding FASTQ_ENCODING_SOLEXA(quality_encoding::solexa);

} // namespace adapterremoval
