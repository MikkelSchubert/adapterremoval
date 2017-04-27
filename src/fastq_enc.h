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
#ifndef FASTQ_ENC_H
#define FASTQ_ENC_H

#include <string>

namespace ar
{

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
const int MAX_PHRED_SCORE = 93;


//! Minimum Solexa score allowed; encodes to ';' with an offset of 64
const int MIN_SOLEXA_SCORE = -5;
//! Maximum Solexa score allowed; encodes to 'h' with an offset of 64
const int MAX_SOLEXA_SCORE = 40;

//! Default character used to separate mate number
const char MATE_SEPARATOR = '/';


/** Exception raised for FASTQ parsing and validation errors. */
class fastq_error : public std::exception
{
public:
    fastq_error(const std::string& message);
    fastq_error(const fastq_error& error);

    virtual ~fastq_error() noexcept;

    /** Returns error message; string is owned by exception. */
    virtual const char* what() const noexcept;

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
    fastq_encoding(char offset = PHRED_OFFSET_33,
                   char max_score = MAX_PHRED_SCORE_DEFAULT);

    virtual ~fastq_encoding();

    /** Encodes a string of Phred+33/66 quality-scores in-place. */
    virtual void encode_string(std::string::iterator it, const std::string::iterator& end) const;
    /** Decodes a string of ASCII values in-place. */
    virtual void decode_string(std::string::iterator it, const std::string::iterator& end) const;

    /** Returns the standard name for this encoding. */
    virtual std::string name() const;

    /**
     * Returns the maximum allowed quality score for input, and the range to
     * range to which output scores are truncated.
     */
    virtual size_t max_score() const;

protected:
    //! Character offset for Phred encoded scores (33 or 64)
    const char m_offset;
    //! Maximum allowed score; used for checking input / truncating output
    const char m_max_score;

private:
    //! Not implemented
    fastq_encoding(const fastq_encoding&);
    //! Not implemented
    fastq_encoding& operator=(const fastq_encoding&);
};


/**
 * Solexa scores encoding by adding '@'; max score is 40.
 *
 * Solexa scores are defined as Q = -10 * log10(p / (1 - p)), and differ from
 * Phred scores for values less than 13. Lossless conversion between the
 * formats is not possible, and since the fastq class stores quality scores as
 * Phred+33 internally, this means that reading Solexa scores is a lossy
 * operation, even if the output is written using Solexa scores.
 */
class fastq_encoding_solexa : public fastq_encoding
{
public:
    /**
     * Create FASTQ Solexa encoding with offset 64, allowing for quality-scores
     * up to a given value (0 - N). Input with higher scores is rejected, and
     * output is truncated to this score.
     */
    fastq_encoding_solexa(unsigned max_score = MAX_PHRED_SCORE_DEFAULT);

    /** Encodes a string of Phred+33/66 quality-scores in-place. */
    virtual void encode_string(std::string::iterator it, const std::string::iterator& end) const;
    /** Decodes a string of ASCII values in-place. */
    virtual void decode_string(std::string::iterator it, const std::string::iterator& end) const;

    /** Returns the standard name for this encoding. */
    std::string name() const;
};


static const fastq_encoding FASTQ_ENCODING_33(PHRED_OFFSET_33);
static const fastq_encoding FASTQ_ENCODING_64(PHRED_OFFSET_64);
static const fastq_encoding FASTQ_ENCODING_SAM(PHRED_OFFSET_33, MAX_PHRED_SCORE);
static const fastq_encoding_solexa FASTQ_ENCODING_SOLEXA;

} // namespace ar

#endif
