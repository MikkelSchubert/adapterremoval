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
#pragma once

#include <stddef.h> // for size_t
#include <string>   // for string
#include <utility>  // for pair
#include <vector>   // for vector

#include "fastq_enc.hpp" // for FASTQ_ENCODING_33, MATE_SEPARATOR

class line_reader_base;
struct mate_info;

/**
 * Represents a FASTQ record with Phred (offset=33) encoded quality scores.
 */
class fastq
{
public:
  /** Constructs a dummy FASTQ record for which all fields are empty. **/
  fastq();

  /**
   * Create a new FASTQ record.
   *
   * @param header FASTQ header, including read name and meta information.
   * @param sequence nucleotide sequence containing the letters "acgtnACGTN."
   * @param qualities phred encoded quality scores
   * @param encoding the encoding used for the quality scores.
   *
   * Nucleotides are converted to uppercase, and dots are replaced with N.
   * Phred scores are converted to to Phred+33 scores, if not already Phred+33.
   *
   * The quality scores are expected to be in the range of 0 .. 40, unless
   * the format is Phred+33, in which case the range 0 .. 41 is accepted.
   */
  fastq(const std::string& header,
        const std::string& sequence,
        const std::string& qualities,
        const fastq_encoding& encoding = FASTQ_ENCODING_33);

  /**
   * Create FASTQ record from a sequence alone.
   *
   * @param header FASTQ header, including read name and meta information.
   * @param sequence nucleotide sequence containing the letters "acgtnACGTN."
   *
   * Works like the full constructor, except that qualities are all 0 ('!').
   */
  fastq(const std::string& header, const std::string& sequence);

  /** Returns true IFF all fields are identical. **/
  bool operator==(const fastq& other) const;

  /** Returns the header (excluding the @) of the record. **/
  const std::string& header() const;
  /** Returns the nucleotide sequence (ACGTN only) of the record. **/
  const std::string& sequence() const;
  /** Returns the Phred+33 encoded scores (0 .. 41) for each base. **/
  const std::string& qualities() const;

  /** Returns the name (excluding the @ and other fields) of the header. **/
  std::string name() const;

  /** Returns the length of the sequence. */
  size_t length() const;

  /** Returns the number of ambiguous nucleotides in the sequence (N). **/
  size_t count_ns() const;

  /** Returns a measure of sequence complexity in the range [0; 1]. **/
  double complexity() const;

  /** The number of bases trimmmed from the 5p and 3p end respectively. **/
  typedef std::pair<size_t, size_t> ntrimmed;

  /**
   * Trims consecutive low-quality bases from the 5'/3' ends of the sequence.
   *
   * @param trim_ns If true, ambiguous bases ('N') are trimmed.
   * @param low_quality Trim bases with a quality score at or below this value.
   * @param preserve5p Only trim from the 3p end if true.
   * @return A pair containing the number of 5' and 3' bases trimmed.
   */
  ntrimmed trim_trailing_bases(const bool trim_ns = true,
                               char low_quality = -1,
                               const bool preserve5p = false);

  /**
   * Trims low-quality bases using a sliding window approach.
   *
   * @param trim_ns If true, ambiguous bases ('N') are trimmed.
   * @param low_quality Trim bases with a quality score at or below this value.
   * @param winlen The length of the sliding window.
   * @param preserve5p Only trim from the 3p end if true.
   * @return A pair containing the number of 5' and 3' bases trimmed.
   */
  ntrimmed trim_windowed_bases(const bool trim_ns = true,
                               char low_quality = -1,
                               const double window_size = 0.1,
                               const bool preserve5p = false);

  /**
   * Truncates the record in place.
   *
   * This function behaves like std::string::substr, except that the
   * substrings (sequence / qualities) are re-assigned to the record itself.
   */
  void truncate(size_t pos = 0, size_t len = std::string::npos);

  /** Reverse complements the record in place. */
  void reverse_complement();

  /**
   * Reads a FASTQ record from a list of lines (without newlines).
   *
   * If a malformed or invalid FASTQ record is encountered, the fastq_error
   * exception is raised. Note that 'this' record is only valid if read
   * returned true. Unlike the constructor, this function does not accept
   * empty headers, or sequences / qualities, as this typically indicates
   * a problem with the source file.
   */
  bool read(line_reader_base& reader,
            const fastq_encoding& encoding = FASTQ_ENCODING_33);

  /**
   * Reads a FASTQ record from a list of lines (without newlines).
   *
   * This is equivalent to `read` except that only minimal validation is
   * performed. `post_process` *must* be called before using bases or quality
   * scores.
   */
  bool read_unsafe(line_reader_base& reader);

  /**
   * Finalizes read, validates sequence and transforms qualities. This function
   * *must* be called for all reads produced by calling `read_unsafe`.
   */
  void post_process(const fastq_encoding& encoding);

  /**
   * Converts a FASTQ record to a string ending with a newline.
   *
   * Only the phred_33 and phred_64 encodings are supported. For phred_64,
   * quality bases are truncated to 0 .. 40, while phred_33 supports quality
   * scores in the range 0 .. 41.
   */
  void into_string(std::string& dst) const;

  /** Converts an error-probability to a Phred+33 encoded quality score. **/
  static char p_to_phred_33(double p);

  /**
   * Attempt to infer the mate separator from a set of paired reads.
   *
   * A separator is considered valid if the reads contain this separator and
   * yield mate 1 for `reads_1` and mate 2 reads for `reads_2`. Reads without
   * mate numbers will be assumed to use '/'. Currently considers '/', '.',
   * and ':' as possible candidate separators.
   *
   * Returns 0 if the separator could not be guessed.
   */
  static char guess_mate_separator(const std::vector<fastq>& reads_1,
                                   const std::vector<fastq>& reads_2);

  /**
   * Validates a pair and normalizes the mate separator.
   *
   * The mate separator character is the character expected as the second-to-
   * last character, if the last character (either '1' or '2') specify the
   * mate number (must be 1 and 2, in that order).
   **/
  static void normalize_paired_reads(fastq& mate1,
                                     fastq& mate2,
                                     char mate_separator = MATE_SEPARATOR);

private:
  /**
   * Converting lower-case nucleotides to uppercase.
   *
   * If the sequence contains letters other than "acgtnACGTN.", a fastq_error
   * is thrown.
   **/
  static void clean_sequence(std::string& sequence);

  /**
   * Trims the read to the specified bases, and returns a pair specifying the
   * number of 5' and 3' bases removed.
   */
  ntrimmed trim_sequence_and_qualities(const size_t left_inclusive,
                                       const size_t right_exclusive);

  //! Header excluding the @ sigil, but (possibly) including meta-info
  std::string m_header;
  //! Nucleotide sequence; contains only uppercase letters "ACGTN"
  std::string m_sequence;
  //! Phred+33 encoded quality scores
  std::string m_qualities;

  //! Needs access to merge sequence/qualities to in-place
  friend class sequence_merger;
};

///////////////////////////////////////////////////////////////////////////////

/**
 * Simple hashing function for nucleotides 'A', 'C', 'G', 'T', returning
 * numbers in the range 0-3. Passing characters other than "ACGT" (uppercase
 * only) will result in hash collisions.
 */
inline size_t
ACGT_TO_IDX(char nt)
{
  return (nt >> 1) & 0x3;
}

/**
 * Inverse of ACGT_TO_IDX. Only values in the range 0 to 3 are allowed.
 */
inline char
IDX_TO_ACGT(size_t idx)
{
  return "ACTG"[idx];
}

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<fastq, fastq> fastq_pair;
typedef std::vector<fastq_pair> fastq_pair_vec;

///////////////////////////////////////////////////////////////////////////////
inline const std::string&
fastq::header() const
{
  return m_header;
}

inline std::string
fastq::name() const
{
  const size_t pos = m_header.find_first_of(' ');
  if (pos != std::string::npos) {
    return m_header.substr(0, pos);
  }

  return m_header;
}

inline const std::string&
fastq::sequence() const
{
  return m_sequence;
}

inline const std::string&
fastq::qualities() const
{
  return m_qualities;
}

inline size_t
fastq::length() const
{
  return m_sequence.length();
}
