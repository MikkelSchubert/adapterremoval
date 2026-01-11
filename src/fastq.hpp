// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef>     // for size_t
#include <iosfwd>      // for ostream
#include <string>      // for string
#include <string_view> // for string_view
#include <utility>     // for pair
#include <vector>      // for vector

namespace adapterremoval {

class buffer;
class fastq_encoding;
class line_reader_base;

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
   *
   * Nucleotides are converted to uppercase, and dots are replaced with N.
   * Phred scores are converted to to Phred+33 scores, if not already Phred+33.
   *
   * The quality scores are expected to be in the range of 0 .. 40, unless
   * the format is Phred+33, in which case the range 0 .. 41 is accepted.
   */
  fastq(std::string_view header, std::string sequence, std::string qualities);

  /**  Create a new FASTQ record (as above) using the quality encoding */
  fastq(std::string_view header,
        std::string sequence,
        std::string qualities,
        const fastq_encoding& encoding);

  /**
   * Create FASTQ record from a sequence alone.
   *
   * @param header FASTQ header, including read name and meta information.
   * @param sequence nucleotide sequence containing the letters "acgtnACGTN."
   *
   * Works like the full constructor, except that qualities are all 0 ('!').
   */
  fastq(std::string_view header, std::string sequence);

  /** Returns true IFF all fields are identical. **/
  bool operator==(const fastq& other) const;

  /** Returns the header (including the @) of the record. **/
  const std::string& header() const { return m_header; }

  /** Returns the nucleotide sequence (ACGTN only) of the record. **/
  const std::string& sequence() const { return m_sequence; }

  /** Returns the Phred+33 encoded scores (0 .. 41) for each base. **/
  const std::string& qualities() const { return m_qualities; }

  /** Returns the name (excluding the @ and other fields) of the header. **/
  [[nodiscard]] std::string_view name(char mate_separator = '\0') const;

  /** Struct collecting basic information contained in a FASTQ header */
  struct fastq_header
  {
    //! Read name without '@' and mate number (if separator is specified)
    std::string_view name{};
    //! Mate number, if separator is specified
    char mate = '\0';
    //! Meta information in the header
    std::string_view meta{};
  };

  /** Returns the name (excluding the @ and other fields) of the header. **/
  [[nodiscard]] fastq_header parse_header(char mate_separator = '\0') const;

  /** Returns the length of the sequence. */
  size_t length() const { return m_sequence.length(); }

  /** Returns the number of ambiguous nucleotides in the sequence (N). **/
  size_t count_ns() const;

  /** Returns a measure of sequence complexity in the range [0; 1]. **/
  double complexity() const;

  /** Returns the mean Phred score of bases in a (non-empty) sequence */
  [[nodiscard]] double mean_quality() const;

  /**
   * Returns true if the read name corresponds to the output from known
   * two-color sequences technologies
   */
  [[nodiscard]] bool is_two_color() const;

  /** The number of bases trimmed from the 5p and 3p end respectively. **/
  using ntrimmed = std::pair<size_t, size_t>;

  /**
   * Trims consecutive low-quality bases from the 5'/3' ends of the sequence.
   *
   * @param trim_ns If true, ambiguous bases ('N') are trimmed.
   * @param low_quality Trim bases with a quality score at or below this value.
   * @param preserve5p Only trim from the 3p end if true.
   * @return A pair containing the number of 5' and 3' bases trimmed.
   */
  ntrimmed trim_trailing_bases(bool trim_ns = true,
                               char low_quality = -1,
                               bool preserve5p = false);

  /**
   * Trims low-quality bases using a sliding window approach.
   *
   * @param trim_ns If true, ambiguous bases ('N') are trimmed.
   * @param low_quality Trim bases with a quality score at or below this value.
   * @param window_size The length of the sliding window.
   * @param preserve5p Only trim from the 3p end if true.
   * @return A pair containing the number of 5' and 3' bases trimmed.
   */
  ntrimmed trim_windowed_bases(bool trim_ns = true,
                               char low_quality = -1,
                               double window_size = 0.1,
                               bool preserve5p = false);

  /**
   * Performs quality based trimming using the modified Mott's algorithm.
   *
   * @param error_limit Error rate limit used for good v. bad quality scores.
   * @param preserve5p Only trim from the 3p end if true.
   * @return A pair containing the number of 5' and 3' bases trimmed.
   */
  ntrimmed mott_trimming(double error_limit, bool preserve5p = false);

  /**
   * Trims the longest poly-X tail, where X is one of the specified nucleotides.
   *
   * @param nucleotides Nucleotides to consider ('A', 'G', 'C', and/or 'T').
   * @param min_length The minimum length of of a poly-G tail.
   * @return The number of 3' bases trimmed for the best matching nucleotide.
   */
  std::pair<char, size_t> poly_x_trimming(std::string_view nucleotides,
                                          size_t min_length);

  /**
   * Truncates the record in place.
   *
   * This function behaves like std::string::substr, except that the
   * substrings (sequence / qualities) are re-assigned to the record itself.
   */
  void truncate(size_t pos = 0, size_t len = std::string::npos);

  /** Reverse complements the record in place. */
  void reverse_complement();

  /** Adds a prefix to the name. */
  void add_prefix_to_name(std::string_view prefix);

  /**
   * Reads a FASTQ record from a list of lines (without newlines).
   *
   * If a malformed or invalid FASTQ record is encountered, the fastq_error
   * exception is raised. Note that 'this' record is only valid if read
   * returned true. Unlike the constructor, this function does not accept
   * empty headers, or sequences / qualities, as this typically indicates
   * a problem with the source file.
   */
  bool read(line_reader_base& reader, const fastq_encoding& encoding);

  /** Like `read`, but post-processing must be manually called afterwards */
  bool read_unsafe(line_reader_base& reader);

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
   * Attempt to infer the mate separator from a single-end reads.
   *
   * A separator is considered valid if the reads contain this separator, and
   * all reads contain a valid mate number (1 or 2). Unsorted and interleaved
   * reads are allowed.
   *
   * Returns 0 if the separator could not be guessed.
   */
  static char guess_mate_separator(const std::vector<fastq>& reads);

  /**
   * Validates a pair and normalizes the mate separator.
   *
   * The mate separator character is the character expected as the second-to-
   * last character, if the last character (either '1' or '2') specify the
   * mate number (must be 1 and 2, in that order).
   **/
  static void validate_paired_reads(const fastq& mate1,
                                    const fastq& mate2,
                                    char mate_separator);

  /**
   * Finalizes read, validates sequence and transforms qualities. This function
   * *must* be called for all reads produced by calling `read_unsafe`.
   */
  void post_process(const fastq_encoding& encoding);

private:
  /**
   * Trims the read to the specified bases, and returns a pair specifying the
   * number of 5' and 3' bases removed.
   */
  ntrimmed trim_sequence_and_qualities(size_t left_inclusive,
                                       size_t right_exclusive);

  //! Header including the @ sigil, but (possibly) including meta-info
  std::string m_header;
  //! Nucleotide sequence; contains only uppercase letters "ACGTN"
  std::string m_sequence;
  //! Phred+33 encoded quality scores
  std::string m_qualities;

  //! Needs access to merge sequence/qualities to in-place
  friend class sequence_merger;
};

using fastq_pair = std::pair<fastq, fastq>;
using fastq_pair_vec = std::vector<fastq_pair>;

/** Stream operator for debugging output */
std::ostream&
operator<<(std::ostream& os, const fastq& value);

} // namespace adapterremoval
