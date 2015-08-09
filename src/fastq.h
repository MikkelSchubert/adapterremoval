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
#ifndef FASTQ_H
#define FASTQ_H

#include <string>

#include "commontypes.h"


const int MIN_PHRED_SCORE =  0;
const int MAX_PHRED_SCORE = 41;
const int MIN_SOLEXA_SCORE = -5;
const int MAX_SOLEXA_SCORE = 41;

// Maximum Phred score allowed by the BAM format
const int MAX_EXTENDED_PHRED_SCORE = '~' - '!';

const int PHRED_OFFSET_33 = 33;
const int PHRED_OFFSET_64 = 64;


class fastq_error : public std::exception
{
public:
    fastq_error(const std::string& message);
    ~fastq_error() throw();

    virtual const char* what() const throw();

private:
    std::string m_message;
};


/**
 * Represents a FASTQ record with Phred (offset=33) encoded quality scores.
 */
class fastq
{
public:
    enum quality_format
    {
        //! Phred scores with offset = 33
        phred_33 = 0,
        //! Phred scores with offset = 64
        phred_64,
        //! Solexa scores (offset = 64); lossily converted to Phred
        solexa,
        //! Quality scores are ignored (set to 0)
        ignored
    };


    /** Constructs a dummy FASTQ record for which all fields are empty. **/
    fastq();

    /**
     *
     * @param qname FASTQ header, including read name and meta information.
     * @param sequence nucleotide sequence containing the letters "acgtnACGTN."
     * @param qualities phred or solexa encoded quality scores
     * @param encoding the encoding used for the quality scores.
     *
     * Nucleotides are converted to uppercase, and dots are replaced with N.
     * Phred scores are calculated from the encoded scores, in the case of
     * Solexa encoding re-encoded under this scheme, otherwise simply decoded
     * to Phred+33 scores, if not already in this format.
     *
     * The quality scores are expected to be in the range of 0 .. 40, unless
     * the format is Phred+33, in which case the range 0 .. 41 is accepted.
     */
    fastq(const std::string& header,
          const std::string& sequence,
          const std::string& qualities,
          quality_format encoding = phred_33);


    /** Returns true IFF all fields are identical. **/
    bool operator==(const fastq& other) const;


    /** Returns the header (excluding the @) of the record. **/
    const std::string& header() const;
    /** Returns the nucleotide sequence (ACGTN only) of the record. **/
    const std::string& sequence() const;
    /** Returns the Phred+33 encoded scores (0 .. 41) for each base. **/
    const std::string& qualities() const;

    /** Returns the length of the sequence. */
    size_t length() const;

    /** Returns the number of ambiguous nucleotides in the sequence (N). **/
    size_t count_ns() const;


    /** The number of bases trimmmed from the 5p and 3p end respectively. **/
    typedef std::pair<size_t, size_t> ntrimmed;

    /**
     * Trims consequtive low-quality bases from the 5'/3' ends of the sequence.
     *
     * @param trim_ns If true, ambiguous bases ('N') are trimmed.
     * @param low_quality Trim bases with a quality score at or below this value.
     * @return A pair containing hte number of bases trimmed from either end.
     */
	ntrimmed trim_low_quality_bases(bool trim_ns = true, char low_quality = -1);

    /**
     * Truncates the record in place.
     *
     * This function behaves like std::string::substr, except that the
     * substrings (sequence / qualities) are re-assigned to the record itself.
     */
    void truncate(size_t pos = 0, size_t len = std::string::npos);

    /** Reverse complements the record in place. */
    void reverse_complement();

    /** Adds a prefix to the header. */
    void add_prefix_to_header(const std::string& prefix);


	/**
	 * Reads a FASTQ record from a list of lines (without newlines).
     *
	 * If a malformed or invalid FASTQ record is encountered, the fastq_error
     * exception is raised. Note that 'this' record is only valid if read
     * returned true. Unlike the constructor, this function does not accept
     * empty headers, or sequences / qualities, as this typically indicates
     * a problem with the source file.
	 */
    bool read(string_vec_citer& begin, const string_vec_citer& end, quality_format encoding = phred_33);

    /**
     * Converts a FASTQ record to a string ending with a newline.
     *
     * Only the phred_33 and phred_64 encodings are supported. For phred_64,
     * quality bases are truncated to 0 .. 40, while phred_33 supports quality
     * scores in the range 0 .. 41.
     */
    std::string to_str(quality_format encoding = phred_33) const;


    /**
     * Converting lower-case nucleotides to uppercase, '.' to N.
     *
     * If the sequence contains letters other than "acgtnACGTN.", a fastq_error
     * is thrown.
     **/
    static void clean_sequence(std::string& sequence);

    /** Converts an error-probability to a Phred+33 encoded quality score. **/
    static char p_to_phred_33(double p);

    /** Validate that two reads form a valid pair. */
    static void validate_paired_reads(const fastq& mate1, const fastq& mate2);

private:
    /** Initializes record; used by constructor and read function. **/
    void process_record(quality_format encoding);

    //! Header excluding the @ sigil, but (possibly) including meta-info
	std::string m_header;
    //! Nucleotide sequence; contains only uppercase letters "ACGTN"
	std::string m_sequence;
    //! Phred+33 encoded quality scores
	std::string m_qualities;
};




///////////////////////////////////////////////////////////////////////////////
inline const std::string& fastq::header() const
{
    return m_header;
}


inline const std::string& fastq::sequence() const
{
    return m_sequence;
}


inline const std::string& fastq::qualities() const
{
    return m_qualities;
}

#endif
