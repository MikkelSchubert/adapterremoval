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

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include "alignment.hpp"
#include "alignment_tables.hpp"
#include "debug.hpp"
#include "fastq.hpp"

#if defined(__SSE__) && defined(__SSE2__)
#include <bitset>
#include <emmintrin.h>

//! Mask of all Ns
const __m128i N_MASK_128 = _mm_set1_epi8('N');

/** Counts the number of masked bytes **/
inline size_t
COUNT_MASKED(__m128i value)
{
  // Generate 16 bit mask from most significant bits of each byte and count bits
  return std::bitset<16>(_mm_movemask_epi8(value)).count();
}
#else
#warning SEE optimizations disabled!
#endif

/**
 * Compares two subsequences in an alignment to a previous (best) alignment.
 *
 * @param best The currently best alignment, used for evaluating this alignment
 * @param current The current alignment to be evaluated (counts are assumed to
 * be zero'd!)
 * @param seq_1_ptr Pointer to the first base in the first sequence in the
 * alignment.
 * @param seq_2_ptr Pointer to the first base in the second sequence in the
 * alignment.
 * @return True if the current alignment is at least as good as the best
 * alignment, false otherwise.
 *
 * If the function returns false, the current alignment cannot be assumed to
 * have been completely evaluated (due to early termination), and hence counts
 * and scores are not reliable. The function assumes uppercase nucleotides.
 */
bool
compare_subsequences(const alignment_info& best,
                     alignment_info& current,
                     const char* seq_1_ptr,
                     const char* seq_2_ptr)
{
  int remaining_bases = current.score = current.length;

#if defined(__SSE__) && defined(__SSE2__)
  while (remaining_bases >= 16 && current.score >= best.score) {
    const __m128i s1 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_1_ptr));
    const __m128i s2 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_2_ptr));

    // Sets 0xFF for every byte where one or both nts is N
    const __m128i ns_mask = _mm_or_si128(_mm_cmpeq_epi8(s1, N_MASK_128),
                                         _mm_cmpeq_epi8(s2, N_MASK_128));

    // Sets 0xFF for every byte where bytes are equal or N
    const __m128i eq_mask = _mm_or_si128(_mm_cmpeq_epi8(s1, s2), ns_mask);

    current.n_ambiguous += COUNT_MASKED(ns_mask);
    current.n_mismatches += 16 - COUNT_MASKED(eq_mask);

    // Matches count for 1, Ns for 0, and mismatches for -1
    current.score =
      current.length - current.n_ambiguous - (current.n_mismatches * 2);

    seq_1_ptr += 16;
    seq_2_ptr += 16;
    remaining_bases -= 16;
  }
#endif

  for (; remaining_bases && current.score >= best.score; --remaining_bases) {
    const char nt_1 = *seq_1_ptr++;
    const char nt_2 = *seq_2_ptr++;

    if (nt_1 == 'N' || nt_2 == 'N') {
      current.n_ambiguous++;
      current.score--;
    } else if (nt_1 != nt_2) {
      current.n_mismatches++;
      current.score -= 2;
    }
  }

  return current.is_better_than(best);
}

alignment_info
pairwise_align_sequences(const alignment_info& best_alignment,
                         const std::string& seq1,
                         const std::string& seq2,
                         int min_offset = std::numeric_limits<int>::min(),
                         int max_offset = std::numeric_limits<int>::max())
{
  const int start_offset =
    std::max<int>(min_offset, -static_cast<int>(seq2.length()) + 1);
  const int end_offset =
    std::min<int>(max_offset, static_cast<int>(seq1.length()) - 1);

  alignment_info best = best_alignment;
  for (int offset = start_offset; offset <= end_offset; ++offset) {
    const size_t initial_seq1_offset = std::max<int>(0, offset);
    const size_t initial_seq2_offset = std::max<int>(0, -offset);
    const size_t length = std::min(seq1.length() - initial_seq1_offset,
                                   seq2.length() - initial_seq2_offset);

    if (static_cast<int>(length) >= best.score) {
      alignment_info current;
      current.offset = offset;
      current.length = length;

      const char* seq_1_ptr = seq1.data() + initial_seq1_offset;
      const char* seq_2_ptr = seq2.data() + initial_seq2_offset;

      if (compare_subsequences(best, current, seq_1_ptr, seq_2_ptr)) {
        best = current;
      }
    }
  }

  return best;
}

struct phred_scores
{
  explicit phred_scores(size_t index)
    : identical_nts(IDENTICAL_NTS[index])
    , different_nts(DIFFERENT_NTS[index])
  {}

  //! Phred score to assign if the two nucleotides are identical
  char identical_nts;
  //! Phred score to assign if the two nucleotides differ
  char different_nts;
};

phred_scores
get_updated_phred_scores(char qual_1, char qual_2)
{
  AR_DEBUG_ASSERT(qual_1 >= qual_2);

  const size_t phred_1 = static_cast<size_t>(qual_1 - PHRED_OFFSET_33);
  const size_t phred_2 = static_cast<size_t>(qual_2 - PHRED_OFFSET_33);
  const size_t index = (phred_1 * (MAX_PHRED_SCORE + 1)) + phred_2;

  AR_DEBUG_ASSERT(index < PHRED_TABLE_SIZE);

  return phred_scores(index);
}

///////////////////////////////////////////////////////////////////////////////
// Public functions

alignment_info::alignment_info()
  : score(0)
  , offset(0)
  , length(0)
  , n_mismatches(0)
  , n_ambiguous(0)
  , adapter_id(-1)
{}

bool
alignment_info::is_better_than(const alignment_info& other) const
{
  if (score > other.score) {
    return true;
  } else if (score == other.score) {
    if (length > other.length) {
      return true;
    } else if (length == other.length) {
      return n_ambiguous < other.n_ambiguous;
    }
  }

  return false;
}

alignment_info
align_single_ended_sequence(const fastq& read,
                            const fastq_pair_vec& adapters,
                            int max_shift)
{
  size_t adapter_id = 0;
  alignment_info best_alignment;
  for (const auto& adapter_pair : adapters) {
    const fastq& adapter = adapter_pair.first;
    const alignment_info alignment =
      pairwise_align_sequences(best_alignment,
                               read.sequence(),
                               adapter.sequence(),
                               -max_shift,
                               std::numeric_limits<int>::max());

    if (alignment.is_better_than(best_alignment)) {
      best_alignment = alignment;
      best_alignment.adapter_id = adapter_id;
    }

    ++adapter_id;
  }

  return best_alignment;
}

alignment_info
align_paired_ended_sequences(const fastq& read1,
                             const fastq& read2,
                             const fastq_pair_vec& adapters,
                             int max_shift)
{
  size_t adapter_id = 0;
  alignment_info best_alignment;
  for (const auto& adapter_pair : adapters) {
    const fastq& adapter1 = adapter_pair.first;
    const fastq& adapter2 = adapter_pair.second;

    const std::string sequence1 = adapter2.sequence() + read1.sequence();
    const std::string sequence2 = read2.sequence() + adapter1.sequence();

    // Only consider alignments where at least one nucleotide from each read
    // is aligned against the other, included shifted alignments to account
    // for missing bases at the 5' ends of the reads.
    const int min_offset = adapter2.length() - read2.length() - max_shift;
    alignment_info alignment =
      pairwise_align_sequences(best_alignment,
                               sequence1,
                               sequence2,
                               min_offset,
                               std::numeric_limits<int>::max());

    if (alignment.is_better_than(best_alignment)) {
      best_alignment = alignment;
      best_alignment.adapter_id = adapter_id;
      // Convert the alignment into an alignment between read 1 & 2 only
      best_alignment.offset -= adapter2.length();
    }

    ++adapter_id;
  }

  return best_alignment;
}

void
truncate_single_ended_sequence(const alignment_info& alignment, fastq& read)
{
  // Given a shift, the alignment of the adapter may start one or more
  // bases before the start of the sequence, leading to a negative offset
  const size_t len = std::max<int>(0, alignment.offset);

  return read.truncate(0, len);
}

size_t
truncate_paired_ended_sequences(const alignment_info& alignment,
                                fastq& read1,
                                fastq& read2)
{
  size_t had_adapter = 0;
  const int template_length =
    std::max<int>(0, static_cast<int>(read2.length()) + alignment.offset);
  if (alignment.offset > static_cast<int>(read1.length())) {
    throw std::invalid_argument("invalid offset");
  } else if (alignment.offset >= 0) {
    // Read1 can potentially extend past read2, but by definition read2
    // cannot extend past read1 when the offset is not negative, so there
    // is no need to edit read2.
    had_adapter += static_cast<size_t>(template_length) < read1.length();
    read1.truncate(0, static_cast<size_t>(template_length));
  } else {
    had_adapter += static_cast<size_t>(template_length) < read1.length();
    had_adapter += static_cast<size_t>(template_length) < read2.length();

    read1.truncate(0, static_cast<size_t>(template_length));
    read2.truncate(
      static_cast<size_t>(static_cast<int>(read2.length()) - template_length));
  }

  return had_adapter;
}

std::string
strip_mate_info(const std::string& header, const char mate_sep)
{
  size_t pos = header.find_first_of(' ');
  if (pos == std::string::npos) {
    pos = header.length();
  }

  if (pos >= 2 && header.at(pos - 2) == mate_sep) {
    const char digit = header.at(pos - 1);

    if (digit == '1' || digit == '2') {
      std::string new_header = header;
      return new_header.erase(pos - 2, 2);
    }
  }

  return header;
}

sequence_merger::sequence_merger()
  : m_mate_sep(MATE_SEPARATOR)
  , m_conservative(false)
{}

void
sequence_merger::set_mate_separator(char sep)
{
  m_mate_sep = sep;
}

void
sequence_merger::set_conservative(bool enabled)
{
  m_conservative = enabled;
}

fastq
sequence_merger::merge(const alignment_info& alignment,
                       const fastq& read1,
                       const fastq& read2)
{
  if (alignment.offset > static_cast<int>(read1.length())) {
    // Gap between the two reads is not allowed
    throw std::invalid_argument("invalid offset");
  }

  // Offset to the first base overlapping read 2
  const size_t read_1_offset =
    static_cast<size_t>(std::max(0, alignment.offset));
  // Offset to the last base overlapping read 1
  const size_t read_2_offset =
    static_cast<int>(read1.length()) - std::max(0, alignment.offset);

  AR_DEBUG_ASSERT(read1.length() - read_1_offset == read_2_offset);

  // Produce draft by merging seq 1 and the parts of seq 2 that extend past seq
  // 1
  std::string merged_seq = read1.sequence();
  merged_seq.append(read2.sequence(), read_2_offset, std::string::npos);
  std::string merged_qual = read1.qualities();
  merged_qual.append(read2.qualities(), read_2_offset, std::string::npos);

  // Pick the best bases for the overlapping part of the reads
  for (size_t i = 0; i < read_2_offset; ++i) {
    char& nt_1 = merged_seq.at(i + read_1_offset);
    char& qual_1 = merged_qual.at(i + read_1_offset);
    const char nt_2 = read2.sequence().at(i);
    const char qual_2 = read2.qualities().at(i);

    if (m_conservative) {
      conservative_merge(nt_1, qual_1, nt_2, qual_2);
    } else {
      original_merge(nt_1, qual_1, nt_2, qual_2);
    }
  }

  // Remove mate number from read, if present, when building new record
  return fastq(m_mate_sep ? strip_mate_info(read1.header(), m_mate_sep)
                          : read1.header(),
               merged_seq,
               merged_qual,
               FASTQ_ENCODING_SAM);
}

void
sequence_merger::original_merge(char& nt_1,
                                char& qual_1,
                                char nt_2,
                                char qual_2)
{
  if (nt_1 == 'N' || nt_2 == 'N') {
    // If one of the bases are N, then we suppose that we just have (at
    // most) a single read at that site and choose that.
    if (nt_1 == 'N' && nt_2 == 'N') {
      qual_1 = PHRED_OFFSET_33;
    } else if (nt_1 == 'N') {
      nt_1 = nt_2;
      qual_1 = qual_2;
    }
  } else if (nt_1 != nt_2 && qual_1 == qual_2) {
    nt_1 = 'N';
    qual_1 = PHRED_OFFSET_33;
  } else {
    // Ensure that nt_1 / qual_1 always contains the preferred nt / score
    // This is an assumption of the g_updated_phred_scores cache.
    if (qual_1 < qual_2) {
      std::swap(nt_1, nt_2);
      std::swap(qual_1, qual_2);
    }

    const phred_scores& new_scores = get_updated_phred_scores(qual_1, qual_2);

    qual_1 =
      (nt_1 == nt_2) ? new_scores.identical_nts : new_scores.different_nts;
  }
}

void
sequence_merger::conservative_merge(char& nt_1,
                                    char& qual_1,
                                    char nt_2,
                                    char qual_2)
{

  if (nt_2 == 'N' || nt_1 == 'N') {
    if (nt_1 == 'N' && nt_2 == 'N') {
      qual_1 = PHRED_OFFSET_33;
    } else if (nt_1 == 'N') {
      nt_1 = nt_2;
      qual_1 = qual_2;
    }
  } else if (nt_1 == nt_2) {
    qual_1 = std::max(qual_1, qual_2);
  } else {
    if (qual_1 < qual_2) {
      nt_1 = nt_2;
      qual_1 = qual_2 - qual_1 + PHRED_OFFSET_33;
    } else if (qual_1 > qual_2) {
      qual_1 = qual_1 - qual_2 + PHRED_OFFSET_33;
    } else {
      // No way to reasonably pick a base
      nt_1 = 'N';
      qual_1 = PHRED_OFFSET_33;
    }
  }
}

bool
extract_adapter_sequences(const alignment_info& alignment,
                          fastq& read1,
                          fastq& read2)
{
  const int template_length =
    std::max(0, static_cast<int>(read2.length()) + alignment.offset);
  if (alignment.offset > static_cast<int>(read1.length())) {
    throw std::invalid_argument("invalid offset");
  }

  read1.truncate(std::min<size_t>(read1.length(), template_length));
  read2.truncate(
    0, std::max<int>(0, static_cast<int>(read2.length()) - template_length));

  return read1.sequence().length() || read2.sequence().length();
}
