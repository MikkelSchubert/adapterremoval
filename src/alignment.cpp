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

#include <algorithm> // for max, min
#include <bitset>    // for bitset
#include <limits>    // for numeric_limits
#include <string>    // for string, operator+
#include <utility>   // for swap, pair

#include "alignment.hpp" // declarations
#include "debug.hpp"     // for AR_REQUIRE
#include "fastq.hpp"     // for fastq, fastq_pair_vec

namespace adapterremoval {

bool
supports_avx2()
{
  return __builtin_cpu_supports("avx2");
}

bool
supports_sse2()
{
  return __builtin_cpu_supports("sse2");
}

bool
compare_subsequences_sse2(const alignment_info& /* best */,
                          alignment_info& /* current */,
                          const char*& /* seq_1_ptr */,
                          const char*& /* seq_2_ptr */,
                          double /* mismatch_threshold */,
                          int& /* remaining_bases */);

bool
compare_subsequences_avx2(const alignment_info& /* best */,
                          alignment_info& /* current */,
                          const char*& /* seq_1_ptr */,
                          const char*& /* seq_2_ptr */,
                          double /* mismatch_threshold */,
                          int& /* remaining_bases */);

/**
 * Compares two subsequences in an alignment to a previous (best) alignment.
 *
 * @param best The currently best alignment, used for evaluating this alignment
 * @param current The current alignment to be evaluated (assumed to be zero'd).
 * @param seq_1_ptr Pointer to the first sequence in the alignment.
 * @param seq_2_ptr Pointer to the second sequence in the alignment.
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
                     const char* seq_2_ptr,
                     const double mismatch_threshold = 1.0)
{
  int remaining_bases = current.score = current.length;

  // Compare 32 bp at a time (if supported)
  static const bool avx2_accelerated = supports_avx2();
  if (avx2_accelerated && !compare_subsequences_avx2(best,
                                                     current,
                                                     seq_1_ptr,
                                                     seq_2_ptr,
                                                     mismatch_threshold,
                                                     remaining_bases)) {
    return false;
  }

  // Compare 16 bp at a time (if supported)
  static const bool sse2_accelerated = supports_sse2();
  if (sse2_accelerated && !compare_subsequences_sse2(best,
                                                     current,
                                                     seq_1_ptr,
                                                     seq_2_ptr,
                                                     mismatch_threshold,
                                                     remaining_bases)) {
    return false;
  }

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
sequence_aligner::pairwise_align_sequences(const alignment_info& best_alignment,
                                           const std::string& seq1,
                                           const std::string& seq2,
                                           int min_offset) const
{
  const int start_offset =
    std::max<int>(min_offset, -static_cast<int>(seq2.length()) + 1);
  const int end_offset = static_cast<int>(seq1.length()) - 1;

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

      if (compare_subsequences(
            best, current, seq_1_ptr, seq_2_ptr, m_mismatch_threshold)) {
        best = current;
      }
    }
  }

  return best;
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
{
}

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

void
alignment_info::truncate_single_end(fastq& read) const
{
  // Given a shift, the alignment of the adapter may start one or more
  // bases before the start of the sequence, leading to a negative offset
  const size_t len = std::max<int>(0, offset);

  return read.truncate(0, len);
}

size_t
alignment_info::truncate_paired_end(fastq& read1, fastq& read2) const
{
  size_t had_adapter = 0;
  const size_t isize = insert_size(read1, read2);

  if (offset >= 0) {
    // Read1 can potentially extend past read2, but by definition read2
    // cannot extend past read1 when the offset is not negative, so there
    // is no need to edit read2.
    had_adapter += isize < read1.length();
    read1.truncate(0, isize);
  } else {
    had_adapter += isize < read1.length();
    had_adapter += isize < read2.length();

    read1.truncate(0, isize);
    read2.truncate(read2.length() - isize);
  }

  return had_adapter;
}

size_t
alignment_info::insert_size(const fastq& read1, const fastq& read2) const
{
  AR_REQUIRE(offset <= static_cast<int>(read1.length()));

  return std::max<int>(0, static_cast<int>(read2.length()) + offset);
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `sequence_aligner`

sequence_aligner::sequence_aligner(const fastq_pair_vec& adapters)
  : m_adapters(adapters)
  , m_mismatch_threshold(1.0)
{
}

void
sequence_aligner::set_mismatch_threshold(double mm)
{
  m_mismatch_threshold = mm;
}

alignment_info
sequence_aligner::align_single_end(const fastq& read, int max_shift) const
{
  size_t adapter_id = 0;
  alignment_info best_alignment;
  for (const auto& adapter_pair : m_adapters) {
    const fastq& adapter = adapter_pair.first;
    const alignment_info alignment = pairwise_align_sequences(
      best_alignment, read.sequence(), adapter.sequence(), -max_shift);

    if (alignment.is_better_than(best_alignment)) {
      best_alignment = alignment;
      best_alignment.adapter_id = adapter_id;
    }

    ++adapter_id;
  }

  return best_alignment;
}

alignment_info
sequence_aligner::align_paired_end(const fastq& read1,
                                   const fastq& read2,
                                   int max_shift) const
{
  size_t adapter_id = 0;
  alignment_info best_alignment;
  for (const auto& adapter_pair : m_adapters) {
    const fastq& adapter1 = adapter_pair.first;
    const fastq& adapter2 = adapter_pair.second;

    const std::string sequence1 = adapter2.sequence() + read1.sequence();
    const std::string sequence2 = read2.sequence() + adapter1.sequence();

    // Only consider alignments where at least one nucleotide from each read
    // is aligned against the other, included shifted alignments to account
    // for missing bases at the 5' ends of the reads.
    const int min_offset = adapter2.length() - read2.length() - max_shift;
    const alignment_info alignment = pairwise_align_sequences(
      best_alignment, sequence1, sequence2, min_offset);

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
strip_mate_info(std::string& header, const char mate_sep)
{
  size_t pos = header.find_first_of(' ');
  if (pos == std::string::npos) {
    pos = header.length();
  }

  if (pos >= 2 && header.at(pos - 2) == mate_sep) {
    const char digit = header.at(pos - 1);

    if (digit == '1' || digit == '2') {
      header.erase(pos - 2, 2);
    }
  }
}

sequence_merger::sequence_merger()
  : m_mate_sep(MATE_SEPARATOR)
  , m_conservative(false)
  , m_max_score(PHRED_OFFSET_MAX)
{
}

void
sequence_merger::set_mate_separator(char sep)
{
  m_mate_sep = sep;
}

void
sequence_merger::merge(const alignment_info& alignment,
                       fastq& read1,
                       const fastq& read2)
{
  // Gap between the two reads is not allowed
  AR_REQUIRE(alignment.offset <= static_cast<int>(read1.length()));

  // Offset to the first base overlapping read 2
  const size_t read_1_offset =
    static_cast<size_t>(std::max(0, alignment.offset));
  // Offset to the last base overlapping read 1
  const size_t read_2_offset =
    static_cast<int>(read1.length()) - std::max(0, alignment.offset);

  AR_REQUIRE(read1.length() - read_1_offset == read_2_offset);

  // Produce draft by merging r1 and the parts of r2 that extend past r1
  read1.m_sequence.append(read2.sequence(), read_2_offset, std::string::npos);
  read1.m_qualities.append(read2.qualities(), read_2_offset, std::string::npos);
  AR_REQUIRE(read1.m_sequence.length() == read1.m_qualities.length());

  // Pick the best bases for the overlapping part of the reads
  for (size_t i = 0; i < read_2_offset; ++i) {
    char& nt_1 = read1.m_sequence.at(i + read_1_offset);
    char& qual_1 = read1.m_qualities.at(i + read_1_offset);
    const char nt_2 = read2.sequence().at(i);
    const char qual_2 = read2.qualities().at(i);

    if (nt_2 == 'N' || nt_1 == 'N') {
      if (nt_1 == 'N') {
        if (nt_2 != 'N') {
          nt_1 = nt_2;
          qual_1 = qual_2;
        } else {
          qual_1 = PHRED_OFFSET_MIN;
        }
      }
    } else if (nt_1 == nt_2) {
      qual_1 = std::max(qual_1, qual_2);
    } else {
      if (qual_1 < qual_2) {
        nt_1 = nt_2;
        qual_1 = qual_2 - qual_1 + PHRED_OFFSET_MIN;
      } else if (qual_1 > qual_2) {
        qual_1 = qual_1 - qual_2 + PHRED_OFFSET_MIN;
      } else {
        // No way to reasonably pick a base
        nt_1 = 'N';
        qual_1 = PHRED_OFFSET_MIN;
      }
    }
  }

  // Remove mate number from read, if present
  if (m_mate_sep) {
    strip_mate_info(read1.m_header, m_mate_sep);
  }
}

bool
extract_adapter_sequences(const alignment_info& alignment,
                          fastq& read1,
                          fastq& read2)
{
  AR_REQUIRE(alignment.offset <= static_cast<int>(read1.length()));
  const int template_length =
    std::max(0, static_cast<int>(read2.length()) + alignment.offset);

  read1.truncate(std::min<size_t>(read1.length(), template_length));
  read2.truncate(
    0, std::max<int>(0, static_cast<int>(read2.length()) - template_length));

  return read1.sequence().length() || read2.sequence().length();
}

} // namespace adapterremoval
