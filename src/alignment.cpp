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
#include <algorithm> // for max, min
#include <string>    // for string, operator+
#include <utility>   // for swap, pair

#include "alignment.hpp"        // declarations
#include "alignment_tables.hpp" // for DIFFERENT_NTS, IDENTICAL_NTS, PHRED_...
#include "debug.hpp"            // for AR_REQUIRE
#include "fastq.hpp"            // for fastq, fastq_pair_vec
#include "simd.hpp"             // for select_compare_function, ...

namespace adapterremoval {

bool
sequence_aligner::pairwise_align_sequences(alignment_info& alignment,
                                           const std::string& seq1,
                                           const std::string& seq2,
                                           const int min_offset) const
{
  int offset =
    // The alignment must involve at least one base from seq2,
    std::max(std::max(min_offset, -static_cast<int>(seq2.length()) + 1),
             // but there's no point aligning pairs too short to matter. This
             // currently only applies to --adapter-list mode.
             alignment.score() - static_cast<int>(seq2.length()));
  int end_offset =
    static_cast<int>(seq1.length()) - std::max(1, alignment.score());

  bool alignment_found = false;
  for (; offset <= end_offset; ++offset) {
    size_t initial_seq1_offset;
    size_t initial_seq2_offset;

    if (offset < 0) {
      initial_seq1_offset = 0;
      initial_seq2_offset = -offset;
    } else {
      initial_seq1_offset = offset;
      initial_seq2_offset = 0;
    }

    const auto length = std::min(seq1.length() - initial_seq1_offset,
                                 seq2.length() - initial_seq2_offset);

    alignment_info current;
    current.offset = offset;
    current.length = length;

    AR_REQUIRE(static_cast<int>(length) >= alignment.score());
    if (m_compare_func(current.n_mismatches,
                       current.n_ambiguous,
                       seq1.data() + initial_seq1_offset,
                       seq2.data() + initial_seq2_offset,
                       length,
                       length - alignment.score()) &&
        current.is_better_than(alignment)) {
      alignment = current;
      alignment_found = true;

      // Alignments involving fewer than `score` bases are not interesting
      end_offset = static_cast<int>(seq1.length()) - alignment.score();
    }
  }

  return alignment_found;
}

struct phred_scores
{
  explicit phred_scores(size_t index)
    : identical_nts(IDENTICAL_NTS.at(index))
    , different_nts(DIFFERENT_NTS.at(index))
  {
  }

  //! Phred score to assign if the two nucleotides are identical
  char identical_nts;
  //! Phred score to assign if the two nucleotides differ
  char different_nts;
};

phred_scores
get_updated_phred_scores(char qual_1, char qual_2)
{
  AR_REQUIRE(qual_1 >= qual_2);

  const auto phred_1 = static_cast<size_t>(qual_1 - PHRED_OFFSET_MIN);
  const auto phred_2 = static_cast<size_t>(qual_2 - PHRED_OFFSET_MIN);
  const size_t index = (phred_1 * (PHRED_SCORE_MAX + 1)) + phred_2;

  return phred_scores(index);
}

///////////////////////////////////////////////////////////////////////////////
// Public functions

bool
alignment_info::is_better_than(const alignment_info& other) const
{
  if (score() > other.score()) {
    return true;
  } else if (score() == other.score()) {
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

sequence_aligner::sequence_aligner(const fastq_pair_vec& adapters,
                                   simd::instruction_set is)
  : m_adapters(adapters)
  , m_compare_func(simd::get_compare_subsequences_func(is))
{
}

alignment_info
sequence_aligner::align_single_end(const fastq& read, int max_shift) const
{
  int adapter_id = 0;
  alignment_info alignment;
  for (const auto& adapter_pair : m_adapters) {
    const auto& adapter = adapter_pair.first.sequence();
    const auto& sequence = read.sequence();
    if (pairwise_align_sequences(alignment, sequence, adapter, -max_shift)) {
      alignment.adapter_id = adapter_id;
    }

    ++adapter_id;
  }

  return alignment;
}

alignment_info
sequence_aligner::align_paired_end(const fastq& read1,
                                   const fastq& read2,
                                   int max_shift) const
{
  size_t adapter_id = 0;
  alignment_info alignment;
  for (const auto& adapter_pair : m_adapters) {
    const fastq& adapter1 = adapter_pair.first;
    const fastq& adapter2 = adapter_pair.second;

    const std::string sequence1 = adapter2.sequence() + read1.sequence();
    const std::string sequence2 = read2.sequence() + adapter1.sequence();

    // Only consider alignments where at least one nucleotide from each read
    // is aligned against the other, included shifted alignments to account
    // for missing bases at the 5' ends of the reads.
    const int min_offset = adapter2.length() - read2.length() - max_shift;

    if (pairwise_align_sequences(alignment, sequence1, sequence2, min_offset)) {
      alignment.adapter_id = adapter_id;
      // Convert the alignment into an alignment between read 1 & 2 only
      alignment.offset -= adapter2.length();
    }

    ++adapter_id;
  }

  return alignment;
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
  , m_merge_strategy(merge_strategy::conservative)
  , m_max_score(PHRED_OFFSET_MAX)
  , m_rng()
{
}

void
sequence_merger::set_mate_separator(char sep)
{
  m_mate_sep = sep;
}

void
sequence_merger::set_merge_strategy(merge_strategy strategy)
{
  m_merge_strategy = strategy;
}

void
sequence_merger::set_max_recalculated_score(char max)
{
  m_max_score = max + '!';
}

void
sequence_merger::set_rng(std::mt19937* rng)
{
  m_rng = rng;
}

void
sequence_merger::merge(const alignment_info& alignment,
                       fastq& read1,
                       const fastq& read2) const
{
  AR_REQUIRE(m_merge_strategy != merge_strategy::none);
  AR_REQUIRE((m_merge_strategy == merge_strategy::original) == !!m_rng);

  // Gap between the two reads is not allowed
  AR_REQUIRE(alignment.offset <= static_cast<int>(read1.length()));

  // Offset to the first base overlapping read 2
  const auto read_1_offset = static_cast<size_t>(std::max(0, alignment.offset));
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

    if (m_merge_strategy == merge_strategy::conservative) {
      conservative_merge(nt_1, qual_1, nt_2, qual_2);
    } else {
      original_merge(nt_1, qual_1, nt_2, qual_2);
    }
  }

  // Remove mate number from read, if present
  if (m_mate_sep) {
    strip_mate_info(read1.m_header, m_mate_sep);
  }
}

void
sequence_merger::original_merge(char& nt_1,
                                char& qual_1,
                                char nt_2,
                                char qual_2) const
{
  if (nt_1 == 'N' || nt_2 == 'N') {
    // If one of the bases are N, then we suppose that we just have (at
    // most) a single read at that site and choose that.
    if (nt_1 == 'N' && nt_2 == 'N') {
      qual_1 = PHRED_OFFSET_MIN;
    } else if (nt_1 == 'N') {
      nt_1 = nt_2;
      qual_1 = qual_2;
    }
  } else if (nt_1 != nt_2 && qual_1 == qual_2) {
    if (m_rng) {
      nt_1 = ((*m_rng)() & 1) ? nt_1 : nt_2;

      const phred_scores& new_scores = get_updated_phred_scores(qual_1, qual_2);
      qual_1 = new_scores.different_nts;
    } else {
      nt_1 = 'N';
      qual_1 = PHRED_OFFSET_MIN;
    }
  } else {
    // Ensure that nt_1 / qual_1 always contains the preferred nt / score
    // This is an assumption of the g_updated_phred_scores cache.
    if (qual_1 < qual_2) {
      std::swap(nt_1, nt_2);
      std::swap(qual_1, qual_2);
    }

    const phred_scores& new_scores = get_updated_phred_scores(qual_1, qual_2);

    qual_1 = std::min<char>(m_max_score,
                            (nt_1 == nt_2) ? new_scores.identical_nts
                                           : new_scores.different_nts);
  }
}

void
sequence_merger::conservative_merge(char& nt_1,
                                    char& qual_1,
                                    char nt_2,
                                    char qual_2) const
{

  if (nt_2 == 'N' || nt_1 == 'N') {
    if (nt_1 == 'N' && nt_2 == 'N') {
      qual_1 = PHRED_OFFSET_MIN;
    } else if (nt_1 == 'N') {
      nt_1 = nt_2;
      qual_1 = qual_2;
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
