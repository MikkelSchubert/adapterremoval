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
#include "alignment.hpp"     // declarations
#include "commontypes.hpp"   // for merge_strategy
#include "debug.hpp"         // for AR_REQUIRE
#include "fastq.hpp"         // for fastq
#include "sequence_sets.hpp" // for adapter_set
#include "simd.hpp"          // for size_t, get_compare_subsequences_func
#include <algorithm>         // for max, min
#include <string>            // for string, operator+
#include <utility>           // for swap, pair

namespace adapterremoval {

bool
sequence_aligner::pairwise_align_sequences(alignment_info& alignment,
                                           const char* seq1,
                                           const size_t seq1_len,
                                           const char* seq2,
                                           const size_t seq2_len,
                                           const int min_offset) const
{
  int offset =
    // The alignment must involve at least one base from seq2,
    std::max(std::max(min_offset, -static_cast<int>(seq2_len) + 1),
             // but there's no point aligning pairs too short to matter. This
             // currently only applies to --adapter-list mode.
             alignment.score() - static_cast<int>(seq2_len));
  int end_offset = static_cast<int>(seq1_len) - std::max(1, alignment.score());

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

    const auto length =
      std::min(seq1_len - initial_seq1_offset, seq2_len - initial_seq2_offset);

    alignment_info current;
    current.offset = offset;
    current.length = length;

    AR_REQUIRE(static_cast<int>(length) >= alignment.score());
    if (m_compare_func(current.n_mismatches,
                       current.n_ambiguous,
                       seq1 + initial_seq1_offset,
                       seq2 + initial_seq2_offset,
                       length,
                       length - alignment.score()) &&
        current.is_better_than(alignment)) {
      alignment = current;
      alignment_found = true;

      // Alignments involving fewer than `score` bases are not interesting
      end_offset = static_cast<int>(seq1_len) - alignment.score();
    }
  }

  return alignment_found;
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

sequence_aligner::sequence_aligner(const adapter_set& adapters,
                                   simd::instruction_set is)
  : m_adapters(adapters)
  , m_compare_func(simd::get_compare_subsequences_func(is))
  , m_padding(simd::padding(is))
{
}

alignment_info
sequence_aligner::align_single_end(const fastq& read, int max_shift)
{
  int adapter_id = 0;
  alignment_info alignment;

  for (const auto& adapter_pair : m_adapters) {
    const std::string_view adapter = adapter_pair.first;

    m_buffer.clear();
    m_buffer += read.sequence();
    m_buffer.append(m_padding, 'N');
    m_buffer += adapter;
    m_buffer.append(m_padding, 'N');

    const char* read_data = m_buffer.data();
    const char* adapter_data = m_buffer.data() + read.length() + m_padding;

    if (pairwise_align_sequences(alignment,
                                 read_data,
                                 read.length(),
                                 adapter_data,
                                 adapter.length(),
                                 -max_shift)) {
      alignment.adapter_id = adapter_id;
    }

    ++adapter_id;
  }

  return alignment;
}

alignment_info
sequence_aligner::align_paired_end(const fastq& read1,
                                   const fastq& read2,
                                   int max_shift)
{
  int adapter_id = 0;
  alignment_info alignment;

  for (const auto& adapter_pair : m_adapters) {
    const std::string_view adapter1 = adapter_pair.first;
    const std::string_view adapter2 = adapter_pair.second;

    m_buffer.clear();
    m_buffer += adapter2;
    m_buffer += read1.sequence();
    m_buffer.append(m_padding, 'N');
    m_buffer += read2.sequence();
    m_buffer += adapter1;
    m_buffer.append(m_padding, 'N');

    const char* sequence1 = m_buffer.data();
    const size_t sequence1_len = adapter2.length() + read1.length();
    const char* sequence2 = m_buffer.data() + sequence1_len + m_padding;
    const size_t sequence2_len = adapter1.length() + read2.length();

    // Only consider alignments where at least one nucleotide from each read
    // is aligned against the other, included shifted alignments to account
    // for missing bases at the 5' ends of the reads.
    const int min_offset = adapter2.length() - read2.length() - max_shift;

    if (pairwise_align_sequences(alignment,
                                 sequence1,
                                 sequence1_len,
                                 sequence2,
                                 sequence2_len,
                                 min_offset)) {
      alignment.adapter_id = adapter_id;
      // Convert the alignment into an alignment between read 1 & 2 only
      alignment.offset -= adapter2.length();
    }

    adapter_id++;
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
  : m_merge_strategy(merge_strategy::maximum)
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
sequence_merger::merge(const alignment_info& alignment,
                       fastq& read1,
                       const fastq& read2)
{
  AR_REQUIRE(m_merge_strategy != merge_strategy::none);

  // Gap between the two reads is not allowed
  AR_REQUIRE(alignment.offset <= static_cast<int>(read1.length()));

  // Offset to the first base overlapping read 2
  const auto read_1_offset = static_cast<size_t>(std::max(0, alignment.offset));
  // Offset to the last base overlapping read 1
  const size_t read_2_offset =
    static_cast<int>(read1.length()) - std::max(0, alignment.offset);

  AR_REQUIRE(read1.length() - read_1_offset == read_2_offset);

  // Produce draft by merging r1 and the parts of r2 that extend past r1
  read1.m_sequence.append(read2.sequence(), read_2_offset);
  read1.m_qualities.append(read2.qualities(), read_2_offset);
  AR_REQUIRE(read1.m_sequence.length() == read1.m_qualities.length());

  // Pick the best bases for the overlapping part of the reads
  for (size_t i = 0; i < read_2_offset; ++i) {
    char& nt_1 = read1.m_sequence.at(i + read_1_offset);
    char& qual_1 = read1.m_qualities.at(i + read_1_offset);
    const char nt_2 = read2.sequence().at(i);
    const char qual_2 = read2.qualities().at(i);

    if (nt_2 == 'N' || nt_1 == 'N') {
      if (nt_1 == 'N' && nt_2 == 'N') {
        qual_1 = PHRED_OFFSET_MIN;
      } else if (nt_1 == 'N') {
        nt_1 = nt_2;
        qual_1 = qual_2;
      }
    } else if (nt_1 == nt_2) {
      if (m_merge_strategy == merge_strategy::maximum) {
        qual_1 = std::max(qual_1, qual_2);
      } else {
        qual_1 =
          std::min<char>(m_max_score, qual_1 + qual_2 - PHRED_OFFSET_MIN);
      }
    } else {
      if (qual_1 < qual_2) {
        nt_1 = nt_2;
        qual_1 = qual_2 - qual_1 + PHRED_OFFSET_MIN;
        m_mismatches_resolved++;
      } else if (qual_1 > qual_2) {
        qual_1 = qual_1 - qual_2 + PHRED_OFFSET_MIN;
        m_mismatches_resolved++;
      } else {
        // No way to reasonably pick a base
        nt_1 = 'N';
        qual_1 = PHRED_OFFSET_MIN;
        m_mismatches_unresolved++;
      }
    }
  }

  m_reads_merged += 2;
  m_bases_merged += read_2_offset;

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
