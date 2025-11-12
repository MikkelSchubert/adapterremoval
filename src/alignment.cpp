// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "alignment.hpp"     // declarations
#include "commontypes.hpp"   // for merge_strategy
#include "debug.hpp"         // for AR_REQUIRE
#include "fastq.hpp"         // for fastq
#include "fastq_enc.hpp"     // for PHRED_OFFSET_MIN
#include "sequence_sets.hpp" // for adapter_set
#include "simd.hpp"          // for size_t, get_compare_subsequences_func
#include <algorithm>         // for max, min
#include <limits>            // for numeric_limits
#include <ostream>           // for ostream
#include <string>            // for string, operator+
#include <string_view>       // for string_view
#include <utility>           // for swap, pair

namespace adapterremoval {

bool
sequence_aligner::pairwise_align_sequences(alignment_info& alignment,
                                           const char* seq1,
                                           const size_t seq1_len,
                                           const char* seq2,
                                           const size_t seq2_len,
                                           const int min_offset,
                                           const int max_offset) const
{
  int offset =
    // The alignment must involve at least one base from seq2,
    std::max(std::max(min_offset, -static_cast<int>(seq2_len) + 1),
             // but there's no point aligning pairs too short to matter. This
             // currently only applies to --adapter-table mode.
             alignment.score() - static_cast<int>(seq2_len));
  // Alignments involving fewer than `score` bases are not interesting, since
  // the maximum possible score is `length` when all bases match
  int end_offset =
    std::min(max_offset,
             static_cast<int>(seq1_len) - std::max(1, alignment.score()));

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
    current.m_offset = offset;
    current.m_length = length;

    if (m_compare_func(current.m_n_mismatches,
                       current.m_n_ambiguous,
                       seq1 + initial_seq1_offset,
                       seq2 + initial_seq2_offset,
                       length,
                       length - alignment.score()) &&
        current.is_better_than(alignment)) {
      alignment = current;
      alignment_found = true;

      end_offset =
        std::min(end_offset,
                 static_cast<int>(seq1_len) - std::max(1, alignment.score()));
    }
  }

  return alignment_found;
}

void
sequence_aligner::update_index(alignment_info& alignment, const int max_offset)
{
  // Convert prioritized adapter index to user-supplied index
  int index = alignment.adapter_id();
  if (index >= 0) {
    AR_REQUIRE(index < static_cast<int>(m_adapters.size()));
    alignment.m_adapter_id = m_adapters.at(index).adapter_id;

    // Prioritize alignments if they involve adapter sequence
    if (alignment.m_offset < max_offset) {
      m_adapters.at(index).hits++;

      while (index > 0 &&
             m_adapters.at(index).hits > m_adapters.at(index - 1).hits) {
        std::swap(m_adapters.at(index), m_adapters.at(index - 1));
        index--;
      }
    }
  }
}

void
sequence_aligner::update_flags(alignment_info& alignment,
                               const bool paired_end) const
{
  AR_REQUIRE(alignment.m_length >=
             alignment.m_n_ambiguous + alignment.m_n_mismatches);

  alignment.m_type = alignment_type::bad;

  // Only pairs of called bases are considered part of the alignment
  const size_t n_aligned = alignment.m_length - alignment.m_n_ambiguous;

  if (alignment.score() > 0 && (paired_end || n_aligned >= m_min_se_overlap)) {
    auto mm_threshold = static_cast<size_t>(m_mismatch_threshold * n_aligned);
    if (n_aligned < 6) {
      mm_threshold = 0;
    } else if (n_aligned < 10) {
      // Allow at most 1 mismatch, possibly set to 0 by the user
      mm_threshold = std::min<size_t>(1, mm_threshold);
    }

    if (alignment.m_n_mismatches <= mm_threshold) {
      alignment.m_type = (paired_end && (n_aligned >= m_merge_threshold))
                           ? alignment_type::mergeable
                           : alignment_type::good;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Public functions

bool
alignment_info::is_better_than(const alignment_info& other) const
{
  if (score() > other.score()) {
    return true;
  } else if (score() == other.score()) {
    if (m_length > other.m_length) {
      return true;
    } else if (m_length == other.m_length) {
      return m_n_ambiguous < other.m_n_ambiguous;
    }
  }

  return false;
}

void
alignment_info::truncate_single_end(fastq& read) const
{
  // Given a shift, the alignment of the adapter may start one or more
  // bases before the start of the sequence, leading to a negative offset
  const size_t len = std::max<int>(0, m_offset);

  return read.truncate(0, len);
}

size_t
alignment_info::truncate_paired_end(fastq& read1, fastq& read2) const
{
  size_t had_adapter = 0;
  const size_t isize = insert_size(read1, read2);

  if (m_offset >= 0) {
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

bool
alignment_info::extract_adapter_sequences(fastq& read1, fastq& read2) const
{
  AR_REQUIRE(m_offset <= static_cast<int>(read1.length()));
  const int template_length =
    std::max(0, static_cast<int>(read2.length()) + m_offset);

  read1.truncate(std::min<size_t>(read1.length(), template_length));
  read2.truncate(
    0,
    std::max<int>(0, static_cast<int>(read2.length()) - template_length));

  return read1.sequence().length() || read2.sequence().length();
}

size_t
alignment_info::insert_size(const fastq& read1, const fastq& read2) const
{
  AR_REQUIRE(m_offset <= static_cast<int>(read1.length()));

  return std::max<int>(0, static_cast<int>(read2.length()) + m_offset);
}

alignment_info
alignment_info::add(const alignment_info& other) const
{
  alignment_info merged;
  merged.m_length = m_length + other.m_length;
  merged.m_n_mismatches = m_n_mismatches + other.m_n_mismatches;
  merged.m_n_ambiguous = m_n_ambiguous + other.m_n_ambiguous;

  return merged;
}

bool
alignment_info::operator==(const alignment_info& other) const
{
  return (m_offset == other.m_offset) && (m_length == other.m_length) &&
         (m_n_mismatches == other.m_n_mismatches) &&
         (m_n_ambiguous == other.m_n_ambiguous) &&
         (m_adapter_id == other.m_adapter_id) && (m_type == other.m_type);
}

////////////////////////////////////////////////////////////////////////////////

alignment_info
simple_alignment(const dna_sequence& seq1, const dna_sequence& seq2)
{
  alignment_info alignment;
  alignment.m_length = std::min(seq1.length(), seq2.length());

  // SIMD acceleration is not used. This removes the need to pad the sequences
  const auto func =
    simd::get_compare_subsequences_func(simd::instruction_set::none);

  func(alignment.m_n_mismatches,
       alignment.m_n_ambiguous,
       seq1.as_string().data(),
       seq2.as_string().data(),
       alignment.m_length,
       // prevent early termination of the alignment
       std::numeric_limits<size_t>::max());

  return alignment;
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `sequence_aligner`

sequence_aligner::sequence_aligner(const adapter_set& adapters,
                                   simd::instruction_set is,
                                   double mismatch_threshold)
  : m_compare_func(simd::get_compare_subsequences_func(is))
  , m_padding(simd::padding(is))
  , m_mismatch_threshold(mismatch_threshold)
  , m_merge_threshold(std::numeric_limits<decltype(m_merge_threshold)>::max())
{
  if (adapters.empty()) {
    // Assumes that reads do not container adapter sequences
    m_adapters.emplace_back();
    m_assume_pretrimmed = true;
  } else {
    for (const auto& it : adapters) {
      m_adapters.push_back({ static_cast<int>(m_adapters.size()),
                             0,
                             it.first.as_string(),
                             it.second.as_string() });
    }
  }
}

alignment_info
sequence_aligner::align_single_end(const fastq& read, unsigned max_shift)
{
  int adapter_id = 0;
  alignment_info alignment;

  for (const auto& it : m_adapters) {
    const std::string_view adapter = it.adapter1;
    if (adapter.empty() ||
        std::min<int>(read.length(), adapter.length()) < alignment.score()) {
      ++adapter_id;
      continue;
    }

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
                                 -max_shift,
                                 read.length())) {
      alignment.m_adapter_id = adapter_id;
    }

    ++adapter_id;
  }

  // All single-end alignments involves adapter sequence
  update_index(alignment, std::numeric_limits<int>::max());
  update_flags(alignment, false);

  return alignment;
}

alignment_info
sequence_aligner::align_paired_end(const fastq& read1,
                                   const fastq& read2,
                                   unsigned max_shift)
{
  int adapter_id = 0;
  int max_adapter_offset = std::numeric_limits<int>::min();
  alignment_info alignment;

  for (const auto& it : m_adapters) {
    const std::string_view adapter1 = it.adapter1;
    const std::string_view adapter2 = it.adapter2;

    const auto max_score = std::min<int>(read1.length(), adapter1.length()) +
                           std::min<int>(read2.length(), adapter2.length());
    // Shifting may align additional bases, depending on read lengths
    if (max_score + static_cast<int>(max_shift) < alignment.score()) {
      ++adapter_id;
      continue;
    }

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
    const int min_offset =
      m_assume_pretrimmed
        ? read1.length() - std::min(read1.length(), read2.length())
        : adapter2.length() - read2.length() - max_shift;

    // Only check for end/end overlaps during the first alignment; subsequent
    // alignments need only consider offsets involving the current adapters
    const int max_offset = (adapter_id ? adapter2.length() : sequence1_len) - 1;

    if (pairwise_align_sequences(alignment,
                                 sequence1,
                                 sequence1_len,
                                 sequence2,
                                 sequence2_len,
                                 min_offset,
                                 max_offset)) {
      alignment.m_adapter_id = adapter_id;
      // Convert the alignment into an alignment between read 1 & 2 only
      alignment.m_offset -= adapter2.length();
      max_adapter_offset = adapter2.length();
    }

    adapter_id++;
  }

  update_index(alignment, max_adapter_offset);
  update_flags(alignment, true);

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
  AR_REQUIRE(alignment.offset() <= static_cast<int>(read1.length()));

  // Offset to the first base overlapping read 2
  const auto read_1_offset =
    static_cast<size_t>(std::max(0, alignment.offset()));
  // Offset to the last base overlapping read 1
  const size_t read_2_offset =
    static_cast<int>(read1.length()) - std::max(0, alignment.offset());

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

std::ostream&
operator<<(std::ostream& os, const alignment_type& value)
{
  switch (value) {
    case alignment_type::bad:
      return os << "alignment_type::bad";
    case alignment_type::good:
      return os << "alignment_type::good";
    case alignment_type::mergeable:
      return os << "alignment_type::mergeable";
    default:                             // GCOVR_EXCL_LINE
      AR_FAIL("invalid alignment_type"); // GCOVR_EXCL_LINE
  }
}

std::ostream&
operator<<(std::ostream& os, const alignment_info& value)
{
  return os << "alignment_info{score=" << value.score()
            << ", adapter_id=" << value.m_adapter_id
            << ", offset=" << value.m_offset << ", length=" << value.m_length
            << ", n_mismatches=" << value.m_n_mismatches
            << ", n_ambiguous=" << value.m_n_ambiguous
            << ", m_type=" << value.type() << "}";
}

} // namespace adapterremoval
