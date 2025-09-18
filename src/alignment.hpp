// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "fastq_enc.hpp" // for MATE_SEPARATOR
#include "simd.hpp"      // for size_t, compare_subsequences_func, instru...
#include <iosfwd>        // for ostream
#include <string>        // for string
#include <string_view>   // for string_view
#include <vector>        // for vector

namespace adapterremoval {

enum class merge_strategy;
class adapter_set;
class fastq;

/** Alignment evaluation from user-specified error-rate and merge threshold */
enum class alignment_type
{
  //! Either unaligned or a poor alignment
  bad = 0,
  //! Good alignment, but cannot be merged (if PE)
  good,
  //! Good alignment that can be merged (if PE)
  mergeable,
};

/**
 * Summarizes an alignment.
 *
 * A single offset value is used to represent the alignment between two
 * sequences, with values ranging from -inf to seq1.len() - 1. The alignment
 * of the first base in each sequence against each other is defined as having
 * the offset 0, with each other offset defined as the relative position of
 * seq 2 base 0 to seq 1 base 0:
 *
 * Seq 1:           aaaaaaaaa
 * seq 2:           bbbbbbbbbbb
 * Offset: 0
 *
 * Seq 1:            aaaaaaaaa
 * Seq 2: bbbbbbbbbbb
 * Offset: -11
 *
 * Seq 1:           aaaaaaaaa
 * Seq 2:                   bbbbbbbbbbb
 * Offset: 8
 *
 * The meaning of the offset is slightly different in SE and PE mode; in SE
 * mode seq2 is the adapter sequence, and the offset therefore unambiguously
 * shows the starting position of the adapter, regardless of the size of the
 * adapter sequence.
 *
 * In PE mode, while the alignment is between seq1+PCR2 and PCR1+seq2, the
 * offset returned is relative to to seq1 and seq2 only. Thus, if '1'
 * represents PCR1 and '2' represents PCR2, the following alignment results
 * from PE mode (ie. offset = -9 rather than 3):
 *
 * Seq 1:    22222222222aaaaaaaaa
 * Seq 2:      bbbbbbbbbbb1111111111
 * Offset: -9
 *
 * Note that the offset can never be greater than len(read 1) - 1, but can be
 * less than -len(seq 2) + 1, if a positive shift is set in PE mode. In PE
 * mode, an offset <= -len(seq 2) indicates that no non-adapter sequence was
 * found, while an offset <= 0 indicates the same for SE mode. Offsets less
 * than -len(seq 2) for PE or less than 0 for SE indicates that bases have been
 * skipped during sequencing, and are discoverable if a shift is set:
 *
 * Read 1:           ...22222222222aaaaaaaaa
 * Read 2:             bbbbbbbbbbb1111111111...
 * Offset: -12
 *
 */
class alignment_info
{
public:
  alignment_info() = default;

  /**
   * Returns true if this is a better alignment than other.
   *
   * When selecting among multiple alignments, the follow criteria are used:
   * 1. The alignment with the highest score is preferred.
   * 2. If score is equal, the longest alignment is preferred.
   * 3. If score and length is equal, the alignment with fewest Ns is preferred.
   */
  [[nodiscard]] bool is_better_than(const alignment_info& other) const;

  /**
   * Truncates a SE read according to the alignment, such that the second read
   * used in the alignment (assumed to represent adapter sequence) is excluded
   * from the read passed to this function.
   */
  void truncate_single_end(fastq& read) const;

  /**
   * Truncate a pair of PE reads, such that any adapter sequence inferred from
   * the alignment is excluded from both mates.
   *
   * @return The number of sequences (0 .. 2) which contained adapter sequence.
   */
  size_t truncate_paired_end(fastq& read1, fastq& read2) const;

  /**
   * Truncates reads such that only adapter sequence remains, returning true if
   * either read contained any adapter sequence.
   */
  bool extract_adapter_sequences(fastq& read1, fastq& read2) const;

  /** Calculates the insert size given a pair of un-truncated reads. */
  [[nodiscard]] size_t insert_size(const fastq& read1,
                                   const fastq& read2) const;

  /** Returns the score used to compare alignments */
  [[nodiscard]] int score() const
  {
    return static_cast<int>(m_length) -
           static_cast<int>(m_n_ambiguous + (2 * m_n_mismatches));
  }

  /** Returns the 0-based ID of best matching adapter or a negative value  */
  [[nodiscard]] int adapter_id() const { return m_adapter_id; }

  /** Returns the offset of the sequence in the pairwise alignment  */
  [[nodiscard]] int offset() const { return m_offset; }

  /** Returns true if the alignment meets minimum requirements */
  [[nodiscard]] alignment_type type() const { return m_type; }

  /** Compares two alignments for equality. This also compares quality flags */
  bool operator==(const alignment_info& other) const;

private:
  //! 0-based ID of best matching adapter or a negative value if not set.
  int m_adapter_id = -1;

  //! Zero based id of the adapter which offered the best alignment. Is less
  //! than zero if no alignment was found.
  int m_offset = 0;

  //! The number of base-pairs included in the alignment. This number
  //! includes both bases aligned between the two mates (in PE mode) and the
  //! number of bases aligned between mates and adapter sequences.
  size_t m_length = 0;

  //! Number of positions in the alignment in which the two sequences were
  //! both called (not N) but differed
  size_t m_n_mismatches = 0;
  //! Number of positions in the alignment where one or both bases were N.
  size_t m_n_ambiguous = 0;

  //! Specifies if the alignment is considered a "good" alignment
  alignment_type m_type = alignment_type::bad;

  /** Creates debug representation of an adapter set */
  friend std::ostream& operator<<(std::ostream&, const alignment_info&);

  /** Class responsible for setting up `alignment_info` instances */
  friend class sequence_aligner;
  /** Required for testing */
  friend struct ALN;
};

class sequence_aligner
{
public:
  explicit sequence_aligner(const adapter_set& adapters,
                            simd::instruction_set is,
                            double mismatch_threshold);

  /** Sets the minimum required overlap for "good" SE alignments */
  void set_min_se_overlap(size_t n) { m_min_se_overlap = n; }

  /** Sets the minimum number of non-N aligned bases required to merge */
  void set_merge_threshold(size_t n) { m_merge_threshold = n; }

  /**
   * Attempts to align adapters sequences against a SE read.
   *
   * @param read A read potentially containing adapter sequences
   * @param max_shift Allow up to this number of missing bases at the 5' end of
   *                  the read, when aligning the adapter.
   * @return The best alignment, or a length 0 alignment if not aligned.
   *
   * The best alignment is selected using alignment_info::is_better_than.
   */
  alignment_info align_single_end(const fastq& read, unsigned max_shift);

  /**
   * Attempts to align PE mates, along with any adapter pairs.
   *
   * @param read1 A mate 1 read potentially containing adapter sequences
   * @param read2 A mate 2 read potentially containing adapter sequences
   * @param max_shift Allow up to this number of missing bases at the 5' end of
   *                  both mate reads.
   * @return The best alignment, or a length 0 alignment if not aligned.
   *
   * The alignment is carried out following the concatenation of adapter2 and
   * read1, and the concatenation of read2 and adapter1, resulting in this
   * alignment:
   *
   *    adapter2-read1
   *    read2-adapter1
   *
   * Note the returned offset is relative read1, not to adapter2 + read1,
   * and can be used to directly infer the alignment between read1 and read2.
   */
  alignment_info align_paired_end(const fastq& read1,
                                  const fastq& read2,
                                  unsigned max_shift);

private:
  /**
   * Perform pairwise alignment between two sequences.
   *
   * @param alignment Current best alignment.
   * @param seq1 First sequence to align (mate 1).
   * @param seq1_len The (unpadded) length of seq1.
   * @param seq2 Second sequence to align (mate 2 or adapter).
   * @param seq2_len The (unpadded) length of seq2.
   * @param min_offset Search for alignments from this offset.
   * @param max_offset Search for alignments until and including this offset.
   * @return true if a better alignment was found
   */
  bool pairwise_align_sequences(alignment_info& alignment,
                                const char* seq1,
                                size_t seq1_len,
                                const char* seq2,
                                size_t seq2_len,
                                int min_offset,
                                int max_offset) const;

  /** Sets the user-facing adapter ID and prioritizes adapter sequences */
  void update_index(alignment_info& alignment, int max_offset);
  /** Sets "is_good" and "can_merge" based on thresholds */
  void update_flags(alignment_info& alignment, bool paired_end) const;

  //! SIMD instruction set to use for sequence comparisons
  const simd::compare_subsequences_func m_compare_func;
  //! Padding required by chosen SIMD instructions
  const size_t m_padding;
  //! Max error rate allowed for "good" alignments
  const double m_mismatch_threshold;
  //! Number of aligned (non-N, including adapters) bases required for merging
  size_t m_merge_threshold;
  //! Minimum alignment length for "good" SE alignments, including Ns
  size_t m_min_se_overlap = 1;

  //! Internal buffer used to combine adapters and reads
  std::string m_buffer{};

  struct adapter_pair
  {
    //! The original adapter ID
    int adapter_id = 0;
    //! The number of best alignments that included this adapter
    int hits = 0;
    //! The forward adapter sequence
    std::string adapter1{};
    //! The reverse adapter sequence
    std::string adapter2{};
  };

  //! Vector of adapter IDs and number of matching alignments with that adapter
  std::vector<adapter_pair> m_adapters{};
};

/**
 * Class for merging two sequence fragments into a single sequence, either
 * picking the highest quality base and its associated quality score, or
 * recalculating the quality score of matching/mismatching bases using the
 * original AR methodology.
 */
class sequence_merger
{
public:
  sequence_merger();

  /**
   * Sets the expected mate separator. This is used to trim mate numbers from
   * merged reads.
   */
  void set_mate_separator(char sep = MATE_SEPARATOR);

  /** Set the strategy used when merging bases. */
  void set_merge_strategy(merge_strategy strategy);
  /** Sets the maximum base quality score for recalculated scores. */
  void set_max_recalculated_score(char max);

  /**
   * Merges two overlapping, trimmed reads into a single sequence. Bases and
   * quality scores are assigned based on the merge_strategy chosen.
   * The sequences are assumed to have been trimmed using the given alignment.
   * This function will produce undefined results if that is not the case!
   */
  void merge(const alignment_info& alignment, fastq& read1, const fastq& read2);

  /* Returns number of reads merged */
  size_t reads_merged() const { return m_reads_merged; }

  /* Returns number of bases merged */
  size_t bases_merged() const { return m_bases_merged; }

  /* Returns number of mismatches where a higher quality base was selected */
  size_t mismatches_resolved() const { return m_mismatches_resolved; }

  /* Returns number of mismatches where there was no higher quality base */
  size_t mismatches_unresolved() const { return m_mismatches_unresolved; }

private:
  //! Mate separator used in read names
  char m_mate_sep = MATE_SEPARATOR;
  //! Strategy used when merging reads
  merge_strategy m_merge_strategy;
  //! Maximum score when recalculating qualities in non-conservative mode
  char m_max_score = PHRED_OFFSET_MAX;

  //! Total number of reads merged
  size_t m_reads_merged = 0;
  //! The  number of bases merged
  size_t m_bases_merged = 0;
  //! The number of mismatches where a higher quality base was selected
  size_t m_mismatches_resolved = 0;
  //! The number of mismatches where there was no higher quality base
  size_t m_mismatches_unresolved = 0;
};

/** Stream operator for debugging output */
std::ostream&
operator<<(std::ostream& os, const alignment_info& value);

} // namespace adapterremoval
