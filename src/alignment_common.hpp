/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <cstddef> // for size_t

namespace adapterremoval {

class fastq;

using std::size_t;

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
struct alignment_info
{
  /** Defaults to unaligned (len = 0), for adapter_id -1. **/
  alignment_info();

  /**
   * Returns true if this is a better alignment than other.
   *
   * When selecting among multiple alignments, the follow criteria are used:
   * 1. The alignment with the highest score is preferred.
   * 2. If score is equal, the longest alignment is preferred.
   * 3. If score and length is equal, the alignment with fewest Ns is preferred.
   */
  bool is_better_than(const alignment_info& other) const;

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

  /** Calculates the insert size given a pair of un-truncated reads. */
  size_t insert_size(const fastq& read1, const fastq& read2) const;

  /** Returns the score used to compare alignments */
  inline int score() const
  {
    return static_cast<int>(length) -
           static_cast<int>(n_ambiguous + 2 * n_mismatches);
  }

  //! Zero based id of the adapter which offered the best alignment. Is less
  //! than zero if no alignment was found.
  int offset;
  //! The number of base-pairs included in the alignment. This number
  //! includes both bases aligned between the two mates (in PE mode) and the
  //! number of bases aligned between mates and adapter sequences.
  size_t length;
  //! Number of positions in the alignment in which the two sequences were
  //! both called (not N) but differed
  size_t n_mismatches;
  //! Number of positions in the alignment where one or both bases were N.
  size_t n_ambiguous;
  //! 0-based ID of best matching adapter or a negative value if not set.
  int adapter_id;
};

/** Returns true if the current CPU supports SSE2 instructions */
bool
supports_sse2();

/** Returns true if the current CPU supports AVX2 instructions */
bool
supports_avx2();

/**
 * Compares two subsequences in an alignment to a previous (best) alignment.
 *
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
compare_subsequences(alignment_info& current,
                     const char* seq_1_ptr,
                     const char* seq_2_ptr,
                     const double mismatch_threshold = 1.0);

bool
compare_subsequences_std(alignment_info& current,
                         const char* seq_1_ptr,
                         const char* seq_2_ptr,
                         size_t max_mismatches,
                         int remaining_bases);

bool
compare_subsequences_sse2(alignment_info& current,
                          const char* seq_1_ptr,
                          const char* seq_2_ptr,
                          size_t max_mismatches,
                          int remaining_bases);

bool
compare_subsequences_avx2(alignment_info& current,
                          const char* seq_1_ptr,
                          const char* seq_2_ptr,
                          size_t max_mismatches,
                          int remaining_bases);

} // namespace adapterremoval
