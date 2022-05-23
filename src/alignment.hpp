/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
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

#include "fastq.hpp"     // for fastq_pair_vec, fastq
#include "fastq_enc.hpp" // for MATE_SEPARATOR

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

  //! Alignment score; equal to length - n_ambiguous - 2 * n_mismatches;
  int score;
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

class sequence_aligner
{
public:
  explicit sequence_aligner(const fastq_pair_vec& adapters);

  /** Set mismatch threshold for alignments returned by the aligner. */
  void set_mismatch_threshold(double mm);

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
  alignment_info align_single_end(const fastq& read, int max_shift) const;

  /**
   * Attempts to align PE mates, along with any adapter pairs.
   *
   * @param read1 A mate 1 read potentially containing adapter sequences
   * @param read2 A mate 2 read potentially containing adapter sequences
   * @param max_shift Allow up to this number of missing bases at the 5' end of
   *                  both mate reads.
   * @return The best alignment, or a length 0 alignment if not aligned.
   *
   * The alignment is carried out following the concatenation of pcr2 and read1,
   * and the concatenation of read2 and pcr1, resulting in this alignment:
   *
   *                pcr2-read1
   *                read2-pcr1
   *
   * Note the returned offset is relative read1, not to adapter2 + read1,
   * and can be used to directly infer the alignment between read1 and read2.
   */
  alignment_info align_paired_end(const fastq& read1,
                                  const fastq& read2,
                                  int max_shift) const;

private:
  /**
   * Perform pairwise alignment between two sequences.
   *
   * @param best_alignment Do not return alignments worse than this alignment.
   * @param seq1 First sequence to align (mate 1).
   * @param seq2 Second sequence to align (mate 2 or adapter).
   * @param min_offset Search for alignments from this offset.
   */
  alignment_info pairwise_align_sequences(const alignment_info& best_alignment,
                                          const std::string& seq1,
                                          const std::string& seq2,
                                          int min_offset) const;

  //! Adapter sequences against which to align the sequences
  const fastq_pair_vec& m_adapters;
  //! Maximum acceptable error rate
  double m_mismatch_threshold;
};

/**
 * Class for merging two sequence fragments into a single sequence, either
 * picking the highest quality base and its assosiated quality score, or
 * recalulating the quality score of matching/mismatching bases using the
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

  /**
   * If enabled, the sequence merger performs a more conservative merging
   * inspired by fastq-join and NGmerge, with the added caveat that mismatches
   * without a higher quality choice are assigned 'N'.
   */
  void set_conservative(bool enabled = false);

  /** Sets the maximum base quality score for recaculated scores. */
  void set_max_recalculated_score(char max);

  /**
   * Merges two overlapping reads into a single sequence, recalculating the
   * quality in one of two ways. If `conservative` mode is enabled, the
   * highest quality score of the two bases is used for matches, and the
   * difference is used for matches. Otherwise an updated score is caculated
   * bases on the original quality scores using the original algorithm
   * implemented in AdapterRemoval.
   *
   * Note that the sequences are assumed to have been trimmed using the
   * function, and this function will produce undefined results if this is not
   * the case!
   */
  void merge(const alignment_info& alignment, fastq& read1, const fastq& read2);

private:
  /** The original merging algorithm implemented in AdapterRemoval. */
  void original_merge(char& nt_1, char& qual_1, char nt_2, char qual_2);

  /** Alternative merging algorithm added in 2.4.0. */
  void conservative_merge(char& nt_1, char& qual_1, char nt_2, char qual_2);

  //! Mate separator used in read names
  char m_mate_sep;
  //! Whether to recalculate scores using the conservative or the standard mode
  bool m_conservative;
  //! Maximum score when recalculating qualities in non-conservative mode
  char m_max_score;
};

/**
 * Truncates reads such that only adapter sequence remains.
 *
 * @return True if either or both reads contained adapter sequence.
 *
 * Reads that do not contain any adapter sequence are completely truncated,
 * such no bases remain of the original sequence.
 */
bool
extract_adapter_sequences(const alignment_info& alignment,
                          fastq& pcr1,
                          fastq& pcr2);
