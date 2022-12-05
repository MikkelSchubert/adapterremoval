/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
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
#pragma once

#include <random>   // for mt19937
#include <stddef.h> // for size_t

#include "alignment_common.hpp" // for alignment_info
#include "commontypes.hpp"      // for merge_strategy
#include "fastq.hpp"            // for fastq_pair_vec, fastq
#include "fastq_enc.hpp"        // for MATE_SEPARATOR

namespace adapterremoval {

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
   * @param alignment Current best alignment.
   * @param seq1 First sequence to align (mate 1).
   * @param seq2 Second sequence to align (mate 2 or adapter).
   * @param min_offset Search for alignments from this offset.
   * @return true if a better alignment was found
   */
  bool pairwise_align_sequences(alignment_info& alignment,
                                const std::string& seq1,
                                const std::string& seq2,
                                const int min_offset) const;

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

  /** Set the strategy used when merging bases. */
  void set_merge_strategy(merge_strategy strategy);
  /** Sets the maximum base quality score for recaculated scores. */
  void set_max_recalculated_score(char max);
  /* Set the RNG used for when performing "original" merging. */
  void set_rng(std::mt19937* rng = nullptr);

  /**
   * Merges two overlapping, trimmed reads into a single sequence. Bases and
   * quality scores are assigned based on the merge_strategy chosen.
   * The sequences are assumed to have been trimmed using the given alignment.
   * This function will produce undefined results if that is not the case!
   *
   * TODO: Simplify to take a "size_t overlap" param instead of "alignment".
   */
  void merge(const alignment_info& alignment,
             fastq& read1,
             const fastq& read2) const;

private:
  /** The original merging algorithm implemented in AdapterRemoval. */
  void original_merge(char& nt_1, char& qual_1, char nt_2, char qual_2) const;

  /** Alternative merging algorithm added in 2.4.0. */
  void conservative_merge(char& nt_1,
                          char& qual_1,
                          char nt_2,
                          char qual_2) const;

  //! Mate separator used in read names
  char m_mate_sep;
  //! Strategy used when merging reads
  merge_strategy m_merge_strategy;
  //! Maximum score when recalculating qualities in non-conservative mode
  char m_max_score;
  //! Optional RNG for picking bases at random for mismatches with the same
  //! quality
  std::mt19937* m_rng;
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

} // namespace adapterremoval
