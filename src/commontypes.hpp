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
#pragma once

#include <string>
#include <vector>

namespace adapterremoval {

class fastq;

using string_vec = std::vector<std::string>;
using string_vec_citer = string_vec::const_iterator;
using string_pair = std::pair<std::string, std::string>;
using string_pair_vec = std::vector<string_pair>;

using fastq_vec = std::vector<fastq>;

/** Different file-types read / generated by AdapterRemoval. */
enum class read_type : size_t
{
  /** Mate 1 reads, either read or written by AR. */
  mate_1 = 0,
  /** Mate 2 reads, either read or written by AR. */
  mate_2,
  /** Overlapping PE reads merged into a single sequence. */
  merged,
  /** PE reads for which the mate has been discarded. */
  singleton,
  /** Discarded reads; e.g. too short reads. */
  discarded,

  //! End value; not to be used as an argument.
  max,
};

/** Enum describing the user-requested output format for processed reads */
enum class output_format
{
  //! Uncompressed FASTQ reads
  fastq,
  //! Gzip compressed FASTQ reads
  fastq_gzip,
  //! Unaligned SAM (Sequence Alignment Map) records
  sam,
  //! Gzip compressed unaligned SAM (Sequence Alignment Map) records
  sam_gzip,

  // TODO:
  //! Unaligned BAM (Binary Alignment Map) records
  // bam,
  //! Uncompressed, unaligned BAM (Binary Alignment Map) records
  // ubam,
};

/** Strategy used when merging reads */
enum class merge_strategy
{
  /** Merging disabled */
  none,
  /** Assign max(Q1,Q2) to matches, abs(Q1-Q2) to mismatches, N to ambiguous */
  maximum,
  /** Assign Q1+Q2 to matches, abs(Q1-Q2) to mismatches, N to ambiguous */
  additive,
};

enum class trimming_strategy
{
  /** Quality trimming disabled */
  none,
  //! Quality trimming using the modified Mott algorithm
  mott,
  //! Sliding window based quality trimming
  window,
  //! The original quality trimming algorithms
  per_base,
};

} // namespace adapterremoval
