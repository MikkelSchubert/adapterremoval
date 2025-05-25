// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <string>      // for string
#include <string_view> // for string_view
#include <vector>      // for vector

namespace adapterremoval {

using string_vec = std::vector<std::string>;
using string_vec_citer = string_vec::const_iterator;
using string_pair = std::pair<std::string, std::string>;
using string_pair_vec = std::vector<string_pair>;

/** Different read files read or written by AdapterRemoval */
enum class read_file : size_t
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

/** Represents types of reads read or written;  */
enum class read_type
{
  //! SE read
  se,
  //! SE read that failed QC
  se_fail,
  //! PE mate 1 read
  pe_1,
  //! PE mate 1 read that failed QC
  pe_1_fail,
  //! PE mate 1 read
  pe_2,
  //! PE mate 1 read that failed QC
  pe_2_fail,
  //! PE mate 1 read for which the mate 2 read failed QC
  singleton_1,
  //! PE mate 2 read for which the mate 1 read failed QC
  singleton_2,
  //! Merged PE reads
  merged,
  //! Merged PE reads that failed QC
  merged_fail,
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
  //! Unaligned BAM (Binary Alignment Map) records
  bam,
  //! Uncompressed, unaligned BAM (Binary Alignment Map) records
  ubam,
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

/** The user-specified or derived orientation of barcodes */
enum class barcode_orientation
{
  //! The orientation of barcode is unspecified
  unspecified = 0,
  //! The barcode is in forward orientation
  forward,
  //! The barcode is in reverse orientation
  reverse,
};

/** The orientation of barcodes in user-provided tables */
enum class barcode_table_orientation
{
  //! The orientation of barcodes in the table is unspecified
  unspecified = 0,
  //! All barcodes in the table are in forward orientation
  forward,
  //! All barcodes in the table are in reverse orientation
  reverse,
  //! The user has provided an oritentation per barcode
  explicit_,
};

//! Path used to indicate that a file is not needed
const std::string_view DEV_NULL = "/dev/null";
//! Path used to indicate that data should be read from STDIN
const std::string_view DEV_STDIN = "/dev/stdout";
//! Path used to indicate that data should be written to STDOUT
const std::string_view DEV_STDOUT = "/dev/stdout";
//! Path used to indicate that data should be written to STDOUT
const std::string_view DEV_STDERR = "/dev/stderr";
//! Filename indicating that data should be read from stdin/written to stdout
const std::string_view DEV_PIPE = "-";

} // namespace adapterremoval
