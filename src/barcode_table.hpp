// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "sequence.hpp" // for sequence_pair_vec
#include <array>        // for array
#include <cstddef>      // for size_t
#include <cstdint>      // for int32_t
#include <vector>       // for vector

namespace adapterremoval {

class sample_set;
class fastq;
struct barcode_match;
struct next_subsequence;

/** The sample/barcode IDs associated with a set of sequences */
struct barcode_key
{
  constexpr static const int32_t unidentified = -1;
  constexpr static const int32_t ambiguous = -2;

  //! Sample identified by 1 or more barcodes
  int32_t sample = unidentified;
  //! The specific barcode(s) used to identify this sample
  int32_t barcode = unidentified;

  bool operator==(const barcode_key& other) const
  {
    return sample == other.sample && barcode == other.barcode;
  }

  bool operator<(const barcode_key& other) const
  {
    return sample < other.sample ||
           (sample == other.sample && barcode < other.barcode);
  }
};

using barcode_pair = std::pair<sequence_pair, barcode_key>;
using barcode_vec = std::vector<barcode_pair>;

/**
 * Struct representing node in quad-tree; children are referenced using the
 * corresponding index in the vector representing the tree; -1 is used to
 * represent unassigned children.
 */
struct barcode_node
{
  barcode_node();

  std::array<int32_t, 4> children;
  barcode_key key;
};

/**
 *
 */
class barcode_table
{
public:
  barcode_table(const sample_set& samples,
                size_t max_mm,
                size_t max_mm_r1,
                size_t max_mm_r2);

  [[nodiscard]] barcode_key identify(const fastq& read_r1) const;
  [[nodiscard]] barcode_key identify(const fastq& read_r1,
                                     const fastq& read_r2) const;

  /** Returns the length of barcodes in read 1  */
  [[nodiscard]] size_t length_1() const { return m_barcode_1_len; }

  /** Returns the length of barcodes in read 2  */
  [[nodiscard]] size_t length_2() const { return m_barcode_2_len; }

private:
  barcode_match lookup(const char* seq,
                       int32_t parent,
                       size_t max_global_mismatches,
                       const char* next) const;

  barcode_match lookup_with_mm(const char* seq,
                               int32_t parent,
                               size_t max_global_mismatches,
                               size_t max_local_mismatches,
                               const char* next) const;

  std::vector<barcode_node> m_nodes{};
  size_t m_max_mismatches = 0;
  size_t m_max_mismatches_r1 = 0;
  size_t m_max_mismatches_r2 = 0;
  size_t m_barcode_1_len = 0;
  size_t m_barcode_2_len = 0;
};

} // namespace adapterremoval
