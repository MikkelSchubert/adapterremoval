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

#include "sequence.hpp" // for sequence_pair_vec
#include <array>        // for array
#include <cstddef>      // for size_t
#include <cstdint>      // for int32_t
#include <vector>       // for vector

namespace adapterremoval {

class sample_set;
class fastq;
struct next_subsequence;

/**
 * Struct representing node in quad-tree; children are referenced using the
 * corresponding index in the vector representing the tree; -1 is used to
 * represent unassigned children.
 */
struct demultiplexer_node
{
  demultiplexer_node();

  std::array<int32_t, 4> children;
  int value;
};

using demux_node_vec = std::vector<demultiplexer_node>;

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

  barcode_table(const sequence_pair_vec& barcodes,
                size_t max_mm,
                size_t max_mm_r1,
                size_t max_mm_r2);

  int identify(const fastq& read_r1) const;
  int identify(const fastq& read_r1, const fastq& read_r2) const;

  [[nodiscard]] size_t length_1() const { return m_barcode_1_len; }

  [[nodiscard]] size_t length_2() const { return m_barcode_2_len; }

  constexpr static const int no_match = -1;
  constexpr static const int ambiguous = -2;

private:
  struct candidate
  {
    explicit candidate(int barcode = barcode_table::no_match,
                       size_t mismatches = -1);

    int barcode;
    size_t mismatches;
  };

  candidate lookup(const char* seq,
                   int parent,
                   size_t max_global_mismatches,
                   const next_subsequence* next) const;

  candidate lookup_with_mm(const char* seq,
                           int parent,
                           size_t max_global_mismatches,
                           size_t max_local_mismatches,
                           const next_subsequence* next) const;

  demux_node_vec m_nodes{};
  size_t m_max_mismatches = 0;
  size_t m_max_mismatches_r1 = 0;
  size_t m_max_mismatches_r2 = 0;
  size_t m_barcode_1_len = 0;
  size_t m_barcode_2_len = 0;
};

} // namespace adapterremoval
