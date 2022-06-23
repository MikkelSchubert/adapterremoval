/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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

#include <algorithm> // for max
#include <array>     // for array
#include <exception> // for exception
#include <stddef.h>  // for size_t
#include <string>    // for string
#include <vector>    // for vector

#include "fastq.hpp" // for fastq_pair_vec

namespace adapterremoval {

struct next_subsequence;

/** Exception raised for FASTQ parsing and validation errors. */
class barcode_error : public std::exception
{
public:
  barcode_error(const std::string& message);
  barcode_error(const barcode_error& error);

  virtual ~barcode_error() override;

  /** Returns error message; string is owned by exception. */
  virtual const char* what() const noexcept override;

private:
  //! Error message associated with exception.
  std::string m_message;
};

/**
 * Struct representing node in quad-tree; children are referenced using the
 * corresponding indice in the vector representing the tree; -1 is used to
 * represent unassigned children.
 */
struct demultiplexer_node
{
  demultiplexer_node();

  std::array<int, 4> children;
  int value;
};

typedef std::vector<demultiplexer_node> demux_node_vec;

/**
 *
 */
class barcode_table
{
public:
  barcode_table(const fastq_pair_vec& barcodes,
                size_t max_mm,
                size_t max_mm_r1,
                size_t max_mm_r2);

  int identify(const fastq& read_r1) const;
  int identify(const fastq& read_r1, const fastq& read_r2) const;

  static const int no_match = -1;
  static const int ambigious = -2;

protected:
  struct candidate
  {
    explicit candidate(int barcode = barcode_table::no_match,
                       size_t mismatches = -1);

    int barcode;
    size_t mismatches;
  };

  candidate lookup(const char* seq,
                   int parent,
                   const size_t max_global_mismatches,
                   const next_subsequence* next) const;

  candidate lookup_with_mm(const char* seq,
                           int parent,
                           const size_t max_global_mismatches,
                           const size_t max_local_mismatches,
                           const next_subsequence* next) const;

  demux_node_vec m_nodes;
  size_t m_max_mismatches;
  size_t m_max_mismatches_r1;
  size_t m_max_mismatches_r2;
  size_t m_barcode_1_len;
  size_t m_barcode_2_len;
};

} // namespace adapterremoval
