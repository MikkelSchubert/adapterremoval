/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
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

#include <cstdlib>
#include <vector>

#include "commontypes.hpp"
#include "vecutils.hpp"

class demultiplex_reads;
class userconfig;
class reads_processor;

/** Object used to collect summary statistics for trimming and other tasks. */
struct statistics
{
  statistics();

  //! Number of collapsed reads
  size_t number_of_collapsed;
  //! Total number of nucleotides left after trimming, collapsing, filtering
  size_t total_number_of_nucleotides;
  //! Total number of reads left after trimming, collapsing, filtering
  size_t total_number_of_good_reads;

  //! Number of reads / pairs with adapters trimmed
  std::vector<size_t> number_of_reads_with_adapter;

  //! Number of unaligned reads; not enough bases, too many mismatches, etc.
  size_t unaligned_reads;
  //! Number of alignments with enough aligned bases, not too many mismatches
  size_t well_aligned_reads;
  //! Number of alignments with a zero or lower score
  size_t poorly_aligned_reads;
  //! Number of retained mate 1 reads
  size_t keep1;
  //! Number of discarded mate 1 reads
  size_t discard1;
  //! Number of retained mate 2 reads
  size_t keep2;
  //! Number of discarded mate 2 reads
  size_t discard2;
  //! Total number of reads / pairs processed
  size_t records;

  /** Increment the number of reads with of a given type / length. */
  void inc_length_count(read_type type, size_t length);

  //! Per read-type length distributions of reads
  std::vector<std::vector<size_t>> read_lengths;

  /** Combine statistics objects, e.g. those used by different threads. */
  statistics& operator+=(const statistics& other);
};

/** Object used to collect summary statistics for demultiplexing. */
struct demux_statistics
{
  demux_statistics();

  bool empty() const;
  void resize(size_t n);

  size_t total() const;

  //! Number of reads / pairs identified for a given barcode / pair of
  //! barcodes
  std::vector<size_t> barcodes;
  //! Number of reads / pairs with no hits
  size_t unidentified;
  //! Number of reads / pairs with no single best hit
  size_t ambiguous;
};
