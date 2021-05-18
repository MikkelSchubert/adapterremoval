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
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>

#include "statistics.hpp"

namespace ar {

statistics::statistics()
  : number_of_collapsed(0)
  , total_number_of_nucleotides(0)
  , total_number_of_good_reads(0)
  , number_of_reads_with_adapter()
  , unaligned_reads(0)
  , well_aligned_reads(0)
  , poorly_aligned_reads(0)
  , keep1(0)
  , discard1(0)
  , keep2(0)
  , discard2(0)
  , records(0)
  , read_lengths()
{}

void
statistics::inc_length_count(read_type type, size_t length)
{
  if (length >= read_lengths.size()) {
    read_lengths.resize(
      length + 1, std::vector<size_t>(static_cast<size_t>(read_type::max)));
  }

  ++read_lengths.at(length).at(static_cast<size_t>(type));
}

statistics&
statistics::operator+=(const statistics& other)
{
  number_of_collapsed += other.number_of_collapsed;
  total_number_of_nucleotides += other.total_number_of_nucleotides;
  total_number_of_good_reads += other.total_number_of_good_reads;

  unaligned_reads += other.unaligned_reads;
  well_aligned_reads += other.well_aligned_reads;
  poorly_aligned_reads += other.poorly_aligned_reads;
  keep1 += other.keep1;
  discard1 += other.discard1;
  keep2 += other.keep2;
  discard2 += other.discard2;

  records += other.records;

  merge_vectors(number_of_reads_with_adapter,
                other.number_of_reads_with_adapter);
  merge_sub_vectors(read_lengths, other.read_lengths);

  return *this;
}

demux_statistics::demux_statistics(const size_t n_barcodes)
  : barcodes(n_barcodes)
  , unidentified(0)
  , ambiguous(0)
{}

size_t
demux_statistics::total() const
{
  size_t total = unidentified + ambiguous;
  for (size_t i = 0; i < barcodes.size(); ++i) {
    total += barcodes.at(i);
  }

  return total;
}

} // namespace ar
