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

#include "debug.hpp"
#include "fastq.hpp"
#include "statistics.hpp"

fastq_statistics::fastq_statistics()
  : m_length_dist()
  , m_quality_dist()
  , m_uncalled_pos()
  , m_uncalled_quality_pos()
  , m_called_pos(4)
  , m_quality_pos(4)
{}

void
fastq_statistics::process(const fastq& read)
{
  m_length_dist.inc(read.length());

  const std::string& sequence = read.sequence();
  const std::string& qualities = read.qualities();

  for (size_t i = 0; i < sequence.length(); ++i) {
    const auto nuc = sequence.at(i);
    if (nuc == 'N') {
      m_uncalled_pos.inc(i);
      m_uncalled_quality_pos.inc(i, qualities.at(i));
    } else {
      const auto nuc_i = ACGT_TO_IDX(nuc);

      m_called_pos.at(nuc_i).inc(i);
      m_quality_pos.at(nuc_i).inc(i, qualities.at(i) - PHRED_OFFSET_33);
    }

    m_quality_dist.inc(qualities.at(i) - PHRED_OFFSET_33);
  }
}

fastq_statistics&
fastq_statistics::operator+=(const fastq_statistics& other)
{
  m_length_dist += other.m_length_dist;
  m_quality_dist += other.m_quality_dist;
  m_uncalled_pos += other.m_uncalled_pos;
  m_uncalled_quality_pos += other.m_uncalled_quality_pos;

  for (size_t i = 0; i < 4; ++i) {
    m_called_pos.at(i) += other.m_called_pos.at(i);
    m_quality_pos.at(i) += other.m_quality_pos.at(i);
  }

  return *this;
}

trimming_statistics::trimming_statistics()
  : read_1()
  , read_2()
  , merged()
  , discarded()
  , number_of_merged_reads()
  , number_of_reads_with_adapter()
{}

trimming_statistics&
trimming_statistics::operator+=(const trimming_statistics& other)
{
  read_1 += other.read_1;
  read_2 += other.read_2;
  merged += other.merged;
  discarded += other.discarded;

  number_of_merged_reads += other.number_of_merged_reads;
  number_of_reads_with_adapter += other.number_of_reads_with_adapter;

  return *this;
}

demultiplexing_statistics::demultiplexing_statistics()
  : barcodes()
  , unidentified(0)
  , ambiguous(0)
{}

bool
demultiplexing_statistics::empty() const
{
  return barcodes.empty();
}

void
demultiplexing_statistics::resize(size_t n)
{
  barcodes.resize(n);
}

size_t
demultiplexing_statistics::total() const
{
  size_t total = unidentified + ambiguous;
  for (size_t i = 0; i < barcodes.size(); ++i) {
    total += barcodes.at(i);
  }

  return total;
}
