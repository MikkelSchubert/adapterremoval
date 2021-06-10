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
#include "counts.hpp"
#include "vecutils.hpp"

class demultiplex_reads;
class fastq;
class reads_processor;
class userconfig;

/** Class used to collect statistics about pre/post-processed FASTQ reads. */
class fastq_statistics
{
public:
  fastq_statistics();

  void process(const fastq& read);

  inline const counts& length_dist() const { return m_length_dist; }
  inline const counts& quality_dist() const { return m_quality_dist; }
  inline const counts& uncalled_pos() const { return m_uncalled_pos; }
  inline const counts& uncalled_quality_pos() const
  {
    return m_uncalled_quality_pos;
  }
  inline const counts& nucleotides_pos(char nuc) const
  {
    return m_called_pos.at(ACGT_TO_IDX(nuc));
  }
  inline const counts& qualities_pos(char nuc) const
  {
    return m_quality_pos.at(ACGT_TO_IDX(nuc));
  }

  /** Sum statistics, e.g. those used by different threads. */
  fastq_statistics& operator+=(const fastq_statistics& other);

private:
  /** Length distribution. */
  counts m_length_dist;
  /** Quality distribution. */
  counts m_quality_dist;
  /** Count of uncalled bases per position. */
  counts m_uncalled_pos;
  /** Sum of qualities of Ns per position. */
  counts m_uncalled_quality_pos;
  /** Count of A/C/G/T per position; indexed using ACGT_TO_IDX. */
  std::vector<counts> m_called_pos;
  /** Sum of qualities of A/C/G/Ts per position; indexed using ACGT_TO_IDX. */
  std::vector<counts> m_quality_pos;
};

/** Object used to collect summary statistics for trimming. */
struct trimming_statistics
{
  trimming_statistics();

  //! Statistics for first reads
  fastq_statistics read_1;
  //! Statistics for second reads
  fastq_statistics read_2;
  //! Statistics for discarded reads
  fastq_statistics merged;
  //! Statistics for discarded reads
  fastq_statistics discarded;

  //! Number of reads merged (incl. discarded reads)
  size_t number_of_merged_reads;
  //! Number of reads / pairs with adapters trimmed
  counts number_of_reads_with_adapter;

  /** Combine statistics objects, e.g. those used by different threads. */
  trimming_statistics& operator+=(const trimming_statistics& other);
};

/** Object used to collect summary statistics for demultiplexing. */
struct demultiplexing_statistics
{
  demultiplexing_statistics();

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

// FIXME: Rename to ... something better
struct ar_statistics
{
  inline ar_statistics()
    : input_1()
    , input_2()
    , demultiplexing()
    , trimming()
  {}

  fastq_statistics input_1;
  fastq_statistics input_2;

  demultiplexing_statistics demultiplexing;
  std::vector<trimming_statistics> trimming;
};
