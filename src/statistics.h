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
#ifndef STATISTICS_H
#define STATISTICS_H

#include <cstdlib>
#include <vector>


struct statistics
{
    statistics()
      : number_of_full_length_collapsed(0)
      , number_of_truncated_collapsed(0)
      , total_number_of_nucleotides(0)
      , total_number_of_good_reads(0)
      , number_of_reads_with_adapter()
      , number_of_barcodes_trimmed()
      , unaligned_reads(0)
      , well_aligned_reads(0)
      , poorly_aligned_reads(0)
      , keep1(0)
      , discard1(0)
      , keep2(0)
      , discard2(0)
      , records(0)
    {
    }

    size_t number_of_full_length_collapsed;
    size_t number_of_truncated_collapsed;
    size_t total_number_of_nucleotides;
    size_t total_number_of_good_reads;

    std::vector<size_t> number_of_reads_with_adapter;
    std::vector<size_t> number_of_barcodes_trimmed;

    size_t unaligned_reads;
    size_t well_aligned_reads;
    size_t poorly_aligned_reads;
    size_t keep1;
    size_t discard1;
    size_t keep2;
    size_t discard2;

    size_t records;

    enum read_type {
        rt_mate_1 = 0,
        rt_mate_2,
        rt_singleton,
        rt_collapsed,
        rt_collapsed_truncated,
        rt_discarded
    };

    void inc_length_count(read_type type, size_t length) {
        if (length >= read_lengths.size()) {
            const size_t nfields = static_cast<size_t>(rt_discarded) + 1;
            read_lengths.resize(length + 1, std::vector<size_t>(nfields));
        }

        ++read_lengths.at(length).at(static_cast<size_t>(type));
    }

    std::vector<std::vector<size_t> > read_lengths;
};

#endif
