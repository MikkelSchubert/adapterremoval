/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "simd.hpp" // for size_t, compare_subsequences_std
#include <cstddef>  // for size_t

namespace adapterremoval {

namespace simd {

namespace {

//! Number of nucleotides to compare per (unrolled) loop
const size_t LOOP_COUNT = 8;

} // namespace

bool
compare_subsequences_std(size_t& out_mismatches,
                         size_t& out_ambiguous,
                         const char* seq_1,
                         const char* seq_2,
                         size_t length,
                         size_t max_penalty)
{
  // Local accumulators are used since references prevent loop unrolling
  size_t n_mismatches = out_mismatches;
  size_t n_ambiguous = out_ambiguous;

  while (length >= LOOP_COUNT) {
    // A constant loop count enables loop-unrolling for increased throughput
    for (size_t i = 0; i < LOOP_COUNT; ++i) {
      const char nt_1 = *seq_1++;
      const char nt_2 = *seq_2++;

      if (nt_1 == 'N' || nt_2 == 'N') {
        n_ambiguous++;
      } else if (nt_1 != nt_2) {
        n_mismatches++;
      }
    }

    if (2 * n_mismatches + n_ambiguous > max_penalty) {
      return false;
    }

    length -= LOOP_COUNT;
  }

  for (; length; --length) {
    const char nt_1 = *seq_1++;
    const char nt_2 = *seq_2++;

    if (nt_1 == 'N' || nt_2 == 'N') {
      n_ambiguous++;
    } else if (nt_1 != nt_2) {
      n_mismatches++;
    }
  }

  out_ambiguous = n_ambiguous;
  out_mismatches = n_mismatches;

  return 2 * n_mismatches + n_ambiguous <= max_penalty;
}

} // namespace simd

} // namespace adapterremoval
