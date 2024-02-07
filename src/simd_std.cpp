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

bool
compare_subsequences_std(size_t& n_mismatches,
                         size_t& n_ambiguous,
                         const char* seq_1,
                         const char* seq_2,
                         size_t length,
                         size_t max_penalty)
{
  size_t penalty = 2 * n_mismatches + n_ambiguous;
  for (; length; --length) {
    const char nt_1 = *seq_1++;
    const char nt_2 = *seq_2++;

    if (nt_1 == 'N' || nt_2 == 'N') {
      n_ambiguous++;
      penalty++;
    } else if (nt_1 != nt_2) {
      n_mismatches++;
      penalty += 2;
    }

    if (penalty > max_penalty) {
      return false;
    }
  }

  return true;
}

} // namespace simd

} // namespace adapterremoval
