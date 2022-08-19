/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "alignment.hpp" // header
#include <bitset>        // for bitset
#include <immintrin.h>   // for _mm256_set1_epi8, __m256i, _mm256_lo...

namespace adapterremoval {

/** Counts the number of masked bytes **/
inline size_t
count_masked_avx2(const __m256i& value)
{
  // Generate 16 bit mask from most significant bits of each byte and count bits
  return std::bitset<32>(_mm256_movemask_epi8(value)).count();
}

bool
compare_subsequences_avx2(const alignment_info& best,
                          alignment_info& current,
                          const char*& seq_1_ptr,
                          const char*& seq_2_ptr,
                          double mismatch_threshold,
                          int& remaining_bases)
{
  //! Mask of all Ns
  const __m256i n_mask = _mm256_set1_epi8('N');

  while (remaining_bases >= 32 && current.score >= best.score) {
    const __m256i s1 =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(seq_1_ptr));
    const __m256i s2 =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(seq_2_ptr));

    // Sets 0xFF for every byte where one or both nts is N
    const __m256i ns_mask = _mm256_or_si256(_mm256_cmpeq_epi8(s1, n_mask),
                                            _mm256_cmpeq_epi8(s2, n_mask));

    // Sets 0xFF for every byte where bytes are equal or N
    const __m256i eq_mask = _mm256_or_si256(_mm256_cmpeq_epi8(s1, s2), ns_mask);

    current.n_ambiguous += count_masked_avx2(ns_mask);
    current.n_mismatches += 32 - count_masked_avx2(eq_mask);
    if (current.n_mismatches >
        (current.length - current.n_ambiguous) * mismatch_threshold) {
      return false;
    }

    // Matches count for 1, Ns for 0, and mismatches for -1
    current.score =
      current.length - current.n_ambiguous - (current.n_mismatches * 2);

    seq_1_ptr += 32;
    seq_2_ptr += 32;
    remaining_bases -= 32;
  }

  return true;
}

} // namespace adapterremoval
