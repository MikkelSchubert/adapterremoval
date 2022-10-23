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
#include "alignment.hpp" // header
#include <bitset>        // for bitset
#include <emmintrin.h>   // for _mm_cmpeq_epi8, __m128i, _mm_loadu_s...

namespace adapterremoval {

/** Counts the number of masked bytes **/
inline size_t
count_masked_sse2(const __m128i& value)
{
  // Generate 16 bit mask from most significant bits of each byte and count bits
  return std::bitset<16>(_mm_movemask_epi8(value)).count();
}

bool
compare_subsequences_sse2(const alignment_info& best,
                          alignment_info& current,
                          const char*& seq_1_ptr,
                          const char*& seq_2_ptr,
                          double mismatch_threshold,
                          int& remaining_bases)
{
  //! Mask of all Ns
  const __m128i n_mask = _mm_set1_epi8('N');

  while (remaining_bases >= 16 && current.score >= best.score) {
    const __m128i s1 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_1_ptr));
    const __m128i s2 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_2_ptr));

    // Sets 0xFF for every byte where one or both nts is N
    const __m128i ns_mask =
      _mm_or_si128(_mm_cmpeq_epi8(s1, n_mask), _mm_cmpeq_epi8(s2, n_mask));

    // Sets 0xFF for every byte where bytes are equal or N
    const __m128i eq_mask = _mm_or_si128(_mm_cmpeq_epi8(s1, s2), ns_mask);

    current.n_ambiguous += count_masked_sse2(ns_mask);
    current.n_mismatches += 16 - count_masked_sse2(eq_mask);
    if (current.n_mismatches >
        (current.length - current.n_ambiguous) * mismatch_threshold) {
      return false;
    }

    // Matches count for 1, Ns for 0, and mismatches for -1
    current.score =
      current.length - current.n_ambiguous - (current.n_mismatches * 2);

    seq_1_ptr += 16;
    seq_2_ptr += 16;
    remaining_bases -= 16;
  }

  return true;
}

} // namespace adapterremoval
