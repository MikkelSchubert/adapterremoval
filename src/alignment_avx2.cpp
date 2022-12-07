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
#include "alignment_common.hpp" // header
#include <bitset>               // for bitset
#include <immintrin.h>          // for _mm256_set1_epi8, __m256i, _mm256_lo...

namespace adapterremoval {

namespace {

/** Counts the number of masked bytes **/
inline auto
count_masked_avx2(__m256i value)
{
  // Generate 32 bit mask from most significant bits of each byte and count bits
  return std::bitset<32>(_mm256_movemask_epi8(value)).count();
}

} // namespace

bool
compare_subsequences_avx2(alignment_info& current,
                          const char* seq_1,
                          const char* seq_2,
                          const size_t max_mismatches,
                          int length)
{
  //! Mask of all Ns
  const __m256i n_mask = _mm256_set1_epi8('N');

  while (length >= 32) {
    const auto s1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(seq_1));
    const auto s2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(seq_2));

    // Sets 0xFF for every byte where one or both nts is N
    const __m256i ns_mask = _mm256_or_si256(_mm256_cmpeq_epi8(s1, n_mask),
                                            _mm256_cmpeq_epi8(s2, n_mask));

    // Sets 0xFF for every byte where bytes are equal or N
    const __m256i eq_mask = _mm256_or_si256(_mm256_cmpeq_epi8(s1, s2), ns_mask);

    current.n_ambiguous += count_masked_avx2(ns_mask);
    current.n_mismatches += 32 - count_masked_avx2(eq_mask);
    if (current.n_mismatches > max_mismatches) {
      return false;
    }

    seq_1 += 32;
    seq_2 += 32;
    length -= 32;
  }

  return compare_subsequences_sse2(
    current, seq_1, seq_2, max_mismatches, length);
}

} // namespace adapterremoval
