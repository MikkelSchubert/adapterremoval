/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "simd.hpp"    // declarations
#include <bitset>      // for bitset
#include <cstddef>     // for size_t
#include <immintrin.h> // for _mm512_set1_epi8, __m512i, _mm512_lo...

#ifdef AR_SUPPORTS_AVX512

namespace adapterremoval {

namespace simd {

namespace {

/** Counts the number of masked bytes **/
inline auto
count_masked_avx512(__mmask64 value)
{
  // Generate 64 bit mask from most significant bits of each byte and count bits
  return std::bitset<64>(value).count();
}

} // namespace

bool
compare_subsequences_avx512(size_t& n_mismatches,
                            size_t& n_ambiguous,
                            const char* seq_1,
                            const char* seq_2,
                            size_t length,
                            size_t max_penalty)
{
  //! Mask of all Ns
  const auto n_mask = _mm512_set1_epi8('N');

  while (length) {
    __m512i s1;
    __m512i s2;

    if (length < 64) {
      const auto mask = 0xFFFFFFFFFFFFFFFFLLU >> (64 - length);
      s1 = _mm512_maskz_loadu_epi8(mask, seq_1);
      s2 = _mm512_maskz_loadu_epi8(mask, seq_2);

      length = 0;
    } else {
      s1 = _mm512_loadu_epi8(seq_1);
      s2 = _mm512_loadu_epi8(seq_2);

      seq_1 += 64;
      seq_2 += 64;
      length -= 64;
    }

    // Sets 0xFF for every byte where one or both nts is N
    const auto ns_mask =
      _mm512_cmpeq_epu8_mask(s1, n_mask) | _mm512_cmpeq_epu8_mask(s2, n_mask);

    // Sets 0xFF for every byte where bytes are equal or N
    const auto eq_mask = _mm512_cmpeq_epu8_mask(s1, s2) | ns_mask;

    n_mismatches += 64 - count_masked_avx512(eq_mask);
    if (2 * n_mismatches + n_ambiguous > max_penalty) {
      return false;
    }

    // Early termination is almost always due to mismatches, so updating the
    // number of Ns after the above check saves time in the common case.
    n_ambiguous += count_masked_avx512(ns_mask);
  }

  return true;
}

} // namespace simd

} // namespace adapterremoval

#endif
