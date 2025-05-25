// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "simd.hpp"    // declarations
#include <algorithm>   // for min
#include <bitset>      // for bitset
#include <cstddef>     // for size_t
#include <emmintrin.h> // for _mm_cmpeq_epi8, __m128i, _mm_loadu_s...

namespace adapterremoval {

namespace simd {

namespace {

/** Counts the number of masked bytes **/
inline auto
count_masked_sse2(__m128i value)
{
  // Generate 16 bit mask from most significant bits of each byte and count bits
  return std::bitset<16>(_mm_movemask_epi8(value)).count();
}

} // namespace

bool
compare_subsequences_sse2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty)
{
  //! Mask of all Ns
  const __m128i n_mask = _mm_set1_epi8('N');

  // The sequence is assumed to be padded with 15 'N's not included in `length`.
  while (length) {
    const auto s1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_1));
    const auto s2 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_2));

    // Sets 0xFF for every byte where one or both nts is N
    const auto ns_mask =
      _mm_or_si128(_mm_cmpeq_epi8(s1, n_mask), _mm_cmpeq_epi8(s2, n_mask));

    // Sets 0xFF for every byte where bytes are equal or N
    const auto eq_mask = _mm_or_si128(_mm_cmpeq_epi8(s1, s2), ns_mask);

    n_mismatches += 16 - count_masked_sse2(eq_mask);
    if (2 * n_mismatches + n_ambiguous > max_penalty) {
      return false;
    }

    // Fragment length without 'N' padding
    size_t unpadded_length = std::min<size_t>(16, length);

    // Early termination is almost always due to mismatches, so updating the
    // number of Ns after the above check saves time in the common case.
    n_ambiguous += count_masked_sse2(ns_mask) - (16 - unpadded_length);

    seq_1 += 16;
    seq_2 += 16;
    length -= unpadded_length;
  }

  return true;
}

} // namespace simd

} // namespace adapterremoval
