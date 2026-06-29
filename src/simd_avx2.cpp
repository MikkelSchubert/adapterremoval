// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "simd.hpp"    // declarations
#include <algorithm>   // for min
#include <bitset>      // for bitset
#include <cstddef>     // for size_t
#include <immintrin.h> // for _mm256_set1_epi8, __m256i, _mm256_lo...

namespace adapterremoval::simd {

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
compare_subsequences_avx2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty)
{
  //! Mask of all Ns
  const __m256i n_mask = _mm256_set1_epi8('N');

  // The sequence is assumed to be padded with 31 'N's not included in `length`.
  while (length) {
    const auto s1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(seq_1));
    const auto s2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(seq_2));

    // Sets 0xFF for every byte where one or both nts is N
    const auto ns_mask = _mm256_or_si256(_mm256_cmpeq_epi8(s1, n_mask),
                                         _mm256_cmpeq_epi8(s2, n_mask));
    const auto n_count = count_masked_avx2(ns_mask);

    if (n_count) {
      // Sets 0xFF for every byte where bytes are equal or N
      const __m256i eq_mask =
        _mm256_or_si256(_mm256_cmpeq_epi8(s1, s2), ns_mask);

      // Early termination is almost always due to mismatches, so updating the
      // number of Ns after the penalty check saves time in the common case.
      n_mismatches += 32 - count_masked_avx2(eq_mask);
      if ((2 * n_mismatches) + n_ambiguous > max_penalty) {
        return false;
      }

      const size_t unpadded_length = std::min<size_t>(32, length);
      n_ambiguous += n_count - (32 - unpadded_length);

      seq_1 += unpadded_length;
      seq_2 += unpadded_length;
      length -= unpadded_length;
    } else {
      n_mismatches += 32 - count_masked_avx2(_mm256_cmpeq_epi8(s1, s2));
      if ((2 * n_mismatches) + n_ambiguous > max_penalty) {
        return false;
      }

      seq_1 += 32;
      seq_2 += 32;
      length -= 32;
    }
  }

  return (2 * n_mismatches) + n_ambiguous <= max_penalty;
}

} // namespace adapterremoval::simd
