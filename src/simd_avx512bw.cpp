// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "debug.hpp"   // for AR_UNLIKELY
#include "simd.hpp"    // declarations
#include <bitset>      // for bitset
#include <cstddef>     // for size_t
#include <immintrin.h> // for _mm512_set1_epi8, __m512i, _mm512_lo...

namespace adapterremoval::simd {

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
    const __m512i s1 = _mm512_loadu_epi8(seq_1);
    const __m512i s2 = _mm512_loadu_epi8(seq_2);

    // Sets 1 for every bit where one or both nucleotides is N
    const auto ns_mask =
      _mm512_cmpeq_epu8_mask(s1, n_mask) | _mm512_cmpeq_epu8_mask(s2, n_mask);

    // Sets 1 for every bit where nucleotides are equal or N
    n_mismatches +=
      count_masked_avx512(~(_mm512_cmpeq_epu8_mask(s1, s2) | ns_mask));

    // Early termination is almost always due to mismatches, so updating the
    // number of Ns after the penalty check saves time in the common case.
    if ((2 * n_mismatches) + n_ambiguous > max_penalty) {
      return false;
    }

    if (ns_mask) {
      // Either read contains an 'N' or we've reached the padding
      const size_t unpadded_length = std::min<size_t>(64, length);

      n_ambiguous += count_masked_avx512(ns_mask) - (64 - unpadded_length);

      seq_1 += unpadded_length;
      seq_2 += unpadded_length;
      length -= unpadded_length;
    } else {
      seq_1 += 64;
      seq_2 += 64;
      length -= 64;
    }
  }

  return (2 * n_mismatches) + n_ambiguous <= max_penalty;
}

} // namespace adapterremoval::simd
