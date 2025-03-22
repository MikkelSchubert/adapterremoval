// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
// #include "simd.hpp"    // declarations
#include <algorithm>  // for min
#include <arm_neon.h> // for vdupq_n_u8, vld1q_u8, vorrq_u8, ...
#include <cstddef>    // for size_t

namespace adapterremoval {

namespace simd {

namespace {

/** Counts the number of masked bytes **/
inline auto
count_masked_neon(uint8x16_t value)
{
  // Extract one bit per uint8_t and sum across vector
  return vaddvq_u8(vandq_u8(vdupq_n_u8(0x01), value));
}

} // namespace

bool
compare_subsequences_neon(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty)
{
  //! Mask of all Ns
  const auto n_mask = vdupq_n_u8('N');

  // The sequence is assumed to be padded with 15 'N's not included in `length`.
  while (length) {
    auto s1 = vld1q_u8(reinterpret_cast<const uint8_t*>(seq_1));
    auto s2 = vld1q_u8(reinterpret_cast<const uint8_t*>(seq_2));

    // Sets 0xFF for every byte where one or both nts is N
    const auto ns_mask = vorrq_u8(vceqq_u8(s1, n_mask), vceqq_u8(s2, n_mask));

    // Sets 0xFF for every byte where bytes are equal or N
    const auto eq_mask = vorrq_u8(vceqq_u8(s1, s2), ns_mask);

    n_mismatches += 16 - count_masked_neon(eq_mask);
    if (2 * n_mismatches + n_ambiguous > max_penalty) {
      return false;
    }

    // Fragment length without 'N' padding
    size_t unpadded_length = std::min<size_t>(16, length);

    // Early termination is almost always due to mismatches, so updating the
    // number of Ns after the above check saves time in the common case.
    n_ambiguous += count_masked_neon(ns_mask) - (16 - unpadded_length);

    seq_1 += 16;
    seq_2 += 16;
    length -= unpadded_length;
  }

  return true;
}

} // namespace simd

} // namespace adapterremoval
