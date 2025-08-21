// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "simd.hpp"   // declarations
#include <algorithm>  //
#include <arm_neon.h> // for vdupq_n_u8, vld1q_u8, vorrq_u8, ...
#include <cstddef>    // for size_t

namespace adapterremoval::simd {

namespace {

//! The number of bytes in a NEON uint8x16_t
const size_t NEON_BLOCK_SIZE = 16;
//! The number of uint8x16_t to process between before summing up stats
const size_t NEON_BATCH_SIZE = NEON_BLOCK_SIZE * 8;

} // namespace

bool
compare_subsequences_neon(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          const size_t length,
                          size_t max_penalty)
{
  //! Mask of all Ns
  const auto n_mask = vdupq_n_u8('N');

  // The sequence is assumed to be padded with 15 'N's not included in `length`.
  auto i_length = static_cast<int64_t>(length);
  while (i_length > 0) {
    auto a_matches = vdupq_n_u8(0);
    auto a_ambiguous = vdupq_n_u8(0);

    size_t chunk = 0;
    for (; chunk < NEON_BATCH_SIZE && i_length > 0; chunk += NEON_BLOCK_SIZE) {
      auto s1 = vld1q_u8(reinterpret_cast<const uint8_t*>(seq_1));
      auto s2 = vld1q_u8(reinterpret_cast<const uint8_t*>(seq_2));

      // Sets 0xFF for every byte where one or both nts is N
      const auto ns_mask = vorrq_u8(vceqq_u8(s1, n_mask), vceqq_u8(s2, n_mask));

      // Sets 0xFF for every byte where bytes are equal or N
      const auto eq_mask = vorrq_u8(vceqq_u8(s1, s2), ns_mask);

      a_matches = vaddq_u8(a_matches, vshrq_n_u8(eq_mask, 7));
      a_ambiguous = vaddq_u8(a_ambiguous, vshrq_n_u8(ns_mask, 7));

      seq_1 += NEON_BLOCK_SIZE;
      seq_2 += NEON_BLOCK_SIZE;
      i_length -= NEON_BLOCK_SIZE;
    }

    n_mismatches += chunk - vaddvq_u8(a_matches);
    n_ambiguous += vaddvq_u8(a_ambiguous) + std::min<int64_t>(0, i_length);
    if (2 * n_mismatches + n_ambiguous > max_penalty) {
      return false;
    }
  }

  return true;
}

} // namespace adapterremoval::simd
