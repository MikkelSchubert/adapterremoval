// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "pragmas.hpp" // for AR_UNROLL
#include "simd.hpp"    // for size_t, compare_subsequences_std
#include <cstddef>     // for size_t

namespace adapterremoval::simd {

bool
compare_subsequences_none(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty)
{
  // Local accumulators are used since references prevent loop unrolling
  size_t tmp_mismatches = n_mismatches;
  size_t tmp_ambiguous = n_ambiguous;

  AR_UNROLL(8)
  for (; length; --length) {
    const char nt_1 = *seq_1++;
    const char nt_2 = *seq_2++;

    if (nt_1 == 'N' || nt_2 == 'N') {
      tmp_ambiguous++;
    } else if (nt_1 != nt_2) {
      tmp_mismatches++;
    }
  }

  n_ambiguous = tmp_ambiguous;
  n_mismatches = tmp_mismatches;

  return (2 * n_mismatches) + n_ambiguous <= max_penalty;
}

} // namespace adapterremoval::simd
