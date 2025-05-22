// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

namespace adapterremoval {

namespace simd {

enum class instruction_set
{
  none,
  sse2,
  avx2,
  avx512,
  neon,
};

/** Returns vector of supports instruction sets */
std::vector<instruction_set>
supported();

/** Returns human-readable name of instruction set */
const char*
name(instruction_set value);

/** Returns the amount of padding expected when using this instruction set */
size_t
padding(instruction_set value);

/***/
using compare_subsequences_func = bool (*)(size_t& n_mismatches,
                                           size_t& n_ambiguous,
                                           const char* seq_1,
                                           const char* seq_2,
                                           size_t length,
                                           size_t max_penalty);

/**
 * Returns a pair-wise alignment function for a given instruction set.
 *
 * Does not check if the instruction set is supported (see `supported`). The
 * alignment functions expect that sequences are 'N' padded with the amount
 * specified by calling the `padding` function with the same `instruction_set`.
 */
compare_subsequences_func
get_compare_subsequences_func(instruction_set is);

} // namespace simd

} // namespace adapterremoval
