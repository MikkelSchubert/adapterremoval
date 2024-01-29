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
#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

namespace adapterremoval {

using std::size_t;

namespace simd {

enum class instruction_set
{
  none,
  sse2,
  avx2,
  avx512,
};

/** Returns vector of supports instruction sets */
std::vector<instruction_set>
supported();

/** Returns human-readable name of instruction set */
const char*
name(instruction_set value);

bool
compare_subsequences_std(size_t& n_mismatches,
                         size_t& n_ambiguous,
                         const char* seq_1,
                         const char* seq_2,
                         size_t length,
                         size_t max_penalty);

bool
compare_subsequences_sse2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty);

bool
compare_subsequences_avx2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty);

bool
compare_subsequences_avx512(size_t& n_mismatches,
                            size_t& n_ambiguous,
                            const char* seq_1,
                            const char* seq_2,
                            size_t length,
                            size_t max_penalty);

using compare_subsequences_func = decltype(&compare_subsequences_std);

compare_subsequences_func
get_compare_subsequences_func(instruction_set is);

} // namespace simd

} // namespace adapterremoval
