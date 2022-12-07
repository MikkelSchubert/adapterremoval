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
#include <cstddef> // for size_t

namespace adapterremoval {

class fastq;

using std::size_t;

/** Returns true if the current CPU supports SSE2 instructions */
bool
supports_sse2();

/** Returns true if the current CPU supports AVX2 instructions */
bool
supports_avx2();

bool
compare_subsequences_std(size_t& n_mismatches,
                         size_t& n_ambiguous,
                         const char* seq_1_ptr,
                         const char* seq_2_ptr,
                         size_t max_mismatches,
                         int remaining_bases);

bool
compare_subsequences_sse2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1_ptr,
                          const char* seq_2_ptr,
                          size_t max_mismatches,
                          int remaining_bases);

bool
compare_subsequences_avx2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1_ptr,
                          const char* seq_2_ptr,
                          size_t max_mismatches,
                          int remaining_bases);

using compare_subsequences_func = decltype(&compare_subsequences_std);

compare_subsequences_func
select_compare_subsequences_func();

} // namespace adapterremoval
