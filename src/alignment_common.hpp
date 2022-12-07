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

/**
 * Compares two subsequences in an alignment to a previous (best) alignment.
 *
 * @param n_mismatches The current number of mismathces.
 * @param n_ambiguous The current number of ambiguous bases.
 * @param seq_1_ptr Pointer to the first sequence in the alignment.
 * @param seq_2_ptr Pointer to the second sequence in the alignment.
 * @return True if the current alignment is at least as good as the best
 * alignment, false otherwise.
 *
 * If the function returns false, the current alignment cannot be assumed to
 * have been completely evaluated (due to early termination), and hence counts
 * and scores are not reliable. The function assumes uppercase nucleotides.
 */
bool
compare_subsequences(size_t& n_mismatches,
                     size_t& n_ambiguous,
                     const char* seq_1_ptr,
                     const char* seq_2_ptr,
                     double mismatch_threshold,
                     int remaining_bases);

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

} // namespace adapterremoval
