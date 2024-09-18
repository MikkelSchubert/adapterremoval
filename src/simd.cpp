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
#include "simd.hpp"  // declarations
#include "debug.hpp" // for AR_FAIL

#ifdef MESON
#include "config-ar3.hpp" // for HAVE_SSE2, HAVE_AVX2, HAVE_AVX512...
#else
#define HAVE_SSE2 1
#define HAVE_AVX2 1

#if __GNUC__ >= 11 || (defined(__clang_major__) && __clang_major__ >= 8)
#define HAVE_AVX512 1
#else
#pragma GCC warning "AVX512 support requires GCC >= 11.0 or Clang >= 8.0"
#define HAVE_AVX512 0
#endif
#endif

namespace adapterremoval {

namespace simd {

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

std::vector<instruction_set>
supported()
{
  std::vector<instruction_set> choices = { instruction_set::none };

#if HAVE_SSE2
  if (__builtin_cpu_supports("sse2")) {
    choices.push_back(instruction_set::sse2);
  }
#endif

#if HAVE_AVX2
  if (__builtin_cpu_supports("avx2")) {
    choices.push_back(instruction_set::avx2);
  }
#endif

#if HAVE_AVX512
  if (__builtin_cpu_supports("avx512bw")) {
    choices.push_back(instruction_set::avx512);
  }
#endif

  return choices;
}

const char*
name(instruction_set value)
{
  switch (value) {
    case instruction_set::none:
      return "none";
    case instruction_set::sse2:
      return "SSE2";
    case instruction_set::avx2:
      return "AVX2";
    case instruction_set::avx512:
      return "AVX512";
    default:
      AR_FAIL("SIMD function not implemented!");
  }
}

size_t
padding(instruction_set value)
{
  switch (value) {
    case instruction_set::none:
      return 0;
    case instruction_set::sse2:
#if HAVE_SSE2
      return 15;
#else
      AR_FAIL("SSE2 not supported!");
#endif
    case instruction_set::avx2:
#if HAVE_AVX2
      return 31;
#else
      AR_FAIL("AVX2 not supported!");
#endif
    case instruction_set::avx512:
#if HAVE_AVX512
      return 0;
#else
      AR_FAIL("AVX512 not supported!");
#endif
    default:
      AR_FAIL("SIMD function not implemented!");
  }
}

compare_subsequences_func
get_compare_subsequences_func(instruction_set is)
{
  switch (is) {
    case instruction_set::none:
      return &compare_subsequences_std;
    case instruction_set::sse2:
#if HAVE_SSE2
      return &compare_subsequences_sse2;
#else
      AR_FAIL("SSE2 not supported!");
#endif
    case instruction_set::avx2:
#if HAVE_AVX2
      return &compare_subsequences_avx2;
#else
      AR_FAIL("AVX2 not supported!");
#endif
    case instruction_set::avx512:
#if HAVE_AVX512
      return &compare_subsequences_avx512;
#else
      AR_FAIL("AVX512 not supported!");
#endif
    default:
      AR_FAIL("SIMD function not implemented!");
  }
}

} // namespace simd

} // namespace adapterremoval
