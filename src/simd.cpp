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
#include "simd.hpp"
#include "debug.hpp"
#include <string> // for basic_string

namespace adapterremoval {

namespace simd {

std::vector<instruction_set>
supported()
{
  std::vector<instruction_set> choices = { instruction_set::none };

  if (__builtin_cpu_supports("sse2")) {
    choices.push_back(instruction_set::sse2);
  }

  if (__builtin_cpu_supports("avx2")) {
    choices.push_back(instruction_set::avx2);
  }

  if (__builtin_cpu_supports("avx512bw")) {
    choices.push_back(instruction_set::avx512);
  }

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

compare_subsequences_func
get_compare_subsequences_func(instruction_set is)
{
  switch (is) {
    case instruction_set::none:
      return &compare_subsequences_std;
    case instruction_set::sse2:
      return &compare_subsequences_sse2;
    case instruction_set::avx2:
      return &compare_subsequences_avx2;
    case instruction_set::avx512:
      return &compare_subsequences_avx512;
    default:
      AR_FAIL("SIMD function not implemented!");
  }
}

} // namespace simd

} // namespace adapterremoval
