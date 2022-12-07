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

namespace adapterremoval {

bool
supports_avx2()
{
  return __builtin_cpu_supports("avx2");
}

bool
supports_sse2()
{
  return __builtin_cpu_supports("sse2");
}

compare_subsequences_func
select_compare_subsequences_func()
{
  if (supports_avx2()) {
    return &compare_subsequences_avx2;
  } else if (supports_sse2()) {
    return &compare_subsequences_sse2;
  } else {
    return &compare_subsequences_std;
  }
}

} // namespace adapterremoval
