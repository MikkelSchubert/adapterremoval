/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2021 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <chrono> // for duration, high_resolution_clock, system_clock::time..
#include <mutex>  // for mutex
#include <random> // for mt19937

#include "utilities.hpp"

namespace adapterremoval {

uint32_t
prng_seed()
{
  static std::mutex lock;
  std::unique_lock<std::mutex> guard(lock);

  using clock = std::chrono::high_resolution_clock;
  static std::mt19937 rng(clock::now().time_since_epoch().count());

  return rng();
}

} // namespace adapterremoval
