// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2021 Mikkel Schubert <mikkelsch@gmail.com>
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
