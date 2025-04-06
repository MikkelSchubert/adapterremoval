// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2021 Mikkel Schubert <mikkelsch@gmail.com>
#include <mutex>  // for mutex
#include <random> // for random_device

#include "utilities.hpp"

namespace adapterremoval {

unsigned int
prng_seed()
{
  static std::mutex lock;
  static std::random_device source;

  std::unique_lock<std::mutex> guard(lock);
  return source();
}

} // namespace adapterremoval
