// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "threading.hpp" // declarations
#include "debug.hpp"     // for AR_REQUIRE
#include "logging.hpp"   // for log
#include <string>        // for string

namespace adapterremoval {

void
finalize_reusable_objects(const std::string& description,
                          const size_t n_allocated,
                          const size_t n_released,
                          const size_t n_acquired)
{
  const auto level =
    (n_allocated == n_released) ? log::level::debug : log::level::error;

  log::log_stream(level) << "Cache statistics for " << description
                         << ": objects created = " << n_allocated
                         << ", objects released = " << n_released
                         << ", allocations saved = "
                         << n_acquired - n_allocated;

  // too many objects being released is caught in `release()`
  AR_REQUIRE(n_allocated == n_released,
             "not all objects returned to " + description + " cache");
}

} // namespace adapterremoval
