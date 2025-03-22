// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <chrono> // for steady_clock, time_point

namespace adapterremoval {

/**
 * Steady timer for measuring durations.
 */
class monotonic_timer
{
  using steady_clock = std::chrono::steady_clock;
  using time_point = std::chrono::time_point<steady_clock>;

public:
  /** Constructor. */
  monotonic_timer();

  /** Returns the duration in seconds since the timer was created. */
  double duration() const;

private:
  //! Starting time of the timer.
  time_point m_start_time;
};

} // namespace adapterremoval
