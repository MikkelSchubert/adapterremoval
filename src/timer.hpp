/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
