// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "timer.hpp" // declarations

namespace adapterremoval {

monotonic_timer::monotonic_timer()
  : m_start_time(steady_clock::now())
{
}

double
monotonic_timer::duration() const
{
  const auto end = steady_clock::now();
  const std::chrono::duration<double> diff = end - m_start_time;

  return diff.count();
}

} // namespace adapterremoval
