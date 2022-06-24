/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <iomanip> // for operator<<, setfill, setw
#include <locale>  // for numpunct, use_facet, locale
#include <sstream> // for ostringstream

#include "logging.hpp"  // for log
#include "progress.hpp" // declarations
#include "strutils.hpp" // for
#include "threads.hpp"  // for print_locker

namespace adapterremoval {

//! Print progress report every N items
const size_t REPORT_EVERY = 1e6;

std::string
format_time(double seconds)
{
  std::ostringstream stream;
  stream.precision(1);
  stream << std::setfill('0');

  if (seconds > 60 * 60) {
    stream << static_cast<size_t>(seconds) / (60 * 60) << ":";
    stream << std::setw(2);
  }

  if (seconds > 60) {
    stream << (static_cast<size_t>(seconds) % (60 * 60)) / 60 << ":";
    stream << std::setw(4) << std::setfill('0');
  }

  stream << std::fixed << (static_cast<size_t>(seconds * 100) % 6000) / 100.0
         << "s";

  return stream.str();
}

progress_timer::progress_timer(const std::string& what)
  : m_what(what)
  , m_total(0)
  , m_current(0)
  , m_timer()
  , m_last_time()
{
}

void
progress_timer::increment(size_t inc)
{
  m_total += inc;
  m_current += inc;

  if (m_current >= REPORT_EVERY) {
    const double current_time = m_timer.duration();

    do_print(m_current, current_time - m_last_time);

    m_current = 0;
    m_last_time = current_time;
  }
}

void
progress_timer::finalize() const
{
  do_print(m_total, m_timer.duration(), true);
}

void
progress_timer::do_print(size_t items, double seconds, bool finalize) const
{
  size_t rate = static_cast<size_t>(items / seconds);

  if (finalize) {
    log::info() << "Processed " << format_thousand_sep(m_total) << " " << m_what
                << " in " << format_time(m_timer.duration()) << "; "
                << format_rough_number(rate) << " " << m_what
                << "/s on average";
  } else {
    log::info() << "Processed " << static_cast<size_t>(m_total / REPORT_EVERY)
                << " M " << m_what << " in " << format_time(m_timer.duration())
                << "; " << format_rough_number(rate) << " " << m_what << "/s";
  }
}

} // namespace adapterremoval
