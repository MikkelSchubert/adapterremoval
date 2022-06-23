/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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

#include "logging.hpp" // for log
#include "threads.hpp" // for print_locker
#include "timer.hpp"   // declarations

namespace adapterremoval {

//! Print progress report every N items
const size_t REPORT_EVERY = 1e6;

std::string
thousands_sep(size_t number)
{
  if (!number) {
    return std::string(1, '0');
  }

  static const std::locale locale;
  static const std::numpunct<char>& facet =
    std::use_facet<std::numpunct<char>>(locale);

  std::string str;
  for (size_t i = 1; number; ++i) {
    str.append(1, '0' + number % 10);
    number /= 10;

    if (number && (i % 3 == 0)) {
      str.append(1, facet.thousands_sep());
    }
  }

  return std::string(str.rbegin(), str.rend());
}

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
  if (rate >= 10000) {
    rate = (rate / 1000) * 1000;
  }

  if (finalize) {
    log::info() << "Processed " << thousands_sep(m_total) << " " << m_what
                << " in " << format_time(m_timer.duration()) << "; "
                << thousands_sep(rate) << " " << m_what << "/s on average";
  } else {
    log::info() << "Processed " << thousands_sep(m_total) << " " << m_what
                << " in " << format_time(m_timer.duration()) << "; "
                << thousands_sep(rate) << " " << m_what << "/s";
  }
}

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
