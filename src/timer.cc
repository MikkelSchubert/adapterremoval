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
#include <locale>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sys/time.h>

#include "timer.h"
#include "threads.h"

namespace ar
{

//! Print progress report every N items
const size_t REPORT_EVERY = 1e6;
//! Number of blocks to store for calculating mean rate
const size_t AVG_BLOCKS = 10;


double get_current_time()
{
    struct timeval timestamp;
    gettimeofday(&timestamp, NULL);

    return timestamp.tv_sec + timestamp.tv_usec / 1e6;
}


std::string thousands_sep(size_t number)
{
    if (!number) {
        return std::string(1, '0');
    }

    static const std::locale locale;
    static const std::numpunct<char>& facet = std::use_facet<std::numpunct<char> >(locale);

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


std::string format_time(double seconds)
{
    std::stringstream stream;
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

    stream << std::fixed
           << (static_cast<size_t>(seconds * 100) % 6000) / 100.0 << "s";

    return stream.str();
}


timer::timer(const std::string& what)
  : m_what(what)
  , m_total(0)
  , m_first_time(get_current_time())
  , m_counts()
{
    m_counts.push_back(time_count_pair(get_current_time(), 0));
}


void timer::increment(size_t inc)
{
    m_total += inc;
    m_counts.back().second += inc;

    if (m_counts.back().second >= REPORT_EVERY) {
        const double current_time = get_current_time();
        // Number of seconds since oldest block was created
        const double seconds = current_time - m_counts.front().first;

        size_t current_total = 0;
        for (time_count_deque::iterator it = m_counts.begin(); it != m_counts.end(); ++it) {
            current_total += it->second;
        }

        do_print(static_cast<size_t>(current_total / seconds), current_time);

        m_counts.push_back(time_count_pair(current_time, 0));
        while (m_counts.size() > AVG_BLOCKS) {
            m_counts.pop_front();
        }
    }
}


void timer::finalize() const
{
    const double current_time = get_current_time();
    const double seconds = current_time - m_first_time;

    do_print(static_cast<size_t>(m_total / seconds), current_time, true);
}


void timer::do_print(size_t rate, double current_time, bool finalize) const
{
    print_locker lock;

    if (finalize) {
        std::cerr << "\rProcessed a total of ";
    } else {
        std::cerr << "\rProcessed ";
    }

    if (rate > 10000) {
        rate = (rate / 1000) * 1000;
    }

    std::cerr << thousands_sep(m_total) << " " << m_what << " in "
              << format_time(current_time - m_first_time) << "; "
              << thousands_sep(rate) << " " << m_what << " per second ";

    if (finalize) {
        std::cerr << "on average ..." << std::endl;
    } else {
        std::cerr << "...";
        std::cerr.flush();
    }
}

} // namespace ar
