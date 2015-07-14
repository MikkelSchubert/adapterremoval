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


const size_t REPORT_EVERY = 1e6;


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


timer::timer(const std::string& what, bool muted)
  : m_what(what)
  , m_counter(0)
  , m_first_time(get_current_time())
  , m_muted(muted)
{
}


void timer::increment(size_t inc)
{
    m_counter += inc;
    if (!m_muted && (m_counter % REPORT_EVERY == 0)) {
        const double last_time = get_current_time();
        const double rate = (last_time - m_first_time) / (m_counter / REPORT_EVERY);

        std::stringstream stream;
        stream << "\rProcessed " << thousands_sep(m_counter) << " " << m_what
               << " in " << format_time(last_time - m_first_time) << "; "
               << format_time(rate) << " per " << thousands_sep(REPORT_EVERY)
               << " " << m_what << " ...";

        std::cerr << stream.str();
        std::cerr.flush();
    }
}


void timer::finalize() const
{
    if (!m_muted) {
        const double current_time = get_current_time();
        if (m_counter >= REPORT_EVERY) {
            std::cerr << "\n";
        }

        std::cerr << "Processed a total of " << thousands_sep(m_counter) << " "
                  << m_what << " in " << format_time(current_time - m_first_time)
                  << " ..." << std::endl;
    }
}
