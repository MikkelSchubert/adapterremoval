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
#ifndef TIMER_H
#define TIMER_H

#include <string>


/**
 * Simply class for reporting current progress of a run.
 *
 * Every 1 million records / reads / etc processed, when the counter is
 * incremented using the 'increment' function, a progress report is printed:
 *   "Processed last 1,000,000 pairs in 14.1s; 2,000,000 pairs in 28.9s ..."
 *
 * A final summary is printed using the 'finalize' function:
 *   "Processed a total of 4,000,000 reads in 31.9s ..."
 */
class timer
{
public:
    /* Constructor.
     *
     * @param what Short name of what is being processed, for use in reports.
     * @param muted If true, no messages are printed.
     *
     * The 'muted' parameter is included to simply handling of optional
     * progress reports, without having to add conditional logic around calls
     * to the functions 'increment' and 'finalize'.
     */
    timer(const std::string& what, bool muted = false);

    /** Increment the progress, and (possibly) print a status report. */
    void increment();

    /** Print final summary based on the number of increments. */
    void finalize() const;

private:
    //! Description of what is being processed.
    std::string m_what;
    //! Number of items processed (total)
    size_t m_counter;
    //! Starting time (in seconds) of the timer.
    double m_first_time;
    //! Last time (in seconds) that a report was given.
    double m_last_time;
    //! If true, no output is produced by calling 'increment' or 'finalize'.
    bool m_muted;
};


#endif
