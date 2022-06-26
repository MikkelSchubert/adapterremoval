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
#pragma once

#include <deque>
#include <mutex>
#include <string>
#include <thread>

#include "timer.hpp"

namespace adapterremoval {

enum class progress_type
{
  //! No progress reports
  none,
  //! Messages written to the log
  simple,
  //! An animated spinner plus messages
  spinner,
};

/**
 * Timer for reporting current progress of a run.
 *
 * The timer takes two forms:
 *  1. A log based timer that prints the current progress every 1M reads.
 *  2. An animated timer (spinner) that updates the current number of processed
 *     reads once every second. This rate is intentionaly lower than the speed
 *     of the since it feels noisy to have the statistics constantly updating.
 *
 * A final summary is printed using the 'finalize' function for both modes
 */
class progress_timer
{
public:
  /** Initializes a progress timer and starts the animation if enabled */
  progress_timer(progress_type type);

  /** Destructor. Silently ends the animation if not finalized. */
  ~progress_timer();

  /** Increment the progress, and (possibly) print a status report. */
  void increment(size_t inc = 1);

  /** Print final summary of the reads processed and ends the animation. */
  void finalize();

private:
  /** Starts the loop (if enabled) and hides the cursor */
  void start();
  /** Loop function for spinning thread */
  void loop();
  /** Stop the looping animation */
  void stop();

  //! The kind of progress messages to be printed
  const progress_type m_type;
  //! Total number of items processed
  size_t m_total;
  //! Number of items processed since last update
  size_t m_current;
  //! Starting time (in seconds) of the timer.
  monotonic_timer m_timer;
  //! Time (in seconds) of the last update
  double m_last_time;

  //! Thread used for animated spinner
  std::thread m_spinner;
  //! Mutex protecting the following member variables when using a spinner
  std::mutex m_lock;
  //! Indicates if the spinner is active
  bool m_spinning;
};

} // namespace adapterremoval
