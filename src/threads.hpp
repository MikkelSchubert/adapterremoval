/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include <exception> // for exception
#include <mutex>     // for lock_guard, mutex
#include <string>    // for string

namespace adapterremoval {

/**
 * Exception thrown for threading related errors, including errors with
 * threads, mutexes, and conditionals.
 */
class thread_error : public std::exception
{
public:
  /** Constructor; takes an error-message. */
  thread_error(const std::string& message);
  /** Copy-constructor; takes an exiting error. */
  thread_error(const thread_error& error);
  /** Destructor; does nothing. */
  virtual ~thread_error() override;

  /** Returns error message; lifetime is the same as the object. */
  virtual const char* what() const noexcept override;

private:
  //! User provided error message
  std::string m_message;
};

/**
 * This exception may be thrown by a task to abort the thread; error-messages
 * are assumed to have already been printed by the thrower, and no further
 * messages are printed.
 */
class thread_abort : public thread_error
{
public:
  thread_abort();
  thread_abort(const thread_abort& other) = default;

  virtual ~thread_abort() override;
};

} // namespace adapterremoval
