/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
  explicit thread_error(const std::string& message);

  /** Returns error message; lifetime is the same as the object. */
  const char* what() const noexcept override;

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
};

} // namespace adapterremoval
