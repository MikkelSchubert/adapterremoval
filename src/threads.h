/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#ifndef THREADS_H
#define THREADS_H

#include <string>
#include <thread>
#include <mutex>


namespace ar
{

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
    ~thread_error() noexcept;

    /** Returns error message; lifetime is the same as the object. */
    virtual const char* what() const noexcept;

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


/**
 * Locker for using stdout / stderr.
 *
 * Any useage of stdout and / or stderr should be preceded by creating a
 * print_locker object. This ensures that output from different threads is
 * not interleaved, regardless of the destination of these pipes.
 */
class print_locker
{
public:
    /*
     * Locks the mutex (blocking). If flush_stderr is true, and
     * partial_stderr_output has been called, then a newline is first
     * written to stderr.
     */
    print_locker(bool flush_stderr=true);

    //! Unlocks the mutex
    ~print_locker();

    //! Call to indicate that a partial line has been written to STDERR.
    void partial_stderr_output();

private:
    //! Not implemented
    print_locker(const print_locker&);
    //! Not implemented
    print_locker& operator=(const print_locker&);

    std::lock_guard<std::mutex> m_lock;
};


} // namespace ar

#endif
