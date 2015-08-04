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

#ifdef AR_PTHREAD_SUPPORT
#include <pthread.h>
#endif


/**
 * Exception thrown for threading related errors, including errors with
 * threads, mutexes, and conditionals.
 */
class thread_error : public std::exception
{
public:
    /** Constructor; takes an error-message. */
    thread_error(const std::string& message);
    /** Destructor; does nothing. */
    ~thread_error() throw();

    /** Returns error message; lifetime is the same as the object. */
    virtual const char* what() const throw();

private:
    //! User provided error message
    std::string m_message;
};


/**
 * Wrapper around PThread Mutex.
 *
 * Locking is performed using the mutex_locker object.
 */
class mutex
{
public:
    /** Constructor. Creates unlocked mutex. */
    mutex();
    /** Destructor. Mutex must be unlocked when destructing; exits on error. */
    ~mutex();

private:
    //! Not implemented
    mutex(const mutex&);
    //! Not implemented
    mutex& operator=(const mutex&);

    friend class mutex_locker;
    friend class conditional;

#ifdef AR_PTHREAD_SUPPORT
    pthread_mutex_t m_mutex;
#endif
};


/**
 * Simple conditional.
 *
 *
 */
class conditional
{
public:
    /** Constructor; initializes conditional. */
    conditional();
    /** Destructor; destroys conditional, exits on error. */
    ~conditional();

    /** Wait for and consume signal; non-blocking if signals are queued. */
    void wait();

    /** Signal a waiting thread, or queue signal if no threads are waiting. */
    void signal();


#ifdef AR_PTHREAD_SUPPORT
private:
    //! Mutex assosiated with conditional; not exposed.
    mutex m_mutex;
    //! Raw conditional; not exposed.
    pthread_cond_t m_cond;
    //! Number of queued signals.
    volatile unsigned m_count;
#endif
};


/** Simple class for automatic locking / unlocking of mutexes. **/
class mutex_locker
{
public:
    //! Locks the mutex (blocking)
    mutex_locker(mutex& to_lock);

    //! Unlocks the mutex
    virtual ~mutex_locker();

private:
    //! Not implemented
    mutex_locker(const mutex_locker&);
    //! Not implemented
    mutex_locker& operator=(const mutex_locker&);

#ifdef AR_PTHREAD_SUPPORT
    pthread_mutex_t& m_mutex;
#endif
};


/**
 * Locker for using stdout / stderr.
 *
 * Any useage of stdout and / or stderr should be preceeded by creating a
 * print_locker object. This ensures that output from different threads is
 * not interleaved, regardless of the destination of these pipes.
 */
class print_locker : public mutex_locker
{
public:
    //! Locks the mutex (blocking)
    print_locker();

    //! Unlocks the mutex
    ~print_locker();

private:
    //! Not implemented
    print_locker(const mutex_locker&);
    //! Not implemented
    print_locker& operator=(const mutex_locker&);
};


/**
 * Basic atomic counter.
 */
class atomic_counter
{
public:
    /** Initialize counter to the given value. */
    atomic_counter(size_t init = 0);

    /** Returns the current value (which may have changed already). */
    size_t current() const;

    /** Increment the current value. */
    size_t increment();

    /** Decrement the current value. */
    size_t decrement();

private:
    //! Mutex used to control access to the counter/
    mutable mutex m_lock;
    //! Raw counter value; access controlled using m_lock;
    size_t m_count;
};


#endif
