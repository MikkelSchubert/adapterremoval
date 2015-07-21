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

#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <unistd.h>

// #include <mach/mach_time.h>

#include "threads.h"

#include <ctime>

///////////////////////////////////////////////////////////////////////////////
// exceptions

thread_error::thread_error(const std::string& message)
    : std::exception()
    , m_message(message)
{
}


thread_error::~thread_error() throw()
{
}


const char* thread_error::what() const throw()
{
    return m_message.c_str();
}


///////////////////////////////////////////////////////////////////////////////
// mutex

mutex::mutex()
    : m_mutex()
{
    switch (pthread_mutex_init(&m_mutex, NULL)) {
        case 0:
            break;

        case EAGAIN:
            throw thread_error("mutex::mutex: out of resouces");
            break;

        case ENOMEM:
            throw thread_error("mutex::mutex: out of memory");
            break;

        case EPERM:
            throw thread_error("mutex::mutex: insufficient permissions");
            break;

        case EBUSY:
            throw thread_error("mutex::mutex: mutex already initialied");
            break;

        case EINVAL:
            throw thread_error("mutex::mutex: invalid attributes");
            break;

        default:
            throw thread_error("Unknown error in mutex::mutex");
    }
}


mutex::~mutex()
{
    const int error = pthread_mutex_destroy(&m_mutex);
    if (error) {
        print_locker lock;

        switch (error) {
            case EBUSY:
                std::cerr << "mutex::~mutex: mutex busy" << std::endl;
                break;

            case EINVAL:
                std::cerr << "mutex::~mutex: invalid mutex" << std::endl;
                break;

            default:
                std::cerr << "Unknown error in mutex::~mutex: " << error << std::endl;
                break;
        }

        std::exit(1);
    }
}


///////////////////////////////////////////////////////////////////////////////
// conditional

conditional::conditional()
  : m_mutex()
  , m_cond()
  , m_count(0)
{
    switch (pthread_cond_init(&m_cond, NULL)) {
        case 0:
            break;

        case EAGAIN:
            throw thread_error("conditional::conditional: insufficient resouces to init conditional");

        case ENOMEM:
            throw thread_error("conditional::conditional: insufficient memory to init conditional");

        case EBUSY:
            throw thread_error("conditional::conditional: conditional already initialized");

        case EINVAL:
            throw thread_error("conditional::conditional: invalid attributes");

        default:
            throw thread_error("Unknown error in conditional::conditional");
    }
}


conditional::~conditional()
{
    const int error = pthread_cond_destroy(&m_cond);
    if (error) {
        print_locker lock;

        switch (error) {
            case EBUSY:
                std::cerr << "conditional::~conditional: conditional busy" << std::endl;
                break;

            case EINVAL:
                std::cerr << "conditional::~conditional: invalid conditional" << std::endl;
                break;

            default:
                std::cerr << "Unknown error in conditional::~conditional: " << error << std::endl;
                break;
        }

        std::exit(1);
    }
}


void conditional::wait()
{
    mutex_locker lock(m_mutex);

    while (!m_count) {
        switch (pthread_cond_wait(&m_cond, &m_mutex.m_mutex)) {
            case 0:
                break;

            case EINVAL:
                throw thread_error("conditional::wait: invalid conditional or mutex");

            case EPERM:
                throw thread_error("conditional::wait: mutex not owned by thread");

            default:
                throw thread_error("Unknown error in conditional::wait");
        }
    }

    --m_count;
}


void conditional::signal()
{
    mutex_locker lock(m_mutex);

    ++m_count;

    switch (pthread_cond_signal(&m_cond)) {
        case 0:
            break;

        case EINVAL:
            throw thread_error("conditional::signal: invalid conditional");

        default:
            throw thread_error("Unknown error in conditional::signal");
    }
}


///////////////////////////////////////////////////////////////////////////////
// mutex_locker

mutex_locker::mutex_locker(mutex& to_lock)
  : m_mutex(to_lock.m_mutex)
{
    switch (pthread_mutex_lock(&m_mutex)) {
        case 0:
            break;

        case EINVAL:
            throw thread_error("mutex::lock: mutex not initialized");

        case EAGAIN:
            throw thread_error("mutex::lock: max number of recursive locks exceeded");

        case EDEADLK:
            throw thread_error("mutex::try_lock: deadlock detected");

        default:
            throw thread_error("Unknown error in mutex::lock");
    }
}


mutex_locker::~mutex_locker()
{
    const int error = pthread_mutex_unlock(&m_mutex);
    if (error) {
        print_locker lock;

        switch (error) {
            case EINVAL:
                std::cerr << "mutex_lock::~mutex_lock: mutex not initialized" << std::endl;
                break;

            case EAGAIN:
                std::cerr << "mutex_lock::~mutex_lock: max number of recursive locks exceeded" << std::endl;
                break;

            case EPERM:
                std::cerr << "mutex_lock::~mutex_lock: the current thread does not own the mutex" << std::endl;
                break;

            default:
                std::cerr << "Unknown error in mutex_lock::~mutex_lock: " << error << std::endl;
                break;
        }

        std::exit(1);
    }
}


///////////////////////////////////////////////////////////////////////////////
// print_locker

//! Shared mutex for STDOUT / STDERR
static mutex s_print_mutex;


print_locker::print_locker()
  : mutex_locker(s_print_mutex)
{
}


print_locker::~print_locker()
{
}
