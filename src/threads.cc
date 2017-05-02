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

#include "threads.h"

#include <ctime>

namespace ar
{

///////////////////////////////////////////////////////////////////////////////
// exceptions

thread_error::thread_error(const std::string& message)
    : std::exception()
    , m_message(message)
{
}


thread_error::thread_error(const thread_error& error)
    : std::exception()
    , m_message(error.m_message)
{
}


thread_error::~thread_error() noexcept
{
}


const char* thread_error::what() const noexcept
{
    return m_message.c_str();
}


thread_abort::thread_abort()
  : thread_error("abort thread")
{
}


///////////////////////////////////////////////////////////////////////////////
// print_locker

//! Shared mutex for STDOUT / STDERR
static std::mutex s_print_mutex;

//! Shared bool indicating if STDERR contains a partial line.
static bool s_stderr_is_incomplete = false;


print_locker::print_locker(bool flush_stderr)
  : m_lock(s_print_mutex)
{
    if (flush_stderr && s_stderr_is_incomplete) {
        s_stderr_is_incomplete = false;
        std::cerr << std::endl;
    }
}


print_locker::~print_locker()
{
}


void print_locker::partial_stderr_output() {
    s_stderr_is_incomplete = true;
}


} // namespace ar
