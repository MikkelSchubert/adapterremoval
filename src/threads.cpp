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
#include "threads.hpp"

namespace adapterremoval {

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

thread_error::~thread_error() {}

const char*
thread_error::what() const noexcept
{
  return m_message.c_str();
}

thread_abort::thread_abort()
  : thread_error("abort thread")
{
}

thread_abort::~thread_abort() {}

} // namespace adapterremoval
