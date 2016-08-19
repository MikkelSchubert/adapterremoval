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
#include "debug.h"

#include <cstdlib>
#include <iostream>
#include <sstream>

namespace ar
{

#ifdef AR_TEST_BUILD
assert_failed::assert_failed(const assert_failed& errror)
    : m_what(errror.m_what)
{
}


assert_failed::assert_failed(const std::string& what)
    : m_what(what)
{
}


assert_failed::~assert_failed() noexcept
{
}

const char* assert_failed::what() const noexcept
{
    return m_what.c_str();
}
#endif


void debug_raise_assert(const char* filename, size_t lineno, const char* what)
{
    std::stringstream message;
    message << "\nFATAL ERROR:\n"
            << "Debug assertion failed in '" << filename << "', line "
            << lineno << ": " << what << "\n\n"
            << "This should not happen! Please file a bug-report at\n    "
            << "https://github.com/MikkelSchubert/adapterremoval/issues/new";

#ifdef AR_TEST_BUILD
    throw assert_failed(message.str());
#else
    std::cerr << message.str() << std::endl;

    std::abort();
#endif
}

} // namespace ar
