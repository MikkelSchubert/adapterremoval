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
#ifndef AR_DEBUG_H
#define AR_DEBUG_H

#include <stdexcept>
#include <string>


/**
 * Aborts after printing the filename, line-number, and message, plus
 * instructions for how to report the problem.
 */
void debug_raise_assert(const char* filename, size_t lineno,
                        const char* what) __attribute__ ((noreturn));


/** Custom assert which prints various information on failure; always enabled. */
#define AR_DEBUG_ASSERT(test) \
    do { \
        if (!(test)) { \
            debug_raise_assert(__FILE__, __LINE__, #test); \
        } \
    } while (0)


/** Raise an assert failure with a user-specified message. */
#define AR_DEBUG_FAIL(msg) \
    debug_raise_assert(__FILE__, __LINE__, msg)


#endif
