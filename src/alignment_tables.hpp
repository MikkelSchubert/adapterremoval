/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
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

#include <array>    // for array
#include <stddef.h> // for size_t

namespace adapterremoval {

const size_t PHRED_TABLE_SIZE = 8836;

/**
 * Table of Phred scores to assigned for identical positions during merging.
 *
 * Position is calculated as phred_1 * (MAX_PHRED_SCORE + 1) + phred_2,
 * assuming that phred_1 >= phred_2.
 */
extern const std::array<signed char, PHRED_TABLE_SIZE> IDENTICAL_NTS;

/**
 * Table of Phred scores to assigned for mismatching positions during merging.
 *
 * Position is calculated as phred_1 * (MAX_PHRED_SCORE + 1) + phred_2,
 * assuming that phred_1 >= phred_2.
 */
extern const std::array<signed char, PHRED_TABLE_SIZE> DIFFERENT_NTS;

} // namespace adapterremoval
