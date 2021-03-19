/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#ifndef ALIGNMENT_TABLES_H
#define ALIGNMENT_TABLES_H

#include <sys/types.h>

namespace ar {

const size_t PHRED_TABLE_SIZE = 8836;

/**
 * Table of Phred scores to assigned for identical positions during collapsing.
 *
 * Position is calculated as phred_1 * (MAX_PHRED_SCORE + 1) + phred_2,
 * assuming that phred_1 >= phred_2.
 */
extern const signed char IDENTICAL_NTS[PHRED_TABLE_SIZE];

/**
 * Table of Phred scores to assigned for mismatching positions during
 * collapsing.
 *
 * Position is calculated as phred_1 * (MAX_PHRED_SCORE + 1) + phred_2,
 * assuming that phred_1 >= phred_2.
 */
extern const signed char DIFFERENT_NTS[PHRED_TABLE_SIZE];

} // namespace ar

#endif
