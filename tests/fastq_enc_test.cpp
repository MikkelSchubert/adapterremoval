/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
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
#include <limits>
#include <stdexcept>
#include <gtest/gtest.h>

#include "testing.hpp"
#include "debug.hpp"
#include "fastq.hpp"


namespace ar
{

///////////////////////////////////////////////////////////////////////////////
// Names (default objects)

TEST(fastq_enc, global_objects__name)
{
    ASSERT_STREQ("Phred+33", FASTQ_ENCODING_33.name());
    ASSERT_STREQ("Phred+64", FASTQ_ENCODING_64.name());
    ASSERT_STREQ("Phred+33", FASTQ_ENCODING_SAM.name());
    ASSERT_STREQ("Solexa", FASTQ_ENCODING_SOLEXA.name());
}

}
