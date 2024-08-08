/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "buffer.hpp"      // for buffer
#include "commontypes.hpp" // for fastq_vec
#include "fastq.hpp"       // for fastq, fastq::ntrimmed, ACGTN, ACGT
#include "serializer.hpp"  // for fastq_serializer
#include "testing.hpp"     // for catch.hpp, StringMaker

// Ignore nucleotide and quality strings
// spell-checker:ignoreRegExp /"[!-~]+"/g
// Ignore nucleotide comments
// spell-checker:ignoreRegExp /\W[acgtnACGTN]+\W/g

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Writing to stream

TEST_CASE("Writing_to_stream_phred_33", "[fastq::fastq]")
{
  const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
  buffer buf;
  fastq_serializer::record(buf, record, fastq_flags::se, '\0', read_group());

  REQUIRE(
    std::string_view(reinterpret_cast<const char*>(buf.data()), buf.size()) ==
    "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n");
}

TEST_CASE("Writing_to_stream_phred_33_explicit", "[fastq::fastq]")
{
  const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
  buffer buf;
  fastq_serializer::record(buf, record, fastq_flags::se, '\0', read_group());

  REQUIRE(
    std::string_view(reinterpret_cast<const char*>(buf.data()), buf.size()) ==
    "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n");
}

} // namespace adapterremoval
