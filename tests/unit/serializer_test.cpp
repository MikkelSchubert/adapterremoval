// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "buffer.hpp"        // for buffer
#include "commontypes.hpp"   // for fastq_vec
#include "fastq.hpp"         // for fastq, fastq::ntrimmed, ACGTN, ACGT
#include "sequence_sets.hpp" // sample_seequences
#include "serializer.hpp"    // for fastq_serializer
#include "testing.hpp"       // for catch.hpp, StringMaker

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
  fastq_serializer::record(buf,
                           record,
                           sample_sequences{},
                           serializer_settings{});

  REQUIRE(
    std::string_view(reinterpret_cast<const char*>(buf.data()), buf.size()) ==
    "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n");
}

TEST_CASE("Writing_to_stream_phred_33_explicit", "[fastq::fastq]")
{
  const fastq record = fastq("record_1", "ACGTACGATA", "!$#$*68CGJ");
  buffer buf;
  fastq_serializer::record(buf,
                           record,
                           sample_sequences{},
                           serializer_settings{});

  REQUIRE(
    std::string_view(reinterpret_cast<const char*>(buf.data()), buf.size()) ==
    "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n");
}

} // namespace adapterremoval
