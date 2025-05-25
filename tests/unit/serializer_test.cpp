// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "buffer.hpp"        // for buffer
#include "fastq.hpp"         // for fastq, fastq::ntrimmed, ACGTN, ACGT
#include "main.hpp"          // for VERSION
#include "read_group.hpp"    // for read_group
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for sample
#include "serializer.hpp"    // for serializer
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <ostream>           // for basic_ostream, operator<<, ostream
#include <string>            // for string
#include <string_view>       // for string_view

// Ignore nucleotide and quality strings
// spell-checker:ignoreRegExp /"[!-~]+"/g
// Ignore nucleotide comments
// spell-checker:ignoreRegExp /\W[acgtnACGTN]+\W/g

namespace adapterremoval {

namespace {

//! The program version with the leading 'v' removed; e.g. "3.0.1"
const std::string VERSION_NO_V{ VERSION.substr(1) };

constexpr std::string_view EXTREMELY_LONG_NAME =
  "123456789-123456789-123456789-123456789-123456789-123456789-123456789-"
  "123456789-123456789-123456789-123456789-123456789-123456789-123456789-"
  "123456789-123456789-123456789-123456789-123456789-123456789-123456789-"
  "123456789-123456789-123456789-123456789-1234";

// Basic named sample with barcodes, but otherwise no special properties
const sample BASIC_SAMPLE_WITH_BARCODES{ "foo",
                                         dna_sequence{ "ACGT" },
                                         dna_sequence{ "TGCA" },
                                         barcode_orientation::unspecified };

} // namespace

////////////////////////////////////////////////////////////////////////////////
// Implementation of enum debug serialization

std::ostream&
operator<<(std::ostream& os, const read_file& value)
{
  switch (value) {
    case read_file::mate_1:
      return os << "read_file::mate_1";
    case read_file::mate_2:
      return os << "read_file::mate_2";
    case read_file::merged:
      return os << "read_file::merged";
    case read_file::singleton:
      return os << "read_file::singleton";
    case read_file::discarded:
      return os << "read_file::discarded";
    case read_file::max:
    default:
      return os << "read_file{?}";
  }
}

///////////////////////////////////////////////////////////////////////////////
// FASTQ header serialization

TEST_CASE("Writing FASTQ header to buffer", "[serializer::fastq]")
{
  serializer s{ GENERATE(output_format::fastq, output_format::fastq_gzip) };

  // Sample information is only used in the case of --demultiplex-only
  if (GENERATE(true, false)) {
    auto sample{ BASIC_SAMPLE_WITH_BARCODES };

    // Read-group information is not used for FASTQ records
    if (GENERATE(true, false)) {
      sample.set_read_group(read_group{ "SM:my-sample" });
    }

    s.set_sample(sample);
  }

  buffer buf;

  SECTION("header")
  {
    // since FASTQ has no header, arguments should not matter
    s.header(buf, { "adapterremoval3", "--blah" });
    REQUIRE(buf == buffer{});
  }

  fastq record{ "record_1", "ACGTACGATA", "!$#$*68CGJ" };

  SECTION("basic record")
  {
    s.record(buf, record, read_meta{ read_type::se });
    REQUIRE(buf == "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n"_buffer);
  }
}

///////////////////////////////////////////////////////////////////////////////
// FASTQ record serialization

TEST_CASE("Writing FASTQ records when only demultiplexing")
{
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  // Read-group information is not used for FASTQ records
  if (GENERATE(true, false)) {
    sample.set_read_group(read_group{ "SM:my-sample" });
  }

  serializer s{ output_format::fastq };
  s.set_sample(sample);
  s.set_demultiplexing_only(true);

  buffer buf;
  fastq record{ "record_1", "ACGTACGATA", "!$#$*68CGJ" };
  s.record(buf, record, read_meta{ read_type::se });

  // The read header should include the barcodes when only demultiplexing
  REQUIRE(buf == "@record_1 BC:ACGT-TGCA\nACGTACGATA\n+\n!$#$*68CGJ\n"_buffer);
}

TEST_CASE("Writing FASTQ with mate separator")
{
  fastq record{ "record_1/1", "ACGTACGATA", "!$#$*68CGJ" };

  serializer s{ output_format::fastq };
  // This shouldn't matter, as mate separators are not removed for FASTQ reads
  s.set_mate_separator(GENERATE('\0', '/'));

  buffer buf;
  s.record(buf, record, read_meta(read_type::pe_1));
  REQUIRE(buf == "@record_1/1\nACGTACGATA\n+\n!$#$*68CGJ\n"_buffer);
}

TEST_CASE("FASTQ is the same for all read and sub-formats")
{
  const auto type = GENERATE(read_type::se,
                             read_type::se_fail,
                             read_type::pe_1,
                             read_type::pe_1_fail,
                             read_type::pe_2,
                             read_type::pe_2_fail,
                             read_type::singleton_1,
                             read_type::singleton_2,
                             read_type::merged,
                             read_type::merged_fail,
                             read_type::se,
                             read_type::se_fail,
                             read_type::pe_1,
                             read_type::pe_1_fail,
                             read_type::pe_2,
                             read_type::pe_2_fail,
                             read_type::singleton_1,
                             read_type::singleton_2,
                             read_type::merged,
                             read_type::merged_fail);
  const auto format = GENERATE(output_format::fastq, output_format::fastq_gzip);

  buffer buf;
  fastq record{ "record_1", "ACGTACGATA", "!$#$*68CGJ" };

  const read_meta meta{ type };
  const serializer s{ format };
  s.record(buf, record, meta);

  REQUIRE(buf == "@record_1\nACGTACGATA\n+\n!$#$*68CGJ\n"_buffer);
}

TEST_CASE("Writing empty FASTQ to buffer", "[serializer::fastq]")
{
  buffer buf;
  fastq record{ "record_1", "", "" };

  serializer s{ output_format::fastq };
  s.record(buf, record, read_meta(read_type::se));

  REQUIRE(buf == "@record_1\n\n+\n\n"_buffer);
}

TEST_CASE("Writing FASTQ with meta-data to buffer", "[serializer::fastq]")
{
  buffer buf;
  fastq record{ "record_1 length=5", "ACGTA", "68CGJ" };

  serializer s{ output_format::fastq };
  s.record(buf, record, read_meta(read_type::se));

  REQUIRE(buf == "@record_1 length=5\nACGTA\n+\n68CGJ\n"_buffer);
}

///////////////////////////////////////////////////////////////////////////////
// SAM header serialization

TEST_CASE("Writing SAM header to buffer", "[serializer::fastq]")
{
  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };

  if (GENERATE(true, false)) {
    s.set_sample(BASIC_SAMPLE_WITH_BARCODES);
  }

  s.header(buf, { "adapterremoval3", "--blah" });
  REQUIRE(buf == "@HD\tVN:1.6\tSO:unsorted\n@PG\tID:adapterremoval\t"
                 "PN:adapterremoval\tCL:adapterremoval3 --blah\t"
                 "VN:3.0.0-alpha3\n"_buffer);
}

TEST_CASE("Writing SAM read-group header to buffer", "[serializer::fastq]")
{
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.set_read_group({});

  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };
  s.set_sample(sample);

  s.header(buf, { "adapterremoval3", "--blah" });
  REQUIRE(buf == "@HD\tVN:1.6\tSO:unsorted\n@RG\tID:foo\tSM:foo\t"
                 "BC:ACGT-TGCA\n@PG\tID:adapterremoval\t"
                 "PN:adapterremoval\tCL:adapterremoval3 --blah\t"
                 "VN:3.0.0-alpha3\n"_buffer);
}

TEST_CASE("Writing SAM read-group header to buffer with multiple barcodes",
          "[serializer::fastq]")
{
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.add_barcodes(dna_sequence{ "TTGG" },
                      dna_sequence{ "AGTT" },
                      barcode_orientation::unspecified);
  sample.set_read_group({});

  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };
  s.set_sample(sample);

  s.header(buf, { "adapterremoval3", "--blah" });
  REQUIRE(buf == "@HD\tVN:1.6\tSO:unsorted\n@RG\tID:foo.1\tSM:foo\t"
                 "BC:ACGT-TGCA\n@RG\tID:foo.2\tSM:foo\tBC:TTGG-AGTT\n"
                 "@PG\tID:adapterremoval\tPN:adapterremoval\t"
                 "CL:adapterremoval3 --blah\tVN:3.0.0-alpha3\n"_buffer);
}

///////////////////////////////////////////////////////////////////////////////
// SAM record serialization

TEST_CASE("serialize SAM record without sample")
{
  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };
  s.set_demultiplexing_only(GENERATE(true, false));

  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };
  s.record(buf, record, read_meta{ read_type::se });

  REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                 "!$#$*68CGJ\n"_buffer);
}

TEST_CASE("serialize SAM record with sample")
{
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };
  s.set_demultiplexing_only(GENERATE(true, false));
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.set_read_group(read_group{});
  s.set_sample(sample);

  buffer buf;
  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };
  s.record(buf, record, read_meta{ read_type::se });

  REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t!$#$*68CGJ\t"
                 "RG:Z:foo\n"_buffer);
}

TEST_CASE("serialize SAM record with multiple barcodes")
{
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.add_barcodes(dna_sequence{ "TTGG" },
                      dna_sequence{ "AGTT" },
                      barcode_orientation::unspecified);
  sample.set_read_group({});

  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };
  s.set_sample(sample);

  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };

  SECTION("first/default barcodes")
  {
    s.record(buf, record, read_meta{ read_type::se });
    REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t!$#$*68CGJ\t"
                   "RG:Z:foo.1\n"_buffer);
  }

  SECTION("second barcodes")
  {
    auto meta = read_meta{ read_type::se }.barcode(1);
    s.record(buf, record, meta);
    REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t!$#$*68CGJ\t"
                   "RG:Z:foo.2\n"_buffer);
  }
}

TEST_CASE("serialize SAM record with mate separator")
{
  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };

  fastq record{ "record/1", "ACGTACGATA", "!$#$*68CGJ" };

  SECTION("without mate sep")
  {
    s.record(buf, record, read_meta{ read_type::se });
    REQUIRE(buf == "record/1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }

  SECTION("with mate sep")
  {
    s.set_mate_separator('/');
    s.record(buf, record, read_meta{ read_type::se });
    REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }
}

TEST_CASE("serialize read types for SAM")
{
  buffer buf;
  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };

  SECTION("single end")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::se, read_type::merged) });
    REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }

  SECTION("single end failed QC")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::se_fail, read_type::merged_fail) });
    REQUIRE(buf == "record\t516\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }

  SECTION("paired end mate 1, mate 2 passed/failed")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::pe_1, read_type::singleton_1) });
    REQUIRE(buf == "record\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }

  SECTION("paired end mate 1 failed, mate 2 passed/failed")
  {
    s.record(buf, record, read_meta{ read_type::pe_1_fail });
    REQUIRE(buf == "record\t589\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }

  SECTION("paired end mate 2, mate 1 passed/failed")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::pe_2, read_type::singleton_2) });
    REQUIRE(buf == "record\t141\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }

  SECTION("paired end mate 2 failed, mate 1 passed/failed")
  {
    s.record(buf, record, read_meta{ read_type::pe_2_fail });
    REQUIRE(buf == "record\t653\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                   "!$#$*68CGJ\n"_buffer);
  }
}

TEST_CASE("serialize empty SAM record")
{
  fastq record{ "record", "", "" };

  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };
  s.record(buf, record, read_meta{ read_type::se });

  REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"_buffer);
}

TEST_CASE("serialize SAM record from FASTQ with meta-data")
{
  buffer buf;
  serializer s{ GENERATE(output_format::sam, output_format::sam_gzip) };

  fastq record{ "record length=NA", "ACGTACGATA", "!$#$*68CGJ" };
  s.record(buf, record, read_meta{ read_type::se });

  // The meta-data is (currently) not serialized
  REQUIRE(buf == "record\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGATA\t"
                 "!$#$*68CGJ\n"_buffer);
}

///////////////////////////////////////////////////////////////////////////////
// BAM header serialization

TEST_CASE("Writing BAM header to buffer", "[serializer::fastq]")
{
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };

  // Sample information is only written if there is a read group, either
  // directly from the user or automatically assigned when demultiplexing
  if (GENERATE(true, false)) {
    sample sample{ BASIC_SAMPLE_WITH_BARCODES };
    s.set_sample(sample);
  }

  buffer buf;
  s.header(buf, { "adapterremoval3", "--blah" });
  REQUIRE(buf == "BAM\x01i\x00\x00\x00@HD\tVN:1.6\tSO:unsorted\n"
                 "@PG\tID:adapterremoval\tPN:adapterremoval\t"
                 "CL:adapterremoval3 --blah\tVN:3.0.0-alpha3\n"
                 "\x00\x00\x00\x00"_buffer);
}

TEST_CASE("Writing BAM read-group header to buffer", "[serializer::fastq]")
{

  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.set_read_group(read_group{ "LB:lib" });

  serializer s{ GENERATE(output_format::bam, output_format::ubam) };
  s.set_sample(sample);

  buffer buf;
  s.header(buf, { "adapterremoval3", "--blah" });
  REQUIRE(buf == "BAM\x01\x8f\x00\x00\x00@HD\tVN:1.6\tSO:unsorted\n@RG\t"
                 "ID:foo\tLB:lib\tSM:foo\tBC:ACGT-TGCA\n@PG\t"
                 "ID:adapterremoval\tPN:adapterremoval\tCL:"
                 "adapterremoval3 --blah\tVN:3.0.0-alpha3\n"
                 "\x00\x00\x00\x00"_buffer);
}

TEST_CASE("Writing BAM read-group header to buffer with multiple barcodes",
          "[serializer::fastq]")
{
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.add_barcodes(dna_sequence{ "TTGG" },
                      dna_sequence{ "AGTT" },
                      barcode_orientation::unspecified);
  sample.set_read_group({});

  serializer s{ GENERATE(output_format::bam, output_format::ubam) };
  s.set_sample(sample);

  buffer buf;
  s.header(buf, { "adapterremoval3", "--blah" });
  REQUIRE(buf == "BAM\x01\xab\x00\x00\x00@HD\tVN:1.6\tSO:unsorted\n@RG\t"
                 "ID:foo.1\tSM:foo\tBC:ACGT-TGCA\n@RG\tID:foo.2\t"
                 "SM:foo\tBC:TTGG-AGTT\n@PG\tID:adapterremoval\t"
                 "PN:adapterremoval\tCL:adapterremoval3 "
                 "--blah\tVN:3.0.0-alpha3\n\x00\x00\x00\x00"_buffer);
}

////////////////////////////////////////////////////////////////////////////////
// BAM record serialization

TEST_CASE("serialize BAM record without sample")
{
  buffer buf;
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };
  s.set_demultiplexing_only(GENERATE(true, false));

  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };
  s.record(buf, record, read_meta{ read_type::se });

  REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                 "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                 "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                 "\x03\t\x15\x17\"&)"_buffer);
}

TEST_CASE("serialize BAM record with uneven length sequence")
{
  buffer buf;
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };
  s.set_demultiplexing_only(GENERATE(true, false));

  fastq record{ "record", "ACGTACGAT", "!$#$*68CG" };
  s.record(buf, record, read_meta{ read_type::se });

  REQUIRE(buf == "5\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                 "\x00\x00\x04\x00\t\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                 "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x80\x00\x03\x02"
                 "\x03\t\x15\x17\"&"_buffer);
}

TEST_CASE("serialize BAM record with sample")
{
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };
  s.set_demultiplexing_only(GENERATE(true, false));
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.set_read_group(read_group{});
  s.set_sample(sample);

  buffer buf;
  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };
  s.record(buf, record, read_meta{ read_type::se });

  REQUIRE(buf == "=\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                 "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                 "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                 "\x03\t\x15\x17\"&)RGZfoo\x00"_buffer);
}

TEST_CASE("serialize BAM record with multiple barcodes")
{
  sample sample{ BASIC_SAMPLE_WITH_BARCODES };
  sample.add_barcodes(dna_sequence{ "TTGG" },
                      dna_sequence{ "AGTT" },
                      barcode_orientation::unspecified);
  sample.set_read_group({});

  buffer buf;
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };
  s.set_sample(sample);

  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };

  SECTION("first/default barcodes")
  {
    s.record(buf, record, read_meta{ read_type::se });
    REQUIRE(buf == "?\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                   "\x03\t\x15\x17\"&)RGZfoo.1\x00"_buffer);
  }

  SECTION("second barcodes")
  {
    auto meta = read_meta{ read_type::se }.barcode(1);
    s.record(buf, record, meta);
    REQUIRE(buf == "?\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                   "\x03\t\x15\x17\"&)RGZfoo.2\x00"_buffer);
  }
}

TEST_CASE("serialize BAM record with mate separator")
{
  buffer buf;
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };

  fastq record{ "record/1", "ACGTACGATA", "!$#$*68CGJ" };

  SECTION("without mate sep")
  {
    s.record(buf, record, read_meta{ read_type::se });
    REQUIRE(buf == "8\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\t\x00H\x12"
                   "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record/1\x00\x12H\x12\x41\x81\x00\x03"
                   "\x02\x03\t\x15\x17\"&)"_buffer);
  }

  SECTION("with mate sep")
  {
    s.set_mate_separator('/');
    s.record(buf, record, read_meta{ read_type::se });
    REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                   "\x03\t\x15\x17\"&)"_buffer);
  }
}

TEST_CASE("serialize read types for BAM")
{
  buffer buf;
  fastq record{ "record", "ACGTACGATA", "!$#$*68CGJ" };
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };

  SECTION("single end")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::se, read_type::merged) });
    REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                   "\x03\t\x15\x17\"&)"_buffer);
  }

  SECTION("single end failed QC")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::se_fail, read_type::merged_fail) });
    REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00\x04\x02\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                   "\x03\t\x15\x17\"&)"_buffer);
  }

  SECTION("paired end mate 1, mate 2 passed/failed")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::pe_1, read_type::singleton_1) });
    REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00M\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff"
                   "\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02\x03"
                   "\t\x15\x17\"&)"_buffer);
  }

  SECTION("paired end mate 1 failed, mate 2 passed/failed")
  {
    s.record(buf, record, read_meta{ read_type::pe_1_fail });
    REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00M\x02\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff"
                   "\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02\x03"
                   "\t\x15\x17\"&)"_buffer);
  }

  SECTION("paired end mate 2, mate 1 passed/failed")
  {
    s.record(buf,
             record,
             read_meta{ GENERATE(read_type::pe_2, read_type::singleton_2) });
    REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00\x8d\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                   "\x03\t\x15\x17\"&)"_buffer);
  }

  SECTION("paired end mate 2 failed, mate 1 passed/failed")
  {
    s.record(buf, record, read_meta{ read_type::pe_2_fail });
    REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                   "\x00\x00\x8d\x02\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                   "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                   "\x03\t\x15\x17\"&)"_buffer);
  }
}

TEST_CASE("serialize empty BAM record")
{
  fastq record{ "record", "", "" };

  buffer buf;
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };
  s.record(buf, record, read_meta{ read_type::se });

  REQUIRE(buf == "'\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                 "\x00\x00\x04\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                 "\xff\x00\x00\x00\x00record\x00"_buffer);
}

TEST_CASE("serialize BAM record from FASTQ with meta-data")
{
  buffer buf;
  serializer s{ GENERATE(output_format::bam, output_format::ubam) };

  fastq record{ "record length=NA", "ACGTACGATA", "!$#$*68CGJ" };
  s.record(buf, record, read_meta{ read_type::se });

  // The meta-data is (currently) not serialized
  REQUIRE(buf == "6\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x07\x00H\x12"
                 "\x00\x00\x04\x00\n\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff"
                 "\xff\x00\x00\x00\x00record\x00\x12H\x12\x41\x81\x00\x03\x02"
                 "\x03\t\x15\x17\"&)"_buffer);
}

////////////////////////////////////////////////////////////////////////////////
// Common SAM/BAM tests

TEST_CASE("Serializing too long read names")
{
  // The maximum allowed read name length for SAM/BAM
  static_assert(EXTREMELY_LONG_NAME.size() == 254);

  auto name = std ::string{ EXTREMELY_LONG_NAME };
  fastq record_254{ name, "ACGT", "!!!!" };
  fastq record_255{ name + "5", "ACGT", "!!!!" };

  const read_meta meta{ read_type::se };
  buffer buf;

  SECTION("fastq is always valid")
  {
    serializer s{ GENERATE(output_format::fastq, output_format::fastq_gzip) };
    CHECK_NOTHROW(s.record(buf, record_254, meta));
    CHECK_NOTHROW(s.record(buf, record_255, meta));
  }

  SECTION("SAM/BAM allows only 254 characters")
  {
    serializer s{ GENERATE(output_format::sam,
                           output_format::sam_gzip,
                           output_format::bam,
                           output_format::ubam) };

    CHECK_NOTHROW(s.record(buf, record_254, meta));
    REQUIRE_THROWS_WITH(s.record(buf, record_255, meta),
                        Catch::Contains("Cannot encode read as SAM/BAM; read "
                                        "name is longer than 254 characters"));
  }
}

TEST_CASE("invalid read names")
{
  const read_meta meta{ read_type::se };
  fastq record{ "invalid\tname", "ACGT", "!!!!" };
  buffer buf;

  SECTION("FASTQ allows anything")
  {
    serializer s{ GENERATE(output_format::fastq, output_format::fastq_gzip) };
    CHECK_NOTHROW(s.record(buf, record, meta));
  }

  SECTION("SAM/BAM limits valid characters")
  {
    serializer s{ GENERATE(output_format::sam,
                           output_format::sam_gzip,
                           output_format::bam,
                           output_format::ubam) };

    REQUIRE_THROWS_WITH(s.record(buf, record, meta),
                        Catch::Contains("Cannot encode read as SAM/BAM; read "
                                        "name contains characters other than "
                                        "the allowed"));
  }
}

////////////////////////////////////////////////////////////////////////////////
// Tests of read meta data

TEST_CASE("read type to file type mapping")
{
  SECTION("--out-file1")
  {
    read_meta meta{ GENERATE(read_type::se, read_type::pe_1) };
    CHECK(meta.get_file() == read_file::mate_1);
  }

  SECTION("--out-file2")
  {
    read_meta meta{ GENERATE(read_type::pe_2) };
    CHECK(meta.get_file() == read_file::mate_2);
  }

  SECTION("--out-singleton")
  {
    read_meta meta(GENERATE(read_type::singleton_1, read_type::singleton_2));
    CHECK(meta.get_file() == read_file::singleton);
  }

  SECTION("--out-singleton")
  {
    read_meta meta(read_type::merged);
    CHECK(meta.get_file() == read_file::merged);
  }

  SECTION("--out-discarded")
  {
    read_meta meta(GENERATE(read_type::se_fail,
                            read_type::pe_1_fail,
                            read_type::pe_2_fail,
                            read_type::merged_fail));
    CHECK(meta.get_file() == read_file::discarded);
  }
}

} // namespace adapterremoval
