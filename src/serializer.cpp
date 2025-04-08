// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "serializer.hpp"  // for fastq_serializer
#include "buffer.hpp"      // for buffer
#include "commontypes.hpp" // for output_format
#include "debug.hpp"       // for AR_REQUIRE, AR_FAIL
#include "fastq.hpp"       // for fastq
#include "fastq_enc.hpp"   // for PHRED_OFFSET_MIN
#include "main.hpp"        // for VERSION
#include "strutils.hpp"    // for join_text
#include <string_view>     // for string_view

namespace adapterremoval {

class userconfig;

namespace {

//! Standard header for BAM files prior to compression
constexpr std::string_view BAM_HEADER{ "BAM\1", 4 };
//! Standard header for SAM/BAM files
constexpr std::string_view SAM_HEADER = "@HD\tVN:1.6\tSO:unsorted\n";

/**
 * Flags mapping onto SAM/BAM flags
 *
 * 0x1 = read paired
 * 0x4 = read unmapped
 * 0x8 = mate unmapped
 * 0x40 = mate 1
 * 0x80 = mate 2
 * 0x200 = failed QC
 */

std::string_view
flags_to_sam(fastq_flags flags)
{
  switch (flags) {
    case fastq_flags::se:
      return "4";
    case fastq_flags::se_fail:
      return "516";
    case fastq_flags::pe_1:
      return "77";
    case fastq_flags::pe_1_fail:
      return "589";
    case fastq_flags::pe_2:
      return "141";
    case fastq_flags::pe_2_fail:
      return "653";
    default:
      AR_FAIL("invalid fastq flags");
  }
}

uint16_t
flags_to_bam(fastq_flags flags)
{
  switch (flags) {
    case fastq_flags::se:
      return 4;
    case fastq_flags::se_fail:
      return 516;
    case fastq_flags::pe_1:
      return 77;
    case fastq_flags::pe_1_fail:
      return 589;
    case fastq_flags::pe_2:
      return 141;
    case fastq_flags::pe_2_fail:
      return 653;
    default:
      AR_FAIL("invalid fastq flags");
  }
}

void
sequence_to_bam(buffer& buf, const std::string& seq)
{
  const auto size = buf.size();

  uint8_t pair = 0;
  for (size_t i = 0; i < seq.length(); ++i) {
    pair = (pair << 4) | "\0\1\0\2\10\0\20\4"[seq[i] & 0x7];

    if (i % 2) {
      buf.append_u8(pair);
      pair = 0;
    }
  }

  if (seq.length() % 2) {
    buf.append_u8(pair << 4);
  }

  AR_REQUIRE(buf.size() - size == (seq.length() + 1) / 2);
}

void
qualities_to_bam(buffer& buf, const std::string& quals)
{
  for (const auto c : quals) {
    buf.append_u8(c - PHRED_OFFSET_MIN);
  }
}

std::string
create_sam_header(const string_vec& args, const sample& s)
{
  std::string header{ SAM_HEADER };

  // @RG
  for (const auto& it : s) {
    if (it.has_read_group) {
      header.append(it.read_group_.header());
      header.append("\n");
    }
  }

  // @PG
  header.append("@PG\tID:adapterremoval\tPN:adapterremoval\tCL:");
  header.append(join_text(args, " "));
  header.append("\tVN:");
  header.append(VERSION.substr(1)); // version without leading v
  header.append("\n");

  return header;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for `fastq_serializer`

void
fastq_serializer::header(buffer& /* buf */,
                         const string_vec& /* args */,
                         const sample& /* s */)
{
}

void
fastq_serializer::record(buffer& buf,
                         const fastq& record,
                         const sample_sequences& sequences,
                         const serializer_settings& settings)
{
  buf.append(record.header());
  if (settings.demultiplexing_only) {
    buf.append(" BC:");
    buf.append(sequences.barcode_1);
    if (sequences.barcode_2.length()) {
      buf.append_u8('-');
      buf.append(sequences.barcode_2);
    }
  }
  buf.append_u8('\n');
  buf.append(record.sequence());
  buf.append("\n+\n", 3);
  buf.append(record.qualities());
  buf.append_u8('\n');
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `sam_serializer`

void
sam_serializer::header(buffer& buf, const string_vec& args, const sample& s)
{
  buf.append(create_sam_header(args, s));
}

void
sam_serializer::record(buffer& buf,
                       const fastq& record,
                       const sample_sequences& sequences,
                       const serializer_settings& settings)
{
  buf.append(record.name(settings.mate_separator)); // 1. QNAME
  buf.append_u8('\t');
  buf.append(flags_to_sam(settings.flags)); // 2. FLAG
  buf.append("\t"
             "*\t" // 3. RNAME
             "0\t" // 4. POS
             "0\t" // 5. MAPQ
             "*\t" // 6. CIGAR
             "*\t" // 7. RNEXT
             "0\t" // 8. PNEXT
             "0\t" // 9. TLEN
  );
  if (record.length()) {
    buf.append(record.sequence()); // 10. SEQ
    buf.append_u8('\t');
    buf.append(record.qualities()); // 11. QUAL
  } else {
    buf.append("*\t" // 10. SEQ
               "*"   // 11. QUAL
    );
  }

  if (sequences.has_read_group) {
    buf.append("\tRG:Z:");
    buf.append(sequences.read_group_.id());
  }

  buf.append("\n");
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `bam_serializer`

void
bam_serializer::header(buffer& buf, const string_vec& args, const sample& s)
{
  const auto sam_header = create_sam_header(args, s);

  buf.append(BAM_HEADER);            // magic
  buf.append_u32(sam_header.size()); // l_text
  buf.append(sam_header);            // terminating NUL not required
  buf.append_u32(0);                 // n_ref
}

void
bam_serializer::record(buffer& buf,
                       const fastq& record,
                       const sample_sequences& sequences,
                       const serializer_settings& settings)
{
  const size_t block_size_pos = buf.size();
  buf.append_u32(0);  // block size (preliminary)
  buf.append_i32(-1); // refID
  buf.append_i32(-1); // pos

  const auto name = record.name(settings.mate_separator).substr(0, 255);
  buf.append_u8(name.length() + 1); // l_read_name
  buf.append_u8(0);                 // mapq
  buf.append_u16(4680);             // bin (c.f. specification 4.2.1)
  buf.append_u16(0);                // n_cigar
  buf.append_u16(flags_to_bam(settings.flags)); // flags

  buf.append_u32(record.length()); // l_seq
  buf.append_i32(-1);              // next_refID
  buf.append_i32(-1);              // next_pos
  buf.append_i32(0);               // tlen

  buf.append(name); // read_name + NUL terminator
  buf.append_u8(0);
  // no cigar operations
  sequence_to_bam(buf, record.sequence());
  qualities_to_bam(buf, record.qualities());

  if (sequences.has_read_group) {
    // RG:Z:${ID} tag
    buf.append("RGZ");
    buf.append(sequences.read_group_.id());
    buf.append_u8(0); // NUL
  }

  const size_t block_size = buf.size() - block_size_pos - 4;
  buf.put_u32(block_size_pos, block_size); // block size (final)
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `serializer`

serializer::serializer(output_format format)
{
  switch (format) {
    case output_format::fastq:
    case output_format::fastq_gzip:
      m_header = fastq_serializer::header;
      m_record = fastq_serializer::record;
      break;
    case output_format::sam:
    case output_format::sam_gzip:
      m_header = sam_serializer::header;
      m_record = sam_serializer::record;
      break;
    case output_format::bam:
    case output_format::ubam:
      m_header = bam_serializer::header;
      m_record = bam_serializer::record;
      break;
    default:
      AR_FAIL("invalid output format");
  }
}

void
serializer::header(buffer& buf, const string_vec& args) const
{
  m_header(buf, args, m_sample);
}

void
serializer::record(buffer& buf,
                   const fastq& record,
                   const fastq_flags flags,
                   const size_t barcode) const
{
  m_record(
    buf,
    record,
    m_sample.at(barcode),
    serializer_settings{ flags, m_mate_separator, m_demultiplexing_only });
}

} // namespace adapterremoval
