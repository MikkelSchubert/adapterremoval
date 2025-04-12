// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "serializer.hpp"  // declarations
#include "buffer.hpp"      // for buffer
#include "commontypes.hpp" // for output_format
#include "debug.hpp"       // for AR_REQUIRE, AR_FAIL
#include "errors.hpp"      // for serializing_error
#include "fastq.hpp"       // for fastq
#include "fastq_enc.hpp"   // for PHRED_OFFSET_MIN
#include "main.hpp"        // for VERSION
#include "strutils.hpp"    // for join_text
#include <sstream>         // for ostringstream
#include <string_view>     // for string_view

namespace adapterremoval {

using namespace std::literals;

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
read_type_to_sam(read_type flags)
{
  switch (flags) {
    case read_type::se:
    case read_type::merged:
      return "4";
    case read_type::se_fail:
    case read_type::merged_fail:
      return "516";
    case read_type::pe_1:
    case read_type::singleton_1:
      return "77";
    case read_type::pe_1_fail:
      return "589";
    case read_type::pe_2:
    case read_type::singleton_2:
      return "141";
    case read_type::pe_2_fail:
      return "653";
    default:                          // GCOVR_EXCL_LINE
      AR_FAIL("invalid fastq flags"); // GCOVR_EXCL_LINE
  }
}

uint16_t
read_type_to_bam(read_type flags)
{
  switch (flags) {
    case read_type::se:
    case read_type::merged:
      return 4;
    case read_type::se_fail:
    case read_type::merged_fail:
      return 516;
    case read_type::pe_1:
    case read_type::singleton_1:
      return 77;
    case read_type::pe_1_fail:
      return 589;
    case read_type::pe_2:
    case read_type::singleton_2:
      return 141;
    case read_type::pe_2_fail:
      return 653;
    default:                          // GCOVR_EXCL_LINE
      AR_FAIL("invalid fastq flags"); // GCOVR_EXCL_LINE
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

std::string_view
prepare_name(const fastq& record, const char mate_separator)
{
  auto name = record.name(mate_separator);
  if (name.length() > 254) {
    std::ostringstream os;
    os << "Cannot encode read as SAM/BAM; read name is longer than 254 "
       << "characters: len(" << log_escape(name) << ") == " << name.length();

    throw serializing_error(os.str());
  }

  for (const char c : name) {
    if ((c < '!' || c > '?') && (c < 'A' || c > '~')) {
      std::ostringstream os;
      os << "Cannot encode read as SAM/BAM; read name contains characters "
         << "other than the allowed [!-?A-~]: " << log_escape(name);

      throw serializing_error(os.str());
    }
  }

  return name;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for `fastq_serializer`

void
serializer::fastq_record(buffer& buf,
                         const fastq& record,
                         const read_meta& /* meta */,
                         const sample_sequences& sequences) const
{
  buf.append(record.header());
  if (m_demultiplexing_only && !sequences.barcode_1.empty()) {
    buf.append(" BC:");
    buf.append(sequences.barcode_1);
    if (sequences.barcode_2.length()) {
      buf.append_u8('-');
      buf.append(sequences.barcode_2);
    }
  }
  buf.append_u8('\n');
  buf.append(record.sequence());
  buf.append("\n+\n"sv);
  buf.append(record.qualities());
  buf.append_u8('\n');
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `sam_serializer`

void
serializer::sam_header(buffer& buf, const string_vec& args, const sample& s)
{
  buf.append(create_sam_header(args, s));
}

void
serializer::sam_record(buffer& buf,
                       const fastq& record,
                       const read_meta& meta,
                       const sample_sequences& sequences) const
{
  buf.append(prepare_name(record, m_mate_separator)); // 1. QNAME
  buf.append_u8('\t');
  buf.append(read_type_to_sam(meta.m_type)); // 2. FLAG
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
serializer::bam_header(buffer& buf, const string_vec& args, const sample& s)
{
  const auto sam_header = create_sam_header(args, s);

  buf.append(BAM_HEADER);            // magic
  buf.append_u32(sam_header.size()); // l_text
  buf.append(sam_header);            // terminating NUL not required
  buf.append_u32(0);                 // n_ref
}

void
serializer::bam_record(buffer& buf,
                       const fastq& record,
                       const read_meta& meta,
                       const sample_sequences& sequences) const
{
  const size_t block_size_pos = buf.size();
  buf.append_u32(0);  // block size (preliminary)
  buf.append_i32(-1); // refID
  buf.append_i32(-1); // pos

  const auto name = prepare_name(record, m_mate_separator);
  buf.append_u8(name.length() + 1); // l_read_name
  buf.append_u8(0);                 // mapq
  buf.append_u16(4680);             // bin (c.f. specification 4.2.1)
  buf.append_u16(0);                // n_cigar
  buf.append_u16(read_type_to_bam(meta.m_type)); // flags

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
// Implementations for `read_meta`

read_file
read_meta::get_file() const noexcept
{
  switch (m_type) {
    case read_type::se:
    case read_type::pe_1:
      return read_file::mate_1;
    case read_type::pe_2:
      return read_file::mate_2;
    case read_type::se_fail:
    case read_type::pe_1_fail:
    case read_type::pe_2_fail:
    case read_type::merged_fail:
      return read_file::discarded;
    case read_type::singleton_1:
    case read_type::singleton_2:
      return read_file::singleton;
    case read_type::merged:
      return read_file::merged;
    default:                         // GCOVR_EXCL_LINE
      AR_FAIL("invalid read flags"); // GCOVR_EXCL_LINE
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `serializer`

serializer::serializer(output_format format)
  : m_format(format)
{
}

void
serializer::header(buffer& buf, const string_vec& args) const
{
  switch (m_format) {
    case output_format::fastq:
    case output_format::fastq_gzip:
      // No header
      break;
    case output_format::sam:
    case output_format::sam_gzip:
      sam_header(buf, args, m_sample);
      break;
    case output_format::bam:
    case output_format::ubam:
      bam_header(buf, args, m_sample);
      break;
    default:                            // GCOVR_EXCL_LINE
      AR_FAIL("invalid output format"); // GCOVR_EXCL_LINE
  }
}

void
serializer::record(buffer& buf, fastq&& record, read_meta meta) const
{
  const auto& sequences = m_sample.at(meta.m_barcode);
  switch (m_format) {
    case output_format::fastq:
    case output_format::fastq_gzip:
      fastq_record(buf, record, meta, sequences);
      break;
    case output_format::sam:
    case output_format::sam_gzip:
      sam_record(buf, record, meta, sequences);
      break;
    case output_format::bam:
    case output_format::ubam:
      bam_record(buf, record, meta, sequences);
      break;
    default:                            // GCOVR_EXCL_LINE
      AR_FAIL("invalid output format"); // GCOVR_EXCL_LINE
  }
}

} // namespace adapterremoval
