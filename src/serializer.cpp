/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "serializer.hpp"  // for fastq_serialzier
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
create_sam_header(const string_vec& args, const read_group& rg)
{
  std::string header{ SAM_HEADER };

  // @RG
  header.append(rg.header());
  header.append("\n");

  // @PG
  header.append("@PG\tID:adapterremoval\tPN:adapterremoval\tCL:");
  header.append(join_text(args, " "));
  header.append("\tVN:");
  header.append(VERSION.substr(1)); // version without leading v
  header.append("\n");

  return header;
}

/** Unescapes the escape sequences "\\" and "\t" in a read-group string */
std::string
unescape_read_group(std::string_view value)
{
  std::string result;

  bool in_escape = false;
  for (auto c : value) {
    if (in_escape) {
      if (c == '\\') {
        result.push_back('\\');
      } else if (c == 't') {
        result.push_back('\t');
      } else {
        throw std::invalid_argument("invalid escape sequence " +
                                    log_escape(std::string("\\") + c));
      }

      in_escape = false;
    } else if (c == '\\') {
      in_escape = true;
    } else {
      result.push_back(c);
    }
  }

  if (in_escape) {
    throw std::invalid_argument("incomplete escape sequence at end of string");
  }

  return result;
}

constexpr bool
is_ascii_alpha(const char c)
{
  return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

constexpr bool
is_ascii_alphanum(const char c)
{
  return is_ascii_alpha(c) || (c >= '0' && c <= '9');
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for `read_group`

read_group::read_group()
  : m_header{ "@RG\tID:1\tPG:adapterremoval" }
  , m_id{ "ID:1" }
{
}

read_group::read_group(std::string_view value_)
  : m_header{ "@RG" }
{
  using invalid = std::invalid_argument;
  auto value = unescape_read_group(value_);

  // It's not unreasonable to except users to try to specify a full @RG line
  if (starts_with(value, "@RG\t")) {
    value = value.substr(4);
  }

  std::string program{ "PG:adapterremoval" };

  if (!value.empty()) {
    for (const auto& field : split_text(value, '\t')) {
      for (const auto c : field) {
        if (c < ' ' || c > '~') {
          throw invalid("only characters in the range ' ' to '~' are allowed");
        }
      }

      if (field.size() < 4) {
        throw invalid("tags must be at least 4 characters long");
      } else if (!is_ascii_alpha(field.at(0))) {
        throw invalid("first character in tag name must be a letter");
      } else if (!is_ascii_alphanum(field.at(1))) {
        throw invalid(
          "second character in tag name must be a letter or number");
      } else if (field.at(2) != ':') {
        throw invalid("third character in tag must be a colon");
      } else if (starts_with(field, "ID:")) {
        if (!m_id.empty()) {
          throw invalid("multiple ID tags found");
        }

        m_id = field.substr(3);
      } else if (starts_with(field, "PG:")) {
        // Allow the user to override the PG field, just in case
        program = field;
      }
    }
  }

  // Generic read-group key if the user did not supply one
  if (m_id.empty()) {
    m_id = "1";
    m_header += "\tID:1";
  }

  if (!value.empty()) {
    m_header += "\t";
    m_header += value;
  }

  m_header.append("\t");
  m_header.append(program);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `fastq_serializer`

void
fastq_serializer::header(buffer& /* buf */,
                         const string_vec& /* args */,
                         const read_group& /* rg */)
{
}

void
fastq_serializer::record(buffer& buf,
                         const fastq& record,
                         fastq_flags /* flags */,
                         char /* mate_separator */,
                         const read_group& /* rg */)
{
  buf.append_u8('@');
  buf.append(record.header());
  buf.append_u8('\n');
  buf.append(record.sequence());
  buf.append("\n+\n", 3);
  buf.append(record.qualities());
  buf.append_u8('\n');
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `sam_serializer`

void
sam_serializer::header(buffer& buf,
                       const string_vec& args,
                       const read_group& rg)
{
  buf.append(create_sam_header(args, rg));
}

void
sam_serializer::record(buffer& buf,
                       const fastq& record,
                       fastq_flags flags,
                       char mate_separator,
                       const read_group& rg)
{
  buf.append(record.name(mate_separator)); // 1. QNAME
  buf.append_u8('\t');
  buf.append(flags_to_sam(flags)); // 2. FLAG
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

  buf.append("\tRG:Z:");
  buf.append(rg.id());
  buf.append("\tPG:Z:adapterremoval\n");
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `bam_serializer`

void
bam_serializer::header(buffer& buf,
                       const string_vec& args,
                       const read_group& rg)
{
  const auto sam_header = create_sam_header(args, rg);

  buf.append(BAM_HEADER);            // magic
  buf.append_u32(sam_header.size()); // l_text
  buf.append(sam_header);            // terminating NUL not required
  buf.append_u32(0);                 // n_ref
}

void
bam_serializer::record(buffer& buf,
                       const fastq& record,
                       fastq_flags flags,
                       char mate_separator,
                       const read_group& rg)
{
  const size_t block_size_pos = buf.size();
  buf.append_u32(0);  // block size (preliminary)
  buf.append_i32(-1); // refID
  buf.append_i32(-1); // pos

  const auto name = record.name(mate_separator).substr(0, 255);
  buf.append_u8(name.length() + 1);    // l_read_name
  buf.append_u8(0xFF);                 // mapq
  buf.append_u16(4680);                // bin (c.f. specification 4.2.1)
  buf.append_u16(0);                   // n_cigar
  buf.append_u16(flags_to_bam(flags)); // flags

  buf.append_u32(record.length()); // l_seq
  buf.append_i32(-1);              // next_refID
  buf.append_i32(-1);              // next_pos
  buf.append_i32(0);               // tlen

  buf.append(name); // read_name + NUL terminator
  buf.append_u8(0);
  // no cigar operations
  sequence_to_bam(buf, record.sequence());
  qualities_to_bam(buf, record.qualities());

  // PG:Z:adapterremoval tag
  buf.append("RGZ");
  buf.append(rg.id());
  buf.append_u8(0); // NUL

  // PG:Z:adapterremoval tag
  buf.append("PGZadapterremoval");
  buf.append_u8(0); // NUL

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
  m_header(buf, args, m_read_group);
}

void
serializer::record(buffer& buf,
                   const fastq& record,
                   const fastq_flags flags) const
{
  m_record(buf, record, flags, m_mate_separator, m_read_group);
}

} // namespace adapterremoval
