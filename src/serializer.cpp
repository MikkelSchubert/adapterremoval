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

namespace adapterremoval {

class userconfig;

namespace {

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

void
write_fastq_record(buffer& buf, const fastq& record)
{
  buf.append_u8('@');
  buf.append(record.header());
  buf.append_u8('\n');
  buf.append(record.sequence());
  buf.append("\n+\n", 3);
  buf.append(record.qualities());
  buf.append_u8('\n');
}

void
write_sam_record(buffer& buf,
                 const fastq& record,
                 const fastq_flags flags,
                 char mate_separator)
{
  if (mate_separator) {
    switch (flags) {
      case fastq_flags::se:
      case fastq_flags::se_fail:
        mate_separator = '\0';
        break;
      case fastq_flags::pe_1:
      case fastq_flags::pe_1_fail:
      case fastq_flags::pe_2:
      case fastq_flags::pe_2_fail:
        break;
      default:
        AR_FAIL("invalid fastq flags");
    }
  }

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
    buf.append_u8('\n');
  } else {
    buf.append("*\t" // 10. SEQ
               "*\n" // 11. QUAL
    );
  }
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for `fastq_serializer`

void
fastq_serializer::header(buffer& buf, const output_format format)
{
  switch (format) {
    case output_format::fastq:
    case output_format::fastq_gzip:
      break;
    case output_format::sam:
    case output_format::sam_gzip:
      buf.append("@HD\tVN:1.6\tSO:unsorted\n");
      break;
    default:
      AR_FAIL("invalid output format");
  }
}

void
fastq_serializer::record(buffer& buf,
                         const fastq& record,
                         const output_format format,
                         const fastq_flags flags,
                         const char mate_separator)
{
  switch (format) {
    case output_format::fastq:
    case output_format::fastq_gzip:
      return write_fastq_record(buf, record);
    case output_format::sam:
    case output_format::sam_gzip:
      return write_sam_record(buf, record, flags, mate_separator);
    default:
      AR_FAIL("invalid output format");
  }
}

} // namespace adapterremoval
