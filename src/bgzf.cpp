// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "bgzf.hpp"    // declarations
#include "buffer.hpp"  // for buffer
#include "debug.hpp"   // for AR_REQUIRE
#include "deflate.hpp" // for crc32, deflate_block_uncompressed, ...

namespace adapterremoval {

void
bgzf_block_compress(const buffer& src,
                    buffer& dst,
                    const unsigned compression_level)
{
  AR_REQUIRE(src.size() <= BGZF_MAX_INPUT_SIZE);

  dst.resize(0);
  dst.reserve(BGZF_MAX_BLOCK_SIZE);
  dst.append(BGZF_HEADER);
  constexpr auto offset = BGZF_HEADER.size();

  // Leave room for CRC32 and ISIZE appended below
  dst.resize(BGZF_MAX_BLOCK_SIZE - 8);

  switch (compression_level) {
    case 0:
      deflate_block_uncompressed(src, dst, offset);
      break;
    case ISAL_COMPRESSION_LEVEL:
      deflate_block_isal_1(src, dst, offset, true);
      break;
    default:
      deflate_block_libdeflate(src, dst, offset, compression_level);
      break;
  }

  dst.append_u32(crc32(src)); // CRC32
  dst.append_u32(src.size()); // ISIZE

  AR_REQUIRE(dst.size() <= BGZF_MAX_BLOCK_SIZE);

  // Write the final block size; -1 to fit 1..65536 in 16 bit
  dst.put_u16(16, dst.size() - 1);
}

} // namespace adapterremoval
