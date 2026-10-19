// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "bgzf.hpp"     // declarations
#include "buffer.hpp"   // for buffer
#include "debug.hpp"    // for AR_REQUIRE
#include "deflate.hpp"  // for calculate_crc32, deflate_block_uncompressed, ...
#include "errors.hpp"   // for gzip_error
#include "strutils.hpp" // for stringify
#include <cstdint>      // for uint8_t
#include <isa-l.h>      // for isal_inflate_init, isal_inflate_stateless, ...
#include <memory>       // for make_unique
#include <sstream>      // for ostringstream
#include <string_view>  // for string_view

namespace adapterremoval {

namespace {

[[noreturn]] void
throw_isal_inflate_error(const std::string_view what,
                         const int error,
                         const inflate_state* state,
                         const size_t isize)
{
  AR_REQUIRE(state);
  std::ostringstream os;
  os << "error " << what << ": ";

  if (error == ISAL_END_INPUT) {
    os << "BGZF block is truncated";
  } else if (error == ISAL_OUT_OVERFLOW) {
    os << "BGZF block has too small ISIZE; content is likely corrupt";
  } else if (error == ISAL_INCORRECT_CHECKSUM) {
    if (state->block_state != ISAL_BLOCK_FINISH) {
      os << "incorrect header checksum; BGZF block is corrupt";
    } else if (isize == state->total_out) {
      os << "incorrect checksum; content is corrupt";
    } else {
      os << "BGZF block has too large ISIZE; content is likely corrupt";
    }
  } else {
    // for most part, the exact error is not useful to the user
    os << "invalid BGZF block";
  }

  os << " (error = " << error << ")";
  throw gzip_error(os.str());
}

} // namespace

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

  dst.append_u32(calculate_crc32(src)); // CRC32
  dst.append_u32(src.size());           // ISIZE

  AR_REQUIRE(dst.size() <= BGZF_MAX_BLOCK_SIZE);

  // Write the final block size; -1 to fit 1..65536 in 16 bit
  dst.put_u16(16, dst.size() - 1);
}

void
bgzf_block_decompress(const buffer& src, buffer& dst, size_t dst_offset)
{
  AR_REQUIRE(dst_offset <= dst.size());
  if (src.size() < BGZF_EOF.size()) {
    throw gzip_error("malformed BGZF block; block is too small");
  }

  const auto block_size = bgzf_block_size(src);
  if (block_size == 0) {
    throw gzip_error("not a BGZF block");
  } else if (src.size() != block_size) {
    throw gzip_error("BGZF block has mismatching BSIZE");
  }

  const auto isize = src.get_u32(block_size - 4);
  if (isize > BGZF_MAX_BLOCK_SIZE) {
    throw gzip_error("BGZF block has invalid ISIZE = " + stringify(isize));
  }

  dst.resize(isize + dst_offset);

  // avoid zeroing giant struct on every call due to hardening flags
  thread_local auto state = std::make_unique<inflate_state>();

  isal_inflate_init(state.get());
  state->crc_flag = ISAL_GZIP;
  state->next_in = const_cast<uint8_t*>(src.data());
  state->avail_in = block_size;
  state->next_out = dst.data() + dst_offset;
  state->avail_out = dst.size() - dst_offset;

  auto result = isal_inflate_stateless(state.get());
  if (result != ISAL_DECOMP_OK) {
    throw_isal_inflate_error("decompressing BGZF block",
                             result,
                             state.get(),
                             isize);
  } else if (state->avail_in) {
    throw gzip_error("BGZF block contains " + stringify(state->avail_in) +
                     " bytes of trailing garbage");
  }
}

size_t
bgzf_block_size(const buffer& src)
{
  const auto* data = src.data();

  // Required values in BGZF:
  if (src.size() >= 18   //
      && data[0] == 0x1f // ID1
      && data[1] == 0x8b // ID2
      && data[2] == 0x8  // Compression method
      && data[3] & 0x4   // Flag FEXTRA must be set
      // MTIME (4 bytes), XFL (1 byte), and OS (1 byte) are ignored
      && src.get_u16(10) >= 6 // XLEN
      && data[12] == 0x42     // Subfield identifier 1
      && data[13] == 0x43     // Subfield identifier 2
      && data[14] == 2        // Subfield len [0]
      && data[15] == 0        // Subfield len [1]
  ) {
    // This will misidentify BGZF files with leading subfields, which is
    // permitted by the standard, but matches the behavior of the reference
    // `bgzip` implementation as of the time of writing:
    //   https://github.com/samtools/htslib/blob/7c5e3e7ebcdf90c8f96afd8a06d75ffa5603e417/bgzf.c#L414
    return static_cast<size_t>(src.get_u16(16)) + 1;
  }

  return 0;
}

} // namespace adapterremoval
