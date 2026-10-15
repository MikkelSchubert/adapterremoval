// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef>     // for size_t
#include <string_view> // for string_view

namespace adapterremoval {

class buffer;

//! Maximum size of uncompressed input to BGZF blocks; same as htslib
constexpr size_t BGZF_MAX_INPUT_SIZE = 0xff00;
//! Maximum allowed size of BGZF blocks including header and tail
constexpr size_t BGZF_MAX_BLOCK_SIZE = 0x10000;

//! Default bgzf header, as described in the SAM spec. v1.6 section 4.1.
//! Includes 2 trailing placeholder bytes for total block size (BSIZE)
constexpr std::string_view BGZF_HEADER = {
  "\37\213\10\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0",
  18
};
//! EOF marker for bgzf files; see SAM spec. v1.6 section 4.1.2
constexpr std::string_view BGZF_EOF = {
  "\37\213\10\4\0\0\0\0\0\377\6\0\102\103\2\0\33\0\3\0\0\0\0\0\0\0\0\0",
  28,
};

/**
 * Creates bgzf block using libdeflate, isa-l, or using no compression
 *
 * @param src Buffer containing data written as bgzf block
 * @param dst Write bgzf block to this buffer, replacing current content
 * @param compression_level Compression level in the range 0 to 12
 */
void
bgzf_block_compress(const buffer& src, buffer& dst, unsigned compression_level);

} // namespace adapterremoval
