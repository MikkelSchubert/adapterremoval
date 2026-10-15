// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "bgzf.hpp"     // for BGZF_MAX_BLOCK_SIZE, BGZF_MAX_INPUT_SIZE, ...
#include "buffer.hpp"   // for buffer
#include "deflate.hpp"  // for MAX_COMPRESSION_LEVEL
#include "errors.hpp"   // for assert_failed
#include "testing.hpp"  // for TEST_CASE, REQUIRE, CHECK, ...
#include <cstddef>      // for size_t
#include <cstdint>      // uint8_t
#include <libdeflate.h> // for libdeflate_crc32
#include <memory>       // for unique_ptr

#define PARAMETERIZE_COMPRESSION GENERATE(0, 1, 2, 3, 4, MAX_COMPRESSION_LEVEL)

using decompressor_ptr =
  std::unique_ptr<libdeflate_decompressor,
                  decltype(&libdeflate_free_decompressor)>;

namespace adapterremoval {

TEST_CASE("bgzf_block_compress gzip compresses with BSIZE, ISIZE, CRC32")
{
  buffer src;

  // Plain data or worst case (incompressible) data
  if (GENERATE(true, false)) {
    src = "AAAAAAAAABBBBBBBBBBBBCCCCCCCCCCCCDDDDDDDDDDEEEEEEEEEEE"_buffer;
  } else {
    src = random_data_buffer(BGZF_MAX_INPUT_SIZE);
  }

  buffer dst;
  bgzf_block_compress(src, dst, PARAMETERIZE_COMPRESSION);

  // Expected values per the SAM/BAM specification (12 Aug 2025)
  CHECK(dst.at(0) == 31);                   // ID1
  CHECK(dst.at(1) == 139);                  // ID2
  CHECK(dst.at(2) == 8);                    // CM
  CHECK(dst.at(3) == 4);                    // FLG
  CHECK(dst.get_u16(10) == 6);              // XLEN; assuming one extra subfield
  CHECK(dst.at(12) == 66);                  // SI1
  CHECK(dst.at(13) == 67);                  // SI2
  CHECK(dst.get_u16(14) == 2);              // SLEN
  CHECK(dst.get_u16(16) == dst.size() - 1); // BSIZE

  // CRC32
  CHECK(dst.get_u32(dst.size() - 8) ==
        libdeflate_crc32(0, src.data(), src.size()));
  // ISIZE
  CHECK(dst.get_u32(dst.size() - 4) == src.size());

  decompressor_ptr decompressor{ libdeflate_alloc_decompressor(),
                                 libdeflate_free_decompressor };
  REQUIRE(decompressor.get());

  buffer roundtrip{ BGZF_MAX_BLOCK_SIZE };

  size_t actual_out_nbytes_ret = 0;
  auto result = libdeflate_gzip_decompress(decompressor.get(),
                                           dst.data(),
                                           dst.size(),
                                           roundtrip.data(),
                                           roundtrip.size(),
                                           &actual_out_nbytes_ret);
  REQUIRE(result == LIBDEFLATE_SUCCESS);
  REQUIRE(actual_out_nbytes_ret == src.size());
  roundtrip.resize(actual_out_nbytes_ret);
  REQUIRE(src == roundtrip);
}

TEST_CASE("bgzf_block_compress rejects input larger than BGZF_MAX_INPUT_SIZE")
{
  auto src = random_data_buffer(BGZF_MAX_INPUT_SIZE + 1);
  buffer dst;

  REQUIRE_THROWS_AS(bgzf_block_compress(src, dst, PARAMETERIZE_COMPRESSION),
                    assert_failed);
}

TEST_CASE("bgzf_block_compress overwrites dst")
{
  auto src = random_data_buffer(1024);
  auto full_dst = random_data_buffer(1024);
  buffer empty_dst;

  const auto compression_level = PARAMETERIZE_COMPRESSION;
  REQUIRE_NOTHROW(bgzf_block_compress(src, full_dst, compression_level));
  REQUIRE_NOTHROW(bgzf_block_compress(src, empty_dst, compression_level));
  REQUIRE(full_dst == empty_dst);
}

TEST_CASE("bgzf_block_compress supports uncompressed blocks")
{
  auto src = "AAAAAAAAABBBBBBBBBBBBCCCCCCCCCCCCDDDDDDDDDDEEEEEEEEEEE"_buffer;
  buffer dst;

  bgzf_block_compress(src, dst, 0);
  REQUIRE(dst.at(18) == 1); // // BFINAL=1, BTYPE=00 (uncompressed)
  REQUIRE(dst.get_u16(19) == src.size());                         // LEN
  REQUIRE(dst.get_u16(21) == static_cast<uint16_t>(~src.size())); // NLEN

  buffer cdata;
  cdata.append(dst.data() + 23, dst.size() - 31);
  REQUIRE(src == cdata);
}

TEST_CASE("bgzf_block_compress rejects compression levels above max")
{
  auto src = "AAAAAAAAABBBBBBBBBBBBCCCCCCCCCCCCDDDDDDDDDDEEEEEEEEEEE"_buffer;
  buffer dst;

  REQUIRE_THROWS_AS(bgzf_block_compress(src, dst, MAX_COMPRESSION_LEVEL + 1),
                    assert_failed);
}

} // namespace adapterremoval
