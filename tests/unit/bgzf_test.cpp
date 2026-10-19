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
#include <limits>       // numeric_limits
#include <memory>       // for unique_ptr

#define PARAMETERIZE_COMPRESSION GENERATE(0, 1, 2, 3, 4, MAX_COMPRESSION_LEVEL)

using decompressor_ptr =
  std::unique_ptr<libdeflate_decompressor,
                  decltype(&libdeflate_free_decompressor)>;

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// bgzf_block_compress

TEST_CASE("bgzf_block_compress gzip compresses with BSIZE, ISIZE, CRC32")
{
  // Run the below for each compression level
  const auto compression_level = PARAMETERIZE_COMPRESSION;

  buffer src;

  SECTION("empty buffer")
  {
    REQUIRE(src.size() == 0); // silence `--warn NoAssertions`
  }

  SECTION("small compressible buffer")
  {
    src = random_text_buffer(1024);
    REQUIRE(src.size() == 1024); // silence `--warn NoAssertions`
  }

  SECTION("large compressible buffer")
  {
    src = random_text_buffer(BGZF_MAX_INPUT_SIZE);
    REQUIRE(src.size() == BGZF_MAX_INPUT_SIZE); // silence `--warn NoAssertions`
  }

  SECTION("large incompressible buffer")
  {
    src = random_data_buffer(BGZF_MAX_INPUT_SIZE);
    REQUIRE(src.size() == BGZF_MAX_INPUT_SIZE); // silence `--warn NoAssertions`
  }

  buffer dst;
  bgzf_block_compress(src, dst, compression_level);

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
  auto src = random_text_buffer(BGZF_MAX_INPUT_SIZE + 1);
  buffer dst;

  REQUIRE_THROWS_AS(bgzf_block_compress(src, dst, PARAMETERIZE_COMPRESSION),
                    assert_failed);
}

TEST_CASE("bgzf_block_compress overwrites dst")
{
  auto src = random_text_buffer(1024);
  auto full_dst = random_text_buffer(1024);
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

////////////////////////////////////////////////////////////////////////////////
// bgzf_block_size

TEST_CASE("bgzf_block_size accepts template header")
{
  buffer buf;
  buf.append(BGZF_HEADER);
  REQUIRE(bgzf_block_size(buf) == 1); // minimum BGZF block size
}

TEST_CASE("bgzf_block_size handles max block sizes")
{
  buffer buf;
  buf.append(BGZF_HEADER);
  buf.put_u16(16, 65535);

  REQUIRE(bgzf_block_size(buf) == 65536);
}

TEST_CASE("bgzf_block_size fails on short buffer")
{
  buffer buf;
  buf.append(BGZF_HEADER);
  REQUIRE(buf.size() == 18);
  buf.resize(GENERATE(range(0, 18)));

  REQUIRE(bgzf_block_size(buf) == 0);
}

TEST_CASE("bgzf_block_size fails on non-gzip")
{
  buffer buf;
  buf.append(BGZF_HEADER);
  REQUIRE(bgzf_block_size(buf) == 1);
  const auto offset = GENERATE(0, 1, 2);
  buf.put_u8(offset, buf.at(offset) ^ 1); // ID1, ID2, CM
  REQUIRE(bgzf_block_size(buf) == 0);
}

TEST_CASE("bgzf_block_size fails on missing extra flag")
{
  buffer src;
  src.append(BGZF_HEADER);
  REQUIRE(bgzf_block_size(src) == 1);
  REQUIRE(src.at(3) & 4); // FEXTRA
  src.put_u8(3, src.at(3) ^ 4);
  REQUIRE(bgzf_block_size(src) == 0);
}

TEST_CASE("bgzf_block_size fails on no/short extra fields")
{
  buffer src;
  src.append(BGZF_HEADER);
  REQUIRE(bgzf_block_size(src) == 1);
  REQUIRE(src.get_u16(10) == 6); // XLEN
  src.put_u16(10, GENERATE(range(0, 6)));
  REQUIRE(bgzf_block_size(src) == 0);
}

TEST_CASE("bgzf_block_size fails on wrong extra field values")
{
  buffer src;
  src.append(BGZF_HEADER);
  REQUIRE(bgzf_block_size(src) == 1);
  REQUIRE(src.get_u16(10) == 6); // XLEN

  SECTION("wrong SI1")
  {
    src.put_u8(12, src.at(12) + GENERATE(-1, 1));
    REQUIRE(true); // silence `--warn NoAssertions`
  }

  SECTION("wrong SI2")
  {
    src.put_u8(13, src.at(13) + GENERATE(-1, 1));
    REQUIRE(true); // silence `--warn NoAssertions`
  }

  SECTION("wrong subfield length 1")
  {
    src.put_u8(14, src.at(14) + GENERATE(-1, 1));
    REQUIRE(true); // silence `--warn NoAssertions`
  }

  SECTION("wrong subfield length 2")
  {
    src.put_u8(15, src.at(15) + GENERATE(-1, 1));
    REQUIRE(true); // silence `--warn NoAssertions`
  }

  REQUIRE(bgzf_block_size(src) == 0);
}

TEST_CASE("bgzf_block_size allows extra flags")
{
  buffer src;
  src.append(BGZF_HEADER);
  REQUIRE(bgzf_block_size(src) == 1);
  src.put_u8(3, 0xF);
  REQUIRE(bgzf_block_size(src) == 1);
}

TEST_CASE("bgzf_block_size ignores MTIME")
{
  buffer buf;
  buf.append(BGZF_HEADER);
  buf.put_u32(4, GENERATE(0, std::numeric_limits<uint32_t>::max()));
  REQUIRE(bgzf_block_size(buf) == 1);
}

TEST_CASE("bgzf_block_size ignores XFL")
{
  buffer buf;
  buf.append(BGZF_HEADER);
  buf.put_u8(8, GENERATE(0, 2, 4));
  REQUIRE(bgzf_block_size(buf) == 1);
}

TEST_CASE("bgzf_block_size ignores OS")
{
  buffer buf;
  buf.append(BGZF_HEADER);
  buf.put_u8(9, GENERATE(range(0, 14), 255));
  REQUIRE(bgzf_block_size(buf) == 1);
}

TEST_CASE("bgzf_block_size returns size of bgzf_block_compress'd block")
{
  buffer src;
  src.append(random_text_buffer(1024));
  buffer dst;
  bgzf_block_compress(src, dst, PARAMETERIZE_COMPRESSION);
  REQUIRE(bgzf_block_size(dst) == dst.size());
}

////////////////////////////////////////////////////////////////////////////////
// bgzf_block_decompress

TEST_CASE("bgzf_block_decompress fails on too short buffer")
{
  buffer buf;
  buf.append(BGZF_EOF); // minimal block size
  REQUIRE(buf.size() == 28);
  const auto size = GENERATE(range(0, 28));
  buf.put_u16(16, size - 1);
  buf.resize(size);

  buffer dst;
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(buf, dst),
                          gzip_error,
                          "malformed BGZF block; block is too small");
}

TEST_CASE("bgzf_block_decompress fails on non-bgzf buffer")
{
  buffer buf;
  buf.append(BGZF_EOF); // minimal block size
  buf.put_u8(12, '?');  // SI1

  buffer dst;
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(buf, dst),
                          gzip_error,
                          "not a BGZF block");
}

TEST_CASE("bgzf_block_decompress fails on BSIZE mismatch")
{
  buffer buf;
  buf.append(BGZF_EOF); // minimal block size
  buf.put_u16(16, buf.size() + GENERATE(-2, 0));

  buffer dst;
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(buf, dst),
                          gzip_error,
                          "BGZF block has mismatching BSIZE");
}

TEST_CASE("bgzf_block_decompress fails on CRC32 mismatch")
{
  const auto data = random_text_buffer(GENERATE(0, 1024));
  const auto c32 = calculate_crc32(data);
  buffer compressed;
  bgzf_block_compress(data, compressed, PARAMETERIZE_COMPRESSION);
  REQUIRE(compressed.get_u32(compressed.size() - 8) == c32);
  buffer uncompressed;
  REQUIRE_NOTHROW(bgzf_block_decompress(compressed, uncompressed));
  compressed.put_u32(compressed.size() - 8, c32 + 1);
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(compressed, uncompressed),
                          gzip_error,
                          "incorrect checksum; content is corrupt");
}

TEST_CASE("bgzf_block_decompress fails on HCRC mismatch")
{
  const auto data = random_text_buffer(1024);

  buffer bgzf;
  bgzf.append(BGZF_HEADER);
  bgzf.append_u16(0); // Dummy value
  bgzf.resize(BGZF_MAX_BLOCK_SIZE);
  deflate_block_libdeflate(data,
                           bgzf,
                           BGZF_HEADER.size() + 2,
                           PARAMETERIZE_COMPRESSION);
  bgzf.append_u32(calculate_crc32(data));
  bgzf.append_u32(data.size());
  bgzf.put_u16(16, bgzf.size() - 1);
  bgzf.put_u8(3, bgzf.at(3) | 2); // FHCRC

  buffer header;
  header.append(bgzf.data(), BGZF_HEADER.size());
  const auto hcrc = calculate_crc32(header);
  bgzf.put_u16(BGZF_HEADER.size(), hcrc & 0xFFFF);

  buffer dst;

  SECTION("matching HCRC")
  {
    bgzf_block_decompress(bgzf, dst);
    REQUIRE(dst == data);
  }

  SECTION("mismatching HCRC")
  {
    bgzf.put_u8(BGZF_HEADER.size(), bgzf.at(BGZF_HEADER.size()) ^ 1);
    REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(bgzf, dst),
                            gzip_error,
                            "incorrect header checksum");
  }
}

TEST_CASE("bgzf_block_decompress fails on too small ISIZE")
{
  const auto data = random_text_buffer(GENERATE(1, 1024));
  buffer bgzf;
  bgzf_block_compress(data, bgzf, PARAMETERIZE_COMPRESSION);
  REQUIRE(bgzf.get_u32(bgzf.size() - 4) == data.size());
  buffer uncompressed;
  REQUIRE_NOTHROW(bgzf_block_decompress(bgzf, uncompressed));

  bgzf.put_u32(bgzf.size() - 4, data.size() - 1);
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(bgzf, uncompressed),
                          gzip_error,
                          "BGZF block has too small ISIZE");
}

TEST_CASE("bgzf_block_decompress fails on too large ISIZE mismatch")
{
  const auto data = random_text_buffer(GENERATE(0, 1024));
  buffer bgzf;
  bgzf_block_compress(data, bgzf, PARAMETERIZE_COMPRESSION);
  REQUIRE(bgzf.get_u32(bgzf.size() - 4) == data.size());
  buffer uncompressed;
  REQUIRE_NOTHROW(bgzf_block_decompress(bgzf, uncompressed));

  bgzf.put_u32(bgzf.size() - 4, data.size() + 1);
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(bgzf, uncompressed),
                          gzip_error,
                          "BGZF block has too large ISIZE");
}

TEST_CASE("bgzf_block_decompress fails on end of input")
{
  const auto data = random_text_buffer(1024);
  buffer bgzf;
  bgzf_block_compress(data, bgzf, PARAMETERIZE_COMPRESSION);

  // drop crc32
  bgzf.put_u32(bgzf.size() - 8, bgzf.get_u32(bgzf.size() - 4));
  bgzf.resize(bgzf.size() - 4);
  bgzf.put_u16(16, bgzf.size() - 1);

  buffer uncompressed;
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(bgzf, uncompressed),
                          gzip_error,
                          "BGZF block is truncated");
}

TEST_CASE("bgzf_block_decompress fails on too large isize")
{
  const auto data = random_text_buffer(1024);
  buffer bgzf;
  bgzf_block_compress(data, bgzf, PARAMETERIZE_COMPRESSION);
  REQUIRE(bgzf.get_u32(bgzf.size() - 4) == data.size());
  buffer uncompressed;
  REQUIRE_NOTHROW(bgzf_block_decompress(bgzf, uncompressed));
  bgzf.put_u32(bgzf.size() - 4, BGZF_MAX_BLOCK_SIZE + 1);
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(bgzf, uncompressed),
                          gzip_error,
                          "BGZF block has invalid ISIZE");
}

TEST_CASE("bgzf_block_decompress fails on non-final")
{
  const bool final = GENERATE(true, false);
  const auto data = random_text_buffer(1024);
  const auto crc32 = calculate_crc32(data);

  buffer bgzf;
  bgzf.append(BGZF_HEADER);
  bgzf.resize(BGZF_MAX_INPUT_SIZE);
  deflate_block_isal_1(data, bgzf, BGZF_HEADER.size(), final);
  bgzf.append_u32(crc32);
  bgzf.append_u32(data.size());
  bgzf.put_u16(16, bgzf.size() - 1);

  buffer dst;
  if (final) {
    REQUIRE_NOTHROW(bgzf_block_decompress(bgzf, dst));
    REQUIRE(dst == data);
  } else {
    // The exact error depends on the (compressed) content
    REQUIRE_THROWS_AS(bgzf_block_decompress(bgzf, dst), gzip_error);
  }
}

TEST_CASE("bgzf_block_decompress fails on corrupt data")
{
  const auto data = random_text_buffer(1024);
  buffer bgzf;
  bgzf_block_compress(data, bgzf, PARAMETERIZE_COMPRESSION);

  // flip a bit
  AR_REQUIRE(bgzf.size() >= 300);
  bgzf.put_u8(200, bgzf.at(200) ^ 1);

  buffer dst;
  // The exact error depends on compression level/method
  REQUIRE_THROWS_AS(bgzf_block_decompress(bgzf, dst), gzip_error);
}

TEST_CASE("bgzf_block_decompress fails on corrupt block type")
{
  const auto data = random_text_buffer(1024);
  buffer bgzf;
  bgzf_block_compress(data, bgzf, PARAMETERIZE_COMPRESSION);

  // nonsense header
  bgzf.put_u8(19, 0xFF);

  buffer dst;
  // The exact error depends on compression level/method
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(bgzf, dst),
                          gzip_error,
                          "invalid BGZF block");
}

TEST_CASE("bgzf_block_decompress fails on trailing data")
{
  buffer compressed;
  bgzf_block_compress(random_text_buffer(1024),
                      compressed,
                      PARAMETERIZE_COMPRESSION);

  // this probably cannot happen outside of specifically crafted input
  const auto tail_start = compressed.size() - 8;
  compressed.append_u32(0xDEADBEEF);
  compressed.append_u32(compressed.get_u32(tail_start));
  compressed.append_u32(compressed.get_u32(tail_start + 4));
  compressed.put_u16(16, compressed.size() - 1);

  buffer uncompressed;
  REQUIRE_THROWS_CONTAINS(bgzf_block_decompress(compressed, uncompressed),
                          gzip_error,
                          "BGZF block contains 12 bytes of trailing garbage");
}

TEST_CASE("bgzf_block_decompress may preserve existing data")
{
  const auto src = random_text_buffer(GENERATE(0, 1024, BGZF_MAX_INPUT_SIZE));

  buffer bgzf;
  bgzf_block_compress(src, bgzf, PARAMETERIZE_COMPRESSION);

  buffer dst;

  SECTION("appending to end")
  {
    dst.append_u32(0xDEADBEEF);
    REQUIRE(dst.size() == 4);
  }

  SECTION("overwriting tail")
  {
    dst.append_u32(0xDEADBEEF);
    dst.append_u32(0xBEEFDEAD);
    REQUIRE(dst.size() == 8);
  }

  bgzf_block_decompress(bgzf, dst, 4);

  buffer expected;
  expected.append_u32(0xDEADBEEF);
  expected.append(src);

  REQUIRE(dst == expected);
}

TEST_CASE("bgzf_block_decompress overwrites if offset not set")
{
  const auto src = random_text_buffer(GENERATE(0, 1024, BGZF_MAX_INPUT_SIZE));

  buffer bgzf;
  bgzf_block_compress(src, bgzf, PARAMETERIZE_COMPRESSION);

  buffer dst;
  dst.append_u32(0xDEADBEEF);
  bgzf_block_decompress(bgzf, dst);

  REQUIRE(dst == src);
}

TEST_CASE("bgzf_block_decompress requires offset inside buffer")
{
  buffer bgzf;
  bgzf_block_compress(random_text_buffer(1024), bgzf, PARAMETERIZE_COMPRESSION);

  buffer dst;
  dst.reserve(GENERATE(0, 1024)); // capacity must not matter

  SECTION("past EOF for empty buffer")
  {
    REQUIRE_THROWS_AS(bgzf_block_decompress(bgzf, dst, 1), assert_failed);
  }

  SECTION("past EOF for non-empty buffer")
  {
    dst.append_u32(0xDEADBEEF);
    REQUIRE_THROWS_AS(bgzf_block_decompress(bgzf, dst, 5), assert_failed);
  }
}

////////////////////////////////////////////////////////////////////////////////
// bgzf_block_decompress for extra fields

TEST_CASE("bgzf_block_decompress accepts extra flags")
{
  buffer bgzf;
  bgzf.append(BGZF_HEADER);

  if (GENERATE(true, false)) {
    bgzf.put_u8(3, bgzf.at(3) | 1); // FTEXT
  }

  if (GENERATE(true, false)) {
    // Trailing FEXTRA fields
    bgzf.append_u8(':');
    bgzf.append_u8(')');
    bgzf.append_u16(4);
    bgzf.append_u32(0xDEADBEEF);
    bgzf.put_u16(10, bgzf.get_u16(10) + 8);
  }

  if (GENERATE(true, false)) {
    bgzf.put_u8(3, bgzf.at(3) | 8); // FNAME
    bgzf.append("some filename.ext");
    bgzf.append_u8(0);
  }

  if (GENERATE(true, false)) {
    bgzf.put_u8(3, bgzf.at(3) | 16); // FCOMMENT
    bgzf.append("some comment");
    bgzf.append_u8(0);
  }

  size_t hcrc_offset = 0;
  if (GENERATE(true, false)) {
    bgzf.put_u8(3, bgzf.at(3) | 2); // FHCRC
    hcrc_offset = bgzf.size();
    bgzf.append_u16(0); // Dummy value
  }

  const auto data = random_text_buffer(1024);
  const auto crc32 = calculate_crc32(data);
  const auto offset = bgzf.size();

  bgzf.resize(BGZF_MAX_BLOCK_SIZE);
  deflate_block_libdeflate(data, bgzf, offset, PARAMETERIZE_COMPRESSION);
  bgzf.append_u32(crc32);
  bgzf.append_u32(data.size());
  bgzf.put_u16(16, bgzf.size() - 1);

  if (hcrc_offset) {
    buffer header;
    header.append(bgzf.data(), hcrc_offset);
    const auto hcrc32 = calculate_crc32(header);
    bgzf.put_u16(hcrc_offset, hcrc32 & 0xFFFF); // CRC16
  }

  buffer dst;
  bgzf_block_decompress(bgzf, dst);
  REQUIRE(dst == data);
}

////////////////////////////////////////////////////////////////////////////////
// roundtrips

TEST_CASE("bgzf roundtrip for typical data")
{
  // small and typical block sizes
  const auto src = random_text_buffer(GENERATE(1024, BGZF_MAX_INPUT_SIZE));

  buffer bgzf;
  bgzf_block_compress(src, bgzf, PARAMETERIZE_COMPRESSION);
  buffer dst;
  bgzf_block_decompress(bgzf, dst, 0);

  REQUIRE(dst == src);
}

TEST_CASE("bgzf roundtrip for incompressible data")
{
  // small and worst-case (incompressible) blocks
  const auto src = random_data_buffer(GENERATE(1024, BGZF_MAX_INPUT_SIZE));

  buffer bgzf;
  bgzf_block_compress(src, bgzf, PARAMETERIZE_COMPRESSION);
  buffer dst;
  bgzf_block_decompress(bgzf, dst, 0);

  REQUIRE(dst == src);
}

TEST_CASE("bgzf roundtrip for empty buffer")
{
  buffer bgzf;
  bgzf_block_compress(buffer{}, bgzf, PARAMETERIZE_COMPRESSION);
  buffer dst;
  bgzf_block_decompress(bgzf, dst, 0);
  REQUIRE(dst.size() == 0);
}

TEST_CASE("bgzf roundtrip for EOF block")
{
  buffer bgzf;
  bgzf.append(BGZF_EOF);
  buffer dst;
  bgzf_block_decompress(bgzf, dst, 0);
  REQUIRE(dst.size() == 0);
}

TEST_CASE("bgzf roundtrip for maximum possible input size")
{
  // bgzip limits input size to 0xFF00, less than the format allows, to ensure
  // that there is room for it, but other bgzf writers may be less conservative
  const auto data = random_text_buffer(BGZF_MAX_BLOCK_SIZE);
  const auto crc32 = calculate_crc32(data);

  buffer bgzf;
  bgzf.append(BGZF_HEADER);
  bgzf.resize(BGZF_MAX_BLOCK_SIZE);
  deflate_block_isal_1(data, bgzf, BGZF_HEADER.size(), true);
  bgzf.append_u32(crc32);
  bgzf.append_u32(data.size());
  bgzf.put_u16(16, bgzf.size() - 1);

  buffer dst;
  bgzf_block_decompress(bgzf, dst);

  REQUIRE(dst == data);
}

TEST_CASE("bgzf roundtrip for bgzip output")
{
  // via `bgzip --binary | xxd -i`
  constexpr std::array<uint8_t, 114> data = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
    0x42, 0x43, 0x02, 0x00, 0x71, 0x00, 0x25, 0x88, 0x41, 0x0a, 0x80, 0x30,
    0x10, 0x03, 0xef, 0xfd, 0x8a, 0xa0, 0x56, 0x10, 0x11, 0x3d, 0x74, 0x4d,
    0x6b, 0xf6, 0x01, 0x7b, 0x17, 0x41, 0x3f, 0xe0, 0xff, 0x0f, 0x6e, 0x75,
    0x02, 0x21, 0x99, 0xf4, 0xdc, 0xe7, 0x75, 0xc4, 0xb1, 0x3f, 0x62, 0x17,
    0x03, 0xcc, 0x48, 0x9a, 0xc0, 0x20, 0x84, 0x3f, 0xf9, 0x2f, 0x04, 0xe2,
    0x5b, 0xe8, 0xc6, 0xaa, 0xfe, 0x52, 0x2b, 0x34, 0x41, 0x55, 0x77, 0x67,
    0xdb, 0x32, 0x2b, 0x5a, 0x21, 0x77, 0x64, 0x65, 0x29, 0x6b, 0x5a, 0x9a,
    0x84, 0x69, 0x1e, 0xa4, 0x45, 0x29, 0xae, 0x9d, 0x1c, 0x5e, 0xfa, 0xe2,
    0x3d, 0x87, 0x76, 0x00, 0x00, 0x00,
  };

  buffer bgzf;
  bgzf.append(data.data(), data.size());
  buffer dst;
  bgzf_block_decompress(bgzf, dst, 0);
  REQUIRE(dst == "@read_150_1/1\n"
                 "CTTGGGTACTCAGCCTTAGGGTACCACATAGAGTTATCTTATATATTATA\n"
                 "+\n"
                 "HHHFFFFBBDGGGGGHHHHHGGFCDHGEE<@;+@C792A.CEEGGFFFFD\n"_buffer);
}

} // namespace adapterremoval
