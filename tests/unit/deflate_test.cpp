// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "bgzf.hpp"     // for BGZF_MAX_BLOCK_SIZE, BGZF_MAX_INPUT_SIZE
#include "buffer.hpp"   // for buffer
#include "deflate.hpp"  // for deflate_block_uncompressed, ...
#include "errors.hpp"   // for assert_failed
#include "testing.hpp"  // for TEST_CASE, REQUIRE, CHECK, ...
#include <libdeflate.h> // for libdeflate_crc32

#define PARAMETERIZE_COMPRESSION GENERATE(0, 1, 2, 3, 4, MAX_COMPRESSION_LEVEL)

namespace adapterremoval {

namespace {

using decompressor_ptr =
  std::unique_ptr<libdeflate_decompressor,
                  decltype(&libdeflate_free_decompressor)>;

/** Decompress a buffer containing a (final) deflate block starting at offset */
buffer
inflate_buffer(const buffer& src, size_t offset = 0)
{
  REQUIRE(offset <= src.size());
  decompressor_ptr decompressor{ libdeflate_alloc_decompressor(),
                                 libdeflate_free_decompressor };
  REQUIRE(decompressor.get());

  buffer dst{ 1024 };
  size_t actual_out_nbytes_ret = 0;
  const auto result = libdeflate_deflate_decompress(decompressor.get(),
                                                    src.data() + offset,
                                                    src.size() - offset,
                                                    dst.data(),
                                                    dst.size(),
                                                    &actual_out_nbytes_ret);
  REQUIRE(result == 0);
  dst.resize(actual_out_nbytes_ret);
  return dst;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// isal stream head and tail

TEST_CASE("gzip_stream_header_isal minimal checks", "[deflate]")
{
  const auto buf = gzip_stream_header_isal();
  CHECK(buf.size() == 10);
  CHECK(buf.at(0) == 0x1f);
  CHECK(buf.at(1) == 0x8b);
}

TEST_CASE("gzip_stream_tail_isal minimal checks", "[deflate]")
{
  const auto buf = gzip_stream_tail_isal(0xDEADBEEF, 1234567890);
  CHECK(buf.size() == 8);
  CHECK(buf.get_u32(0) == 0xDEADBEEF);
  CHECK(buf.get_u32(4) == 1234567890);
}

////////////////////////////////////////////////////////////////////////////////
// deflate_block_isal_1

TEST_CASE("deflate_block_isal_1 roundtrip")
{
  buffer src;
  src.append(GENERATE("", "Holy hole in a donut!"));

  SECTION("empty destination")
  {
    buffer comp{ 1024 };
    deflate_block_isal_1(src, comp, 0, true);
    buffer dst = inflate_buffer(comp, 0);
    REQUIRE(src == dst);
  }
  SECTION("non-empty destination")
  {
    buffer comp{ 1024 };
    comp.put_u32(0, 0xDEADBEEF);
    deflate_block_isal_1(src, comp, 8, true);
    REQUIRE(comp.get_u32(0) == 0xDEADBEEF);
    buffer dst = inflate_buffer(comp, 8);
    REQUIRE(src == dst);
  }
}

TEST_CASE("deflate_block_isal_1 must fit in output buffer")
{
  const auto src = "RANDOM BULLSHIT GO!!!!"_buffer;
  const bool eof = GENERATE(true, false);

  SECTION("OK")
  {
    buffer dst{ 32 };
    REQUIRE_NOTHROW(deflate_block_isal_1(src, dst, 0, eof));
  }

  SECTION("OK; extra size for offset")
  {
    buffer dst{ 32 };
    dst.put_u16(0, 0xBEEF);
    REQUIRE_NOTHROW(deflate_block_isal_1(src, dst, 2, eof));
    REQUIRE(dst.get_u16(0) == 0xBEEF);
  }

  SECTION("Not OK")
  {
    buffer dst{ 32 };
    REQUIRE_THROWS_AS(deflate_block_isal_1(src, dst, 16, eof), assert_failed);
  }
}

TEST_CASE("deflate_block_isal_1 offset must be inside output buffer")
{
  const auto src = "But for me, it was Tuesday"_buffer;
  const auto eof = GENERATE(true, false);

  buffer dst{ GENERATE(0, 32) };
  const auto offset = dst.size() + GENERATE(0, 1);
  REQUIRE_THROWS_AS(deflate_block_isal_1(src, dst, offset, eof), assert_failed);
}

TEST_CASE("deflate_block_isal_1 worst case for stream compression")
{
  auto src = random_data_buffer(BGZF_MAX_INPUT_SIZE);
  buffer dst{ BGZF_MAX_BLOCK_SIZE };
  const bool eof = GENERATE(true, false);
  REQUIRE_NOTHROW(deflate_block_isal_1(src, dst, 0, eof));
}

////////////////////////////////////////////////////////////////////////////////
// isal gzip stream

TEST_CASE("isal gzip stream with head and tail roundtrip")
{
  const auto src = "Excuse me, but can I distract you for a second?"_buffer;
  auto comp = gzip_stream_header_isal();
  size_t offset = comp.size();
  comp.resize(1024);
  deflate_block_isal_1("Excuse me, but can I "_buffer, comp, offset, false);
  offset = comp.size();
  comp.resize(1024);
  deflate_block_isal_1("distract you for a second?"_buffer, comp, offset, true);
  comp.append(gzip_stream_tail_isal(calculate_crc32(src), src.size()));

  decompressor_ptr decompressor{ libdeflate_alloc_decompressor(),
                                 libdeflate_free_decompressor };
  REQUIRE(decompressor.get());

  buffer dst{ 1024 };
  size_t actual_out_nbytes_ret = 0;
  const auto result = libdeflate_gzip_decompress(decompressor.get(),
                                                 comp.data(),
                                                 comp.size(),
                                                 dst.data(),
                                                 dst.size(),
                                                 &actual_out_nbytes_ret);
  REQUIRE(result == 0);
  dst.resize(actual_out_nbytes_ret);
  REQUIRE(src == dst);
}

////////////////////////////////////////////////////////////////////////////////
// deflate_block_libdeflate

TEST_CASE("deflate_block_libdeflate roundtrip")
{
  buffer src;
  src.append(GENERATE("", "Holy hole in a donut!"));
  const auto compression_level = PARAMETERIZE_COMPRESSION;

  SECTION("empty destination")
  {
    buffer comp{ 1024 };
    deflate_block_libdeflate(src, comp, 0, compression_level);
    buffer dst = inflate_buffer(comp, 0);
    REQUIRE(src == dst);
  }

  SECTION("non-empty destination")
  {
    buffer comp{ 1024 };
    comp.put_u32(0, 0xDEADBEEF);
    deflate_block_libdeflate(src, comp, 8, compression_level);
    REQUIRE(comp.get_u32(0) == 0xDEADBEEF);
    buffer dst = inflate_buffer(comp, 8);
    REQUIRE(src == dst);
  }
}

TEST_CASE("deflate_block_libdeflate must fit in output buffer")
{
  const auto src = "RANDOM BULLSHIT GO!!!!"_buffer;
  const auto compression_level = PARAMETERIZE_COMPRESSION;

  SECTION("OK")
  {
    buffer dst{ 32 };
    REQUIRE_NOTHROW(deflate_block_libdeflate(src, dst, 0, compression_level));
  }

  SECTION("OK; extra size for offset")
  {
    buffer dst{ 32 };
    dst.put_u16(0, 0xBEEF);
    REQUIRE_NOTHROW(deflate_block_libdeflate(src, dst, 2, compression_level));
    REQUIRE(dst.get_u16(0) == 0xBEEF);
  }

  SECTION("Not OK")
  {
    buffer dst{ 32 };
    REQUIRE_THROWS_AS(deflate_block_libdeflate(src, dst, 16, compression_level),
                      assert_failed);
  }
}

TEST_CASE("deflate_block_libdeflate offset must be inside output buffer")
{
  const auto src = "But for me, it was Tuesday"_buffer;
  const auto compression_level = PARAMETERIZE_COMPRESSION;

  buffer dst{ GENERATE(0, 32) };
  const auto offset = dst.size() + GENERATE(0, 1);
  REQUIRE_THROWS_AS(
    deflate_block_libdeflate(src, dst, offset, compression_level),
    assert_failed);
}

////////////////////////////////////////////////////////////////////////////////
// deflate_block_uncompressed

TEST_CASE("deflate_block_uncompressed roundtrip")
{
  buffer src;
  src.append(GENERATE("", "Sweet great googly moogly!"));

  SECTION("empty destination")
  {
    buffer comp{ 1024 };
    deflate_block_uncompressed(src, comp, 0);
    buffer dst = inflate_buffer(comp, 0);
    REQUIRE(src == dst);
  }

  SECTION("non-empty destination")
  {
    buffer comp{ 1024 };
    comp.put_u32(0, 0xDEADBEEF);
    deflate_block_uncompressed(src, comp, 8);
    REQUIRE(comp.get_u32(0) == 0xDEADBEEF);
    buffer dst = inflate_buffer(comp, 8);
    REQUIRE(src == dst);
  }
}

TEST_CASE("deflate_block_uncompressed must fit in output buffer")
{
  const auto src = "RANDOM BULLSHIT GO!!!!"_buffer;

  SECTION("OK")
  {
    buffer dst{ src.size() + 5 }; // flag + 2 * LEN
    REQUIRE_NOTHROW(deflate_block_uncompressed(src, dst, 0));
  }

  SECTION("OK; extra size for offset")
  {
    buffer dst{ src.size() + 7 }; // 2 + flag + 2 * LEN
    dst.put_u16(0, 0xBEEF);
    REQUIRE_NOTHROW(deflate_block_uncompressed(src, dst, 2));
    REQUIRE(dst.get_u16(0) == 0xBEEF);
  }

  SECTION("Not OK")
  {
    buffer dst{ src.size() + 5 }; // flag + 2 * LEN
    REQUIRE_THROWS_AS(deflate_block_uncompressed(src, dst, GENERATE(1, 2)),
                      assert_failed);
  }
}

TEST_CASE("deflate_block_uncompressed offset must be inside output buffer")
{
  const auto src = "But for me, it was Tuesday"_buffer;

  buffer dst{ GENERATE(0, 32) };
  const auto offset = dst.size() + GENERATE(0, 1);
  REQUIRE_THROWS_AS(deflate_block_uncompressed(src, dst, offset),
                    assert_failed);
}

TEST_CASE("deflate_block_uncompressed size is max 16 bit")
{
  buffer dst{ 0x20000 };

  SECTION("OK")
  {
    auto src = random_data_buffer(0xFFFF);
    REQUIRE_NOTHROW(deflate_block_uncompressed(src, dst, 0));
  }

  SECTION("too big")
  {
    auto src = random_data_buffer(0xFFFF + 1);
    REQUIRE_THROWS_AS(deflate_block_uncompressed(src, dst, 0), assert_failed);
  }
}

////////////////////////////////////////////////////////////////////////////////
// calculate_crc32

TEST_CASE("calculate_crc32 on empty buffer")
{
  REQUIRE(calculate_crc32(buffer{}, 0) == 0);
  REQUIRE(calculate_crc32(buffer{}, 0xDEADBEEF) == 0xDEADBEEF);
}

TEST_CASE("calculate_crc32 across buffers")
{
  uint32_t crc32 = 0;
  crc32 = calculate_crc32("Lorem ipsum dolor sit amet, "_buffer, crc32);
  crc32 = calculate_crc32(buffer{}, crc32);
  crc32 = calculate_crc32("consectetur  adipiscing elit"_buffer, crc32);

  REQUIRE(crc32 == calculate_crc32("Lorem ipsum dolor sit amet, "
                                   "consectetur  adipiscing elit"_buffer));
}

} // namespace adapterremoval
