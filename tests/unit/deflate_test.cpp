// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "bgzf.hpp"      // for BGZF_MAX_BLOCK_SIZE, BGZF_MAX_INPUT_SIZE
#include "buffer.hpp"    // for buffer
#include "deflate.hpp"   // for deflate_block_uncompressed, ...
#include "errors.hpp"    // for assert_failed
#include "testing.hpp"   // for TEST_CASE, REQUIRE, CHECK, ...
#include "testutils.hpp" // for random_data_buffer
#include <isa-l/igzip_lib.h>
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
    REQUIRE_THROWS_CONTAINS(deflate_block_isal_1(src, dst, 16, eof),
                            assert_failed,
                            "isal_deflate_stateless returned -1");
  }
}

TEST_CASE("deflate_block_isal_1 offset must be inside output buffer")
{
  const auto src = "But for me, it was Tuesday"_buffer;
  const auto eof = GENERATE(true, false);

  buffer dst{ GENERATE(0U, 32U) };
  const auto offset = dst.size() + GENERATE(0, 1);
  REQUIRE_THROWS_CONTAINS(deflate_block_isal_1(src, dst, offset, eof),
                          assert_failed,
                          "~b4c7eead893019a6~");
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
    REQUIRE_THROWS_CONTAINS(
      deflate_block_libdeflate(src, dst, 16, compression_level),
      assert_failed,
      "~90ec9180db72efa4~");
  }
}

TEST_CASE("deflate_block_libdeflate offset must be inside output buffer")
{
  const auto src = "But for me, it was Tuesday"_buffer;
  const auto compression_level = PARAMETERIZE_COMPRESSION;

  buffer dst{ GENERATE(0U, 32U) };
  const auto offset = dst.size() + GENERATE(0, 1);
  REQUIRE_THROWS_CONTAINS(
    deflate_block_libdeflate(src, dst, offset, compression_level),
    assert_failed,
    "~14bf91174bef80b5~");
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
    REQUIRE_THROWS_CONTAINS(
      deflate_block_uncompressed(src, dst, GENERATE(1, 2)),
      assert_failed,
      "~be70eb4c5e2663e0~");
  }
}

TEST_CASE("deflate_block_uncompressed offset must be inside output buffer")
{
  const auto src = "But for me, it was Tuesday"_buffer;

  buffer dst{ GENERATE(0U, 32U) };
  const auto offset = dst.size() + GENERATE(0, 1);
  REQUIRE_THROWS_CONTAINS(deflate_block_uncompressed(src, dst, offset),
                          assert_failed,
                          "~be70eb4c5e2663e0~");
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
    REQUIRE_THROWS_CONTAINS(deflate_block_uncompressed(src, dst, 0),
                            assert_failed,
                            "~8c4967be494147cb~");
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

////////////////////////////////////////////////////////////////////////////////
// checked_isal_inflate / checked_isal_inflate_stateless

namespace {

struct isal_inflator
{
  isal_inflator() = default;

  isal_inflator(const buffer& buf)
    : src(gzip_data_into_buffer(buf))
    , dst(buf.size())
    , isize(buf.size())
  {
  }

  enum method
  {
    // decompression using `checked_isal_inflate`
    stateful,
    // decompression using `checked_isal_inflate_stateless`
    stateless,
  };

  void operator()(method m, std::string_view filename = {})
  {
    inflate_state state;
    isal_inflate_init(&state);
    state.crc_flag = ISAL_GZIP;
    state.next_in = src.data();
    state.avail_in = src.size();
    state.next_out = dst.data();
    state.avail_out = dst.size();

    switch (m) {
      case stateful:
        checked_isal_inflate(&state, eof, filename);
        return;
      case stateless:
        checked_isal_inflate_stateless(&state, isize, filename);
        return;
        break;
      default:                     // GCOVR_EXCL_LINE
        AR_FAIL("invalid method"); // GCOVR_EXCL_LINE
    }
  }

  //! Compressed data
  buffer src{};
  //! Destination for decompressed data; defaults to sufficient capacity for src
  buffer dst{};
  //! ISIZE for stateless decompression
  size_t isize = 0;
  //! Final block for stateful inflation?
  bool eof = true;
};

} // namespace

#define PARAMETERIZE_METHOD                                                    \
  GENERATE(isal_inflator::method::stateful, isal_inflator::method::stateless)

#define REQUIRE_INFLATE_THROWS_FOR(inflator, method, msg)                      \
  SECTION("without filename")                                                  \
  {                                                                            \
    REQUIRE_THROWS_STARTS_WITH(inflator(method), gzip_error, msg);             \
  }                                                                            \
                                                                               \
  SECTION("with filename")                                                     \
  {                                                                            \
    REQUIRE_THROWS_STARTS_WITH(inflator(method, "abcd0123"),                   \
                               gzip_error,                                     \
                               "error while decompressing 'abcd0123': " msg);  \
  }

#define REQUIRE_INFLATE_THROWS(inflator, message)                              \
  REQUIRE_INFLATE_THROWS_FOR(inflator, PARAMETERIZE_METHOD, message)

TEST_CASE("isal_compressor basic check")
{
  auto text = random_text_buffer(1024);
  isal_inflator inflator{ text };
  inflator(PARAMETERIZE_METHOD);

  REQUIRE(inflator.dst == text);
}

TEST_CASE("checked inflate fails on nullptr state")
{
  REQUIRE_THROWS_CONTAINS(checked_isal_inflate(nullptr, true),
                          assert_failed,
                          "~a43067e9e983afc7~");
  REQUIRE_THROWS_CONTAINS(checked_isal_inflate_stateless(nullptr, 1234),
                          assert_failed,
                          "~b235c2b9f5f022a8~");
}

TEST_CASE("checked inflate fails on end of input")
{
  isal_inflator inflator{ random_text_buffer(1024) };
  inflator.src.resize(inflator.src.size() - GENERATE(1, 5, 9, 32));

  REQUIRE_INFLATE_THROWS(inflator, "file is truncated");
}

TEST_CASE("checked inflate stateless fails if output is too large")
{
  isal_inflator inflator{ random_text_buffer(1024) };
  inflator.dst.resize(1023);

  REQUIRE_INFLATE_THROWS_FOR(inflator,
                             isal_inflator::stateless,
                             "content is larger than ISIZE");
}

TEST_CASE("checked inflate fail on HCRC mismatch")
{
  const auto data = random_text_buffer(1024);
  const auto header = "\37\213\10\2\0\0\0\0\0\377"_buffer;

  isal_inflator inflator;
  inflator.src.append(header);
  const size_t offset = inflator.src.size();
  inflator.src.append_u16(calculate_crc32(header) & 0xFFFF); // FHCRC
  inflator.src.resize(BGZF_MAX_BLOCK_SIZE);
  deflate_block_libdeflate(data,
                           inflator.src,
                           offset + 2,
                           PARAMETERIZE_COMPRESSION);
  inflator.src.append_u32(calculate_crc32(data));
  inflator.src.append_u32(data.size());
  inflator.dst.resize(data.size());
  inflator.isize = data.size();

  SECTION("matching HCRC")
  {
    REQUIRE_NOTHROW(inflator(PARAMETERIZE_METHOD));
    REQUIRE(inflator.dst == data);
  }

  inflator.src.put_u8(offset, inflator.src.at(offset) ^ 1);
  REQUIRE_INFLATE_THROWS(inflator, "incorrect header checksum");
}

TEST_CASE("checked inflate stateful fails on CRC32 mismatch")
{
  isal_inflator inflator{ random_text_buffer(1024) };
  const size_t offset = inflator.src.size() - 8;
  inflator.src.put_u32(offset, ~inflator.src.get_u32(offset));

  REQUIRE_INFLATE_THROWS_FOR(inflator,
                             isal_inflator::stateful,
                             "incorrect ISIZE or content checksum");
}

TEST_CASE("checked inflate stateless fails on CRC32 mismatch")
{
  isal_inflator inflator{ random_text_buffer(1024) };
  const size_t offset = inflator.src.size() - 8;
  inflator.src.put_u32(offset, ~inflator.src.get_u32(offset));

  REQUIRE_INFLATE_THROWS_FOR(inflator,
                             isal_inflator::stateless,
                             "incorrect content checksum");
}

TEST_CASE("checked inflate stateful fails on ISIZE mismatch")
{
  isal_inflator inflator{ random_text_buffer(1024) };
  auto offset = inflator.src.size() - 4;
  inflator.src.put_u32(offset, inflator.src.get_u32(offset) + GENERATE(-1, 1));

  REQUIRE_INFLATE_THROWS_FOR(inflator,
                             isal_inflator::stateful,
                             "incorrect ISIZE or content checksum");
}

TEST_CASE("checked inflate stateless fails on ISIZE mismatch")
{
  isal_inflator inflator{ random_text_buffer(1024) };
  auto offset = inflator.src.size() - 4;
  inflator.src.put_u32(offset, inflator.src.get_u32(offset) + GENERATE(-1, 1));
  inflator.isize -= 4;

  REQUIRE_INFLATE_THROWS_FOR(inflator,
                             isal_inflator::stateless,
                             "size of content does not match ISIZE");
}

TEST_CASE("checked inflate fails on corrupt block type")
{
  isal_inflator inflator{ random_text_buffer(1024) };
  inflator.src.put_u8(19, 0xFF); // nonsense header

  REQUIRE_INFLATE_THROWS(inflator, "invalid deflate block");
}

TEST_CASE("checked inflate fails on invalid gzip file")
{
  isal_inflator inflator{};
  inflator.isize = 1024;
  inflator.src = random_data_buffer(inflator.isize);
  inflator.dst.resize(inflator.isize);

  REQUIRE_INFLATE_THROWS(inflator, "not a valid gzip file");
}

////////////////////////////////////////////////////////////////////////////////
// inflate functions with nullptr

TEST_CASE("checked infate fails for non-zero avail_in with nullptr")
{
  buffer dst{ 1024 };
  inflate_state state;
  state.crc_flag = ISAL_GZIP;
  state.next_in = nullptr;
  state.avail_in = GENERATE(1, 1024);
  state.next_out = dst.data();
  state.avail_out = dst.size();

  SECTION("stateful")
  {
    REQUIRE_THROWS_CONTAINS(checked_isal_inflate(&state, true),
                            assert_failed,
                            "~c573155ef8a7fd67~");
  }

  SECTION("stateless")
  {
    REQUIRE_THROWS_CONTAINS(checked_isal_inflate_stateless(&state, 1234),
                            assert_failed,
                            "~34d6fcc3d35446d7~");
  }
}

TEST_CASE("checked inflate fails for non-zero avail_out with nullptr")
{
  buffer src{ 1024 };
  inflate_state state;
  state.crc_flag = ISAL_GZIP;
  state.next_in = src.data();
  state.avail_in = src.size();
  state.next_out = nullptr;
  state.avail_out = GENERATE(1, 1024);

  SECTION("stateful")
  {
    REQUIRE_THROWS_CONTAINS(checked_isal_inflate(&state, true),
                            assert_failed,
                            "~ea9fa8f9702f136f~");
  }

  SECTION("stateless")
  {
    REQUIRE_THROWS_CONTAINS(checked_isal_inflate_stateless(&state, 1234),
                            assert_failed,
                            "~044ddde4f5fc4b72~");
  }
}

TEST_CASE("checked inflate allows null/empty source")
{
  inflate_state state;
  state.crc_flag = ISAL_GZIP;

  buffer src;
  if (GENERATE(true, false)) {
    REQUIRE(src.data() == nullptr);
    src.reserve(1024);
  }

  state.next_in = src.data();
  state.avail_in = src.size();

  buffer dst{ 1024 };
  state.next_out = dst.data();
  state.avail_out = dst.size();

  SECTION("stateful")
  {
    REQUIRE_NOTHROW(checked_isal_inflate(&state, false));
    REQUIRE(state.next_in == src.data());
    REQUIRE(state.avail_in == src.size());
    REQUIRE(state.total_out == 0);
  }

  SECTION("stateless")
  {
    REQUIRE_THROWS_CONTAINS(checked_isal_inflate_stateless(&state, 1024),
                            gzip_error,
                            "file is truncated");
  }
}

TEST_CASE("checked inflate allows null/empty destination")
{
  inflate_state state;
  state.crc_flag = ISAL_GZIP;

  buffer src = gzip_data_into_buffer(random_text_buffer(1024));
  state.next_in = src.data();
  state.avail_in = src.size();

  buffer dst;
  if (GENERATE(true, false)) {
    REQUIRE(dst.data() == nullptr);
    dst.reserve(1024);
  }

  state.next_out = dst.data();
  state.avail_out = dst.size();

  SECTION("stateful")
  {
    REQUIRE_NOTHROW(checked_isal_inflate(&state, false));
    REQUIRE(state.next_out == dst.data());
    REQUIRE(state.avail_out == dst.size());
    REQUIRE(state.total_out == 0);
  }

  SECTION("stateless")
  {
    REQUIRE_THROWS_CONTAINS(checked_isal_inflate_stateless(&state, 1024),
                            gzip_error,
                            "content is larger than ISIZE");
  }
}

} // namespace adapterremoval
