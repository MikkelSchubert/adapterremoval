// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "buffer.hpp"   // declarations
#include "testing.hpp"  // for TEST_CASE, REQUIRE, ...
#include <cstring>      // for memset
#include <strutils.hpp> // for log_escape
#include <utility>      // for move

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// implementations of buffer test functions

/**
 * Creates debug representation of an buffer; should only be used on fully
 * initialized buffers, to avoid reading uninitialized values.
 */
std::ostream&
operator<<(std::ostream& os, const buffer& value)
{
  const std::string_view buf{ reinterpret_cast<const char*>(value.data()),
                              value.size() };

  return os << "buffer{data=" << log_escape(buf)
            << ", capacity=" << value.capacity() << ", size=" << value.size()
            << "}";
}

/** Initialize a buffer using a string literal */
buffer
operator""_buffer(const char* s, size_t l)
{
  buffer buf;
  buf.append(s, l);

  return buf;
}

////////////////////////////////////////////////////////////////////////////////
// testing helper functions

TEST_CASE("buffer to string")
{
  buffer buf;
  buf.append("test string");
  buf.append_u32(0xDEADBEEF);

  CHECK(Catch::fallbackStringifier(buf) ==
        "buffer{data='test string\\xef\\xbe\\xad\\xde', capacity=16, size=15}");
}

TEST_CASE("buffer literals")
{
  CHECK(buffer{} == ""_buffer);
  CHECK_FALSE(buffer{} == "test string"_buffer);

  SECTION("string")
  {
    buffer buf;
    buf.append("test string");
    CHECK(buf == "test string"_buffer);
  }

  SECTION("unprintable")
  {
    buffer buf;
    buf.append_u32(0xDEADBEEF);
    CHECK(buf == "\xef\xbe\xad\xde"_buffer);
  }
}

TEST_CASE("buffer equality operator")
{
  CHECK(buffer{} == buffer{});

  SECTION("reserved space is not considered")
  {
    buffer buff;
    buff.reserve(16);
    CHECK(buff == buffer{});
  }

  SECTION("content is considered")
  {
    buffer buf1;
    buffer buf2;
    REQUIRE(buf1 == buf2);

    buf1.append("foo bar");
    REQUIRE_FALSE(buf1 == buf2);

    buf2.append("foo bar");
    REQUIRE(buf1 == buf2);

    buf2.append_u32(0xDEADBEEF);
    REQUIRE_FALSE(buf1 == buf2);

    buf1.append_u32(0xDEADBEEF);
    REQUIRE(buf1 == buf2);
  }
}

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("buffer defaults to empty/nullptr")
{
  buffer buf;

  REQUIRE(buf.size() == 0);
  REQUIRE(buf.capacity() == 0);
  REQUIRE(buf.data() == nullptr);
}

TEST_CASE("buffer with size allocates")
{
  buffer buf(7);

  REQUIRE(buf.size() == 7);
  REQUIRE(buf.capacity() == 7);
  REQUIRE(buf.data() != nullptr);
  REQUIRE(buf.data() == static_cast<const buffer&>(buf).data());
}

TEST_CASE("buffer resize capacity when needed")
{
  buffer buf(7);

  buf.resize(3);
  REQUIRE(buf.size() == 3);
  REQUIRE(buf.capacity() == 7);

  buf.resize(8);
  REQUIRE(buf.size() == 8);
  REQUIRE(buf.capacity() == 8);
}

TEST_CASE("buffer resizing up is allowed")
{
  buffer buf(7);
  buf.resize(3);
  buf.resize(7);
  REQUIRE(buf.size() == 7);
  REQUIRE(buf.capacity() == 7);
}

TEST_CASE("buffer move constructor clears source")
{
  buffer src(7);
  auto* src_ptr = src.data();
  buffer dst(std::move(src));

  REQUIRE(src.size() == 0);
  REQUIRE(src.capacity() == 0);
  REQUIRE(src.data() == nullptr);

  REQUIRE(dst.size() == 7);
  REQUIRE(dst.capacity() == 7);
  REQUIRE(dst.data() == src_ptr);
}

TEST_CASE("buffer move assignment clears source")
{
  buffer src(7);
  auto* src_ptr = src.data();
  buffer dst(7);

  dst = std::move(src);

  REQUIRE(src.size() == 0);
  REQUIRE(src.capacity() == 0);
  REQUIRE(src.data() == nullptr);

  REQUIRE(dst.size() == 7);
  REQUIRE(dst.capacity() == 7);
  REQUIRE(dst.data() == src_ptr);
}

TEST_CASE("buffer append u8")
{
  buffer buf;
  buf.append_u8(0xDE);
  buf.append_u8(0xAD);
  buf.append_u8(0xBE);
  buf.append_u8(0xEF);

  REQUIRE(buf.size() == 4);
  REQUIRE(buf.at(0) == 0xDE);
  REQUIRE(buf.at(1) == 0xAD);
  REQUIRE(buf.at(2) == 0xBE);
  REQUIRE(buf.at(3) == 0xEF);
}

TEST_CASE("buffer append u16")
{
  buffer buf;
  buf.append_u16(0xDEAD);
  buf.append_u16(0xBEEF);

  REQUIRE(buf.size() == 4);
  REQUIRE(buf.at(0) == 0xAD);
  REQUIRE(buf.at(1) == 0xDE);
  REQUIRE(buf.at(2) == 0xEF);
  REQUIRE(buf.at(3) == 0xBE);
}

TEST_CASE("buffer append_u32 is le")
{
  buffer buf;
  buf.append_u32(0xDEADBEEF);

  REQUIRE(buf.at(3) == 0xDE);
  REQUIRE(buf.at(2) == 0xAD);
  REQUIRE(buf.at(1) == 0xBE);
  REQUIRE(buf.at(0) == 0xEF);
}

TEST_CASE("buffer append_i32 is le")
{
  buffer buf;
  buf.append_i32(0xDEADBEEF);

  REQUIRE(buf.at(3) == 0xDE);
  REQUIRE(buf.at(2) == 0xAD);
  REQUIRE(buf.at(1) == 0xBE);
  REQUIRE(buf.at(0) == 0xEF);
}

TEST_CASE("buffer append_u32 with offset")
{
  buffer buf;
  buf.append_u8(0xFF);
  buf.append_u32(0xDEADBEEF);
  buf.append_u8(0xFF);

  REQUIRE(buf.size() == 6);

  REQUIRE(buf.at(5) == 0xFF);
  REQUIRE(buf.at(4) == 0xDE);
  REQUIRE(buf.at(3) == 0xAD);
  REQUIRE(buf.at(2) == 0xBE);
  REQUIRE(buf.at(1) == 0xEF);
  REQUIRE(buf.at(0) == 0xFF);
}

TEST_CASE("buffer put u16 is le")
{
  buffer buf;
  buf.append_u32(0xFFFFFFFF);
  buf.append_u32(0xFFFFFFFF);

  buf.put_u16(2, 0xDEAD);
  buf.put_u16(4, 0xBEEF);

  REQUIRE(buf.size() == 8);
  REQUIRE(buf.at(0) == 0xFF);
  REQUIRE(buf.at(1) == 0xFF);
  REQUIRE(buf.at(2) == 0xAD);
  REQUIRE(buf.at(3) == 0xDE);
  REQUIRE(buf.at(4) == 0xEF);
  REQUIRE(buf.at(5) == 0xBE);
  REQUIRE(buf.at(6) == 0xFF);
  REQUIRE(buf.at(7) == 0xFF);
}

TEST_CASE("buffer put u32 is le")
{
  buffer buf;
  buf.append_u32(0xFFFFFFFF);
  buf.append_u32(0xFFFFFFFF);

  buf.put_u32(2, 0xDEADBEEF);

  REQUIRE(buf.size() == 8);
  REQUIRE(buf.at(0) == 0xFF);
  REQUIRE(buf.at(1) == 0xFF);
  REQUIRE(buf.at(2) == 0xEF);
  REQUIRE(buf.at(3) == 0xBE);
  REQUIRE(buf.at(4) == 0xAD);
  REQUIRE(buf.at(5) == 0xDE);
  REQUIRE(buf.at(6) == 0xFF);
  REQUIRE(buf.at(7) == 0xFF);
}

TEST_CASE("append buffer")
{
  buffer buf1;
  buf1.append_u32(0xDEADBEEF);

  buffer buf2;
  buf2.append_u8(0xFF);
  buf2.append(buf1);

  REQUIRE(buf1.size() == 4);
  REQUIRE(buf2.size() == 5);
  REQUIRE(buf2.at(0) == 0xFF);
  REQUIRE(buf2.at(1) == 0xEF);
  REQUIRE(buf2.at(2) == 0xBE);
  REQUIRE(buf2.at(3) == 0xAD);
  REQUIRE(buf2.at(4) == 0xDE);
}

} // namespace adapterremoval
