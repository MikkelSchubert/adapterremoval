/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "buffer.hpp"
#include "testing.hpp"
#include <cstring>

namespace adapterremoval {

TEST_CASE("buffer defaults to empty/nullptr")
{
  buffer buf;

  REQUIRE(buf.size() == 0);
  REQUIRE(buf.capacity() == 0);
  REQUIRE(buf.get() == nullptr);
  REQUIRE(buf.get_signed() == nullptr);
}

TEST_CASE("buffer with size allocates")
{
  buffer buf(7);
  const buffer& cbuf = buf;

  REQUIRE(buf.size() == 7);
  REQUIRE(buf.capacity() == 7);
  REQUIRE(buf.get() != nullptr);
  REQUIRE(buf.get() == cbuf.get());
  REQUIRE(buf.get_signed() != nullptr);
  REQUIRE(buf.get_signed() == cbuf.get_signed());
  REQUIRE(static_cast<void*>(buf.get()) ==
          static_cast<void*>(buf.get_signed()));
}

TEST_CASE("buffer resize changes only apparent size")
{
  buffer buf(7);

  buf.resize(3);
  REQUIRE(buf.size() == 3);
  REQUIRE(buf.capacity() == 7);
}

TEST_CASE("buffer resize must be less than capacity")
{
  buffer buf(7);

  REQUIRE_THROWS_AS(buf.resize(8), assert_failed);
  REQUIRE(buf.size() == 7);
  REQUIRE(buf.capacity() == 7);
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
  auto src_ptr = src.get();
  auto src_signed_ptr = src.get_signed();
  buffer dst(std::move(src));

  REQUIRE(src.size() == 0);
  REQUIRE(src.capacity() == 0);
  REQUIRE(src.get() == nullptr);
  REQUIRE(src.get_signed() == nullptr);

  REQUIRE(dst.size() == 7);
  REQUIRE(dst.capacity() == 7);
  REQUIRE(dst.get() == src_ptr);
  REQUIRE(dst.get_signed() == src_signed_ptr);
}

TEST_CASE("buffer move assignment clears source")
{
  buffer src(7);
  auto src_ptr = src.get();
  auto src_signed_ptr = src.get_signed();
  buffer dst(7);

  dst = std::move(src);

  REQUIRE(src.size() == 0);
  REQUIRE(src.capacity() == 0);
  REQUIRE(src.get() == nullptr);
  REQUIRE(src.get_signed() == nullptr);

  REQUIRE(dst.size() == 7);
  REQUIRE(dst.capacity() == 7);
  REQUIRE(dst.get() == src_ptr);
  REQUIRE(dst.get_signed() == src_signed_ptr);
}

TEST_CASE("buffer move assignment to self is no-op")
{
  buffer buf(7);
  auto buf_ptr = buf.get();
  auto buf_signed_ptr = buf.get_signed();

  buf = std::move(buf);
  REQUIRE(buf.size() == 7);
  REQUIRE(buf.capacity() == 7);
  REQUIRE(buf.get() == buf_ptr);
  REQUIRE(buf.get_signed() == buf_signed_ptr);
}

TEST_CASE("buffer write_u32 is le")
{
  buffer buf(4);
  buf.write_u32(0, 0xDEADBEEF);
  auto ptr = buf.get();

  REQUIRE(ptr[3] == 0xDE);
  REQUIRE(ptr[2] == 0xAD);
  REQUIRE(ptr[1] == 0xBE);
  REQUIRE(ptr[0] == 0xEF);
}

TEST_CASE("buffer write_u32 with offset")
{
  buffer buf(6);
  auto ptr = buf.get();
  ::bzero(ptr, buf.size());

  buf.write_u32(1, 0xDEADBEEF);

  REQUIRE(ptr[5] == 0);
  REQUIRE(ptr[4] == 0xDE);
  REQUIRE(ptr[3] == 0xAD);
  REQUIRE(ptr[2] == 0xBE);
  REQUIRE(ptr[1] == 0xEF);
  REQUIRE(ptr[0] == 0);
}

TEST_CASE("buffer write_u32 requires space")
{
  SECTION("capacity must be sufficient")
  {
    buffer buf(3);
    REQUIRE_THROWS_AS(buf.write_u32(0, 0xDEADBEEF), assert_failed);
  }

  SECTION("size must be sufficient")
  {
    buffer buf(4);
    buf.resize(3);
    REQUIRE_THROWS_AS(buf.write_u32(0, 0xDEADBEEF), assert_failed);
  }

  SECTION("size from offset must be sufficient")
  {
    buffer buf(4);
    REQUIRE_THROWS_AS(buf.write_u32(1, 0xDEADBEEF), assert_failed);
  }
}

} // namespace adapterremoval
