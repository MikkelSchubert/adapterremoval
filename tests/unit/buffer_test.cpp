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
#include "buffer.hpp"  // declarations
#include "errors.hpp"  // for assert_failed
#include "testing.hpp" // for catch.hpp, StringMaker
#include <cstring>     // for memset
#include <utility>     // for move

namespace adapterremoval {

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

TEST_CASE("buffer write_u32 is le")
{
  buffer buf;
  buf.append_u32(0xDEADBEEF);
  const auto* ptr = buf.data();

  REQUIRE(ptr[3] == 0xDE);
  REQUIRE(ptr[2] == 0xAD);
  REQUIRE(ptr[1] == 0xBE);
  REQUIRE(ptr[0] == 0xEF);
}

TEST_CASE("buffer write_u32 with offset")
{
  buffer buf;
  buf.append_u8(0);
  buf.append_u32(0xDEADBEEF);
  buf.append_u8(0);

  REQUIRE(buf.size() == 6);

  const auto* ptr = buf.data();
  REQUIRE(ptr[5] == 0);
  REQUIRE(ptr[4] == 0xDE);
  REQUIRE(ptr[3] == 0xAD);
  REQUIRE(ptr[2] == 0xBE);
  REQUIRE(ptr[1] == 0xEF);
  REQUIRE(ptr[0] == 0);
}

} // namespace adapterremoval
