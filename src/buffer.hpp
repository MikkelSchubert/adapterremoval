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
#pragma once

#include "debug.hpp" // for AR_REQUIRE
#include <algorithm> // for max
#include <cstddef>   // for size_t
#include <cstdint>   // for uint32_t
#include <vector>    // for vector

namespace adapterremoval {

class buffer
{
public:
  /** Creates a buffer of the specified size, if non-zero */
  explicit buffer(size_t size = 0);

  /** Returns pointer to the buffer, if any */
  unsigned char* get();
  /** Returns pointer to the buffer, if any */
  const unsigned char* get() const;
  /** Returns pointer to the buffer, if any */
  char* get_signed();
  /** Returns pointer to the buffer, if any */
  const char* get_signed() const;

  /** Returns the current size of the buffer */
  size_t size() const;
  /** Returns the capacity of the buffer */
  size_t capacity() const;

  /** Changes the reported size of the buffer; must be 0 <= x <= capacity. */
  void resize(size_t size);

  /** Copies uint32 into buffer */
  void write_u32(size_t offset, uint32_t value);

  buffer(buffer&& other) noexcept = default;
  buffer& operator=(buffer&& other) noexcept = default;

  buffer(const buffer& other) = delete;
  buffer& operator=(const buffer& other) = delete;

  ~buffer() = default;

private:
  /** Construction intentionally omitted, to prevent default initialization */
  struct no_init
  {
    unsigned char value;
  };

  static_assert(sizeof(no_init) == sizeof(decltype(no_init::value)),
                "unexpected size of buffer::no_init");
  static_assert(alignof(no_init) == alignof(decltype(no_init::value)),
                "unexpected alignment of buffer::no_init");

  //! Backing buffer
  std::vector<no_init> m_buffer;
};

using buffer_vec = std::vector<buffer>;

////////////////////////////////////////////////////////////////////////////////

inline buffer::buffer(size_t size)
  : m_buffer()
{
  m_buffer.resize(size);
}

inline unsigned char*
buffer::get()
{
  return reinterpret_cast<unsigned char*>(m_buffer.data());
}

inline const unsigned char*
buffer::get() const
{
  return reinterpret_cast<const unsigned char*>(m_buffer.data());
}

inline char*
buffer::get_signed()
{
  return reinterpret_cast<char*>(get());
}

inline const char*
buffer::get_signed() const
{
  return reinterpret_cast<const char*>(get());
}

inline size_t
buffer::size() const
{
  return m_buffer.size();
}

inline size_t
buffer::capacity() const
{
  return m_buffer.capacity();
}

inline void
buffer::resize(size_t size)
{
  AR_REQUIRE(size <= capacity());
  m_buffer.resize(size);
}

inline void
buffer::write_u32(size_t offset, uint32_t value)
{
  AR_REQUIRE(offset + 4 <= size());

  unsigned char* dst = get() + offset;
  dst[0] = value & 0xFFu;
  dst[1] = (value >> 8) & 0xFFu;
  dst[2] = (value >> 16) & 0xFFu;
  dst[3] = (value >> 24) & 0xFFu;
}

} // namespace adapterremoval
