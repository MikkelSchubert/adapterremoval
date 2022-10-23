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
#include <memory>    // for unique_ptr
#include <vector>    // for vector

namespace adapterremoval {

class buffer
{
public:
  /** Creates a buffer of the specified size, if non-zero */
  buffer(size_t size = 0);

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

  buffer(buffer&& other) noexcept;
  buffer& operator=(buffer&& other) noexcept;

  buffer(const buffer& other) = delete;
  buffer& operator=(const buffer& other) = delete;

  ~buffer() = default;

private:
  //! Size of the buffer
  size_t m_size;
  //! Capacity of the buffer
  size_t m_capacity;
  //! Backing buffer
  std::unique_ptr<unsigned char[]> m_buffer;
};

using buffer_vec = std::vector<buffer>;

////////////////////////////////////////////////////////////////////////////////

inline buffer::buffer(size_t size)
  : m_size(size)
  , m_capacity(size)
  , m_buffer()
{
  if (m_size != 0) {
    m_buffer = std::make_unique<unsigned char[]>(m_size);
  }
}

inline buffer::buffer(buffer&& other) noexcept
  : buffer(0)
{
  *this = std::move(other);
}

inline buffer&
buffer::operator=(buffer&& other) noexcept
{
  if (this != &other) {
    m_size = other.m_size;
    m_capacity = other.m_capacity;
    m_buffer = std::move(other.m_buffer);

    other.m_size = 0;
    other.m_capacity = 0;
  }

  return *this;
}

inline unsigned char*
buffer::get()
{
  return m_buffer.get();
}

inline const unsigned char*
buffer::get() const
{
  return m_buffer.get();
}

inline char*
buffer::get_signed()
{
  return reinterpret_cast<char*>(m_buffer.get());
}

inline const char*
buffer::get_signed() const
{
  return reinterpret_cast<const char*>(m_buffer.get());
}

inline size_t
buffer::size() const
{
  return m_size;
}

inline size_t
buffer::capacity() const
{
  return m_capacity;
}

inline void
buffer::resize(size_t size)
{
  AR_REQUIRE(size <= m_capacity);
  m_size = size;
}

inline void
buffer::write_u32(size_t offset, uint32_t value)
{
  AR_REQUIRE(offset + 4 <= m_size);

  m_buffer[offset + 0] = value & 0xFF;
  m_buffer[offset + 1] = (value >> 8) & 0xFF;
  m_buffer[offset + 2] = (value >> 16) & 0xFF;
  m_buffer[offset + 3] = (value >> 24) & 0xFF;
}

} // namespace adapterremoval
