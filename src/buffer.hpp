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

#include "debug.hpp"   // for AR_REQUIRE
#include <algorithm>   // for max
#include <cstddef>     // for size_t
#include <cstdint>     // for uint32_t
#include <cstring>     // for memcpy
#include <string_view> // for string_view
#include <vector>      // for vector

namespace adapterremoval {

class buffer
{
  struct no_init;

public:
  /** Creates a buffer of the specified size, if non-zero */
  explicit buffer(size_t size = 0) { resize(size); }

  buffer(buffer&& other) noexcept
  {
    free();
    std::swap(m_data, other.m_data);
    std::swap(m_capacity, other.m_capacity);
    std::swap(m_size, other.m_size);
  }

  buffer& operator=(buffer&& other) noexcept
  {
    free();
    std::swap(m_data, other.m_data);
    std::swap(m_capacity, other.m_capacity);
    std::swap(m_size, other.m_size);
    return *this;
  }

  ~buffer() { free(); }

  /** Returns pointer to the buffer, if any */
  [[nodiscard]] inline uint8_t* data() { return m_data; }

  /** Returns pointer to the buffer, if any */
  [[nodiscard]] inline const uint8_t* data() const { return m_data; }

  /** Returns the nth value in the buffer */
  [[nodiscard]] uint8_t at(size_t n) const
  {
    AR_REQUIRE(n < m_size);
    return m_data[n];
  }

  /** Returns the current size of the buffer */
  [[nodiscard]] inline size_t size() const { return m_size; }

  /** Returns the capacity of the buffer */
  [[nodiscard]] inline size_t capacity() const { return m_capacity; }

  /** Changes the reported size of the buffer */
  inline void resize(size_t n)
  {
    if (n > m_capacity) {
      reserve(n);
    }

    m_size = n;
  }

  /** Reserve additional space in the buffer */
  inline void reserve(size_t n)
  {
    if (n > m_capacity) {
      if (m_data) {
        m_data = reinterpret_cast<uint8_t*>(std::realloc(m_data, n));
      } else {
        m_data = reinterpret_cast<uint8_t*>(std::malloc(n));
      }

      AR_REQUIRE(m_data);
      m_capacity = n;
    }
  }

  /** Append a string to the buffer */
  inline void append(std::string_view data)
  {
    append(data.data(), data.size());
  }

  /** Append a string to the buffer */
  inline void append(const buffer& data) { append(data.data(), data.size()); }

  /** Append data to the buffer */
  inline void append(const char* data, size_t length)
  {
    append(reinterpret_cast<const uint8_t*>(data), length);
  }

  /** Append data to the buffer */
  inline void append(const uint8_t* data, size_t length)
  {
    grow_to_fit(m_size + length);
    std::memcpy(m_data + m_size, data, length);
    m_size += length;
  }

  /** Append a byte to the buffer */
  inline void append_u8(uint8_t value)
  {
    grow_to_fit(m_size + 1);
    m_data[m_size++] = value;
  }

  /** Append 16 bit integer to buffer in LE orientation */
  inline void append_u16(uint16_t value)
  {
    grow_to_fit(m_size + 2);
    m_data[m_size + 0] = value & 0xFFU;
    m_data[m_size + 1] = (value >> 8U) & 0xFFU;
    m_size += 2;
  }

  /** Append 32 bit integer to buffer in LE orientation */
  inline void append_u32(uint32_t value)
  {
    grow_to_fit(m_size + 4);
    m_data[m_size + 0] = value & 0xFFU;
    m_data[m_size + 1] = (value >> 8U) & 0xFFU;
    m_data[m_size + 2] = (value >> 16U) & 0xFFU;
    m_data[m_size + 3] = (value >> 24U) & 0xFFU;
    m_size += 4;
  }

  /** Append signed 32 bit integer to buffer in LE orientation */
  inline void append_i32(int32_t value)
  {
    append_u32(static_cast<uint32_t>(value));
  }

  /** Insert 16 bit integer into buffer in LE orientation at offset */
  inline void put_u16(size_t offset, uint16_t value)
  {
    AR_REQUIRE(offset + 2 <= size());
    m_data[offset + 0] = value & 0xFFU;
    m_data[offset + 1] = (value >> 8U) & 0xFFU;
  }

  /** Insert 32 bit integer into buffer in LE orientation at offset */
  inline void put_u32(size_t offset, uint32_t value)
  {
    AR_REQUIRE(offset + 4 <= size());
    m_data[offset + 0] = value & 0xFFU;
    m_data[offset + 1] = (value >> 8U) & 0xFFU;
    m_data[offset + 2] = (value >> 16U) & 0xFFU;
    m_data[offset + 3] = (value >> 24U) & 0xFFU;
  }

  buffer(const buffer& other) = delete;
  buffer& operator=(const buffer& other) = delete;

private:
  inline void free()
  {
    if (m_data) {
      std::free(m_data);
      m_data = nullptr;
      m_capacity = 0;
      m_size = 0;
    }
  }

  inline void grow_to_fit(size_t n)
  {
    if (n > m_capacity) {
      size_t increase_to = std::max<size_t>(1, m_capacity);
      while (n > increase_to) {
        increase_to *= 2;
      }

      reserve(increase_to);
    }
  }

  //! Backing buffer
  uint8_t* m_data = nullptr;
  //! Capacity of backing buffer
  size_t m_capacity = 0;
  //! Size of backing buffer sized
  size_t m_size = 0;
};

using buffer_vec = std::vector<buffer>;

} // namespace adapterremoval
