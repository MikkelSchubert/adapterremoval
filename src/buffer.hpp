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
#include <string>      // for string
#include <string_view> // for string_view
#include <vector>      // for vector

namespace adapterremoval {

class buffer
{
  struct no_init;

public:
  /** Creates a buffer of the specified size, if non-zero */
  explicit buffer(size_t size = 0) { m_buffer.resize(size); }

  /** Returns pointer to the buffer, if any */
  [[nodiscard]] inline unsigned char* data()
  {
    return reinterpret_cast<unsigned char*>(m_buffer.data());
  }

  /** Returns pointer to the buffer, if any */
  [[nodiscard]] inline const unsigned char* data() const
  {
    return reinterpret_cast<const unsigned char*>(m_buffer.data());
  }

  /** Returns the current size of the buffer */
  [[nodiscard]] inline size_t size() const { return m_buffer.size(); }

  /** Returns the capacity of the buffer */
  [[nodiscard]] inline size_t capacity() const { return m_buffer.capacity(); }

  /** Changes the reported size of the buffer; must be 0 <= x <= capacity. */
  inline void resize(size_t size)
  {
    AR_REQUIRE(size <= capacity());
    m_buffer.resize(size);
  }

  /** Reserve additional space in the buffer */
  inline void reserve(size_t size) { m_buffer.reserve(size); }

  /** Append a string to the buffer */
  inline void append(std::string_view data)
  {
    append(data.data(), data.size());
  }

  /** Append data to the buffer */
  inline void append(const char* data, size_t length)
  {
    append(reinterpret_cast<const unsigned char*>(data), length);
  }

  /** Append data to the buffer */
  inline void append(const unsigned char* data, size_t length)
  {
    const auto* ptr = reinterpret_cast<const no_init*>(data);
    m_buffer.insert(m_buffer.end(), ptr, ptr + length);
  }

  /** Append a byte to the buffer */
  inline void append_u8(uint8_t value) { m_buffer.push_back({ value }); }

  /** Copies uint32 into buffer */
  inline void append_u32(uint32_t value)
  {
    append_u8(value & 0xFFU);
    append_u8((value >> 8U) & 0xFFU);
    append_u8((value >> 16U) & 0xFFU);
    append_u8((value >> 24U) & 0xFFU);
  }

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
  std::vector<no_init> m_buffer{};
};

using buffer_vec = std::vector<buffer>;

} // namespace adapterremoval
