// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "buffer_reader.hpp" // declarations
#include "debug.hpp"         // for AR_REQUIRE
#include "pragmas.hpp"       // for AR_UNLIKELY
#include <cstring>           // for memchr
#include <deque>             // for deque
#include <ostream>           // for ostream

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// buffer_reader

buffer_reader::buffer_reader(std::vector<std::string> filenames)
  : m_filenames(std::move(filenames))
  , m_remaining_eof(m_filenames.remaining_filenames())
{
}

buffer_reader::~buffer_reader()
{
  if (m_release) {
    try {
      for (auto& [buf, _] : m_buffers) {
        if (buf.capacity()) {
          m_release(std::move(buf));
        }
      }
    } catch (...) { // NOLINT(bugprone-empty-catch)
    }
  }
}

void
buffer_reader::add(buffer&& buf, bool eof)
{
  AR_REQUIRE(m_remaining_eof > 0);

  if (buf.size() || eof) {
    m_buffers.emplace_back(std::move(buf), eof);
    if (!m_buffer_ptr) {
      next_buffer();
    }

    if (eof) {
      m_remaining_eof--;
    }
  } else if (m_release && buf.capacity()) {
    m_release(std::move(buf));
  }
}

void
buffer_reader::close()
{
  m_read_len = std::string::npos;
  m_remaining_eof = 0;
  m_buffer_ptr = nullptr;
  m_buffer_end = nullptr;
  m_closed = true;

  for (; !m_buffers.empty(); m_buffers.pop_front()) {
    auto& buf = m_buffers.front().first;
    if (buf.capacity() && m_release) {
      m_release(std::move(buf));
    }
  }
}

buffer_reader::result
buffer_reader::getline(std::string& dst)
{
  // Clear the destination before anything is written to it; for simplicity this
  // is done even if nothing ends up being written to it (cache empty or at eof)
  if (m_read_len == std::string::npos) {
    m_read_len = 0;
    dst.clear();
  } else {
    AR_REQUIRE(m_read_len == dst.size());
  }

  bool at_eof = false;
  while (m_buffer_ptr && !at_eof) {
    const auto len = m_buffer_end - m_buffer_ptr;
    const auto* end = static_cast<const char*>(memchr(m_buffer_ptr, '\n', len));
    if (end) {
      // Excluding terminal \n
      dst.append(m_buffer_ptr, end - m_buffer_ptr);
      // Excluding terminal \r; this may have been added in the loop prior to
      // \n being found, if \n is the first character in the buffer. It is
      // therefore easiest to check the final string.
      if (AR_UNLIKELY(!dst.empty() && dst.back() == '\r')) {
        dst.pop_back();
      }

      m_read_len = std::string::npos;
      m_buffer_ptr = end + 1;

      m_filenames.inc_position();

      return result::complete;
    }

    // Can potentially introduce a \r; this is handled above.
    dst.append(m_buffer_ptr, len);
    m_read_len += len;

    if (m_release && m_buffers.front().first.capacity()) {
      // `release` may throw, so let the next `getline()` call resume from here
      m_buffer_ptr = &m_eof_char;
      m_buffer_end = &m_eof_char;

      m_release(std::move(m_buffers.front().first));
    }

    if (m_buffers.front().second) {
      if (m_read_len) {
        m_filenames.inc_position();
        at_eof = true;
      }

      m_filenames.next_file();
    }

    m_buffers.pop_front();
    next_buffer();
  }

  if (at_eof) {
    m_read_len = std::string::npos;
    return result::complete;
  } else if (eof()) {
    return result::eof;
  } else {
    return result::partial;
  }
}

void
buffer_reader::next_buffer()
{
  if (m_buffers.empty()) {
    m_buffer_ptr = nullptr;
    m_buffer_end = nullptr;
  } else {
    const auto& buf = m_buffers.front();

    if (buf.first.size() == 0) {
      // The buffer may not have allocated, so use a known, valid address
      m_buffer_ptr = &m_eof_char;
      m_buffer_end = &m_eof_char;
    } else {
      m_buffer_ptr = reinterpret_cast<const char*>(buf.first.data());
      m_buffer_end = m_buffer_ptr + buf.first.size();
    }
  }
}

std::ostream&
operator<<(std::ostream& os, buffer_reader::result r)
{
  switch (r) { // GCOVR_EXCL_BR_LINE
    case buffer_reader::result::complete:
      return os << "buffer_reader::result::complete";
    case buffer_reader::result::partial:
      return os << "buffer_reader::result::partial";
    case buffer_reader::result::eof:
      return os << "buffer_reader::result::eof";
    default:                                    // GCOVR_EXCL_LINE
      AR_FAIL("invalid buffer_reader::result"); // GCOVR_EXCL_LINE
  }
}

} // namespace adapterremoval
