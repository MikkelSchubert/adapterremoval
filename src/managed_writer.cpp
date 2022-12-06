/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <cerrno>       // for errno, EMFILE
#include <mutex>        // for mutex, lock_guard
#include <stdexcept>    // for runtime_error
#include <system_error> // for errc, errc::too_many_files_open

#if _POSIX_C_SOURCE >= 200112L
#include <fcntl.h> // for posix_fadvise
#endif

#include "debug.hpp"          // for AR_REQUIRE
#include "logging.hpp"        // for log
#include "managed_writer.hpp" // declarations

namespace adapterremoval {

static std::mutex g_writer_lock;
managed_writer* managed_writer::s_head = nullptr;
managed_writer* managed_writer::s_tail = nullptr;
bool managed_writer::s_warning_printed = false;

managed_writer::managed_writer(const std::string& filename)
  : m_filename(filename)
  , m_stream()
  , m_created(false)
  , m_prev(nullptr)
  , m_next(nullptr)
{
  m_stream.exceptions(std::ofstream::failbit | std::ofstream::badbit);
}

managed_writer::~managed_writer()
{
  close();
}

FILE*
managed_writer::fopen(const std::string& filename, const char* mode)
{
  AR_REQUIRE(mode);

  while (true) {
    FILE* handle = ::fopen(filename.c_str(), mode);

    if (handle) {
#if _POSIX_C_SOURCE >= 200112L
      // Hint that we'll (only) be doing sequential reads
      posix_fadvise(fileno(handle), 0, 0, POSIX_FADV_SEQUENTIAL);
      posix_fadvise(fileno(handle), 0, 0, POSIX_FADV_WILLNEED);
      posix_fadvise(fileno(handle), 0, 0, POSIX_FADV_NOREUSE);
#endif

      return handle;
    } else if (errno == EMFILE) {
      std::lock_guard<std::mutex> lock(g_writer_lock);
      managed_writer::close_tail_writer();
    } else {
      return nullptr;
    }
  }
}

void
managed_writer::write_buffer(const buffer& buf, bool flush)
{
  if (buf.size() || flush) {
    std::lock_guard<std::mutex> lock(g_writer_lock);
    managed_writer::open_writer(this);

    if (buf.size()) {
      m_stream.write(buf.get_signed(), buf.size());
    }

    if (flush) {
      m_stream.flush();
    }
  }
}

void
managed_writer::write_buffers(const buffer_vec& buffers, bool flush)
{
  if (buffers.size() || flush) {
    std::lock_guard<std::mutex> lock(g_writer_lock);
    managed_writer::open_writer(this);

    for (auto& buf : buffers) {
      if (buf.size()) {
        m_stream.write(buf.get_signed(), buf.size());
      }
    }

    if (flush) {
      m_stream.flush();
    }
  }
}

void
managed_writer::write_string(const std::string& str, bool flush)
{
  std::lock_guard<std::mutex> lock(g_writer_lock);
  if (str.size() || flush) {
    managed_writer::open_writer(this);

    m_stream.write(str.data(), str.length());
    if (flush) {
      m_stream.flush();
    }
  }
}

void
managed_writer::close()
{
  std::lock_guard<std::mutex> lock(g_writer_lock);

  managed_writer::remove_writer(this);
  if (m_stream.is_open()) {
    m_stream.close();
  }
}

const std::string&
managed_writer::filename() const
{
  return m_filename;
}

void
managed_writer::open_writer(managed_writer* ptr)
{
  const std::ios_base::openmode mode =
    ptr->m_created ? std::ofstream::app : std::ofstream::trunc;

  while (!ptr->m_stream.is_open()) {
    try {
      ptr->m_stream.open(ptr->m_filename, std::ofstream::binary | mode);
      break;
    } catch (const std::ofstream::failure&) {
      if (errno != int(std::errc::too_many_files_open)) {
        throw;
      }
    }

    managed_writer::close_tail_writer();
  }

  if (ptr != s_head) {
    managed_writer::remove_writer(ptr);
    managed_writer::add_head_writer(ptr);
  }

  ptr->m_created = true;
}

void
managed_writer::remove_writer(managed_writer* ptr)
{
  AR_REQUIRE(ptr);
  AR_REQUIRE(!s_head == !s_tail);
  AR_REQUIRE(!s_head || !s_head->m_prev);
  AR_REQUIRE(!s_tail || !s_tail->m_next);

  if (ptr == s_head) {
    s_head = ptr->m_next;
  }

  if (ptr == s_tail) {
    s_tail = ptr->m_prev;
  }

  AR_REQUIRE(!s_head == !s_tail);

  if (ptr->m_prev) {
    ptr->m_prev->m_next = ptr->m_next;
  }

  if (ptr->m_next) {
    ptr->m_next->m_prev = ptr->m_prev;
  }

  ptr->m_prev = nullptr;
  ptr->m_next = nullptr;

  AR_REQUIRE(ptr != s_head);
  AR_REQUIRE(ptr != s_tail);
  AR_REQUIRE(!ptr->m_prev);
  AR_REQUIRE(!ptr->m_next);
  AR_REQUIRE(!s_head || !s_head->m_prev);
  AR_REQUIRE(!s_tail || !s_tail->m_next);
}

void
managed_writer::add_head_writer(managed_writer* ptr)
{
  AR_REQUIRE(!ptr->m_prev);
  AR_REQUIRE(!ptr->m_next);
  AR_REQUIRE(!s_head == !s_tail);
  if (s_head) {
    ptr->m_next = s_head;
    s_head->m_prev = ptr;
  }

  s_head = ptr;

  if (!s_tail) {
    s_tail = ptr;
  }

  AR_REQUIRE(s_head && s_tail);
  AR_REQUIRE(!s_head->m_prev);
  AR_REQUIRE(!s_tail->m_next);
}

void
managed_writer::close_tail_writer()
{
  AR_REQUIRE(!s_head == !s_tail);
  if (!s_warning_printed) {
    log::warn() << "Number of available file-handles (ulimit -n) is too low. "
                << "AdapterRemoval will dynamically close/re-open files as "
                << "required, but performance may suffer as a result.";

    s_warning_printed = true;
  }

  if (s_tail) {
    AR_REQUIRE(s_tail->m_stream.is_open());

    s_tail->m_stream.close();
    managed_writer::remove_writer(s_tail);
    return;
  }

  throw std::runtime_error(
    "available number of file-handles too low; could not open any files");
}

} // namespace adapterremoval
