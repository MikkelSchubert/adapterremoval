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
#include "managed_io.hpp" // declarations
#include "debug.hpp"      // for AR_REQUIRE
#include "errors.hpp"     // for io_error
#include "logging.hpp"    // for log::warn
#include "strutils.hpp"   // for log_escape
#include <cerrno>         // for EMFILE, errno
#include <cstddef>        // for size_t
#include <cstdio>         // for fopen, fread, fwrite, ...
#include <mutex>          // for mutex, lock_guard
#include <string>         // for string

namespace adapterremoval {

class writer_list
{
public:
  writer_list() = default;

  ~writer_list()
  {
    AR_REQUIRE(m_head == nullptr);
    AR_REQUIRE(m_tail == nullptr);
  }

  FILE* fopen(const std::string& filename, const char* mode)
  {
    while (true) {
      FILE* handle = ::fopen(filename.c_str(), mode);

      if (handle) {
        AR_REQUIRE(!::ferror_unlocked(handle));
        return handle;
      } else if (errno == EMFILE) {
        close_one();
      } else {
        throw io_error("failed to open file", errno);
      }
    }
  }

  /** Adds writer to list of inactive writers */
  void add(managed_writer* writer)
  {
    std::lock_guard<std::mutex> lock(m_lock);

    AR_REQUIRE(!writer->m_prev);
    AR_REQUIRE(!writer->m_next);
    AR_REQUIRE(!m_head == !m_tail);
    if (m_head) {
      writer->m_next = m_head;
      m_head->m_prev = writer;
    }

    m_head = writer;

    if (!m_tail) {
      m_tail = writer;
    }

    AR_REQUIRE(m_head && m_tail);
    AR_REQUIRE(!m_head->m_prev);
    AR_REQUIRE(!m_tail->m_next);
  }

  /* Removes the writer from the list of inactive writers */
  void remove(managed_writer* writer)
  {
    std::lock_guard<std::mutex> lock(m_lock);

    AR_REQUIRE(!m_head == !m_tail);
    AR_REQUIRE(!m_head || !m_head->m_prev);
    AR_REQUIRE(!m_tail || !m_tail->m_next);

    if (writer == m_head) {
      m_head = writer->m_next;
    }

    if (writer == m_tail) {
      m_tail = writer->m_prev;
    }

    AR_REQUIRE(!m_head == !m_tail);

    if (writer->m_prev) {
      writer->m_prev->m_next = writer->m_next;
    }

    if (writer->m_next) {
      writer->m_next->m_prev = writer->m_prev;
    }

    writer->m_prev = nullptr;
    writer->m_next = nullptr;

    AR_REQUIRE(writer != m_head);
    AR_REQUIRE(writer != m_tail);
    AR_REQUIRE(!writer->m_prev);
    AR_REQUIRE(!writer->m_next);
    AR_REQUIRE(!m_head || !m_head->m_prev);
    AR_REQUIRE(!m_tail || !m_tail->m_next);
  }

  writer_list(const writer_list&) = delete;
  writer_list(writer_list&&) = delete;
  writer_list& operator=(const writer_list&) = delete;
  writer_list& operator=(writer_list&&) = delete;

private:
  /** Try to close the least recently used writer */
  void close_one()
  {
    AR_REQUIRE(!m_head == !m_tail);
    if (!m_warning_printed) {
      log::warn() << "Number of available file-handles (ulimit -n) is too low. "
                  << "AdapterRemoval will dynamically close/re-open files as "
                  << "required, but performance may suffer as a result.";

      m_warning_printed = true;
    }

    if (m_tail) {
      if (fclose(m_tail->m_file)) {
        m_tail->m_file = nullptr;
        throw io_error("failed to close file", errno);
      }
      m_tail->m_file = nullptr;

      remove(m_tail);
    } else {
      throw io_error(
        "available number of file-handles too low; could not open any files");
    }
  }

  //! Indicates if a performance warning has been printed
  bool m_warning_printed = false;
  //! Most recently used managed_writer
  managed_writer* m_head = nullptr;
  //! Least recently used managed_writer
  managed_writer* m_tail = nullptr;
  //! Lock used to control access to internal state
  std::mutex m_lock{};
};

namespace {

//! List of (currently) inactive writers
writer_list g_writers;

} // namespace

///////////////////////////////////////////////////////////////////////////////

managed_reader::managed_reader(FILE* handle)
  : m_filename("<unknown file>")
  , m_file(handle)
{
  AR_REQUIRE(handle);
}

managed_reader::managed_reader(std::string filename)
  : m_filename(std::move(filename))
  , m_file(g_writers.fopen(m_filename, "rb"))
{
}

managed_reader::~managed_reader()
{
  if (::fclose(m_file)) {
    AR_FAIL(format_io_error("error closing " + log_escape(m_filename), errno));
  }
}

void
managed_reader::close()
{
  const auto ret = ::fclose(m_file);
  m_file = nullptr;

  if (ret) {
    throw io_error("error closing " + log_escape(m_filename), errno);
  }
}

size_t
managed_reader::read(void* buffer, size_t size)
{
  const auto nread = ::fread_unlocked(buffer, 1, size, m_file);
  if (ferror(m_file)) {
    throw io_error("error reading " + log_escape(m_filename), errno);
  }

  return nread;
}

///////////////////////////////////////////////////////////////////////////////

/** Locker  */
class writer_lock
{
public:
  explicit writer_lock(managed_writer* writer)
    : m_writer(writer)
  {
    AR_REQUIRE(m_writer);

    // Remove from global queue to prevent other threads from manipulating it
    g_writers.remove(m_writer);

    if (!m_writer->m_file) {
      m_writer->m_file = g_writers.fopen(m_writer->filename(),
                                         m_writer->m_created ? "ab" : "wb");
      m_writer->m_created = true;
    }
  };

  ~writer_lock()
  {
    if (m_writer->m_file) {
      // Allow this writer to be closed if we run out of file handles
      g_writers.add(m_writer);
    }
  }

  void write(const void* buffer, size_t size)
  {
    AR_REQUIRE(m_writer && m_writer->m_file);
    const auto ret = ::fwrite_unlocked(buffer, 1, size, m_writer->m_file);
    if (ret != size) {
      throw io_error("error writing to " + log_escape(m_writer->m_filename),
                     errno);
    }
  }

  void flush()
  {
    AR_REQUIRE(m_writer && m_writer->m_file);
    if (::fflush_unlocked(m_writer->m_file)) {
      throw io_error("error flushing file " + log_escape(m_writer->filename()),
                     errno);
    }
  }

  writer_lock(const writer_lock&) = delete;
  writer_lock(writer_lock&&) = delete;
  writer_lock& operator=(const writer_lock&) = delete;
  writer_lock& operator=(writer_lock&&) = delete;

private:
  managed_writer* m_writer;
};

///////////////////////////////////////////////////////////////////////////////

managed_writer::managed_writer(std::string filename)
  : m_filename(std::move(filename))
{
}

managed_writer::~managed_writer()
{
  AR_REQUIRE(!m_file);
}

void
managed_writer::write(const buffer& buf, const flush mode)
{
  if (buf.size() || mode == flush::on) {
    writer_lock writer{ this };

    writer.write(buf.get_signed(), buf.size());

    if (mode == flush::on) {
      writer.flush();
    }
  }
}

void
managed_writer::write(const buffer_vec& buffers, const flush mode)
{
  if (!buffers.empty() || mode == flush::on) {
    writer_lock writer{ this };

    for (const auto& buf : buffers) {
      writer.write(buf.get_signed(), buf.size());
    }

    if (mode == flush::on) {
      writer.flush();
    }
  }
}

void
managed_writer::write(const std::string& buf, const flush mode)
{
  if (!buf.empty() || mode == flush::on) {
    writer_lock writer{ this };

    writer.write(buf.data(), buf.length());

    if (mode == flush::on) {
      writer.flush();
    }
  }
}

void
managed_writer::close()
{
  g_writers.remove(this);
  if (m_file) {
    if (fclose(m_file)) {
      m_file = nullptr;
      throw io_error("failed to close file", errno);
    }

    m_file = nullptr;
  }
}

} // namespace adapterremoval
