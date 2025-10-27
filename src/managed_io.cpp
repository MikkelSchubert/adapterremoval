// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "managed_io.hpp"  // declarations
#include "buffer.hpp"      // for buffer
#include "commontypes.hpp" // for DEV_STDOUT, DEV_STDERR, ...
#include "debug.hpp"       // for AR_REQUIRE
#include "errors.hpp"      // for io_error
#include "logging.hpp"     // for log::warn
#include "strutils.hpp"    // for log_escape
#include <algorithm>       // for any_of
#include <cerrno>          // for EMFILE, errno
#include <cstddef>         // for size_t
#include <cstdio>          // for fopen, fread, fwrite, ...
#include <exception>       // for std::exception
#include <fcntl.h>         // for posix_fadvise
#include <mutex>           // for mutex, lock_guard
#include <string>          // for string
#include <string_view>     // for string_view
#include <sys/stat.h>      // for fstat
#include <utility>         // for move
#include <vector>          // for vector

namespace adapterremoval {

class io_manager
{
public:
  static void open(managed_reader* reader)
  {
    if (reader->m_file) {
      return;
    } else if (reader->filename() == DEV_STDIN) {
      reader->m_file = stdin;
    } else if (reader->filename() != DEV_PIPE) {
      reader->m_file = io_manager::fopen(reader->filename(), "rb");

#if (defined(_XOPEN_SOURCE) && _XOPEN_SOURCE >= 600) ||                        \
  (defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L)
      posix_fadvise(fileno(reader->m_file), 0, 0, POSIX_FADV_WILLNEED);
      posix_fadvise(fileno(reader->m_file), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    } else {
      // Merged I/O depends on filenames being identical
      AR_FAIL("unhandled STDIN marker");
    }
  }

  static void open(managed_writer* writer)
  {
    if (writer->m_file) {
      return;
    } else if (writer->filename() == DEV_STDOUT) {
      writer->m_file = stdout;
    } else if (writer->filename() == DEV_STDERR) {
      // Not sure why anyone would do this, but ¯\_(ツ)_/¯
      writer->m_file = stderr;
    } else if (writer->filename() != DEV_PIPE) {
      writer->m_file =
        io_manager::fopen(writer->filename(), writer->m_created ? "ab" : "wb");
    } else {
      // Merged I/O depends on filenames being identical
      AR_FAIL("unhandled STDOUT marker");
    }

    writer->m_created = true;
    writer->m_stream =
      io_manager::is_stream(writer->filename(), writer->m_file);
  }

  /** Adds writer to list of inactive writers */
  static void add(managed_writer* writer)
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
  static void remove(managed_writer* writer)
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

  io_manager() = delete;
  ~io_manager() = delete;
  io_manager(const io_manager&) = delete;
  io_manager(io_manager&&) = delete;
  io_manager& operator=(const io_manager&) = delete;
  io_manager& operator=(io_manager&&) = delete;

private:
  static FILE* fopen(const std::string& filename, const char* mode)
  {
    while (true) {
      FILE* handle = ::fopen(filename.c_str(), mode);

      if (handle) {
        AR_REQUIRE(!::ferror(handle));
        return handle;
      } else if (errno == EMFILE) {
        close_one();
      } else {
        throw io_error("failed to open file", errno);
      }
    }
  }

  static bool is_stream(const std::string& filename, FILE* handle)
  {
    if (handle == stdin || handle == stdout || handle == stderr) {
      return true;
    }

    struct stat statbuf = {};
    if (fstat(fileno(handle), &statbuf) == 0) {
      return S_ISFIFO(statbuf.st_mode);
    }

    log::warn() << "Could not fstat " << log_escape(filename);

    // Assumed to be a stream to be safe
    return true;
  }

  /** Try to close the least recently used writer */
  static void close_one()
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
  static bool m_warning_printed;
  //! Most recently used managed_writer
  static managed_writer* m_head;
  //! Least recently used managed_writer
  static managed_writer* m_tail;
  //! Lock used to control access to internal state
  static std::mutex m_lock;
};

bool io_manager::m_warning_printed = false;
managed_writer* io_manager::m_head = nullptr;
managed_writer* io_manager::m_tail = nullptr;
std::mutex io_manager::m_lock{};

///////////////////////////////////////////////////////////////////////////////

managed_reader::managed_reader(FILE* handle)
  : m_filename("<unknown file>")
  , m_file(handle)
{
  AR_REQUIRE(handle);
}

managed_reader::managed_reader(std::string filename)
  : m_filename(std::move(filename))
{
  io_manager::open(this);
}

managed_reader::~managed_reader()
{
  if (m_file && ::fclose(m_file) != 0) {
    AR_FAIL(format_io_error("error closing " + log_escape(m_filename), errno));
  }
}

void
managed_reader::close()
{
  if (m_file) {
    if (::fclose(m_file) != 0) {
      m_file = nullptr;
      throw io_error("error closing " + log_escape(m_filename), errno);
    }

    m_file = nullptr;
  }
}

size_t
managed_reader::read(void* buffer, size_t size)
{
  AR_REQUIRE(buffer);
  const auto nread = ::fread(buffer, 1, size, m_file);
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
    io_manager::remove(m_writer);
    io_manager::open(m_writer);
  };

  ~writer_lock()
  {
    // Streams cannot be managed, since they cannot be reopened
    if (m_writer->m_file && !m_writer->m_stream) {
      // Allow this writer to be closed if we run out of file handles
      io_manager::add(m_writer);
    }
  }

  void write(const void* buffer, size_t size)
  {
    AR_REQUIRE(buffer || size == 0);
    AR_REQUIRE(m_writer && m_writer->m_file);
    if (size) {
      const auto ret = ::fwrite(buffer, 1, size, m_writer->m_file);
      if (ret != size) {
        throw io_error("error writing to " + log_escape(m_writer->m_filename),
                       errno);
      }
    }
  }

  void flush()
  {
    AR_REQUIRE(m_writer && m_writer->m_file);
    if (::fflush(m_writer->m_file)) {
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

namespace {

bool
any_nonempty_buffers(const std::vector<buffer>& buffers)
{
  return std::any_of(buffers.begin(), buffers.end(), [](const auto& it) {
    return it.size() != 0;
  });
}

} // namespace

managed_writer::managed_writer(std::string filename)
  : m_filename(std::move(filename))
{
}

managed_writer::~managed_writer()
{
  try {
    close();
  } catch (const std::exception&) {
    AR_FAIL("unhandled exception");
  }
}

void
managed_writer::write(const buffer& buf, const flush mode)
{
  if (buf.size() || mode == flush::on) {
    writer_lock writer{ this };

    writer.write(buf.data(), buf.size());

    if (mode == flush::on) {
      writer.flush();
    }
  }
}

void
managed_writer::write(const std::vector<buffer>& buffers, const flush mode)
{
  if (mode == flush::on || any_nonempty_buffers(buffers)) {
    writer_lock writer{ this };

    for (const auto& buf : buffers) {
      writer.write(buf.data(), buf.size());
    }

    if (mode == flush::on) {
      writer.flush();
    }
  }
}

void
managed_writer::write(std::string_view buf, const flush mode)
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
  io_manager::remove(this);
  if (m_file) {
    if (fclose(m_file)) {
      m_file = nullptr;
      throw io_error("failed to close file", errno);
    }

    m_file = nullptr;
  }
}

} // namespace adapterremoval
