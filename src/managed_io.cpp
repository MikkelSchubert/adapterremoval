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
#include <fcntl.h>         // for posix_fadvise
#include <mutex>           // for mutex, unique_lock
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
    // Merged I/O depends on filenames being identical
    AR_REQUIRE(reader->filename() != DEV_PIPE, "unhandled STDIN marker");

    if (reader->m_file) {
      return;
    } else if (reader->filename() == DEV_STDIN) {
      reader->m_file = stdin;
    } else {
      reader->m_file = io_manager::fopen(reader->filename(), "rb");

#if (defined(_XOPEN_SOURCE) && _XOPEN_SOURCE >= 600) ||                        \
  (defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L)
      ::posix_fadvise(::fileno(reader->m_file), 0, 0, POSIX_FADV_WILLNEED);
      ::posix_fadvise(::fileno(reader->m_file), 0, 0, POSIX_FADV_SEQUENTIAL);
      ::posix_fadvise(::fileno(reader->m_file), 0, 0, POSIX_FADV_NOREUSE);
#endif
    }
  }

  static void open(managed_writer* writer, bool lazy = false)
  {
    // Merged I/O depends on filenames being identical
    AR_REQUIRE(writer->filename() != DEV_PIPE, "unhandled STDOUT marker");

    if (writer->m_file) {
      return;
    } else if (writer->filename() == DEV_STDOUT) {
      writer->m_file = stdout;
      writer->m_state = managed_writer::state::streaming;

#ifdef _WIN32
      ::_setmode(::fileno(stdout), O_BINARY);
#endif
    } else if (writer->filename() == DEV_STDERR) {
      // Not sure why anyone would do this, but ¯\_(ツ)_/¯
      writer->m_file = stderr;
      writer->m_state = managed_writer::state::streaming;

#ifdef _WIN32
      ::_setmode(::fileno(stderr), O_BINARY);
#endif
    } else if (!lazy) {
      try {
        if (writer->m_state == managed_writer::state::uninitialized) {
          writer->m_file = io_manager::fopen(writer->filename(), "wb");
        } else if (writer->m_state == managed_writer::state::writing) {
          writer->m_file = io_manager::fopen(writer->filename(), "ab");
        } else {
          AR_FAIL("invalid managed_writer::state");
        }
      } catch (const io_error&) {
        writer->m_state = managed_writer::state::failed;
        throw;
      }

      if (writer->m_state == managed_writer::state::uninitialized) {
        // Assuming it's a stream is harmless unless we run out of file handles
        writer->m_state = managed_writer::state::streaming;

        struct stat statbuf{};
        if (::fstat(::fileno(writer->m_file), &statbuf) == 0) {
          if (S_ISREG(statbuf.st_mode)) {
            // File can be closed and re-opened to if we run out of file handles
            writer->m_state = managed_writer::state::writing;
          }
        } else {
          log::warn() << "Could not fstat " << log_escape(writer->filename());
        }
      }
    }
  }

  /** Adds writer to list of inactive writers */
  static void add(managed_writer* writer)
  {
    const std::unique_lock lock{ m_lock };
    AR_REQUIRE(writer);
    AR_REQUIRE(writer->m_state == managed_writer::state::writing);

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
    const std::unique_lock lock{ m_lock };

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
        return handle;
      } else if (errno == EMFILE) {
        close_one();
      } else {
        throw io_error("failed to open file", errno);
      }
    }
  }

  /** Try to close the least recently used writer */
  static void close_one()
  {
    const std::unique_lock lock{ m_lock };

    AR_REQUIRE(!m_head == !m_tail);
    if (!m_warning_printed) {
      log::warn() << "Number of available file-handles (ulimit -n) is too low. "
                  << "AdapterRemoval will dynamically close/re-open files as "
                  << "required, but performance may suffer as a result.";

      m_warning_printed = true;
    }

    auto* writer = m_tail;
    if (writer) {
      remove(writer);

      if (::fclose(writer->m_file) != 0) {
        writer->m_file = nullptr;
        writer->m_state = managed_writer::state::failed;
        throw io_error("failed to close file", errno);
      }
      writer->m_file = nullptr;
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
  static std::recursive_mutex m_lock;
};

bool io_manager::m_warning_printed = false;
managed_writer* io_manager::m_head = nullptr;
managed_writer* io_manager::m_tail = nullptr;
std::recursive_mutex io_manager::m_lock{};

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
    const auto error_number = errno;

    // this is a non-fatal error if there were no problems reading the file
    log::warn() << format_io_error("error closing " + log_escape(m_filename),
                                   error_number);
  }
}

void
managed_reader::close()
{
  // This will also close non-regular files, such as STDIN. That is intentional
  // as only a single reader should be created per input file or stream
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
  if (::ferror(m_file) != 0) {
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

    switch (m_writer->m_state) {
      case managed_writer::state::writing:
        // Remove from global queue to prevent other threads from touching it
        io_manager::remove(m_writer);
        [[fallthrough]];
      case managed_writer::state::uninitialized:
        io_manager::open(m_writer);
        [[fallthrough]];
      case managed_writer::state::streaming:
        return; // streams are always open and aren't placed on the queue
      case managed_writer::state::finalized:
        AR_FAIL("attempted to re-open finalized writer");
      case managed_writer::state::failed:
        AR_FAIL("attempted to re-open failed writer");
      default:
        AR_FAIL("invalid managed_writer::state");
    }
  }

  ~writer_lock()
  {
    // Allow regular files to be closed if we run out of file handles. Streams
    // cannot be re-opened, and must therefore not be closed once opened
    if (m_writer->m_state == managed_writer::state::writing) {
      AR_REQUIRE(m_writer->m_file);
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
        m_writer->m_state = managed_writer::state::failed;

        throw io_error("error writing to " + log_escape(m_writer->m_filename),
                       errno);
      }
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
  // Identify standard streams but don't open actual files, since that could
  // cause churn when demultiplexing a large number of samples
  io_manager::open(this, true);
}

managed_writer::~managed_writer()
{
  io_manager::remove(this);

  if (m_file && ::fclose(m_file) != 0) {
    // AdapterRemoval should already have aborted before getting to this point,
    // as handles should normally be closed, and pipeline errors triggers abort
    AR_FAIL(format_io_error("error closing " + log_escape(m_filename), errno));
  }
}

void
managed_writer::write(const buffer& buf)
{
  if (buf.size()) {
    writer_lock writer{ this };

    writer.write(buf.data(), buf.size());
  }
}

void
managed_writer::write(const std::vector<buffer>& buffers)
{
  if (any_nonempty_buffers(buffers)) {
    writer_lock writer{ this };

    for (const auto& buf : buffers) {
      writer.write(buf.data(), buf.size());
    }
  }
}

void
managed_writer::write(std::string_view buf)
{
  if (!buf.empty()) {
    writer_lock writer{ this };

    writer.write(buf.data(), buf.length());
  }
}

void
managed_writer::finalize()
{
  AR_REQUIRE(m_state != state::failed);
  if (m_state == managed_writer::state::writing) {
    io_manager::remove(this);
  } else if (m_state == managed_writer::state::uninitialized) {
    // Ensure creation of the file, if nothing was written to it
    io_manager::open(this);
  }

  // This will also close non-regular files, such as STDOUT. That is intentional
  // as only a single writer should be created per output file or stream
  if (m_file) {
    if (::fclose(m_file) != 0) {
      m_file = nullptr;
      m_state = managed_writer::state::failed;

      throw io_error("failed to close file", errno);
    }

    m_file = nullptr;
  }

  m_state = state::finalized;
}

} // namespace adapterremoval
