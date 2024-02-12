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
#include "linereader.hpp"
#include "debug.hpp"          // for AR_FAIL
#include "logging.hpp"        // for warn, log_stream
#include "managed_writer.hpp" // for managed_writer
#include "strutils.hpp"       // for shell_escape
#include <cerrno>             // for errno
#include <cstring>            // for strerror, memchr
#include <isa-l/igzip_lib.h>  // for inflate_state, isal_gzip_header, isal_...
#include <memory>             // for unique_ptr, shared_ptr, __shared_ptr_a...
#include <sstream>            // for operator<<, basic_ostream
#include <stdint.h>           // for uint8_t

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'io_error'

std::string
format_io_msg(const std::string& message, int error_number)
{
  if (error_number) {
    std::ostringstream stream;
    stream << message << " ('" << std::strerror(error_number) << "')";

    return stream.str();
  } else {
    return message;
  }
}

io_error::io_error(const std::string& message, int error_number)
  : std::ios_base::failure(format_io_msg(message, error_number))
  , m_what(format_io_msg(message, error_number))
{
}

const char*
io_error::what() const noexcept
{
  return m_what.c_str();
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_error'

[[noreturn]] void
throw_gzip_error(const std::string& filename,
                 const char* action,
                 const char* error,
                 const char* diagnosis = "file is likely corrupt")
{
  std::ostringstream stream;
  stream << "Error while " << action << " " << shell_escape(filename) << ": "
         << error << "; " << diagnosis;

  throw gzip_error(stream.str());
}

gzip_error::gzip_error(const std::string& message)
  : io_error(message)
{
}

///////////////////////////////////////////////////////////////////////////////
// Helper functions for isa-l

void
check_isal_return_code(int returncode,
                       const std::string& file,
                       const char* action)
{
  switch (returncode) {
    case ISAL_DECOMP_OK:
      return;

    case ISAL_END_INPUT:
      throw_gzip_error(file, action, "end of input reached");

    case ISAL_OUT_OVERFLOW:
      throw_gzip_error(file, action, "end of output reached");

    case ISAL_NAME_OVERFLOW:
      throw_gzip_error(file, action, "end of gzip name buffer reached");

    case ISAL_COMMENT_OVERFLOW:
      throw_gzip_error(file, action, "end of gzip name buffer reached");

    case ISAL_EXTRA_OVERFLOW:
      throw_gzip_error(file, action, "end of extra buffer reached");

    case ISAL_NEED_DICT:
      throw_gzip_error(file, action, "stream needs dictionary to continue");

    case ISAL_INVALID_BLOCK:
      throw_gzip_error(file, action, "invalid deflate block found");

    case ISAL_INVALID_SYMBOL:
      throw_gzip_error(file, action, "invalid deflate symbol found");

    case ISAL_INVALID_LOOKBACK:
      throw_gzip_error(file, action, "invalid lookback distance found");

    case ISAL_INVALID_WRAPPER:
      throw_gzip_error(file, action, "invalid gzip/zlib wrapper found");

    case ISAL_UNSUPPORTED_METHOD:
      throw_gzip_error(file, action, "unsupported compression method");

    case ISAL_INCORRECT_CHECKSUM:
      throw_gzip_error(file, action, "incorrect checksum found");

    default:
      throw_gzip_error(file, action, "unknown error");
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'vec_reader'

vec_reader::vec_reader(const string_vec& lines)
  : m_lines(lines)
  , m_it(m_lines.begin())
{
}

bool
vec_reader::getline(std::string& dst)
{
  if (m_it == m_lines.end()) {
    return false;
  }

  dst = *m_it++;
  return true;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'line_reader'

line_reader::line_reader(FILE* handle)
  : m_filename("<unnamed file>")
  , m_file(handle)
  , m_gzip_stream(nullptr)
  , m_gzip_header(nullptr)
  , m_buffer(nullptr)
  , m_buffer_ptr(nullptr)
  , m_buffer_end(nullptr)
  , m_raw_buffer(std::make_shared<line_buffer>())
  , m_raw_buffer_end(m_raw_buffer->end())
  , m_eof(false)
{
  if (!m_file) {
    throw io_error("could not open file", errno);
  }

  refill_buffers();
}

line_reader::line_reader(const std::string& fpath)
  : line_reader(managed_writer::fopen(fpath, "rb"))
{
  m_filename = fpath;
}

line_reader::~line_reader()
{
  // Pending input is ignored
  if (fclose(m_file)) {
    AR_FAIL(format_io_msg("error closing input file", errno));
  }
}

bool
line_reader::getline(std::string& dst)
{
  dst.clear();

  while (!m_eof) {
    const size_t length = m_buffer_end - m_buffer_ptr;
    auto* ptr = static_cast<char*>(memchr(m_buffer_ptr, '\n', length));
    if (ptr) {
      // Excluding terminal \n
      dst.append(m_buffer_ptr, ptr - m_buffer_ptr);
      if (!dst.empty() && dst.back() == '\r') {
        // Excluding terminal \r; this may have been added in the loop prior to
        // \n being found, if \n is the first character in the buffer. It is
        // therefore easiest to check dst directly.
        dst.pop_back();
      }

      m_buffer_ptr = ptr + 1;
      return true;
    }

    // Can potentially introduce a \r; this is handled above.
    dst.append(m_buffer_ptr, length);
    refill_buffers();
  }

  return !dst.empty();
}

void
line_reader::refill_buffers()
{
  if (m_buffer) {
    if (m_gzip_stream) {
      refill_buffers_gzip();
    } else {
      refill_raw_buffer();
      refill_buffers_uncompressed();
    }
  } else {
    refill_raw_buffer();

    if (identify_gzip()) {
      initialize_buffers_gzip();
    } else {
      refill_buffers_uncompressed();
    }
  }
}

void
line_reader::refill_buffers_uncompressed()
{
  m_buffer = m_raw_buffer;
  m_buffer_ptr = m_raw_buffer->data();
  m_buffer_end = m_raw_buffer_end;
}

void
line_reader::refill_raw_buffer()
{
  const int nread =
    fread(m_raw_buffer->data(), 1, m_raw_buffer->size(), m_file);

  if (ferror(m_file)) {
    throw io_error("read error while filling buffer", errno);
  } else {
    // EOF set only once all data has been consumed
    m_eof = (nread == 0);
    m_raw_buffer_end = m_raw_buffer->data() + nread;
  }
}

bool
line_reader::identify_gzip() const
{
  if (m_raw_buffer_end - m_raw_buffer->data() < 2) {
    return false;
  } else if (m_raw_buffer->at(0) != '\x1f' || m_raw_buffer->at(1) != '\x8b') {
    return false;
  }

  return true;
}

void
line_reader::initialize_buffers_gzip()
{
  m_buffer = std::make_shared<line_buffer>();
  m_buffer_ptr = m_buffer->end();
  m_buffer_end = m_buffer->end();

  m_gzip_stream = std::make_unique<inflate_state>();
  m_gzip_header = std::make_unique<isal_gzip_header>();

  isal_inflate_init(m_gzip_stream.get());
  m_gzip_stream->crc_flag = ISAL_GZIP_NO_HDR_VER;
  m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer->data();
  m_gzip_stream->next_in = reinterpret_cast<uint8_t*>(m_raw_buffer->data());

  isal_gzip_header_init(m_gzip_header.get());
  auto result = isal_read_gzip_header(m_gzip_stream.get(), m_gzip_header.get());
  check_isal_return_code(result, m_filename, "reading gzip header from");
}

void
line_reader::refill_buffers_gzip()
{
  if (!m_gzip_stream->avail_in) {
    refill_raw_buffer();
    m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer->data();
    m_gzip_stream->next_in = reinterpret_cast<uint8_t*>(m_raw_buffer->data());
  }

  m_gzip_stream->avail_out = m_buffer->size();
  m_gzip_stream->next_out = reinterpret_cast<uint8_t*>(m_buffer->data());

  if (m_gzip_stream->avail_in &&
      m_gzip_stream->block_state == isal_block_state::ISAL_BLOCK_FINISH) {
    if (check_next_gzip_block()) {
      isal_inflate_reset(m_gzip_stream.get());
      // Handle headers in subsequent blocks
      m_gzip_stream->crc_flag = ISAL_GZIP;
    } else {
      return;
    }
  }

  check_isal_return_code(
    isal_inflate(m_gzip_stream.get()), m_filename, "decompressing");

  m_buffer_ptr = m_buffer->data();
  m_buffer_end = m_buffer_ptr + (m_buffer->size() - m_gzip_stream->avail_out);

  if (m_eof && !m_gzip_stream->avail_in &&
      m_gzip_stream->block_state != isal_block_state::ISAL_BLOCK_FINISH) {
    throw_gzip_error(m_filename,
                     "decompressing",
                     "unexpected end of file",
                     "file is likely truncated!");
  }
}

bool
line_reader::check_next_gzip_block()
{
  // This matches the behavior of zcat / gzip / igzip
  if ((m_gzip_stream->avail_in > 0 && m_gzip_stream->next_in[0] != 0x1f) ||
      (m_gzip_stream->avail_in > 1 && m_gzip_stream->next_in[1] != 0x8b)) {
    log::warn() << "Ignoring trailing garbage at the end of "
                << shell_escape(m_filename);

    m_buffer_ptr = m_buffer->data();
    m_buffer_end = m_buffer_ptr;
    m_eof = true;

    return false;
  }

  return true;
}

} // namespace adapterremoval
