/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <cerrno>   // for errno
#include <cstdio>   // for BUFSIZ
#include <cstdlib>  // for exit
#include <cstring>  // for strerror
#include <iostream> // for operator<<, basic_ostream, endl, cerr
#include <sstream>  // for stringstream

#include "debug.hpp" // for AR_DEBUG_ASSSERT
#include "linereader.hpp"
#include "managed_writer.hpp" // for managed_writer
#include "threads.hpp"        // for print_locker

//! Size of compressed and uncompressed buffers.
const int BUF_SIZE = 10 * BUFSIZ;

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'io_error'

std::string
format_io_msg(const std::string& message, int error_number)
{
  if (error_number) {
    std::stringstream stream;
    stream << message << " ('" << std::strerror(error_number) << "')";

    return stream.str();
  } else {
    return message;
  }
}

io_error::io_error(const std::string& message, int error_number)
  : std::ios_base::failure(format_io_msg(message, error_number))
{}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_error'

[[noreturn]] void
throw_gzip_error(const char* func, z_stream* zstream, const char* msg)
{
  std::stringstream stream;
  stream << func << ": " << msg;
  if (zstream && zstream->msg) {
    stream << " ('" << zstream->msg << "')";
  }

  throw gzip_error(stream.str());
}

#define THROW_GZIP_ERROR(msg)                                                  \
  throw_gzip_error(__func__, m_gzip_stream.get(), (msg))

gzip_error::gzip_error(const std::string& message)
  : io_error(message)
{}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'line_reader'

line_reader::line_reader(const std::string& fpath)
  : m_filename(fpath)
  , m_file(managed_writer::fopen(fpath, "rb"))
  , m_gzip_stream(nullptr)
#if defined(USE_LIBISAL)
  , m_gzip_header(nullptr)
#endif
  , m_buffer(nullptr)
  , m_buffer_ptr(nullptr)
  , m_buffer_end(nullptr)
  , m_raw_buffer(new char[BUF_SIZE])
  , m_raw_buffer_end(m_raw_buffer + BUF_SIZE)
  , m_eof(false)
{
  if (!m_file) {
    throw io_error("line_reader::open: failed to open file", errno);
  }
}

line_reader::~line_reader()
{
  try {
    close_buffers_gzip();

    delete[] m_raw_buffer;
    m_raw_buffer = nullptr;

    if (fclose(m_file)) {
      throw io_error("line_reader::close: error closing file", errno);
    }
  } catch (const std::exception& error) {
    print_locker lock;
    std::cerr << "Error closing file: " << error.what() << std::endl;
    std::exit(1);
  }
}

bool
line_reader::getline(std::string& dst)
{
  dst.clear();

  while (!m_eof) {
    const size_t length = m_buffer_end - m_buffer_ptr;
    char* ptr = reinterpret_cast<char*>(memchr(m_buffer_ptr, '\n', length));
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
  m_buffer_ptr = m_raw_buffer;
  m_buffer_end = m_raw_buffer_end;
}

void
line_reader::refill_raw_buffer()
{
  const int nread = fread(m_raw_buffer, 1, BUF_SIZE, m_file);

  if (nread == BUF_SIZE) {
    m_raw_buffer_end = m_raw_buffer + BUF_SIZE;
  } else if (ferror(m_file)) {
    throw io_error("line_reader::refill_buffer: error reading file", errno);
  } else {
    // EOF set only once all data has been consumed
    m_eof = (nread == 0);
    m_raw_buffer_end = m_raw_buffer + nread;
  }
}

bool
line_reader::identify_gzip() const
{
  if (m_raw_buffer_end - m_raw_buffer < 2) {
    return false;
  } else if (m_raw_buffer[0] != '\x1f' || m_raw_buffer[1] != '\x8b') {
    return false;
  }

  return true;
}

void
line_reader::initialize_buffers_gzip()
{
  std::cerr << "initialize_buffers_gzip" << std::endl;

  m_buffer = new char[BUF_SIZE];
  m_buffer_ptr = m_buffer + BUF_SIZE;
  m_buffer_end = m_buffer + BUF_SIZE;

#if defined(USE_LIBISAL)
  m_gzip_stream.reset(new inflate_state());
  m_gzip_header.reset(new isal_gzip_header());

  isal_inflate_init(m_gzip_stream.get());
  m_gzip_stream->crc_flag = ISAL_GZIP_NO_HDR_VER;
  m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
  m_gzip_stream->next_in = reinterpret_cast<uint8_t*>(m_raw_buffer);

  const auto ret =
    isal_read_gzip_header(m_gzip_stream.get(), m_gzip_header.get());
#warning TODO: handle isal_read_gzip_header return value
  AR_DEBUG_ASSERT(ret == ISAL_DECOMP_OK);
#else
  m_gzip_stream.reset(new z_stream());
  m_gzip_stream->zalloc = nullptr;
  m_gzip_stream->zfree = nullptr;
  m_gzip_stream->opaque = nullptr;

  m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
  m_gzip_stream->next_in = reinterpret_cast<Bytef*>(m_raw_buffer);

  switch (inflateInit2(m_gzip_stream.get(), 15 + 16)) {
    case Z_OK:
      break;

    case Z_MEM_ERROR:
      THROW_GZIP_ERROR("insufficient memory");

    case Z_VERSION_ERROR:
      THROW_GZIP_ERROR("incompatible zlib version");

    case Z_STREAM_ERROR:
      THROW_GZIP_ERROR("invalid parameters");

    default:
      THROW_GZIP_ERROR("unknown error");
  }
#endif
}

void
line_reader::refill_buffers_gzip()
{
#if defined(USE_LIBISAL)
  if (!m_gzip_stream->avail_in) {
    refill_raw_buffer();
    m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
    m_gzip_stream->next_in = reinterpret_cast<uint8_t*>(m_raw_buffer);
  }

  m_gzip_stream->avail_out = BUF_SIZE;
  m_gzip_stream->next_out = reinterpret_cast<uint8_t*>(m_buffer);

  const auto ret = isal_inflate(m_gzip_stream.get());
#warning TODO: handle isal_inflate return value
  AR_DEBUG_ASSERT(ret == ISAL_DECOMP_OK);

  switch (m_gzip_stream->block_state) {
#warning TODO: handle isal_inflate block states
    case isal_block_state::ISAL_BLOCK_FINISH:
      isal_inflate_reset(m_gzip_stream.get());
  }

  m_buffer_ptr = m_buffer;
  m_buffer_end = m_buffer + (BUF_SIZE - m_gzip_stream->avail_out);
#else
  if (!m_gzip_stream->avail_in) {
    refill_raw_buffer();
    m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
    m_gzip_stream->next_in = reinterpret_cast<Bytef*>(m_raw_buffer);
  }

  m_gzip_stream->avail_out = BUF_SIZE;
  m_gzip_stream->next_out = reinterpret_cast<Bytef*>(m_buffer);
  switch (inflate(m_gzip_stream.get(), Z_NO_FLUSH)) {
    case Z_OK:
    case Z_BUF_ERROR: /* input buffer empty or output buffer full */
      break;

    case Z_STREAM_END:
      // Handle concatenated streams; causes unnecessary reset at EOF
      if (inflateReset(m_gzip_stream.get()) != Z_OK) {
        THROW_GZIP_ERROR("failed to reset stream");
      }
      break;

    case Z_STREAM_ERROR:
      THROW_GZIP_ERROR("inconsistent stream state");

    default:
      THROW_GZIP_ERROR("unknown error");
  }

  m_buffer_ptr = m_buffer;
  m_buffer_end = m_buffer + (BUF_SIZE - m_gzip_stream->avail_out);
#endif
}

void
line_reader::close_buffers_gzip()
{
  if (m_gzip_stream) {
#if defined(USE_LIBISAL)
#warning TODO
    m_gzip_stream.reset();
    m_gzip_header.reset();
#else
    switch (inflateEnd(m_gzip_stream.get())) {
      case Z_OK:
        break;

      case Z_STREAM_ERROR:
        THROW_GZIP_ERROR("stream error");

      default:
        THROW_GZIP_ERROR("unknown error");
    }

    m_gzip_stream.reset();
#endif

    delete[] m_buffer;
    m_buffer = nullptr;
  }
}
