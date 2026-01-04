// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "linereader.hpp"    // declarations
#include "errors.hpp"        // for io_error
#include "logging.hpp"       // for warn, log_stream
#include "managed_io.hpp"    // for managed_writer
#include "strutils.hpp"      // for shell_escape
#include <cstdint>           // for uint8_t
#include <cstring>           // for strerror, memchr
#include <isa-l/igzip_lib.h> // for inflate_state, isal_gzip_header, isal_...
#include <memory>            // for unique_ptr, shared_ptr, __shared_ptr_a...
#include <sstream>           // for operator<<, basic_ostream
#include <utility>           // for move

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Helper functions

namespace {

[[noreturn]] void
throw_gzip_error(const std::string_view filename,
                 const std::string_view action,
                 const std::string_view error,
                 const std::string_view diagnosis = "file is likely corrupt")
{
  std::ostringstream stream;
  stream << "Error while " << action << " " << shell_escape(filename) << ": "
         << error << "; " << diagnosis;

  throw gzip_error(stream.str());
}

void
check_isal_return_code(int returncode,
                       const std::string_view file,
                       const std::string_view action)
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
      throw_gzip_error(file, action, "end of gzip comment buffer reached");

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

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'vec_reader'

vec_reader::vec_reader(std::initializer_list<std::string_view> lines)
  : m_lines(lines.begin(), lines.end())
  , m_it(m_lines.begin())
{
}

vec_reader::vec_reader(string_vec lines)
  : m_lines(std::move(lines))
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
  : m_reader(handle)
  , m_gzip_stream(nullptr)
  , m_gzip_header(nullptr)
  , m_buffer(nullptr)
  , m_buffer_ptr(nullptr)
  , m_buffer_end(nullptr)
  , m_raw_buffer(std::make_shared<line_buffer>())
  , m_raw_buffer_end(m_raw_buffer->begin())
  , m_eof(false)
{
  refill_buffers();
}

line_reader::line_reader(std::string filename)
  : m_reader(std::move(filename))
  , m_gzip_stream(nullptr)
  , m_gzip_header(nullptr)
  , m_buffer(nullptr)
  , m_buffer_ptr(nullptr)
  , m_buffer_end(nullptr)
  , m_raw_buffer(std::make_shared<line_buffer>())
  , m_raw_buffer_end(m_raw_buffer->begin())
  , m_eof(false)
{
  refill_buffers();
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

    if (is_raw_buffer_gzip()) {
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
line_reader::refill_raw_buffer(size_t avail_in)
{
  if (avail_in) {
    // Move unused (compressed) data to the front of the buffer
    std::memmove(m_raw_buffer->data(), m_raw_buffer_end - avail_in, avail_in);
  }

  const size_t nread = m_reader.read(m_raw_buffer->data() + avail_in,
                                     m_raw_buffer->size() - avail_in);

  // EOF set only once all data has been consumed
  m_eof = (nread + avail_in == 0);
  m_raw_buffer_end = m_raw_buffer->data() + nread + avail_in;
}

bool
line_reader::is_raw_buffer_gzip() const
{
  return m_raw_buffer_end - m_raw_buffer->data() > 1 &&
         m_raw_buffer->at(0) == '\x1f' && m_raw_buffer->at(1) == '\x8b';
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
  check_isal_return_code(result,
                         m_reader.filename(),
                         "reading first gzip header from");
}

void
line_reader::refill_buffers_gzip()
{
  m_gzip_stream->avail_out = m_buffer->size();
  m_gzip_stream->next_out = reinterpret_cast<uint8_t*>(m_buffer->data());

  // Refill the buffer if empty or if we need more data to properlyidentify
  // additional gzip blocks and parse their headers. The number of bytes (64) is
  // arbitrary, but should suffice.
  if (!m_gzip_stream->avail_in ||
      (m_gzip_stream->avail_in < 64 &&
       m_gzip_stream->block_state == isal_block_state::ISAL_BLOCK_FINISH)) {
    refill_raw_buffer(m_gzip_stream->avail_in);
    m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer->data();
    m_gzip_stream->next_in = reinterpret_cast<uint8_t*>(m_raw_buffer->data());
  }

  if (m_gzip_stream->block_state == isal_block_state::ISAL_BLOCK_FINISH) {
    if (m_gzip_stream->avail_in > 1 && m_gzip_stream->next_in[0] == 0x1f &&
        m_gzip_stream->next_in[1] == 0x8b) {
      isal_inflate_reset(m_gzip_stream.get());
      // simplify by letting isa-l handle partly buffered headers
      m_gzip_stream->crc_flag = ISAL_GZIP;
    } else if (m_gzip_stream->avail_in) {
      log::warn() << "Ignoring trailing garbage at the end of "
                  << shell_escape(m_reader.filename());

      m_buffer_ptr = m_buffer->data();
      m_buffer_end = m_buffer_ptr;
      m_eof = true;
      return;
    }
  }

  check_isal_return_code(isal_inflate(m_gzip_stream.get()),
                         m_reader.filename(),
                         "decompressing");

  m_buffer_ptr = m_buffer->data();
  m_buffer_end = m_buffer_ptr + (m_buffer->size() - m_gzip_stream->avail_out);

  if (m_eof && !m_gzip_stream->avail_in &&
      m_gzip_stream->block_state != isal_block_state::ISAL_BLOCK_FINISH) {
    throw_gzip_error(m_reader.filename(),
                     "decompressing",
                     "unexpected end of file",
                     "file is likely truncated!");
  }
}

} // namespace adapterremoval
