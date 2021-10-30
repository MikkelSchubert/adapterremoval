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
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>

#include "linereader.hpp"
#include "managed_writer.hpp"
#include "threads.hpp"

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

std::string
format_gzip_msg(const std::string& message, const char* gzip_msg)
{
  if (gzip_msg) {
    std::stringstream stream;
    stream << message << " ('" << gzip_msg << "')";

    return stream.str();
  } else {
    return message;
  }
}

gzip_error::gzip_error(const std::string& message, const char* gzip_msg)
  : io_error(format_gzip_msg(message, gzip_msg))
{}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'line_reader'

line_reader::line_reader(const std::string& fpath)
  : m_filename(fpath)
  , m_file(managed_writer::gzopen(fpath, "rb"))
  , m_buffer(new char[BUF_SIZE])
  , m_buffer_ptr(m_buffer.get())
  , m_buffer_end(m_buffer.get())
  , m_eof(false)
{}

line_reader::~line_reader()
{
  auto errnum = gzclose(m_file);
  if (errnum != Z_OK) {
    std::string error;

    switch (errnum) {
      case Z_STREAM_ERROR:
        error = "~line_reader: invalid file";
        break;

      case Z_ERRNO:
        error = format_io_msg("~line_reader: error closing file", errno);
        break;

      case Z_MEM_ERROR:
        error = "~line_reader: insufficient memory";
        break;

      case Z_BUF_ERROR:
        error = "~line_reader: read ended in stream";
        break;

      default:
        error = "~line_reader: unknown error";
    }

    print_locker lock;
    std::cerr << "Error closing " << m_filename << ": " << error << std::endl;
    std::exit(1);
  }
}

bool
line_reader::getline(std::string& dst)
{
  dst.clear();

  while (!m_eof) {
    const char* start = m_buffer_ptr;
    char* end = m_buffer_ptr;

    for (; end != m_buffer_end; ++end) {
      if (*end == '\n') {
        // Excluding terminal \n
        dst.append(start, end - start);
        if (!dst.empty() && dst.back() == '\r') {
          // Excluding terminal \r; dst is examined, since the \r may
          // have been added separately, if the \r was the last
          // character in the previous buffer fill (see below).
          dst.pop_back();
        }

        m_buffer_ptr = end + 1;
        return true;
      }
    }

    // Can potentially introduce a \r; this is handled above.
    dst.append(start, end - start);
    refill_buffers();
  }

  return !dst.empty();
}

void
line_reader::refill_buffers()
{
  const int nread = gzread(m_file, m_buffer.get(), BUF_SIZE);
  if (nread == -1) {
    int errnum = 0;
    auto message = gzerror(m_file, &errnum);

    throw gzip_error("line_reader: error reading file", message);
  } else {
    // EOF set only once all data has been consumed
    m_eof = (nread == 0);

    m_buffer_ptr = m_buffer.get();
    m_buffer_end = m_buffer.get() + nread;
  }
}
