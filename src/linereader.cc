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
#include <iostream>

#include "linereader.h"
#include "threads.h"

const size_t BUF_SIZE = 4 * 1024;


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'io_error'

io_error::io_error(const std::string& message)
  : m_message(message)
{
}


io_error::~io_error() throw()
{
}


const char* io_error::what() const throw()
{
    return m_message.c_str();
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_error'

gzip_error::gzip_error(const std::string& message)
  : io_error(message)
{
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'line_reader'

line_reader::line_reader()
  : m_file(NULL)
  , m_buffer(new char[BUF_SIZE])
  , m_buffer_ptr(m_buffer + BUF_SIZE)
  , m_buffer_end(m_buffer + BUF_SIZE)
  , m_eof(false)
{
}



line_reader::line_reader(const std::string& fpath)
  : m_file(NULL)
  , m_buffer(new char[BUF_SIZE])
  , m_buffer_ptr(NULL)
  , m_buffer_end(NULL)
  , m_eof(false)
{
    open(fpath);
}


line_reader::~line_reader()
{
    delete[] m_buffer;
    m_buffer = NULL;

    try {
        close();
    } catch (const std::exception& error) {
        print_locker lock;
        std::cerr << "ERROR: " << error.what() << std::endl;
        std::exit(1);
    }
}


void line_reader::open(const std::string& fpath)
{
    close();

    errno = 0;
    m_file = gzopen(fpath.c_str(), "rb");
    if (!m_file) {
        if (errno) {
            throw io_error(std::string("line_reader::open: ") + std::strerror(errno));
        } else {
            throw gzip_error("line_reader::open: failed to open file");
        }
    }

    // Increase buffer size to significantly improve throughput
#if ((ZLIB_VER_MAJOR > 1) \
     || (ZLIB_VER_MAJOR == 1 \
         && ((ZLIB_VER_MINOR > 2) \
             || (ZLIB_VER_MINOR == 2 \
                 && ZLIB_VER_REVISION >= 4))))
    if (gzbuffer(m_file, 256 * 1024)) {
        throw gzip_error("could not set gzbuffer");
    }
#endif

    m_eof = false;
    m_buffer_end = m_buffer + BUF_SIZE;
    m_buffer_ptr = m_buffer_end;
}


bool line_reader::getline(std::string& dst)
{
    dst.clear();

    while (m_file && !m_eof) {
        char* start = m_buffer_ptr;
        char* end = m_buffer_ptr;

        for (; end != m_buffer_end; ++end) {
            if (*end == '\n') {
                // Excluding terminal \n
                dst.append(start, end - start);

                m_buffer_ptr = end + 1;
                return true;
            }
        }

        dst.append(start, end - start);
        refill_buffer();
    }

    return !dst.empty();
}


void line_reader::close()
{
    if (m_file) {
        const int errorcode = gzclose(m_file);
        switch (errorcode) {
            case Z_STREAM_ERROR:
                throw gzip_error("line_reader::close: attempted to close invalid file");

            case Z_ERRNO:
                throw io_error(std::string("line_reader::close: ") + std::strerror(errno));

            case Z_MEM_ERROR:
                throw gzip_error("line_reader::close: out of memory when closing file");

            case Z_BUF_ERROR:
                throw gzip_error("line_reader::close: file not completely read");

            case Z_OK:
                break;

            default:
                throw gzip_error("unknown error in line_reader::close");
        }

        m_file = NULL;
    }
}


bool line_reader::eof() const
{
    return m_eof;
}


bool line_reader::is_open() const
{
    return m_file;
}


void line_reader::refill_buffer()
{
    m_buffer_ptr = m_buffer;
    const int nread = gzread(m_file, m_buffer, BUF_SIZE);

    if (nread > 0) {
        m_buffer_end = m_buffer + nread;
    } else if (nread == 0) {
        m_eof = true;
    } else {
        int error_number = 0;
        const char* gzerror_msg = gzerror(m_file, &error_number);
        if (gzerror_msg) {
            throw gzip_error(gzerror_msg);
        }

        throw io_error(std::strerror(error_number));
    }
}
