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
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "linereader.h"
#include "threads.h"

namespace ar
{

//! Size of compressed and uncompressed buffers.
const int BUF_SIZE = 10 * BUFSIZ;


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'io_error'

std::string format_io_msg(const std::string& message, int error_number)
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
{
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_error'

std::string format_gzip_msg(const std::string& message, const char* gzip_msg)
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
{
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_error'

bzip2_error::bzip2_error(const std::string& message)
  : io_error(message)
{
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'line_reader'

line_reader::line_reader(const std::string& fpath)
  : m_file(fopen(fpath.c_str(), "rb"))
#ifdef AR_GZIP_SUPPORT
  , m_gzip_stream(NULL)
#endif
#ifdef AR_BZIP2_SUPPORT
  , m_bzip2_stream(NULL)
#endif
  , m_buffer(NULL)
  , m_buffer_ptr(NULL)
  , m_buffer_end(NULL)
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
        close();
    } catch (const std::exception& error) {
        print_locker lock;
        std::cerr << "Error closing file: " << error.what() << std::endl;
        std::exit(1);
    }
}


bool line_reader::getline(std::string& dst)
{
    dst.clear();

    while (m_file && !m_eof) {
        const char* start = m_buffer_ptr;
        char* end = m_buffer_ptr;

        for (; end != m_buffer_end; ++end) {
            if (*end == '\n') {
                // Excluding terminal \n
                dst.append(start, end - start);
                if (!dst.empty() && dst.back() == '\r') {
                    // Excluding terminal \r; dst is examined, since the \r may
                    // have been added seperately, if the \r was the last
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


void line_reader::close()
{
    close_buffers_gzip();
    close_buffers_bzip2();

    delete[] m_raw_buffer;
    m_raw_buffer = NULL;

    if (m_file && fclose(m_file)) {
        throw io_error("line_reader::close: error closing file", errno);
    }

    m_file = NULL;
}


bool line_reader::eof() const
{
    return m_eof;
}


bool line_reader::is_open() const
{
    return m_file;
}


void line_reader::refill_buffers()
{
    if (m_buffer) {
#ifdef AR_GZIP_SUPPORT
        if (m_gzip_stream) {
            refill_buffers_gzip();
        } else
#endif

#ifdef AR_BZIP2_SUPPORT
        if (m_bzip2_stream) {
            refill_buffers_bzip2();
        } else
#endif

        {
            refill_raw_buffer();
            refill_buffers_uncompressed();
        }
    } else {
        refill_raw_buffer();

        if (identify_gzip()) {
            initialize_buffers_gzip();
        } else if (identify_bzip2()) {
            initialize_buffers_bzip2();
        } else {
            refill_buffers_uncompressed();
        }
    }
}


void line_reader::refill_buffers_uncompressed()
{
    m_buffer = m_raw_buffer;
    m_buffer_ptr = m_raw_buffer;
    m_buffer_end = m_raw_buffer_end;
}


void line_reader::refill_raw_buffer()
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


bool line_reader::identify_gzip() const
{
    if (m_raw_buffer_end - m_raw_buffer < 2) {
        return false;
    } else if (m_raw_buffer[0] != '\x1f' || m_raw_buffer[1] != '\x8b') {
        return false;
    }

    return true;
}


void line_reader::initialize_buffers_gzip()
{
#ifdef AR_GZIP_SUPPORT
    m_buffer = new char[BUF_SIZE];
    m_buffer_ptr = m_buffer + BUF_SIZE;
    m_buffer_end = m_buffer + BUF_SIZE;

    m_gzip_stream = new z_stream();
    m_gzip_stream->zalloc = Z_NULL;
    m_gzip_stream->zfree = Z_NULL;
    m_gzip_stream->opaque = Z_NULL;

    m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
    m_gzip_stream->next_in = reinterpret_cast<Bytef*>(m_raw_buffer);

    switch (inflateInit2(m_gzip_stream, 15 + 16)) {
        case Z_OK:
            break;

        case Z_MEM_ERROR:
            throw gzip_error("line_reader::initialize_buffers_gzip: insufficient memory",
                             m_gzip_stream ? m_gzip_stream->msg : NULL);

        case Z_VERSION_ERROR:
            throw gzip_error("line_reader::initialize_buffers_gzip: incompatible zlib version",
                             m_gzip_stream ? m_gzip_stream->msg : NULL);

        case Z_STREAM_ERROR:
            throw gzip_error("line_reader::initialize_buffers_gzip: invalid parameters",
                             m_gzip_stream ? m_gzip_stream->msg : NULL);

        default:
            throw gzip_error("line_reader::initialize_buffers_gzip: unknown error",
                             m_gzip_stream ? m_gzip_stream->msg : NULL);
    }
#else
    throw gzip_error("Attempted to read gzipped file, but gzip"
                     "support was not enabled when AdapterRemoval"
                     "was compiled");
#endif
}


void line_reader::refill_buffers_gzip()
{
#ifdef AR_GZIP_SUPPORT
    if (!m_gzip_stream->avail_in) {
        refill_raw_buffer();
        m_gzip_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
        m_gzip_stream->next_in = reinterpret_cast<Bytef*>(m_raw_buffer);
    }

    m_gzip_stream->avail_out = BUF_SIZE;
    m_gzip_stream->next_out = reinterpret_cast<Bytef*>(m_buffer);
    switch (inflate(m_gzip_stream, Z_NO_FLUSH)) {
        case Z_OK:
        case Z_BUF_ERROR: /* input buffer empty or output buffer full */
            break;

        case Z_STREAM_END:
            // Handle concatenated streams; causes unnessesary reset at EOF
            if (inflateReset(m_gzip_stream) != Z_OK) {
                throw gzip_error("line_reader::refill_buffers_gzip: failed to reset stream",
                             m_gzip_stream ? m_gzip_stream->msg : NULL);
            }
            break;

        case Z_STREAM_ERROR:
            throw gzip_error("line_reader::refill_buffers_gzip: inconsistent stream state",
                             m_gzip_stream ? m_gzip_stream->msg : NULL);

        default:
            throw gzip_error("line_reader::refill_buffers_gzip: unknown error",
                             m_gzip_stream ? m_gzip_stream->msg : NULL);
    }

    m_buffer_ptr = m_buffer;
    m_buffer_end = m_buffer + (BUF_SIZE - m_gzip_stream->avail_out);
#endif
}


void line_reader::close_buffers_gzip()
{
#ifdef AR_GZIP_SUPPORT
    if (m_gzip_stream) {
        switch (inflateEnd(m_gzip_stream)) {
            case Z_OK:
                break;

            case Z_STREAM_ERROR:
                throw gzip_error("line_reader::close: stream error",
                                 m_gzip_stream ? m_gzip_stream->msg : NULL);

            default:
                throw gzip_error("Unknown error in line_reader::close",
                                 m_gzip_stream ? m_gzip_stream->msg : NULL);
        }

        delete m_gzip_stream;
        m_gzip_stream = NULL;

        delete[] m_buffer;
        m_buffer = NULL;
    }
#endif
}


#ifdef AR_BZIP2_SUPPORT

void bzip2_initialize_stream(bz_stream* stream)
{
    switch (BZ2_bzDecompressInit(stream, 0, 0)) {
        case BZ_OK:
            break;

        case BZ_CONFIG_ERROR:
            throw bzip2_error("bzip2_initialize_buffer: "
                              "bzip2 library is miscompiled");

        case BZ_PARAM_ERROR:
            throw bzip2_error("bzip2_initialize_buffer: "
                              "invalid parameters during initialization");

        case BZ_MEM_ERROR:
            throw bzip2_error("bzip2_initialize_buffer: "
                              "insufficient memory to initialize");

        case BZ_SEQUENCE_ERROR:
            throw bzip2_error("bzip2_initialize_buffer: bzip2 sequence error");

        default:
            throw bzip2_error("bzip2_initialize_buffer: unknown error");
    }
}


void bzip2_close_stream(bz_stream* stream)
{
    switch (BZ2_bzDecompressEnd(stream)) {
        case BZ_OK:
            break;

        case BZ_PARAM_ERROR:
            throw bzip2_error("bzip2_close_stream: invalid parameters");

        case BZ_SEQUENCE_ERROR:
            throw bzip2_error("bzip2_close_stream: bzip2 sequence error");

        default:
            throw bzip2_error("bzip2_close_stream: unknown bzip2 error");
    }
}

#endif


bool line_reader::identify_bzip2() const
{
    if (m_raw_buffer_end - m_raw_buffer < 4) {
        return false;
    } else if (m_raw_buffer[0] != 'B' || m_raw_buffer[1] != 'Z') {
        // Fixed magic header "BZ"
        return false;
    } else if (m_raw_buffer[2] != 'h' && m_raw_buffer[2] != '0') {
        // bzip2 or bzip1 (deprecated)
        return false;
    } else if (m_raw_buffer[3] < '1' || m_raw_buffer[3] > '9') {
        // Blocksizes; '1' - '9'
        return false;
    }

    return true;
}


void line_reader::initialize_buffers_bzip2()
{
#ifdef AR_BZIP2_SUPPORT
    m_buffer = new char[BUF_SIZE];
    m_buffer_ptr = m_buffer + BUF_SIZE;
    m_buffer_end = m_buffer + BUF_SIZE;

    m_bzip2_stream = new bz_stream();
    m_bzip2_stream->bzalloc = NULL;
    m_bzip2_stream->bzfree = NULL;
    m_bzip2_stream->opaque = NULL;

    m_bzip2_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
    m_bzip2_stream->next_in = m_raw_buffer;

    bzip2_initialize_stream(m_bzip2_stream);
#else
    throw bzip2_error("Attempted to read bzipped file, but bzip2"
                      "support was not enabled when AdapterRemoval"
                      "was compiled");
#endif
}


void line_reader::refill_buffers_bzip2()
{
#ifdef AR_BZIP2_SUPPORT
    if (!m_bzip2_stream->avail_in) {
        refill_raw_buffer();
        m_bzip2_stream->avail_in = m_raw_buffer_end - m_raw_buffer;
        m_bzip2_stream->next_in = m_raw_buffer;
    }

    m_bzip2_stream->avail_out = BUF_SIZE;
    m_bzip2_stream->next_out = m_buffer;
    if (m_bzip2_stream->avail_in) {
        switch (BZ2_bzDecompress(m_bzip2_stream)) {
            case BZ_OK:
                break;

            case BZ_STREAM_END:
                // Close an restart stream, to handle concatenated files
                bzip2_close_stream(m_bzip2_stream);
                bzip2_initialize_stream(m_bzip2_stream);
                break;

            case BZ_PARAM_ERROR:
                throw bzip2_error("line_reader::refill_buffers_bzip2: "
                                  "inconsistent bzip2 parameters");

            case BZ_DATA_ERROR:
            case BZ_DATA_ERROR_MAGIC:
                throw bzip2_error("line_reader::refill_buffers_bzip2: "
                                  "malformed bzip2 file");

            case BZ_MEM_ERROR:
                throw bzip2_error("line_reader::refill_buffers_bzip2: "
                                  "insufficient memory to deflate bzip2 stream");

            case BZ_SEQUENCE_ERROR:
                throw bzip2_error("line_reader::refill_buffers_bzip2: "
                                  "bzip2 sequence error");

            default:
                throw bzip2_error("line_reader::refill_buffers_bzip2: "
                                  "unknown bzip2 error");
        }
    }

    m_buffer_ptr = m_buffer;
    m_buffer_end = m_buffer + (BUF_SIZE - m_bzip2_stream->avail_out);
#endif
}


void line_reader::close_buffers_bzip2()
{
#ifdef AR_BZIP2_SUPPORT
    if (m_bzip2_stream) {
        bzip2_close_stream(m_bzip2_stream);

        delete m_bzip2_stream;
        m_bzip2_stream = NULL;

        delete[] m_buffer;
        m_buffer = NULL;
    }
#endif
}

} // namespace ar
