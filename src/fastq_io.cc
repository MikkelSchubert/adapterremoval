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
#include <stdexcept>
#include <iostream>
#include <cerrno>
#include <cstring>

#include "debug.h"
#include "fastq_io.h"
#include "userconfig.h"

namespace ar
{


size_t read_fastq_reads(fastq_vec& dst, line_reader& reader, size_t offset,
                      const fastq_encoding& encoding)
{
    dst.reserve(FASTQ_CHUNK_SIZE);

    try {
        fastq record;
        for (size_t i = 0; i < FASTQ_CHUNK_SIZE; ++i) {
            if (record.read(reader, encoding)) {
                dst.push_back(record);
            } else {
                break;
            }
        }
    } catch (const fastq_error& error) {
        print_locker lock;
        std::cerr << "Error reading FASTQ record at line "
                  << offset + dst.size()
                  << "; aborting:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;

        throw thread_abort();
    }

    return dst.size();
}



///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_read_chunk'

fastq_read_chunk::fastq_read_chunk(bool eof_)
  : eof(eof_)
  , reads_1()
  , reads_2()
{
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_output_chunk'

fastq_output_chunk::fastq_output_chunk(bool eof_)
  : eof(eof_)
  , count(0)
  , reads()
  , buffers()
{
    reads.reserve(FASTQ_CHUNK_SIZE);
}


fastq_output_chunk::~fastq_output_chunk()
{
    for (buffer_vec::iterator it = buffers.begin(); it != buffers.end(); ++it) {
        delete[] it->second;
    }
}


void fastq_output_chunk::add(const fastq_encoding& encoding,
                             const fastq& read, size_t count_)
{
    count += count_;
    reads.push_back(read.to_str(encoding));
}



///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_single_fastq'

read_single_fastq::read_single_fastq(const fastq_encoding* encoding,
                                     const std::string& filename,
                                     size_t next_step)
  : analytical_step(analytical_step::ordered, true)
  , m_encoding(encoding)
  , m_line_offset(1)
  , m_io_input(filename)
  , m_next_step(next_step)
  , m_eof(false)
{
}


chunk_vec read_single_fastq::process(analytical_chunk* chunk)
{
    AR_DEBUG_ASSERT(chunk == NULL);
    if (!m_io_input.is_open()) {
        return chunk_vec();
    }

    read_chunk_ptr file_chunk(new fastq_read_chunk());

    const size_t n_read = read_fastq_reads(file_chunk->reads_1, m_io_input,
                                           m_line_offset, *m_encoding);

    if (!n_read) {
        // EOF is detected by failure to read any lines, not line_reader::eof,
        // so that unbalanced files can be caught in all cases.
        m_io_input.close();
        file_chunk->eof = true;
        m_eof = true;
    }

    m_line_offset += n_read;

    chunk_vec chunks;
    chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));

    return chunks;
}


void read_single_fastq::finalize()
{
    if (!m_eof) {
        throw thread_error("read_single_fastq::finalize: terminated before EOF");
    }
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_paired_fastq'

read_paired_fastq::read_paired_fastq(const fastq_encoding* encoding,
                                     const std::string& filename_1,
                                     const std::string& filename_2,
                                     size_t next_step)
  : analytical_step(analytical_step::ordered, true)
  , m_encoding(encoding)
  , m_line_offset(1)
  , m_io_input_1(filename_1)
  , m_io_input_2(filename_2)
  , m_next_step(next_step)
  , m_eof(false)
{
}


chunk_vec read_paired_fastq::process(analytical_chunk* chunk)
{
    AR_DEBUG_ASSERT(chunk == NULL);
    if (!m_io_input_1.is_open()) {
        return chunk_vec();
    }

    read_chunk_ptr file_chunk(new fastq_read_chunk());

    const size_t n_read_1 = read_fastq_reads(file_chunk->reads_1, m_io_input_1,
                                             m_line_offset, *m_encoding);
    const size_t n_read_2 = read_fastq_reads(file_chunk->reads_2, m_io_input_2,
                                             m_line_offset, *m_encoding);

    if (n_read_1 != n_read_2) {
        print_locker lock;
        std::cerr << "ERROR: Input --file1 and --file2 contains different "
                  << "numbers of lines; one or the other file may have been "
                  << "truncated. Please correct before continuing!"
                  << std::endl;

        throw thread_abort();
    } else if (!n_read_1) {
        // EOF is detected by failure to read any lines, not line_reader::eof,
        // so that unbalanced files can be caught in all cases.
        m_io_input_1.close();
        m_io_input_2.close();
        file_chunk->eof = true;
        m_eof = true;
    }

    m_line_offset += n_read_1;

    chunk_vec chunks;
    chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));

    return chunks;
}


void read_paired_fastq::finalize()
{
    if (!m_eof) {
        throw thread_error("read_paired_fastq::finalize: terminated before EOF");
    }
}



///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_interleaved_fastq'

read_interleaved_fastq::read_interleaved_fastq(const fastq_encoding* encoding,
                                          const std::string& filename,
                                          size_t next_step)
  : analytical_step(analytical_step::ordered, true)
  , m_encoding(encoding)
  , m_line_offset(1)
  , m_io_input(filename)
  , m_next_step(next_step)
  , m_eof(false)
{
}


chunk_vec read_interleaved_fastq::process(analytical_chunk* chunk)
{
    AR_DEBUG_ASSERT(chunk == NULL);
    if (!m_io_input.is_open()) {
        return chunk_vec();
    }

    read_chunk_ptr file_chunk(new fastq_read_chunk());

    file_chunk->reads_1.reserve(FASTQ_CHUNK_SIZE);
    file_chunk->reads_2.reserve(FASTQ_CHUNK_SIZE);

    try {
        fastq record;
        for (size_t i = 0; i < FASTQ_CHUNK_SIZE; ++i) {
            // Mate 1 reads
            if (record.read(m_io_input, *m_encoding)) {
                file_chunk->reads_1.push_back(record);
            } else {
                break;
            }

            // Mate 2 reads
            if (record.read(m_io_input, *m_encoding)) {
                file_chunk->reads_2.push_back(record);
            } else {
                break;
            }
        }
    } catch (const fastq_error& error) {
        const size_t offset = m_line_offset
            + file_chunk->reads_1.size() * 4
            + file_chunk->reads_2.size() * 4;

        print_locker lock;
        std::cerr << "Error reading FASTQ record starting at line "
                  << offset << ":\n"
                  << cli_formatter::fmt(error.what()) << std::endl;

        throw thread_abort();
    }

    const size_t n_read_1 = file_chunk->reads_1.size();
    const size_t n_read_2 = file_chunk->reads_2.size();

    if (n_read_1 != n_read_2) {
        print_locker lock;
        std::cerr << "ERROR: Interleaved FASTQ file contains uneven number of "
                  << "reads; file may have been truncated! Please correct "
                  << "before continuing!"
                  << std::endl;

        throw thread_abort();
    } else if (!n_read_1) {
        m_io_input.close();
        file_chunk->eof = true;
        m_eof = true;
    }

    m_line_offset += (n_read_1 + n_read_2) * 4;

    chunk_vec chunks;
    chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));

    return chunks;
}


void read_interleaved_fastq::finalize()
{
    if (!m_eof) {
        throw thread_error("read_interleaved_fastq::finalize: terminated before EOF");
    }
}


///////////////////////////////////////////////////////////////////////////////
// Utility function used by both gzip and bzip compression steps

/**
 * Writes a set of lines into a buffer, and returns the size of the buffer and
 * the buffer itself as a pair. */
std::pair<size_t, unsigned char*> build_input_buffer(const string_vec& lines)
{
    size_t buffer_size = 0;
    for (string_vec::const_iterator it = lines.begin(); it != lines.end(); ++it) {
        buffer_size += it->size();
    }

    unsigned char* input_buffer = new unsigned char[buffer_size];
    unsigned char* input_buffer_ptr = input_buffer;
    for (string_vec::const_iterator it = lines.begin(); it != lines.end(); ++it) {
        std::memcpy(input_buffer_ptr, it->data(), it->size());
        input_buffer_ptr += it->size();
    }

    return std::pair<size_t, unsigned char*>(buffer_size, input_buffer);
}


#ifdef AR_BZIP2_SUPPORT

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_fastq'


bzip2_fastq::bzip2_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordered, false)
  , m_buffered_reads(0)
  , m_next_step(next_step)
  , m_stream()
  , m_eof(false)
{
    m_stream.bzalloc = NULL;
    m_stream.bzfree = NULL;
    m_stream.opaque = NULL;

    const int errorcode = BZ2_bzCompressInit(/* strm          = */ &m_stream,
                                             /* blockSize100k = */ config.bzip2_level,
                                             /* verbosity     = */ 0,
                                             /* workFactor    = */ 0);

    switch (errorcode) {
        case BZ_OK:
            break;

        case BZ_MEM_ERROR:
            throw thread_error("bzip2_fastq: not enough memory");

        case BZ_CONFIG_ERROR:
            throw thread_error("bzip2_fastq: miscompiled bzip2 library");

        case BZ_PARAM_ERROR:
            throw thread_error("bzip2_fastq: invalid parameters");

        default:
            throw thread_error("bzip2_fastq: unknown error");
    }
}


void bzip2_fastq::finalize()
{
    if (!m_eof) {
        throw thread_error("bzip2_fastq::finalize: terminated before EOF");
    }

    const int errorcode = BZ2_bzCompressEnd(&m_stream);
    if (errorcode != BZ_OK) {
        print_locker lock;
        switch (errorcode) {
            case BZ_PARAM_ERROR:
                throw thread_error("bzip2_fastq::finalize: parameter error");

            default:
                throw thread_error("Unknown error in bzip2_fastq::finalize");
        }
    }
}


chunk_vec bzip2_fastq::process(analytical_chunk* chunk)
{
    output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
    buffer_vec& buffers = file_chunk->buffers;

    if (m_eof) {
        throw thread_error("bzip2_fastq::process: received data after EOF");
    }

    m_eof = file_chunk->eof;
    if (file_chunk->reads.empty() && !m_eof) {
        return chunk_vec();
    }

    std::pair<size_t, unsigned char*> input_buffer;
    std::pair<size_t, unsigned char*> output_buffer;
    try {
        input_buffer = build_input_buffer(file_chunk->reads);

        m_stream.avail_in = input_buffer.first;
        m_stream.next_in = reinterpret_cast<char*>(input_buffer.second);

        if (m_stream.avail_in || m_eof) {
            int errorcode = -1;

            do {
                output_buffer.first = FASTQ_COMPRESSED_CHUNK;
                output_buffer.second = new unsigned char[FASTQ_COMPRESSED_CHUNK];

                m_stream.avail_out = output_buffer.first;
                m_stream.next_out = reinterpret_cast<char*>(output_buffer.second);

                errorcode = BZ2_bzCompress(&m_stream, m_eof ? BZ_FINISH : BZ_RUN);
                switch (errorcode) {
                    case BZ_RUN_OK:
                    case BZ_FINISH_OK:
                    case BZ_STREAM_END:
                        break;

                    case BZ_FLUSH_OK:
                        throw thread_error("bzip2_fastq::process: BZ_FLUSH_OK");

                    case BZ_PARAM_ERROR:
                        throw thread_error("bzip2_fastq::process: BZ_PARAM_ERROR");

                    case BZ_SEQUENCE_ERROR:
                        throw thread_error("bzip2_fastq::process: sequence error");

                    default:
                        throw thread_error("bzip2_fastq::process: unknown error");
                }

                output_buffer.first = FASTQ_COMPRESSED_CHUNK - m_stream.avail_out;
                if (output_buffer.first) {
                    buffers.push_back(output_buffer);
                } else {
                    delete[] output_buffer.second;
                }

                output_buffer.second = NULL;
            } while (m_stream.avail_in || errorcode == BZ_FINISH_OK);
        }

        delete[] input_buffer.second;
    } catch (...) {
        delete[] input_buffer.second;
        delete[] output_buffer.second;
        throw;
    }

    chunk_vec chunks;
    if (!file_chunk->buffers.empty() || m_eof) {
        file_chunk->count += m_buffered_reads;
        chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));
        m_buffered_reads = 0;
    } else {
        m_buffered_reads += file_chunk->count;
    }

    return chunks;
}

#endif


#ifdef AR_GZIP_SUPPORT

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_fastq'

gzip_fastq::gzip_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordered, false)
  , m_buffered_reads(0)
  , m_next_step(next_step)
  , m_stream()
  , m_eof(false)
{
    m_stream.zalloc = Z_NULL;
    m_stream.zfree = Z_NULL;
    m_stream.opaque = Z_NULL;

    const int errorcode = deflateInit2(/* strm       = */ &m_stream,
                                       /* level      = */ config.gzip_level,
                                       /* method     = */ Z_DEFLATED,
                                       /* windowBits = */ 15 + 16,
                                       /* memLevel   = */ 8,
                                       /* strategy   = */ Z_DEFAULT_STRATEGY);

    switch (errorcode) {
        case Z_OK:
            break;

        case Z_MEM_ERROR:
            throw thread_error("gzip_fastq: not enough memory");

        case Z_STREAM_ERROR:
            throw thread_error("gzip_fastq: invalid parameters");

        case Z_VERSION_ERROR:
            throw thread_error("gzip_fastq: incompatible zlib version");

        default:
            throw thread_error("gzip_fastq: unknown error");
    }
}


void gzip_fastq::finalize()
{
    if (!m_eof) {
        throw thread_error("gzip_fastq::finalize: terminated before EOF");
    }

    const int errorcode = deflateEnd(&m_stream);
    if (errorcode != Z_OK) {
        print_locker lock;
        switch (errorcode) {
            case Z_STREAM_ERROR:
                throw thread_error("gzip_fastq::finalize: stream error");

            case Z_DATA_ERROR:
                throw thread_error("gzip_fastq::finalize: data error");

            default:
                throw thread_error("Unknown error in gzip_fastq::finalize");
        }
    }
}


chunk_vec gzip_fastq::process(analytical_chunk* chunk)
{
    output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
    buffer_vec& buffers = file_chunk->buffers;

    if (m_eof) {
        throw thread_error("bzip2_fastq::process: received data after EOF");
    }

    m_eof = file_chunk->eof;
    if (file_chunk->reads.empty() && !m_eof) {
        return chunk_vec();
    }

    std::pair<size_t, unsigned char*> input_buffer;
    std::pair<size_t, unsigned char*> output_buffer;
    try {
        input_buffer = build_input_buffer(file_chunk->reads);
        file_chunk->reads.clear();

        if (input_buffer.first || m_eof) {
            m_stream.avail_in = input_buffer.first;
            m_stream.next_in = input_buffer.second;
            int returncode = -1;

            do {
                output_buffer.first = FASTQ_COMPRESSED_CHUNK;
                output_buffer.second = new unsigned char[FASTQ_COMPRESSED_CHUNK];

                m_stream.avail_out = output_buffer.first;
                m_stream.next_out = output_buffer.second;

                returncode = deflate(&m_stream, m_eof ? Z_FINISH : Z_NO_FLUSH);
                switch (returncode) {
                    case Z_OK:
                    case Z_STREAM_END:
                        break;

                    case Z_BUF_ERROR:
                        throw thread_error("gzip_fastq::process: buf error");

                    case Z_STREAM_ERROR:
                        throw thread_error("gzip_fastq::process: stream error");

                    default:
                        throw thread_error("gzip_fastq::process: unknown error");
                }

                output_buffer.first = FASTQ_COMPRESSED_CHUNK - m_stream.avail_out;
                if (output_buffer.first) {
                    buffers.push_back(output_buffer);
                } else {
                    delete[] output_buffer.second;
                }

                output_buffer.second = NULL;
            } while (m_stream.avail_out == 0 || (m_eof && returncode != Z_STREAM_END));
        }

        delete[] input_buffer.second;
    } catch (...) {
        delete[] input_buffer.second;
        delete[] output_buffer.second;
        throw;
    }

    chunk_vec chunks;
    if (!file_chunk->buffers.empty() || m_eof) {
        file_chunk->count += m_buffered_reads;
        chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));
        m_buffered_reads = 0;
    } else {
        m_buffered_reads += file_chunk->count;
    }

    return chunks;
}

#endif


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

//! Mutex used to control access to s_timer and s_finalized;
static std::mutex s_timer_lock;
//! Timer used to track trimming progress; accessed by all instances
static timer s_timer = timer("reads");
//! Indicates if 'timer::finalize' has been called.
static bool s_finalized = false;


write_fastq::write_fastq(const std::string& filename)
  : analytical_step(analytical_step::ordered, true)
  , m_output(filename.c_str(), std::ofstream::out | std::ofstream::binary)
  , m_eof(false)
{
    if (!m_output.is_open()) {
        std::string message = std::string("Failed to open file '") + filename + "': ";
        throw std::ofstream::failure(message + std::strerror(errno));
    }

    m_output.exceptions(std::ofstream::failbit | std::ofstream::badbit);
}


chunk_vec write_fastq::process(analytical_chunk* chunk)
{
    output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
    const string_vec& lines = file_chunk->reads;

    if (m_eof) {
        throw thread_error("write_fastq::process: received data after EOF");
    }

    m_eof = file_chunk->eof;
    if (file_chunk->buffers.empty()) {
        for (string_vec::const_iterator it = lines.begin(); it != lines.end(); ++it) {
            m_output << *it;
        }
    } else {
        buffer_vec& buffers = file_chunk->buffers;
        for (buffer_vec::iterator it = buffers.begin(); it != buffers.end(); ++it) {
            if (it->first) {
                m_output.write(reinterpret_cast<char*>(it->second), it->first);
            }
        }
    }

    if (m_eof) {
        m_output.flush();
    }

    std::lock_guard<std::mutex> lock(s_timer_lock);
    s_timer.increment(file_chunk->count);

    return chunk_vec();
}


void write_fastq::finalize()
{
    std::lock_guard<std::mutex> lock(s_timer_lock);
    if (!s_finalized) {
        s_timer.finalize();
        s_finalized = true;
    }

    if (!m_eof) {
        throw thread_error("write_fastq::finalize: terminated before EOF");
    }

    // Close file to trigger any exceptions due to badbit / failbit
    m_output.close();
}

} // namespace ar
