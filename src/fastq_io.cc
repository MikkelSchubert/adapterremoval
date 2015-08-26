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

#include "fastq_io.h"
#include "userconfig.h"


//! Number of lines to read for each data-chunk
const size_t CHUNK_SIZE = 4 * 1000;

#if defined(AR_GZIP_SUPPORT) || defined(AR_BZIP2_SUPPORT)
//! Size of compressed chunks used to transport compressed data
const size_t COMPRESSED_CHUNK = 40 * 1024;
#endif


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_file_chunk'

fastq_file_chunk::fastq_file_chunk(size_t offset_)
  : eof(false)
  , offset(offset_)
  , mates(2)
{
    typedef std::vector<string_vec>::iterator iter;
    for (iter it = mates.begin(); it != mates.end(); ++it) {
        it->reserve(CHUNK_SIZE);
    }
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_output_chunk'

fastq_output_chunk::fastq_output_chunk(bool eof_)
  : eof(eof_)
  , reads()
{
    reads.reserve(CHUNK_SIZE);
}


fastq_output_chunk::~fastq_output_chunk()
{
    for (buffer_vec::iterator it = buffers.begin(); it != buffers.end(); ++it) {
        delete[] it->second;
    }

}



///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_paired_fastq'

read_paired_fastq::read_paired_fastq(const userconfig& config,
                                     read_type mate,
                                     size_t next_step)
  : analytical_step(analytical_step::ordered, true)
  , m_line_offset(1)
  , m_io_input(get_filename(config, mate))
  , m_type(mate)
  , m_next_step(next_step)
{
}


chunk_list read_paired_fastq::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_file_chunk> file_chunk;
    if (chunk) {
        file_chunk.reset(dynamic_cast<fastq_file_chunk*>(chunk));
    } else {
        file_chunk.reset(new fastq_file_chunk(m_line_offset));
    }

    if (!m_io_input.is_open()) {
        return chunk_list();
    }

    string_vec& lines = file_chunk->mates.at(m_type);
    lines.resize(CHUNK_SIZE);

    for (size_t i = 0; i < CHUNK_SIZE; ++i) {
        if (!m_io_input.getline(lines.at(i))) {
            lines.resize(i);
            break;
        }
    }

    if (lines.empty()) {
        // EOF is detected by failure to read any lines, not line_reader::eof,
        // so that unbalanced files can be caught in all cases.
        m_io_input.close();
        file_chunk->eof = true;
    }

    m_line_offset += lines.size();

    chunk_list chunks;
    chunks.push_back(chunk_pair(m_next_step, file_chunk.release()));

    return chunks;
}


std::string read_paired_fastq::get_filename(const userconfig& config, read_type mate)
{
    if (mate == rt_mate_1) {
        return config.input_file_1;
    } else if (mate == rt_mate_2) {
        return config.input_file_2;
    } else {
        throw std::invalid_argument("mate must be rt_mate_1 or rt_mate_2");
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
// Implementations for 'gzip_paired_fastq'


bzip2_paired_fastq::bzip2_paired_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordered, false)
  , m_next_step(next_step)
  , m_stream()
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
            throw thread_error("bzip2_paired_fastq: not enough memory");

        case BZ_CONFIG_ERROR:
            throw thread_error("bzip2_paired_fastq: miscompiled bzip2 library");

        case BZ_PARAM_ERROR:
            throw thread_error("bzip2_paired_fastq: invalid parameters");

        default:
            throw thread_error("bzip2_paired_fastq: unknown error");
    }
}


bzip2_paired_fastq::~bzip2_paired_fastq()
{
    const int errorcode = BZ2_bzCompressEnd(&m_stream);
    if (errorcode != BZ_OK) {
        print_locker lock;
        switch (errorcode) {
            case BZ_PARAM_ERROR:
                std::cerr << "bzip2_paired_fastq::~bzip2_paired_fastq: parameter error" << std::endl;
                break;

            default:
                std::cerr << "Unknown error in bzip2_paired_fastq::~bzip2_paired_fastq: " << errorcode << std::endl;
                break;
        }

        std::exit(1);
    }
}


chunk_list bzip2_paired_fastq::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_output_chunk> file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
    buffer_vec& buffers = file_chunk->buffers;

    if (file_chunk->reads.empty() && !file_chunk->eof) {
        // The empty chunk must still be forwarded, to ensure that tracking of
        // ordered chunks does not break.
        chunk_list chunks;
        chunks.push_back(chunk_pair(m_next_step, file_chunk.release()));
        return chunks;
    }

    std::pair<size_t, unsigned char*> input_buffer;
    std::pair<size_t, unsigned char*> output_buffer;
    try {
        input_buffer = build_input_buffer(file_chunk->reads);

        m_stream.avail_in = input_buffer.first;
        m_stream.next_in = reinterpret_cast<char*>(input_buffer.second);

        do {
            output_buffer.first = COMPRESSED_CHUNK;
            output_buffer.second = new unsigned char[COMPRESSED_CHUNK];

            m_stream.avail_out = output_buffer.first;
            m_stream.next_out = reinterpret_cast<char*>(output_buffer.second);

            if (m_stream.avail_in || file_chunk->eof) {
                switch (BZ2_bzCompress(&m_stream, file_chunk->eof ? BZ_FINISH : BZ_RUN)) {
                    case BZ_RUN_OK:
                    case BZ_FINISH_OK:
                    case BZ_STREAM_END:
                        break;

                    case BZ_SEQUENCE_ERROR:
                        throw thread_error("bzip2_paired_fastq::process: sequence error");

                    default:
                        throw thread_error("bzip2_paired_fastq::process: unknown error");
                }
            }

            output_buffer.first = COMPRESSED_CHUNK - m_stream.avail_out;
            // A buffer must be sent, even if #bytes == 0.
            buffers.push_back(output_buffer);
            output_buffer.second = NULL;
        } while (m_stream.avail_out == 0);

        delete[] input_buffer.second;
    } catch (...) {
        delete[] input_buffer.second;
        delete[] output_buffer.second;
        throw;
    }

    chunk_list chunks;
    chunks.push_back(chunk_pair(m_next_step, file_chunk.release()));
    return chunks;
}

#endif


#ifdef AR_GZIP_SUPPORT

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_paired_fastq'

gzip_paired_fastq::gzip_paired_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordered, false)
  , m_next_step(next_step)
  , m_stream()
{
    m_stream.zalloc = Z_NULL;
    m_stream.zfree = Z_NULL;
    m_stream.opaque = Z_NULL;

    const int errorcode = deflateInit2(/* strm       = */ &m_stream,
                                       /* level      = */ config.gzip_level,
                                       /* method     = */ Z_DEFLATED,
                                       /* windowBits = */ 15 + 16,
                                       /* memLevel   = */ 9,
                                       /* strategy   = */ Z_DEFAULT_STRATEGY);

    switch (errorcode) {
        case Z_OK:
            break;

        case Z_MEM_ERROR:
            throw thread_error("gzip_paired_fastq: not enough memory");

        case Z_STREAM_ERROR:
            throw thread_error("gzip_paired_fastq: invalid parameters");

        case Z_VERSION_ERROR:
            throw thread_error("gzip_paired_fastq: incompatible zlib version");

        default:
            throw thread_error("gzip_paired_fastq: unknown error");
    }
}


gzip_paired_fastq::~gzip_paired_fastq()
{
    const int errorcode = deflateEnd(&m_stream);
    if (errorcode != Z_OK) {
        print_locker lock;
        switch (errorcode) {
            case Z_STREAM_ERROR:
                std::cerr << "gzip_paired_fastq::~gzip_paired_fastq: stream error" << std::endl;
                break;

            case Z_DATA_ERROR:
                std::cerr << "gzip_paired_fastq::~gzip_paired_fastq: data error" << std::endl;
                break;

            default:
                std::cerr << "Unknown error in gzip_paired_fastq::~gzip_paired_fastq: " << errorcode << std::endl;
                break;
        }

        std::exit(1);
    }
}


chunk_list gzip_paired_fastq::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_output_chunk> file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
    buffer_vec& buffers = file_chunk->buffers;

    if (file_chunk->reads.empty() && !file_chunk->eof) {
        // The empty chunk must still be forwarded, to ensure that tracking of
        // ordered chunks does not break.
        chunk_list chunks;
        chunks.push_back(chunk_pair(m_next_step, file_chunk.release()));
        return chunks;
    }

    std::pair<size_t, unsigned char*> input_buffer;
    std::pair<size_t, unsigned char*> output_buffer;
    try {
        input_buffer = build_input_buffer(file_chunk->reads);

        m_stream.avail_in = input_buffer.first;
        m_stream.next_in = input_buffer.second;

        do {
            output_buffer.first = COMPRESSED_CHUNK;
            output_buffer.second = new unsigned char[COMPRESSED_CHUNK];

            m_stream.avail_out = output_buffer.first;
            m_stream.next_out = output_buffer.second;

            switch (deflate(&m_stream, file_chunk->eof ? Z_FINISH : Z_NO_FLUSH)) {
                case Z_OK:
                case Z_STREAM_END:
                case Z_BUF_ERROR: /* End of out / in buffer reached. */
                    break;

                case Z_STREAM_ERROR:
                    throw thread_error("gzip_paired_fastq::process: stream error");

                default:
                    throw thread_error("gzip_paired_fastq::process: unknown error");
            }

            output_buffer.first = COMPRESSED_CHUNK - m_stream.avail_out;
            // A buffer must be sent, even if #bytes == 0.
            buffers.push_back(output_buffer);
            output_buffer.second = NULL;
        } while (m_stream.avail_out == 0);

        delete[] input_buffer.second;
    } catch (...) {
        delete[] input_buffer.second;
        delete[] output_buffer.second;
        throw;
    }

    chunk_list chunks;
    chunks.push_back(chunk_pair(m_next_step, file_chunk.release()));
    return chunks;
}

#endif


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_paired_fastq'

//! Mutex used to control access to s_timer and s_finalized;
static mutex s_timer_lock;
//! Timer used to track trimming progress; accessed by all instances
static timer s_timer = timer("reads", false);
//! Indicates if 'timer::finalize' has been called.
static bool s_finalized = false;


write_paired_fastq::write_paired_fastq(const userconfig& config,
                                       read_type type)
  : analytical_step(analytical_step::ordered, true)
  , m_type(type)
  , m_output()
  , m_progress(!config.quiet)
  , m_pair_ended(config.paired_ended_mode)
{
    switch (m_type) {
        case rt_mate_1:
            if (config.paired_ended_mode) {
                m_output = config.open_with_default_filename("--output1", ".pair1.truncated");
            } else {
                m_output = config.open_with_default_filename("--output1", ".truncated");
            }

            break;

        case rt_mate_2:
            m_output = config.open_with_default_filename("--output2", ".pair2.truncated");
            break;

        case rt_singleton:
            m_output = config.open_with_default_filename("--singleton", ".singleton.truncated");
            break;

        case rt_collapsed:
            m_output = config.open_with_default_filename("--outputcollapsed", ".collapsed");
            break;

        case rt_collapsed_truncated:
            m_output = config.open_with_default_filename("--outputcollapsedtruncated", ".collapsed.truncated");
            break;

        case rt_discarded:
            m_output = config.open_with_default_filename("--discarded", ".discarded");
            break;

        default:
            throw std::invalid_argument("invalid read-type in write_paired_fastq constructor");
    }
}


write_paired_fastq::~write_paired_fastq()
{
}


chunk_list write_paired_fastq::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_output_chunk> file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
    const string_vec& lines = file_chunk->reads;

    if (file_chunk->buffers.empty()) {
        for (string_vec::const_iterator it = lines.begin(); it != lines.end(); ++it) {
            *m_output << *it;
        }
    } else {
        buffer_vec& buffers = file_chunk->buffers;
        for (buffer_vec::iterator it = buffers.begin(); it != buffers.end(); ++it) {
            if (it->first) {
                m_output->write(reinterpret_cast<char*>(it->second), it->first);
            }
        }
    }

    if (file_chunk->eof) {
        m_output->flush();
    }

    if (m_progress) {
        size_t record_count = lines.size();
        if (m_pair_ended && (m_type == rt_collapsed || m_type == rt_collapsed_truncated)) {
            record_count *= 2;
        }

        mutex_locker lock(s_timer_lock);
        s_timer.increment(record_count);
    }

    return chunk_list();
}


void write_paired_fastq::finalize()
{
    mutex_locker lock(s_timer_lock);
    if (m_progress && !s_finalized) {
        s_timer.finalize();
        s_finalized = true;
    }
}

