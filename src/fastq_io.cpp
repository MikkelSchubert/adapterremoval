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
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "debug.hpp"
#include "fastq_io.hpp"
#include "userconfig.hpp"

size_t
read_fastq_reads(fastq_vec& dst,
                 joined_line_readers& reader,
                 size_t offset,
                 const fastq_encoding& encoding,
                 fastq_statistics* stats)
{
  dst.reserve(FASTQ_CHUNK_SIZE);

  try {
    fastq record;
    for (size_t i = 0; i < FASTQ_CHUNK_SIZE; ++i) {
      if (record.read(reader, encoding)) {
        if (stats) {
          stats->process(record);
        }

        dst.push_back(record);
      } else {
        break;
      }
    }
  } catch (const fastq_error& error) {
    print_locker lock;
    std::cerr << "Error reading FASTQ record at line "
              << offset + dst.size() * 4 << "; aborting:\n"
              << cli_formatter::fmt(error.what()) << std::endl;

    throw thread_abort();
  }

  return dst.size() * 4;
}

///////////////////////////////////////////////////////////////////////////////
// Helper functions for 'zlib'

void
checked_deflate_init2(z_streamp stream, unsigned int level)
{
  AR_DEBUG_ASSERT(stream);
  stream->zalloc = nullptr;
  stream->zfree = nullptr;
  stream->opaque = nullptr;

  const int errorcode = deflateInit2(/* strm       = */ stream,
                                     /* level      = */ level,
                                     /* method     = */ Z_DEFLATED,
                                     /* windowBits = */ 15 + 16,
                                     /* memLevel   = */ 8,
                                     /* strategy   = */ Z_DEFAULT_STRATEGY);

  switch (errorcode) {
    case Z_OK:
      break;

    case Z_MEM_ERROR:
      throw thread_error("gzip deflateInit2: not enough memory");

    case Z_STREAM_ERROR:
      throw thread_error("gzip deflateInit2: invalid parameters");

    case Z_VERSION_ERROR:
      throw thread_error("gzip deflateInit2: incompatible zlib version");

    default:
      throw thread_error("gzip deflateInit2: unknown error");
  }
}

int
checked_deflate(z_streamp stream, int flush)
{
  const int errorcode = deflate(stream, flush);
  switch (errorcode) {
    case Z_OK:
    case Z_STREAM_END:
      break;

    case Z_BUF_ERROR:
      throw thread_error("gzip deflate: buf error");

    case Z_STREAM_ERROR:
      throw thread_error("gzip deflate: stream error");

    default:
      throw thread_error("gzip deflate: unknown error");
  }

  return errorcode;
}

void
checked_deflate_end(z_streamp stream)
{
  switch (deflateEnd(stream)) {
    case Z_OK:
      break;

    case Z_STREAM_ERROR:
      throw thread_error("gzip deflateEnd: stream error");

    case Z_DATA_ERROR:
      throw thread_error("gzip deflateEnd: data error");

    default:
      throw thread_error("gzip deflateEnd: unknown error");
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_read_chunk'

fastq_read_chunk::fastq_read_chunk(bool eof_)
  : eof(eof_)
  , reads_1()
  , reads_2()
{}

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
  for (auto& buffer : buffers) {
    delete[] buffer.second;
  }
}

void
fastq_output_chunk::add(const fastq_encoding& encoding,
                        const fastq& read,
                        size_t count_)
{
  count += count_;
  reads.push_back(read.to_str(encoding));
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_single_fastq'

read_single_fastq::read_single_fastq(const fastq_encoding* encoding,
                                     const string_vec& filenames,
                                     size_t next_step,
                                     fastq_statistics* statistics)
  : analytical_step(analytical_step::ordering::ordered, true)
  , m_encoding(encoding)
  , m_line_offset(1)
  , m_io_input(filenames)
  , m_statistics(statistics)
  , m_next_step(next_step)
  , m_eof(false)
  , m_lock()
{
  AR_DEBUG_ASSERT(!filenames.empty());
}

chunk_vec
read_single_fastq::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(chunk == nullptr);
  if (m_eof) {
    return chunk_vec();
  }

  read_chunk_ptr file_chunk(new fastq_read_chunk());

  const size_t n_read = read_fastq_reads(
    file_chunk->reads_1, m_io_input, m_line_offset, *m_encoding, m_statistics);

  if (!n_read) {
    // EOF is detected by failure to read any lines, not line_reader::eof,
    // so that unbalanced files can be caught in all cases.
    file_chunk->eof = true;
    m_eof = true;
  }

  m_line_offset += n_read;

  chunk_vec chunks;
  chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));

  return chunks;
}

void
read_single_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  if (!m_eof) {
    throw thread_error("read_single_fastq::finalize: terminated before EOF");
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_paired_fastq'

read_paired_fastq::read_paired_fastq(const fastq_encoding* encoding,
                                     const string_vec& filenames_1,
                                     const string_vec& filenames_2,
                                     size_t next_step,
                                     fastq_statistics* statistics_1,
                                     fastq_statistics* statistics_2)
  : analytical_step(analytical_step::ordering::ordered, true)
  , m_encoding(encoding)
  , m_line_offset(1)
  , m_io_input_1(filenames_1)
  , m_io_input_2(filenames_2)
  , m_statistics_1(statistics_1)
  , m_statistics_2(statistics_2)
  , m_next_step(next_step)
  , m_eof(false)
  , m_lock()
{
  AR_DEBUG_ASSERT(filenames_1.size() == filenames_2.size());
}

chunk_vec
read_paired_fastq::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(chunk == nullptr);
  if (m_eof) {
    return chunk_vec();
  }

  read_chunk_ptr file_chunk(new fastq_read_chunk());

  const size_t n_read_1 = read_fastq_reads(file_chunk->reads_1,
                                           m_io_input_1,
                                           m_line_offset,
                                           *m_encoding,
                                           m_statistics_1);
  const size_t n_read_2 = read_fastq_reads(file_chunk->reads_2,
                                           m_io_input_2,
                                           m_line_offset,
                                           *m_encoding,
                                           m_statistics_2);

  if (n_read_1 != n_read_2) {
    print_locker lock;
    std::cerr << "ERROR: Input --file1 and --file2 contains different "
              << "numbers of lines; one or the other file may have been "
              << "truncated. Please correct before continuing!" << std::endl;

    throw thread_abort();
  } else if (!n_read_1) {
    // EOF is detected by failure to read any lines, not line_reader::eof,
    // so that unbalanced files can be caught in all cases.
    file_chunk->eof = true;
    m_eof = true;
  }

  m_line_offset += n_read_1;

  chunk_vec chunks;
  chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));

  return chunks;
}

void
read_paired_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  if (!m_eof) {
    throw thread_error("read_paired_fastq::finalize: terminated before EOF");
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_interleaved_fastq'

read_interleaved_fastq::read_interleaved_fastq(const fastq_encoding* encoding,
                                               const string_vec& filenames,
                                               size_t next_step,
                                               fastq_statistics* statistics_1,
                                               fastq_statistics* statistics_2)
  : analytical_step(analytical_step::ordering::ordered, true)
  , m_encoding(encoding)
  , m_line_offset(1)
  , m_io_input(filenames)
  , m_statistics_1(statistics_1)
  , m_statistics_2(statistics_2)
  , m_next_step(next_step)
  , m_eof(false)
  , m_lock()
{
  AR_DEBUG_ASSERT(!filenames.empty());
}

chunk_vec
read_interleaved_fastq::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(chunk == nullptr);
  if (m_eof) {
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
        if (m_statistics_1) {
          m_statistics_1->process(record);
        }

        file_chunk->reads_1.push_back(record);
      } else {
        break;
      }

      // Mate 2 reads
      if (record.read(m_io_input, *m_encoding)) {
        if (m_statistics_2) {
          m_statistics_2->process(record);
        }

        file_chunk->reads_2.push_back(record);
      } else {
        break;
      }
    }
  } catch (const fastq_error& error) {
    const size_t offset = m_line_offset + file_chunk->reads_1.size() * 4 +
                          file_chunk->reads_2.size() * 4;

    print_locker lock;
    std::cerr << "Error reading FASTQ record starting at line " << offset
              << ":\n"
              << cli_formatter::fmt(error.what()) << std::endl;

    throw thread_abort();
  }

  const size_t n_read_1 = file_chunk->reads_1.size();
  const size_t n_read_2 = file_chunk->reads_2.size();

  if (n_read_1 != n_read_2) {
    print_locker lock;
    std::cerr << "ERROR: Interleaved FASTQ file contains uneven number of "
              << "reads; file may have been truncated! Please correct "
              << "before continuing!" << std::endl;

    throw thread_abort();
  } else if (!n_read_1) {
    file_chunk->eof = true;
    m_eof = true;
  }

  m_line_offset += (n_read_1 + n_read_2) * 4;

  chunk_vec chunks;
  chunks.push_back(chunk_pair(m_next_step, std::move(file_chunk)));

  return chunks;
}

void
read_interleaved_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  if (!m_eof) {
    throw thread_error(
      "read_interleaved_fastq::finalize: terminated before EOF");
  }
}

///////////////////////////////////////////////////////////////////////////////
// Utility function used by both gzip and bzip compression steps

/**
 * Writes a set of lines into a buffer, and returns the size of the buffer and
 * the buffer itself as a pair. */
std::pair<size_t, unsigned char*>
build_input_buffer(const string_vec& lines)
{
  size_t buffer_size = 0;
  for (const auto& line : lines) {
    buffer_size += line.size();
  }

  unsigned char* input_buffer = new unsigned char[buffer_size];
  unsigned char* input_buffer_ptr = input_buffer;
  for (const auto& line : lines) {
    std::memcpy(input_buffer_ptr, line.data(), line.size());
    input_buffer_ptr += line.size();
  }

  return std::pair<size_t, unsigned char*>(buffer_size, input_buffer);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'bzip2_fastq'

bzip2_fastq::bzip2_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordering::ordered, false)
  , m_buffered_reads(0)
  , m_next_step(next_step)
  , m_stream()
  , m_eof(false)
  , m_lock()
{
  m_stream.bzalloc = nullptr;
  m_stream.bzfree = nullptr;
  m_stream.opaque = nullptr;

  const int errorcode =
    BZ2_bzCompressInit(/* strm          = */ &m_stream,
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

void
bzip2_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
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

chunk_vec
bzip2_fastq::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
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

        output_buffer.second = nullptr;
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

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_fastq'

gzip_fastq::gzip_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordering::ordered, false)
  , m_buffered_reads(0)
  , m_next_step(next_step)
  , m_stream()
  , m_eof(false)
  , m_lock()
{
  checked_deflate_init2(&m_stream, config.gzip_level);
}

void
gzip_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  if (!m_eof) {
    throw thread_error("gzip_fastq::finalize: terminated before EOF");
  }

  checked_deflate_end(&m_stream);
}

chunk_vec
gzip_fastq::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
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

        returncode = checked_deflate(&m_stream, m_eof ? Z_FINISH : Z_NO_FLUSH);

        output_buffer.first = FASTQ_COMPRESSED_CHUNK - m_stream.avail_out;
        if (output_buffer.first) {
          buffers.push_back(output_buffer);
        } else {
          delete[] output_buffer.second;
        }

        output_buffer.second = nullptr;
      } while (m_stream.avail_out == 0 ||
               (m_eof && returncode != Z_STREAM_END));
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

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'split_fastq'

const size_t BLOCK_SIZE = 64 * 1024;

split_fastq::split_fastq(size_t next_step)
  : analytical_step(analytical_step::ordering::ordered, false)
  , m_next_step(next_step)
  , m_buffer(new unsigned char[BLOCK_SIZE])
  , m_offset()
  , m_count()
  , m_eof(false)
  , m_lock()
{}

void
split_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  if (!m_eof) {
    throw thread_error("split_fastq::finalize: terminated before EOF");
  }

  AR_DEBUG_ASSERT(!m_buffer);
}

chunk_vec
split_fastq::process(analytical_chunk* chunk)
{
  output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));

  AR_DEBUG_LOCK(m_lock);
  if (m_eof) {
    throw thread_error("split_fastq::process: received data after EOF");
  }

  m_eof = file_chunk->eof;

  chunk_vec chunks;
  for (const auto& line : file_chunk->reads) {
    size_t line_offset = 0;
    do {
      const auto n = std::min(line.size() - line_offset, BLOCK_SIZE - m_offset);
      std::memcpy(m_buffer + m_offset, line.data() + line_offset, n);

      line_offset += n;
      m_offset += n;

      if (m_offset == BLOCK_SIZE) {
        output_chunk_ptr block(new fastq_output_chunk());
        block->buffers.push_back(buffer_pair(BLOCK_SIZE, m_buffer));
        block->count = m_count;

        chunks.push_back(chunk_pair(m_next_step, std::move(block)));

        m_buffer = new unsigned char[BLOCK_SIZE];
        m_offset = 0;
        m_count = 0;
      }
    } while (line_offset < line.size());

    m_count++;
  }

  if (m_eof) {
    output_chunk_ptr block(new fastq_output_chunk(true));
    block->buffers.push_back(buffer_pair(m_offset, m_buffer));
    block->count = m_count;

    chunks.push_back(chunk_pair(m_next_step, std::move(block)));

    m_buffer = NULL;
    m_offset = 0;
    m_count = 0;
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_split_fastq'

gzip_split_fastq::gzip_split_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordering::unordered, false)
  , m_config(config)
  , m_next_step(next_step)
{}

chunk_vec
gzip_split_fastq::process(analytical_chunk* chunk)
{
  output_chunk_ptr input_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
  AR_DEBUG_ASSERT(input_chunk->buffers.size() == 1);

  z_stream stream;
  checked_deflate_init2(&stream, m_config.gzip_level);

  chunk_vec chunks;
  buffer_pair input_buffer = input_chunk->buffers.front();
  buffer_pair output_buffer;
  try {
    stream.avail_in = input_buffer.first;
    stream.next_in = input_buffer.second;

    output_buffer.first = BLOCK_SIZE;
    output_buffer.second = new unsigned char[BLOCK_SIZE];

    stream.avail_out = output_buffer.first;
    stream.next_out = output_buffer.second;

    // The easily compressible input should fit in a single output block
    const int returncode = checked_deflate(&stream, Z_FINISH);
    AR_DEBUG_ASSERT(stream.avail_out && returncode == Z_STREAM_END);

    output_chunk_ptr block(new fastq_output_chunk(input_chunk->eof));
    output_buffer.first = BLOCK_SIZE - stream.avail_out;
    block->buffers.push_back(output_buffer);
    block->count = input_chunk->count;
    output_buffer.second = nullptr;

    chunks.push_back(chunk_pair(m_next_step, std::move(block)));
  } catch (...) {
    delete[] output_buffer.second;
    throw;
  }

  checked_deflate_end(&stream);

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

//! Mutex used to control access to s_timer and s_finalized;
static std::mutex s_timer_lock;
//! Timer used to track trimming progress; accessed by all instances
static progress_timer s_timer = progress_timer("reads");
//! Indicates if 'timer::finalize' has been called.
static bool s_finalized = false;

write_fastq::write_fastq(const std::string& filename)
  : analytical_step(analytical_step::ordering::ordered, true)
  , m_output(filename)
  , m_eof(false)
  , m_lock()
{}

chunk_vec
write_fastq::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
  output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));

  if (m_eof) {
    throw thread_error("write_fastq::process: received data after EOF");
  }

  try {
    m_eof = file_chunk->eof;
    if (file_chunk->buffers.empty()) {
      m_output.write_strings(file_chunk->reads, m_eof);
    } else {
      AR_DEBUG_ASSERT(file_chunk->reads.empty());

      m_output.write_buffers(file_chunk->buffers, m_eof);
    }
  } catch (const std::ios_base::failure&) {
    const std::string message =
      std::string("Error writing FASTQ file '") + m_output.filename() + "': ";
    throw std::ofstream::failure(message + std::strerror(errno));
  }

  std::lock_guard<std::mutex> lock(s_timer_lock);
  s_timer.increment(file_chunk->count);

  return chunk_vec();
}

void
write_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  std::lock_guard<std::mutex> lock(s_timer_lock);
  if (!s_finalized) {
    s_timer.finalize();
    s_finalized = true;
  }

  if (!m_eof) {
    throw thread_error("write_fastq::finalize: terminated before EOF");
  }

  // Close file to trigger any exceptions due to badbit / failbit
  try {
    m_output.close();
  } catch (const std::ios_base::failure&) {
    const std::string message =
      std::string("Error closing FASTQ file '") + m_output.filename() + "': ";
    throw std::ofstream::failure(message + std::strerror(errno));
  }
}
