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

#include <algorithm> // for max, min
#include <cerrno>    // for errno
#include <cstring>   // for size_t, strerror, memcpy
#include <iostream>  // for operator<<, basic_ostream, char_traits, endl
#include <utility>   // for move, swap

#ifdef USE_LIBDEFLATE
#include <libdeflate.h>
#endif

#include "debug.hpp"     // for AR_DEBUG_ASSERT, AR_DEBUG_LOCK
#include "fastq.hpp"     // for fastq
#include "fastq_enc.hpp" // for fastq_error
#include "fastq_io.hpp"
#include "statistics.hpp" // for fastq_statistics
#include "strutils.hpp"   // for cli_formatter
#include "threads.hpp"    // for thread_error, print_locker, thread_abort
#include "userconfig.hpp" // for userconfig

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
  , nucleotides()
  , reads_1()
  , reads_2()
{}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_output_chunk'

fastq_output_chunk::fastq_output_chunk(bool eof_)
  : eof(eof_)
  , nucleotides(0)
  , reads()
  , buffers()
{}

void
fastq_output_chunk::add(const fastq& read)
{
  nucleotides += read.length();
  read.into_string(reads);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_fastq'

bool
read_record(joined_line_readers& reader,
            fastq_vec* chunk,
            fastq& record,
            size_t& n_nucleotides)
{
  try {
    if (record.read_unsafe(reader)) {
      chunk->push_back(record);
      n_nucleotides += record.length();

      return true;
    } else {
      return false;
    }
  } catch (const fastq_error& error) {
    print_locker lock;
    std::cerr << "Error reading FASTQ record from '" << reader.filename()
              << "' at line " << reader.linenumber() << "; aborting:\n"
              << cli_formatter::fmt(error.what()) << std::endl;

    throw thread_abort();
  }
}

read_fastq::read_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordering::ordered, true)
  , m_io_input_1_base(config.input_files_1)
  , m_io_input_2_base(config.input_files_2)
  , m_io_input_1(&m_io_input_1_base)
  , m_io_input_2(&m_io_input_2_base)
  , m_next_step(next_step)
  , m_single_end(false)
  , m_eof(false)
  , m_timer("reads")
  , m_lock()
{
  if (config.interleaved_input) {
    AR_DEBUG_ASSERT(config.input_files_2.empty());

    m_io_input_2 = &m_io_input_1_base;
  } else if (config.input_files_2.empty()) {
    m_io_input_2 = &m_io_input_1_base;
    m_single_end = true;
  } else {
    AR_DEBUG_ASSERT(config.input_files_1.size() == config.input_files_2.size());
  }
}

chunk_vec
read_fastq::process(analytical_chunk* chunk)
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(chunk == nullptr);
  if (m_eof) {
    return chunk_vec();
  }

  read_chunk_ptr file_chunk(new fastq_read_chunk());
  auto reads_1 = &file_chunk->reads_1;
  auto reads_2 = m_single_end ? reads_1 : &file_chunk->reads_2;

  fastq record;
  bool eof = false;
  size_t n_nucleotides = 0;
  while (n_nucleotides < INPUT_BLOCK_SIZE * 4 && !eof) {
    eof = !read_record(*m_io_input_1, reads_1, record, n_nucleotides);

    bool eof_2 = !read_record(*m_io_input_2, reads_2, record, n_nucleotides);

    if (eof && !eof_2) {
      print_locker lock;
      std::cerr << "ERROR: More mate 2 reads than mate 1 reads found in '"
                << m_io_input_1->filename() << "'; file may be truncated. "
                << "Please fix before continuing." << std::endl;

      throw thread_abort();
    } else if (eof_2 && !(eof || m_single_end)) {
      print_locker lock;
      std::cerr << "ERROR: More mate 1 reads than mate 2 reads found in '"
                << m_io_input_2->filename() << "'; file may be truncated. "
                << "Please fix before continuing." << std::endl;

      throw thread_abort();
    }

    eof |= eof_2;
  }

  file_chunk->eof = eof;
  m_eof = eof;

  m_timer.increment(file_chunk->reads_1.size());
  m_timer.increment(file_chunk->reads_2.size());

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(file_chunk));

  return chunks;
}

void
read_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(m_eof);

  m_timer.finalize();
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'post_process_fastq'

post_process_fastq::post_process_fastq(const fastq_encoding& encoding,
                                       size_t next_step,
                                       ar_statistics* statistics)
  : analytical_step(analytical_step::ordering::ordered)
  , m_encoding(encoding)
  , m_statistics_1(statistics ? &statistics->input_1 : nullptr)
  , m_statistics_2(statistics ? &statistics->input_2 : nullptr)
  , m_next_step(next_step)
  , m_eof(false)
  , m_lock()
{}

chunk_vec
post_process_fastq::process(analytical_chunk* chunk)
{
  read_chunk_ptr file_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
  AR_DEBUG_ASSERT(!m_eof);
  AR_DEBUG_LOCK(m_lock);

  m_eof = file_chunk->eof;

  for (auto& read : file_chunk->reads_1) {
    read.post_process(m_encoding);
    if (m_statistics_1) {
      m_statistics_1->process(read);
    }
  }

  for (auto& read : file_chunk->reads_2) {
    read.post_process(m_encoding);
    if (m_statistics_2) {
      m_statistics_2->process(read);
    }
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(file_chunk));

  return chunks;
}

void
post_process_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(m_eof);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_fastq'

gzip_fastq::gzip_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordering::ordered, false)
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
  AR_DEBUG_ASSERT(m_eof);

  checked_deflate_end(&m_stream);
}

chunk_vec
gzip_fastq::process(analytical_chunk* chunk)
{
  output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));

  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(!m_eof);

  m_eof = file_chunk->eof;
  if (file_chunk->reads.empty() && !m_eof) {
    return chunk_vec();
  }

  buffer_pair output_buffer;

  m_stream.avail_in = file_chunk->reads.size();
  m_stream.next_in = reinterpret_cast<unsigned char*>(
    const_cast<char*>(file_chunk->reads.data()));

  int returncode = -1;
  do {
    output_buffer.first = OUTPUT_BLOCK_SIZE;
    output_buffer.second.reset(new unsigned char[OUTPUT_BLOCK_SIZE]);

    m_stream.avail_out = output_buffer.first;
    m_stream.next_out = output_buffer.second.get();

    returncode = checked_deflate(&m_stream, m_eof ? Z_FINISH : Z_NO_FLUSH);

    output_buffer.first = OUTPUT_BLOCK_SIZE - m_stream.avail_out;
    if (output_buffer.first) {
      file_chunk->buffers.push_back(std::move(output_buffer));
    }
  } while (m_stream.avail_out == 0 || (m_eof && returncode != Z_STREAM_END));

  file_chunk->reads.clear();

  chunk_vec chunks;
  if (!file_chunk->buffers.empty() || m_eof) {
    chunks.emplace_back(m_next_step, std::move(file_chunk));
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'split_fastq'

split_fastq::split_fastq(size_t next_step)
  : analytical_step(analytical_step::ordering::ordered, false)
  , m_next_step(next_step)
  , m_buffer(new unsigned char[GZIP_BLOCK_SIZE])
  , m_offset()
  , m_eof(false)
  , m_lock()
{}

void
split_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(m_eof);

  AR_DEBUG_ASSERT(!m_buffer);
}

chunk_vec
split_fastq::process(analytical_chunk* chunk)
{
  output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));

  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(!m_eof);
  m_eof = file_chunk->eof;

  chunk_vec chunks;
  const auto& src = file_chunk->reads;
  for (size_t src_offset = 0; src_offset < src.size();) {
    const auto n =
      std::min(src.size() - src_offset, GZIP_BLOCK_SIZE - m_offset);
    std::memcpy(m_buffer.get() + m_offset, src.data() + src_offset, n);

    src_offset += n;
    m_offset += n;

    if (m_offset == GZIP_BLOCK_SIZE) {
      output_chunk_ptr block(new fastq_output_chunk());
      block->buffers.emplace_back(GZIP_BLOCK_SIZE, std::move(m_buffer));

      chunks.emplace_back(m_next_step, std::move(block));

      m_buffer.reset(new unsigned char[GZIP_BLOCK_SIZE]);
      m_offset = 0;
    }
  }

  if (m_eof) {
    output_chunk_ptr block(new fastq_output_chunk(true));
    block->buffers.emplace_back(m_offset, std::move(m_buffer));
    chunks.emplace_back(m_next_step, std::move(block));

    m_offset = 0;
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_split_fastq'

gzip_split_fastq::gzip_split_fastq(const userconfig& config, size_t next_step)
  : analytical_step(analytical_step::ordering::unordered, false)
  , m_config(config)
  , m_next_step(next_step)
  , m_buffers()
{}

chunk_vec
gzip_split_fastq::process(analytical_chunk* chunk)
{
  output_chunk_ptr input_chunk(dynamic_cast<fastq_output_chunk*>(chunk));
  AR_DEBUG_ASSERT(input_chunk->buffers.size() == 1);

  buffer_pair& input_buffer = input_chunk->buffers.front();
  buffer_pair output_buffer;

  // Try to re-use a previous input-buffer
  output_buffer.first = GZIP_BLOCK_SIZE;
  output_buffer.second = m_buffers.try_acquire();
  if (!output_buffer.second) {
    output_buffer.second.reset(new unsigned char[GZIP_BLOCK_SIZE]);
  }

#ifdef USE_LIBDEFLATE
  auto compressor = libdeflate_alloc_compressor(m_config.gzip_level);
  auto compressed_size = libdeflate_gzip_compress(compressor,
                                                  input_buffer.second.get(),
                                                  input_buffer.first,
                                                  output_buffer.second.get(),
                                                  output_buffer.first);
  libdeflate_free_compressor(compressor);

  // The easily compressible input should fit in a single output block
  AR_DEBUG_ASSERT(compressed_size);
#else
  z_stream stream;
  checked_deflate_init2(&stream, m_config.gzip_level);

  stream.avail_in = input_buffer.first;
  stream.next_in = input_buffer.second.get();
  stream.avail_out = output_buffer.first;
  stream.next_out = output_buffer.second.get();
  // The easily compressible input should fit in a single output block
  const int returncode = checked_deflate(&stream, Z_FINISH);
  AR_DEBUG_ASSERT(stream.avail_out && returncode == Z_STREAM_END);

  const auto compressed_size = GZIP_BLOCK_SIZE - stream.avail_out;
  checked_deflate_end(&stream);
#endif

  // Re-use the input chunk and make its buffer available for re-use above
  output_buffer.first = compressed_size;
  std::swap(input_buffer, output_buffer);
  m_buffers.release(output_buffer.second);

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(input_chunk));

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

const std::string STDOUT = "/dev/stdout";

write_fastq::write_fastq(const std::string& filename)
  : analytical_step(analytical_step::ordering::ordered, filename != STDOUT)
  , m_output(filename)
  , m_eof(false)
  , m_lock()
{}

chunk_vec
write_fastq::process(analytical_chunk* chunk)
{
  output_chunk_ptr file_chunk(dynamic_cast<fastq_output_chunk*>(chunk));

  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(!m_eof);

  try {
    m_eof = file_chunk->eof;
    if (file_chunk->buffers.empty()) {
      m_output.write_string(file_chunk->reads, m_eof);
    } else {
      AR_DEBUG_ASSERT(file_chunk->reads.empty());

      m_output.write_buffers(file_chunk->buffers, m_eof);
    }
  } catch (const std::ios_base::failure&) {
    const std::string message =
      std::string("Error writing FASTQ file '") + m_output.filename() + "': ";
    throw std::ofstream::failure(message + std::strerror(errno));
  }

  return chunk_vec();
}

void
write_fastq::finalize()
{
  AR_DEBUG_LOCK(m_lock);
  AR_DEBUG_ASSERT(m_eof);

  // Close file to trigger any exceptions due to badbit / failbit
  try {
    m_output.close();
  } catch (const std::ios_base::failure&) {
    const std::string message =
      std::string("Error closing FASTQ file '") + m_output.filename() + "': ";
    throw std::ofstream::failure(message + std::strerror(errno));
  }
}
