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

#ifdef USE_LIBDEFLATE
#include <libdeflate.h>
#endif

#include "debug.hpp"
#include "fastq_io.hpp"
#include "userconfig.hpp"

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
  , count(0)
  , nucleotides()
  , reads()
  , buffers()
{}

void
fastq_output_chunk::add(const fastq_encoding& encoding,
                        const fastq& read,
                        size_t count_)
{
  count += count_;
  nucleotides += read.length();
  read.into_string(reads, encoding);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_fastq'

bool
read_record(joined_line_readers& reader,
            const fastq_encoding& encoding,
            fastq_statistics* stats,
            fastq_vec* chunk,
            fastq& record,
            size_t& n_nucleotides)
{
  try {
    if (record.read(reader, encoding)) {
      if (stats) {
        stats->process(record);
      }

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

read_fastq::read_fastq(const fastq_encoding* encoding,
                       const string_vec& filenames_1,
                       const string_vec& filenames_2,
                       size_t next_step,
                       bool interleaved,
                       fastq_statistics* statistics_1,
                       fastq_statistics* statistics_2)
  : analytical_step(analytical_step::ordering::ordered, true)
  , m_encoding(encoding)
  , m_io_input_1_base(filenames_1)
  , m_io_input_2_base(filenames_2)
  , m_io_input_1(&m_io_input_1_base)
  , m_io_input_2(&m_io_input_2_base)
  , m_statistics_1(statistics_1)
  , m_statistics_2(statistics_2)
  , m_next_step(next_step)
  , m_single_end(false)
  , m_eof(false)
  , m_timer("reads")
  , m_lock()
{
  if (interleaved) {
    AR_DEBUG_ASSERT(filenames_2.empty());

    m_io_input_2 = &m_io_input_1_base;
  } else if (filenames_2.empty()) {
    m_io_input_2 = &m_io_input_1_base;
    m_statistics_2 = statistics_1;
    m_single_end = true;
  } else {
    AR_DEBUG_ASSERT(filenames_1.size() == filenames_2.size());
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
    eof = !read_record(*m_io_input_1,
                       *m_encoding,
                       m_statistics_1,
                       reads_1,
                       record,
                       n_nucleotides);

    bool eof_2 = !read_record(*m_io_input_2,
                              *m_encoding,
                              m_statistics_2,
                              reads_2,
                              record,
                              n_nucleotides);

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
    file_chunk->count += m_buffered_reads;
    chunks.emplace_back(m_next_step, std::move(file_chunk));
    m_buffered_reads = 0;
  } else {
    m_buffered_reads += file_chunk->count;
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
  , m_count()
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
      block->count = m_count;

      chunks.emplace_back(m_next_step, std::move(block));

      m_buffer.reset(new unsigned char[GZIP_BLOCK_SIZE]);
      m_offset = 0;
      m_count = 0;
    }
  }

  // The next chunk will contain the last of these reads; incrementing the
  // timer then is easier than trying to figure out how many reads are in each
  // buffer. We ignore the case were m_buffer is currently empty to keep
  // things simple.
  m_count += file_chunk->count;

  if (m_eof) {
    output_chunk_ptr block(new fastq_output_chunk(true));
    block->buffers.emplace_back(m_offset, std::move(m_buffer));
    block->count = m_count;

    chunks.emplace_back(m_next_step, std::move(block));

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

  output_chunk_ptr block(new fastq_output_chunk(input_chunk->eof));
  output_buffer.first = compressed_size;
  block->buffers.emplace_back(std::move(output_buffer));
  block->count = input_chunk->count;
  output_buffer.second = nullptr;

  // Make the input buffer available for re-use
  m_buffers.release(input_buffer.second);

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(block));

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
