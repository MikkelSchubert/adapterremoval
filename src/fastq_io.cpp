/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <memory>    // for make_unique
#include <utility>   // for move, swap

#ifdef USE_LIBDEFLATE
#include <libdeflate.h>
#endif

#include "debug.hpp"      // for AR_REQUIRE, AR_REQUIRE_SINGLE_THREAD
#include "fastq.hpp"      // for fastq
#include "fastq_enc.hpp"  // for fastq_error
#include "fastq_io.hpp"   // declarations
#include "logging.hpp"    // for log
#include "statistics.hpp" // for fastq_statistics
#include "strutils.hpp"   // for cli_formatter
#include "threads.hpp"    // for thread_error, print_locker, thread_abort
#include "userconfig.hpp" // for userconfig

namespace adapterremoval {

///////////////////////////////////////////////////////////////////////////////
// Helper functions for 'zlib'

void
checked_deflate_init2(z_streamp stream, unsigned int level)
{
  AR_REQUIRE(stream);
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
{
}

fastq_read_chunk::~fastq_read_chunk() {}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_output_chunk'

fastq_output_chunk::fastq_output_chunk(bool eof_)
  : eof(eof_)
  , nucleotides(0)
  , reads()
  , buffers()
{
}

void
fastq_output_chunk::add(const fastq& read)
{
  nucleotides += read.length();
  read.into_string(reads);
}

fastq_output_chunk::~fastq_output_chunk() {}

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
    log::error() << "Error reading FASTQ record from '" << reader.filename()
                 << "' at line " << reader.linenumber() << "; aborting:\n"
                 << cli_formatter::fmt(error.what());

    throw thread_abort();
  }
}

read_fastq::read_fastq(const userconfig& config, size_t next_step)
  : analytical_step(processing_order::ordered_io, "read_fastq")
  , m_io_input_1_base(config.input_files_1)
  , m_io_input_2_base(config.input_files_2)
  , m_io_input_1(&m_io_input_1_base)
  , m_io_input_2(&m_io_input_2_base)
  , m_next_step(next_step)
  , m_single_end(false)
  , m_eof(false)
  , m_timer(progress_type::spinner)
  , m_head(config.head)
  , m_lock()
{
  if (config.interleaved_input) {
    AR_REQUIRE(config.input_files_2.empty());

    m_io_input_2 = &m_io_input_1_base;
  } else if (config.input_files_2.empty()) {
    m_io_input_2 = &m_io_input_1_base;
    m_single_end = true;
  } else {
    AR_REQUIRE(config.input_files_1.size() == config.input_files_2.size());
  }
}

chunk_vec
read_fastq::process(chunk_ptr chunk)
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!chunk);
  if (m_eof) {
    return chunk_vec();
  }

  auto file_chunk = std::make_unique<fastq_read_chunk>();
  if (m_single_end) {
    read_single_end(file_chunk);
  } else {
    read_paired_end(file_chunk);
  }

  m_eof |= !m_head;
  file_chunk->eof = m_eof;

  m_timer.increment(file_chunk->reads_1.size() + file_chunk->reads_2.size());

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(file_chunk));

  return chunks;
}

void
read_fastq::read_single_end(read_chunk_ptr& chunk)
{
  fastq record;
  size_t n_nucleotides = 0;
  while (n_nucleotides < INPUT_BLOCK_SIZE && m_head && !m_eof) {
    m_eof = !read_record(*m_io_input_1, &chunk->reads_1, record, n_nucleotides);

    if (m_head != std::numeric_limits<unsigned>::max()) {
      m_head--;
    }
  }
}

void
read_fastq::read_paired_end(read_chunk_ptr& chunk)
{
  auto reads_1 = &chunk->reads_1;
  auto reads_2 = &chunk->reads_2;

  fastq record;
  size_t n_nucleotides = 0;
  while (n_nucleotides < INPUT_BLOCK_SIZE && m_head && !m_eof) {
    m_eof = !read_record(*m_io_input_1, reads_1, record, n_nucleotides);

    bool eof_2 = !read_record(*m_io_input_2, reads_2, record, n_nucleotides);

    if (m_eof && !eof_2) {
      log::error() << "More mate 2 reads than mate 1 reads found in '"
                   << m_io_input_1->filename() << "'; file may be truncated. "
                   << "Please fix before continuing.";

      throw thread_abort();
    } else if (eof_2 && !m_eof) {
      log::error() << "More mate 1 reads than mate 2 reads found in '"
                   << m_io_input_2->filename() << "'; file may be truncated. "
                   << "Please fix before continuing.";

      throw thread_abort();
    }

    if (m_head != std::numeric_limits<unsigned>::max()) {
      m_head--;
    }
  }
}

void
read_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);

  m_timer.finalize();
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'post_process_fastq'

post_process_fastq::post_process_fastq(const userconfig& config,
                                       size_t next_step,
                                       statistics* stats)
  : analytical_step(processing_order::ordered, "post_process_fastq")
  , m_encoding(config.io_encoding)
  , m_mate_separator(config.mate_separator)
  , m_statistics_1(stats ? stats->input_1 : nullptr)
  , m_statistics_2(stats ? stats->input_2 : nullptr)
  , m_next_step(next_step)
  , m_eof(false)
  , m_lock()
{
}

chunk_vec
post_process_fastq::process(chunk_ptr chunk)
{
  auto& file_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  AR_REQUIRE(!m_eof);
  AR_REQUIRE_SINGLE_THREAD(m_lock);

  m_eof = file_chunk.eof;

  auto& reads_1 = file_chunk.reads_1;
  auto& reads_2 = file_chunk.reads_2;

  if (reads_1.size() == reads_2.size()) {
    process_paired_end(reads_1, reads_2);
  } else {
    AR_REQUIRE(reads_2.empty());
    process_single_end(reads_1);
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

void
post_process_fastq::process_single_end(fastq_vec& reads_1)
{
  for (auto& read : reads_1) {
    read.post_process(m_encoding);

    if (m_statistics_1) {
      m_statistics_1->process(read);
    }
  }
}

void
post_process_fastq::process_paired_end(fastq_vec& reads_1, fastq_vec& reads_2)
{
  AR_REQUIRE(reads_1.size() == reads_2.size());

  if (!m_mate_separator) {
    // Attempt to determine the mate separator character
    m_mate_separator = fastq::guess_mate_separator(reads_1, reads_2);

    if (!m_mate_separator) {
      // Fall back to the default so that a human-readable error will be thrown
      // during normalization below.
      m_mate_separator = MATE_SEPARATOR;
    }
  }

  auto it_1 = reads_1.begin();
  auto it_2 = reads_2.begin();
  while (it_1 != reads_1.end()) {
    fastq& read_1 = *it_1++;
    fastq& read_2 = *it_2++;

    // Throws if read-names or mate numbering does not match
    fastq::normalize_paired_reads(read_1, read_2, m_mate_separator);

    if (m_statistics_1) {
      m_statistics_1->process(read_1);
      m_statistics_2->process(read_2);
    }
  }
}

void
post_process_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_fastq'

gzip_fastq::gzip_fastq(const userconfig& config, size_t next_step)
  : analytical_step(processing_order::ordered, "gzip_fastq")
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
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);

  checked_deflate_end(&m_stream);
}

chunk_vec
gzip_fastq::process(chunk_ptr chunk)
{
  auto& file_chunk = dynamic_cast<fastq_output_chunk&>(*chunk);

  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!m_eof);

  m_eof = file_chunk.eof;
  if (file_chunk.reads.empty() && !m_eof) {
    return chunk_vec();
  }

  buffer_pair output_buffer;

  m_stream.avail_in = file_chunk.reads.size();
  m_stream.next_in = reinterpret_cast<unsigned char*>(
    const_cast<char*>(file_chunk.reads.data()));

  int returncode = -1;
  do {
    output_buffer.first = OUTPUT_BLOCK_SIZE;
    output_buffer.second.reset(new unsigned char[OUTPUT_BLOCK_SIZE]);

    m_stream.avail_out = output_buffer.first;
    m_stream.next_out = output_buffer.second.get();

    returncode = checked_deflate(&m_stream, m_eof ? Z_FINISH : Z_NO_FLUSH);

    output_buffer.first = OUTPUT_BLOCK_SIZE - m_stream.avail_out;
    if (output_buffer.first) {
      file_chunk.buffers.push_back(std::move(output_buffer));
    }
  } while (m_stream.avail_out == 0 || (m_eof && returncode != Z_STREAM_END));

  file_chunk.reads.clear();

  chunk_vec chunks;
  if (!file_chunk.buffers.empty() || m_eof) {
    chunks.emplace_back(m_next_step, std::move(chunk));
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'split_fastq'

split_fastq::split_fastq(size_t next_step)
  : analytical_step(processing_order::ordered, "split_fastq")
  , m_next_step(next_step)
  , m_buffer(new unsigned char[GZIP_BLOCK_SIZE])
  , m_offset()
  , m_eof(false)
  , m_lock()
{
}

void
split_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);

  AR_REQUIRE(!m_buffer);
}

chunk_vec
split_fastq::process(chunk_ptr chunk)
{
  auto& file_chunk = dynamic_cast<fastq_output_chunk&>(*chunk);

  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!m_eof);
  m_eof = file_chunk.eof;

  chunk_vec chunks;
  const auto& src = file_chunk.reads;
  for (size_t src_offset = 0; src_offset < src.size();) {
    const auto n =
      std::min(src.size() - src_offset, GZIP_BLOCK_SIZE - m_offset);
    std::memcpy(m_buffer.get() + m_offset, src.data() + src_offset, n);

    src_offset += n;
    m_offset += n;

    if (m_offset == GZIP_BLOCK_SIZE) {
      auto block = std::make_unique<fastq_output_chunk>();
      block->buffers.emplace_back(GZIP_BLOCK_SIZE, std::move(m_buffer));

      chunks.emplace_back(m_next_step, std::move(block));

      m_buffer.reset(new unsigned char[GZIP_BLOCK_SIZE]);
      m_offset = 0;
    }
  }

  if (m_eof) {
    auto block = std::make_unique<fastq_output_chunk>(true);
    block->buffers.emplace_back(m_offset, std::move(m_buffer));
    chunks.emplace_back(m_next_step, std::move(block));

    m_offset = 0;
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_split_fastq'

gzip_split_fastq::gzip_split_fastq(const userconfig& config, size_t next_step)
  : analytical_step(processing_order::unordered, "gzip_split_fastq")
  , m_config(config)
  , m_next_step(next_step)
  , m_buffers()
{
}

chunk_vec
gzip_split_fastq::process(chunk_ptr chunk)
{
  auto& input_chunk = dynamic_cast<fastq_output_chunk&>(*chunk);
  AR_REQUIRE(input_chunk.buffers.size() == 1);

  buffer_pair& input_buffer = input_chunk.buffers.front();
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
  AR_REQUIRE(compressed_size);
#else
  z_stream stream;
  checked_deflate_init2(&stream, m_config.gzip_level);

  stream.avail_in = input_buffer.first;
  stream.next_in = input_buffer.second.get();
  stream.avail_out = output_buffer.first;
  stream.next_out = output_buffer.second.get();
  // The easily compressible input should fit in a single output block
  const int returncode = checked_deflate(&stream, Z_FINISH);
  AR_REQUIRE(stream.avail_out && returncode == Z_STREAM_END);

  const auto compressed_size = GZIP_BLOCK_SIZE - stream.avail_out;
  checked_deflate_end(&stream);
#endif

  // Re-use the input chunk and make its buffer available for re-use above
  output_buffer.first = compressed_size;
  std::swap(input_buffer, output_buffer);
  m_buffers.release(output_buffer.second);

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

const std::string STDOUT = "/dev/stdout";

write_fastq::write_fastq(const std::string& filename)
  // Allow disk IO and writing to STDOUT at the same time
  : analytical_step(filename == STDOUT ? processing_order::ordered
                                       : processing_order::ordered_io,
                    "write_fastq")
  , m_output(filename)
  , m_eof(false)
  , m_lock()
{
}

chunk_vec
write_fastq::process(chunk_ptr chunk)
{
  auto& file_chunk = dynamic_cast<fastq_output_chunk&>(*chunk);

  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!m_eof);

  try {
    m_eof = file_chunk.eof;
    if (file_chunk.buffers.empty()) {
      m_output.write_string(file_chunk.reads, m_eof);
    } else {
      AR_REQUIRE(file_chunk.reads.empty());

      m_output.write_buffers(file_chunk.buffers, m_eof);
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
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);

  // Close file to trigger any exceptions due to badbit / failbit
  try {
    m_output.close();
  } catch (const std::ios_base::failure&) {
    const std::string message =
      std::string("Error closing FASTQ file '") + m_output.filename() + "': ";
    throw std::ofstream::failure(message + std::strerror(errno));
  }
}

} // namespace adapterremoval
