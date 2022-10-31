/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <cstring>   // for size_t, memcpy
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

////////////////////////////////////////////////////////////////////////////////
// Helper function for isa-l

namespace {

//! Supported ISAL compression level; used rather than ISAL_DEF_MAX_LEVEL since
//! the ISAL source code implies that may be increased
const uint32_t MAX_ISAL_LEVEL = 3;

uint32_t
isal_level(uint32_t gzip_level)
{
  return std::min<uint32_t>(gzip_level, MAX_ISAL_LEVEL);
}

size_t
isal_buffer_size(size_t level)
{
  switch (level) {
    default:
    case 3:
      return ISAL_DEF_LVL3_DEFAULT;
    case 2:
      return ISAL_DEF_LVL2_DEFAULT;
    case 1:
      return ISAL_DEF_LVL1_DEFAULT;
    case 0:
      return ISAL_DEF_LVL0_DEFAULT;
  }
}

bool
isal_enabled(const userconfig& config, const std::string& filename)
{
  return (config.gzip || ends_with(tolower(filename), ".gz"))
#ifdef USE_LIBDEFLATE
         && config.gzip_level <= MAX_ISAL_LEVEL;
#endif
  ;
}

} // namespace

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

fastq_output_chunk::fastq_output_chunk(bool eof_, uint32_t crc32_)
  : eof(eof_)
  , crc32(crc32_)
  , nucleotides()
  , reads()
  , buffers()
  , uncompressed_size()
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
            fastq_vec& chunk,
            size_t& n_nucleotides,
            const fastq_encoding& encoding)
{
  // Line numbers change as we attempt to read the record, and potentially
  // points to the next record in the case of invalid qualities/nucleotides
  const auto line_number = reader.linenumber();

  try {
    chunk.emplace_back();
    auto& record = chunk.back();

    if (record.read(reader, encoding)) {
      n_nucleotides += record.length();

      return true;
    } else {
      chunk.pop_back();
      return false;
    }
  } catch (const fastq_error& error) {
    log::error() << "Error reading FASTQ record from '" << reader.filename()
                 << "' at line " << line_number << "; aborting:\n"
                 << indent_lines(error.what());

    throw thread_abort();
  }
}

read_fastq::read_fastq(const userconfig& config, size_t next_step)
  : analytical_step(processing_order::ordered_io, "read_fastq")
  , m_encoding(config.io_encoding)
  , m_io_input_1_base(config.input_files_1)
  , m_io_input_2_base(config.input_files_2)
  , m_io_input_1(&m_io_input_1_base)
  , m_io_input_2(&m_io_input_2_base)
  , m_next_step(next_step)
  , m_mate_separator(config.mate_separator)
  , m_single_end(false)
  , m_eof(false)
  , m_timer(config.log_progress)
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
    m_eof = !read_single_end(file_chunk->reads_1);
  } else {
    m_eof = !read_paired_end(file_chunk->reads_1, file_chunk->reads_2);
  }

  file_chunk->eof = m_eof;

  m_timer.increment(file_chunk->reads_1.size() + file_chunk->reads_2.size());

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(file_chunk));

  return chunks;
}

bool
read_fastq::read_single_end(fastq_vec& reads_1)
{
  bool eof = false;
  size_t n_nucleotides = 0;
  while (n_nucleotides < INPUT_BLOCK_SIZE && m_head && !eof) {
    eof = !read_record(*m_io_input_1, reads_1, n_nucleotides, m_encoding);
    m_head--;
  }

  return !eof && m_head;
}

bool
read_fastq::read_paired_end(fastq_vec& reads_1, fastq_vec& reads_2)
{
  bool eof_1 = false;
  bool eof_2 = false;

  size_t n_nucleotides = 0;
  while (n_nucleotides < INPUT_BLOCK_SIZE && m_head && !eof_1 && !eof_2) {
    eof_1 = !read_record(*m_io_input_1, reads_1, n_nucleotides, m_encoding);
    eof_2 = !read_record(*m_io_input_2, reads_2, n_nucleotides, m_encoding);
    m_head--;
  }

  if (eof_1 && !eof_2) {
    log::error() << "More mate 2 reads than mate 1 reads found in '"
                 << m_io_input_1->filename() << "'; file may be truncated. "
                 << "Please fix before continuing.";

    throw thread_abort();
  } else if (eof_2 && !eof_1) {
    log::error() << "More mate 1 reads than mate 2 reads found in '"
                 << m_io_input_2->filename() << "'; file may be truncated. "
                 << "Please fix before continuing.";

    throw thread_abort();
  }

  if (!m_mate_separator) {
    m_mate_separator = identify_mate_separators(reads_1, reads_2);
  }

  auto it_1 = reads_1.begin();
  auto it_2 = reads_2.begin();
  while (it_1 != reads_1.end()) {
    fastq::normalize_paired_reads(*it_1++, *it_2++, m_mate_separator);
  }

  return !eof_1 && !eof_2 && m_head;
}

char
read_fastq::identify_mate_separators(const fastq_vec& reads_1,
                                     const fastq_vec& reads_2) const
{
  AR_REQUIRE(reads_1.size() == reads_2.size());

  // Attempt to determine the mate separator character
  char mate_separator = fastq::guess_mate_separator(reads_1, reads_2);

  if (!mate_separator) {
    // Fall back to the default so that a human-readable error will be thrown
    // during normalization below.
    mate_separator = MATE_SEPARATOR;
  }

  return mate_separator;
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

post_process_fastq::post_process_fastq(size_t next_step, statistics& stats)
  : analytical_step(processing_order::ordered, "post_process_fastq")
  , m_statistics_1(stats.input_1)
  , m_statistics_2(stats.input_2)
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
    m_statistics_1->process(read);
  }
}

void
post_process_fastq::process_paired_end(fastq_vec& reads_1, fastq_vec& reads_2)
{
  AR_REQUIRE(reads_1.size() == reads_2.size());

  for (const auto& read_1 : reads_1) {
    m_statistics_1->process(read_1);
  }

  for (const auto& read_2 : reads_2) {
    m_statistics_2->process(read_2);
  }
}

void
post_process_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'split_fastq'

split_fastq::split_fastq(const userconfig& config,
                         const std::string& filename,
                         size_t next_step)
  : analytical_step(processing_order::ordered, "split_fastq")
  , m_next_step(next_step)
  , m_buffer(GZIP_BLOCK_SIZE)
  , m_offset()
  , m_isal_enabled(isal_enabled(config, filename))
  , m_isal_crc32()
  , m_eof(false)
  , m_lock()
{
}

void
split_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);

  AR_REQUIRE(m_eof);
  AR_REQUIRE(m_buffer.capacity() == 0);
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
      std::min(src.size() - src_offset, m_buffer.size() - m_offset);
    std::memcpy(m_buffer.get() + m_offset, src.data() + src_offset, n);

    src_offset += n;
    m_offset += n;

    if (m_offset == m_buffer.size()) {
      auto block = std::make_unique<fastq_output_chunk>();

      if (m_isal_enabled) {
        m_isal_crc32 =
          crc32_gzip_refl(m_isal_crc32, m_buffer.get(), m_buffer.size());
      }

      block->uncompressed_size = m_offset;
      block->buffers.emplace_back(std::move(m_buffer));

      chunks.emplace_back(m_next_step, std::move(block));

      m_buffer = buffer(GZIP_BLOCK_SIZE);
      m_offset = 0;
    }
  }

  if (m_eof) {
    auto block = std::make_unique<fastq_output_chunk>(true);

    m_buffer.resize(m_offset);
    if (m_isal_enabled) {
      m_isal_crc32 =
        crc32_gzip_refl(m_isal_crc32, m_buffer.get(), m_buffer.size());
    }

    block->crc32 = m_isal_crc32;
    block->uncompressed_size = m_offset;
    block->buffers.emplace_back(std::move(m_buffer));

    chunks.emplace_back(m_next_step, std::move(block));

    m_offset = 0;
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_split_fastq'

gzip_split_fastq::gzip_split_fastq(const userconfig& config,
                                   const std::string& filename,
                                   size_t next_step)
  : analytical_step(processing_order::unordered, "gzip_split_fastq")
  , m_config(config)
  , m_isal_enabled(isal_enabled(config, filename))
  , m_next_step(next_step)
{
}

chunk_vec
gzip_split_fastq::process(chunk_ptr chunk)
{
  auto& input_chunk = dynamic_cast<fastq_output_chunk&>(*chunk);
  AR_REQUIRE(input_chunk.buffers.size() == 1);

  // Enable re-use of the fastq_output_chunk
  buffer& output_buffer = input_chunk.buffers.front();
  buffer input_buffer(GZIP_BLOCK_SIZE);
  std::swap(input_buffer, output_buffer);

  size_t compressed_size = 0;

  if (m_isal_enabled) {
    isal_zstream stream;
    isal_deflate_stateless_init(&stream);

    stream.flush = FULL_FLUSH;
    stream.end_of_stream = input_chunk.eof;

    stream.level = isal_level(m_config.gzip_level);
    stream.level_buf_size = isal_buffer_size(stream.level);
    buffer level_buffer(stream.level_buf_size);
    stream.level_buf = level_buffer.get();

    stream.avail_in = input_buffer.size();
    stream.next_in = input_buffer.get();
    stream.next_out = output_buffer.get();
    stream.avail_out = output_buffer.size();

    switch (isal_deflate_stateless(&stream)) {
      case COMP_OK:
        break;
      case INVALID_FLUSH:
        throw thread_error("isal_deflate_stateless: invalid flush");
      case ISAL_INVALID_LEVEL:
        throw thread_error("isal_deflate_stateless: invalid level");
      case ISAL_INVALID_LEVEL_BUF:
        throw thread_error("isal_deflate_stateless: invalid buffer size");
      default:
        throw thread_error("isal_deflate_stateless: unexpected error");
    }

    // The easily compressible input should fit in a single output block
    AR_REQUIRE(stream.avail_in == 0);

    compressed_size = stream.total_out;
  }
#ifdef USE_LIBDEFLATE
  else {
    auto compressor = libdeflate_alloc_compressor(m_config.gzip_level);
    compressed_size = libdeflate_gzip_compress(compressor,
                                               input_buffer.get(),
                                               input_buffer.size(),
                                               output_buffer.get(),
                                               output_buffer.size());
    libdeflate_free_compressor(compressor);

    // The easily compressible input should fit in a single output block
    AR_REQUIRE(compressed_size);
  }
#endif

  output_buffer.resize(compressed_size);

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

const std::string STDOUT = "/dev/stdout";

write_fastq::write_fastq(const userconfig& config, const std::string& filename)
  // Allow disk IO and writing to STDOUT at the same time
  : analytical_step(filename == STDOUT ? processing_order::ordered
                                       : processing_order::ordered_io,
                    "write_fastq")
  , m_output(filename)
  , m_isal_enabled(isal_enabled(config, filename))
  , m_uncompressed_bytes()
  , m_eof(false)
  , m_lock()
{
  if (m_isal_enabled) {
    buffer level_buf(isal_buffer_size(config.gzip_level));
    buffer output_buf(OUTPUT_BLOCK_SIZE);

    struct isal_zstream stream;
    struct isal_gzip_header header;

    isal_gzip_header_init(&header);
    isal_deflate_init(&stream);
    stream.avail_in = 0;
    stream.flush = NO_FLUSH;
    stream.level = isal_level(config.gzip_level);
    stream.level_buf = level_buf.get();
    stream.level_buf_size = level_buf.size();
    stream.gzip_flag = IGZIP_GZIP_NO_HDR;
    stream.next_out = output_buf.get();
    stream.avail_out = output_buf.size();

    const auto ret = isal_write_gzip_header(&stream, &header);
    AR_REQUIRE(ret == 0, "buffer was not large enough for gzip header");

    output_buf.resize(stream.total_out);

    m_output.write_buffer(output_buf, false);
  }
}

chunk_vec
write_fastq::process(chunk_ptr chunk)
{
  auto& file_chunk = dynamic_cast<fastq_output_chunk&>(*chunk);

  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!m_eof);

  try {
    m_eof = file_chunk.eof;
    m_uncompressed_bytes += file_chunk.uncompressed_size;

    const bool flush = m_eof && !m_isal_enabled;

    if (file_chunk.buffers.empty()) {
      m_output.write_string(file_chunk.reads, flush);
    } else {
      AR_REQUIRE(file_chunk.reads.empty());

      m_output.write_buffers(file_chunk.buffers, flush);
    }

    if (m_eof && m_isal_enabled) {
      buffer trailer(8);
      trailer.write_u32(0, file_chunk.crc32);
      trailer.write_u32(4, m_uncompressed_bytes);

      m_output.write_buffer(trailer, true);
    }
  } catch (const std::ios_base::failure&) {
    std::ostringstream msg;
    msg << "Error writing to FASTQ file " << shell_escape(m_output.filename());

    throw io_error(msg.str(), errno);
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
    std::ostringstream msg;
    msg << "Error closing FASTQ file " << shell_escape(m_output.filename());

    throw io_error(msg.str(), errno);
  }
}

} // namespace adapterremoval
