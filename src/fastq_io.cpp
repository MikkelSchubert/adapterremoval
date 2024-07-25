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

#include "fastq_io.hpp"      // declarations
#include "debug.hpp"         // for AR_REQUIRE, AR_REQUIRE_SINGLE_THREAD
#include "errors.hpp"        // for io_error, gzip_error, fastq_error
#include "fastq.hpp"         // for fastq
#include "fastq_enc.hpp"     // for MATE_SEPARATOR
#include "simd.hpp"          // for size_t
#include "statistics.hpp"    // for fastq_statistics, fastq_stats_ptr, stat...
#include "strutils.hpp"      // for shell_escape, string_vec, ends_with
#include "userconfig.hpp"    // for userconfig
#include <algorithm>         // for max, min
#include <cerrno>            // for errno
#include <cstring>           // for size_t, memcpy
#include <isa-l/crc.h>       // for crc32_gzip_refl
#include <isa-l/igzip_lib.h> // for isal_zstream, isal_deflate_init, isal_d...
#include <libdeflate.h>      // for libdeflate_alloc_compressor, libdeflate...
#include <memory>            // for unique_ptr, make_unique, __shared_ptr_a...
#include <sstream>           // for basic_ostream, basic_ostringstream, ope...
#include <utility>           // for move, swap

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
  return (config.gzip || ends_with(to_lower(filename), ".gz")) &&
         config.gzip_level <= MAX_ISAL_LEVEL;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_read_chunk'

fastq_read_chunk::fastq_read_chunk(bool eof_)
  : eof(eof_)
{
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_output_chunk'

fastq_output_chunk::fastq_output_chunk(bool eof_, uint32_t crc32_)
  : eof(eof_)
  , crc32(crc32_)
{
}

void
fastq_output_chunk::add(const fastq& read)
{
  nucleotides += read.length();
  read.into_string(reads);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_fastq'

enum class read_fastq::file_type
{
  read_1,
  read_2,
  interleaved
};

namespace {

bool
read_record(joined_line_readers& reader, fastq_vec& chunk, size_t& nucleotides)
{
  // Line numbers change as we attempt to read the record, and potentially
  // points to the next record in the case of invalid qualities/nucleotides
  const auto line_number = reader.linenumber();

  try {
    chunk.emplace_back();
    auto& record = chunk.back();

    if (record.read_unsafe(reader)) {
      nucleotides += record.length();

      return true;
    } else {
      chunk.pop_back();
      return false;
    }
  } catch (const fastq_error& error) {
    std::ostringstream stream;
    stream << "Error reading FASTQ record from '" << reader.filename()
           << "' at line " << line_number << "; aborting:\n"
           << indent_lines(error.what());

    throw fastq_error(stream.str());
  }
}

const string_vec&
select_filenames(const userconfig& config, const read_fastq::file_type mode)
{
  switch (mode) {
    case read_fastq::file_type::read_1:
    case read_fastq::file_type::interleaved:
      return config.input_files_1;
    case read_fastq::file_type::read_2:
      return config.input_files_2;
    default:
      AR_FAIL("invalid read_fastq::file_type value");
  }
}

char
identify_mate_separators(const fastq_vec& reads_1, const fastq_vec& reads_2)
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

} // namespace

read_fastq::read_fastq(const userconfig& config,
                       const size_t next_step,
                       const read_fastq::file_type mode)
  : analytical_step(processing_order::ordered, "read_fastq")
  , m_reader(select_filenames(config, mode))
  , m_next_step(next_step)
  , m_mode(mode)
  , m_head(config.head)
  , m_mate_separator(config.mate_separator)
{
}

void
read_fastq::add_steps(scheduler& sch,
                      const userconfig& config,
                      size_t next_step)
{
  if (config.interleaved_input) {
    sch.add<read_fastq>(config, next_step, read_fastq::file_type::interleaved);
  } else if (config.paired_ended_mode) {
    next_step =
      sch.add<read_fastq>(config, next_step, read_fastq::file_type::read_2);
    sch.add<read_fastq>(config, next_step, read_fastq::file_type::read_1);
  } else {
    sch.add<read_fastq>(config, next_step, read_fastq::file_type::read_1);
  }
}

chunk_vec
read_fastq::process(chunk_ptr chunk)
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);

  if (m_mode == file_type::read_1 || m_mode == file_type::interleaved) {
    // The scheduler only terminates when the first step stops returning chunks
    if (m_eof) {
      return {};
    }

    chunk = std::make_unique<fastq_read_chunk>();
  } else {
    AR_REQUIRE(!m_eof && chunk);
  }

  auto& file_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  auto& reads_1 = file_chunk.reads_1;
  auto& reads_2 = file_chunk.reads_2;

  if (m_mode == file_type::read_1 || m_mode == file_type::interleaved) {
    if (m_mode == file_type::read_1) {
      read_single_end(reads_1);
    } else {
      read_interleaved(reads_1, reads_2);
    }
  } else if (m_mode == file_type::read_2) {
    read_single_end(reads_2);
  } else {
    AR_FAIL("invalid file_type value");
  }

  if (m_mode != file_type::read_1) {
    if (reads_1.size() != reads_2.size()) {
      throw fastq_error("Found unequal number of mate 1 and mate 2 reads; "
                        "input files may be truncated. Please fix before "
                        "continuing.");
    } else if (!m_mate_separator) {
      // Mate separators are identified using the first block, in order to
      // reduce the need for locking in the post-processing step
      m_mate_separator = identify_mate_separators(reads_1, reads_2);
    }
  }

  // Head must be checked after the first loop, to produce at least one chunk
  m_eof |= !m_head;
  file_chunk.eof = m_eof;
  file_chunk.mate_separator = m_mate_separator;

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));
  return chunks;
}

void
read_fastq::read_single_end(fastq_vec& reads)
{
  size_t nucleotides = 0;
  for (; nucleotides < INPUT_BLOCK_SIZE && m_head && !m_eof; m_head--) {
    m_eof = !read_record(m_reader, reads, nucleotides);
  }
}

void
read_fastq::read_interleaved(fastq_vec& reads_1, fastq_vec& reads_2)
{
  size_t nucleotides = 0;
  for (; nucleotides < INPUT_BLOCK_SIZE && m_head && !m_eof; m_head--) {
    m_eof = !read_record(m_reader, reads_1, nucleotides);
    read_record(m_reader, reads_2, nucleotides);
  }
}

void
read_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'post_process_fastq'

post_process_fastq::post_process_fastq(const userconfig& config,
                                       size_t next_step,
                                       statistics& stats)
  : analytical_step(processing_order::unordered, "post_process_fastq")
  , m_statistics_1(stats.input_1)
  , m_statistics_2(stats.input_2)
  , m_next_step(next_step)
  , m_encoding(config.io_encoding)
  , m_timer(config.log_progress)
{
  AR_REQUIRE(m_statistics_1 && m_statistics_2);

  for (size_t i = 0; i < config.max_threads; ++i) {
    // Synchronize sampling of mate 1 and mate 2 reads
    const auto seed = prng_seed();

    m_stats_1.emplace_back(config.report_sample_rate, seed);
    m_stats_2.emplace_back(config.report_sample_rate, seed);
  }
}

chunk_vec
post_process_fastq::process(chunk_ptr chunk)
{
  auto& file_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  auto& reads_1 = file_chunk.reads_1;
  auto& reads_2 = file_chunk.reads_2;

  auto stats_1 = m_stats_1.acquire();
  auto stats_2 = m_stats_2.acquire();

  AR_REQUIRE((reads_1.size() == reads_2.size()) || reads_2.empty());
  if (reads_2.empty()) {
    for (auto& read_1 : reads_1) {
      read_1.post_process(m_encoding);
      stats_1->process(read_1);
    }
  } else {
    auto it_1 = reads_1.begin();
    auto it_2 = reads_2.begin();
    for (; it_1 != reads_1.end(); ++it_1, ++it_2) {
      fastq::normalize_paired_reads(*it_1, *it_2, file_chunk.mate_separator);

      it_1->post_process(m_encoding);
      stats_1->process(*it_1);

      it_2->post_process(m_encoding);
      stats_2->process(*it_2);
    }
  }

  m_stats_1.release(stats_1);
  m_stats_2.release(stats_2);

  {
    std::unique_lock<std::mutex> lock(m_timer_lock);
    m_timer.increment(reads_1.size() + reads_2.size());
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

void
post_process_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_timer_lock);

  while (auto next = m_stats_1.try_acquire()) {
    *m_statistics_1 += *next;
  }

  while (auto next = m_stats_2.try_acquire()) {
    *m_statistics_2 += *next;
  }

  m_timer.finalize();
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'split_fastq'

split_fastq::split_fastq(const userconfig& config,
                         const std::string& filename,
                         size_t next_step)
  : analytical_step(processing_order::ordered, "split_fastq")
  , m_next_step(next_step)
  , m_isal_enabled(isal_enabled(config, filename))
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
    isal_zstream stream{};
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
        throw gzip_error("isal_deflate_stateless: invalid flush");
      case ISAL_INVALID_LEVEL:
        throw gzip_error("isal_deflate_stateless: invalid level");
      case ISAL_INVALID_LEVEL_BUF:
        throw gzip_error("isal_deflate_stateless: invalid buffer size");
      default:
        throw gzip_error("isal_deflate_stateless: unexpected error");
    }

    // The easily compressible input should fit in a single output block
    AR_REQUIRE(stream.avail_in == 0);

    compressed_size = stream.total_out;
  } else {
    auto* compressor = libdeflate_alloc_compressor(m_config.gzip_level);
    compressed_size = libdeflate_gzip_compress(compressor,
                                               input_buffer.get(),
                                               input_buffer.size(),
                                               output_buffer.get(),
                                               output_buffer.size());
    libdeflate_free_compressor(compressor);

    // The easily compressible input should fit in a single output block
    AR_REQUIRE(compressed_size);
  }

  output_buffer.resize(compressed_size);

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

write_fastq::write_fastq(const userconfig& config, const std::string& filename)
  // Allow disk IO and writing to STDOUT at the same time
  : analytical_step(processing_order::ordered_io, "write_fastq")
  , m_output(filename)
  , m_isal_enabled(isal_enabled(config, filename))
{
  if (m_isal_enabled) {
    buffer level_buf(isal_buffer_size(config.gzip_level));
    buffer output_buf(OUTPUT_BLOCK_SIZE);

    struct isal_zstream stream = {};
    struct isal_gzip_header header = {};

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

    m_output.write(output_buf);
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

    const auto mode = (m_eof && !m_isal_enabled) ? flush::on : flush::off;

    if (file_chunk.buffers.empty()) {
      m_output.write(file_chunk.reads, mode);
    } else {
      AR_REQUIRE(file_chunk.reads.empty());

      m_output.write(file_chunk.buffers, mode);
    }

    if (m_eof && m_isal_enabled) {
      buffer trailer(8);
      trailer.write_u32(0, file_chunk.crc32);
      trailer.write_u32(4, m_uncompressed_bytes);

      m_output.write(trailer, flush::on);
    }
  } catch (const std::ios_base::failure&) {
    std::ostringstream msg;
    msg << "Error writing to FASTQ file " << shell_escape(m_output.filename());

    throw io_error(msg.str(), errno);
  }

  return {};
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
