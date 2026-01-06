// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "fastq_io.hpp"    // declarations
#include "commontypes.hpp" // for output_format
#include "debug.hpp"       // for AR_REQUIRE, AR_REQUIRE_SINGLE_THREAD
#include "errors.hpp"      // for io_error, gzip_error, fastq_error
#include "fastq.hpp"       // for fastq
#include "fastq_enc.hpp"   // for MATE_SEPARATOR
#include "output.hpp"      // for output_file
#include "progress.hpp"    // for progress_timer
#include "scheduler.hpp"   // for provides analytical_step, chunk_ptr, ...
#include "statistics.hpp"  // for fastq_statistics, fastq_stats_ptr, stat...
#include "strutils.hpp"    // for shell_escape, string_vec, ends_with
#include "userconfig.hpp"  // for userconfig
#include "utilities.hpp"   // for prng_seed
#include <algorithm>       // for max, min
#include <cerrno>          // for errno
#include <isa-l.h>         // for isal_gzip_header
#include <libdeflate.h>    // for libdeflate_alloc_compressor, libdeflate...
#include <memory>          // for unique_ptr, make_unique, __shared_ptr_a...
#include <sstream>         // for basic_ostream, basic_ostringstream, ope...
#include <string>          // for string<<
#include <string_view>     // for string_view
#include <utility>         // for move, swap

namespace adapterremoval {

// Default bgzip header, as described in the SAM spec. v1.6 section 4.1.
// Includes 2 trailing placeholder bytes for total block size (BSIZE)
constexpr std::string_view BGZF_HEADER = {
  "\37\213\10\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0",
  18
};

// Eof of file marker for bgzip files; see SAM spec. v1.6 section 4.1.2
constexpr std::string_view BGZF_EOF = {
  "\37\213\10\4\0\0\0\0\0\377\6\0\102\103\2\0\33\0\3\0\0\0\0\0\0\0\0\0",
  28,
};

////////////////////////////////////////////////////////////////////////////////
// Helper function for isa-l

namespace {

//! The compression level used for block/stream compression with isa-l
constexpr size_t ISAL_COMPRESSION_LEVEL = 1;
//! The default buffer size for compression at level ISAL_COMPRESSION_LEVEL
constexpr size_t ISAL_BUFFER_SIZE = ISAL_DEF_LVL1_SMALL;

/**
 * ISA-l streaming is enabled only at compression level 1, since little
 * difference was observed between levels 1 to 3. However, it still offers a
 * faster compression (with a lower ratio) than libdeflate level 1.
 */
bool
is_isal_streaming_enabled(const output_file file, unsigned compression_level)
{
  switch (file.format) {
    case output_format::fastq:
    case output_format::sam:
    case output_format::bam:
    case output_format::ubam:
      return false;
    case output_format::fastq_gzip:
    case output_format::sam_gzip:
      return compression_level == ISAL_COMPRESSION_LEVEL;
    default:
      AR_FAIL("invalid output format");
  }
}

} // namespace

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
read_record(joined_line_readers& reader, std::vector<fastq>& chunk)
{
  // Line numbers change as we attempt to read the record, and potentially
  // points to the next record in the case of invalid qualities/nucleotides
  const auto start = reader.linenumber();

  try {
    chunk.emplace_back();
    auto& record = chunk.back();

    if (record.read_unsafe(reader)) {
      return true;
    } else {
      chunk.pop_back();
      return false;
    }
  } catch (const fastq_error& error) {
    std::ostringstream stream;
    stream << "Error reading FASTQ record from "
           << reader.filenames(start, reader.linenumber()) << ": "
           << error.what();

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

/** Estimates an upper bound for the required capacity for gzip compression */
constexpr size_t
estimate_capacity(size_t input_size, bool eof)
{
  return input_size                    //
         + BGZF_HEADER.size()          // Standard bgzip header
         + 1                           // BFINAL | BTYPE
         + 4                           // LEN + NLEN
         + 4                           // CRC32
         + 4                           // ISIZE
         + (eof ? BGZF_EOF.size() : 0) // Standard bgzip tail
    ;
}

} // namespace

read_fastq::read_fastq(const userconfig& config,
                       const size_t next_step,
                       const read_fastq::file_type mode,
                       statistics& stats)
  : analytical_step(processing_order::ordered, "read_fastq")
  , m_reader(select_filenames(config, mode))
  , m_next_step(next_step)
  , m_single_end(!config.paired_ended_mode)
  , m_mode(mode)
  , m_head(config.head)
  , m_mate_separator(config.input_mate_separator)
  , m_mate_separator_identified(config.input_mate_separator)
  , m_duplication_1(stats.duplication_1)
  , m_duplication_2(stats.duplication_2)
{
  AR_REQUIRE(m_duplication_1 && m_duplication_2);
}

void
read_fastq::add_steps(scheduler& sch,
                      const userconfig& config,
                      size_t next_step,
                      statistics& stats)
{
  const auto add_step = [&sch, &config, &stats](auto next, auto type) {
    return sch.add<read_fastq>(config, next, type, stats);
  };

  if (config.interleaved_input) {
    add_step(next_step, read_fastq::file_type::interleaved);
  } else if (config.paired_ended_mode) {
    next_step = add_step(next_step, read_fastq::file_type::read_2);
    add_step(next_step, read_fastq::file_type::read_1);
  } else {
    add_step(next_step, read_fastq::file_type::read_1);
  }
}

chunk_vec
read_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<fastq_chunk>(data);
  AR_REQUIRE_SINGLE_THREAD(m_lock);

  if (m_mode == file_type::read_1 || m_mode == file_type::interleaved) {
    // The scheduler only terminates when the first step stops returning chunks
    if (m_eof) {
      return {};
    }

    chunk = std::make_unique<fastq_chunk>();
  } else {
    AR_REQUIRE(!m_eof && chunk);
  }

  auto& reads_1 = chunk->reads_1;
  auto& reads_2 = chunk->reads_2;

  if (m_mode == file_type::read_1 || m_mode == file_type::interleaved) {
    if (m_mode == file_type::read_1) {
      read_single_end(reads_1, *m_duplication_1);
    } else {
      read_interleaved(reads_1, reads_2);
    }
  } else if (m_mode == file_type::read_2) {
    read_single_end(reads_2, *m_duplication_2);
  } else {
    AR_FAIL("invalid file_type value");
  }

  if (m_single_end) {
    AR_REQUIRE(m_mode == file_type::read_1);

    if (!m_mate_separator_identified) {
      // Attempt to determine the mate separator character
      m_mate_separator = fastq::guess_mate_separator(reads_1);
      m_mate_separator_identified = true;
    }
  } else if (m_mode != file_type::read_1) {
    if (reads_1.size() != reads_2.size()) {
      throw fastq_error("Found unequal number of mate 1 and mate 2 reads; "
                        "input files may be truncated. Please fix before "
                        "continuing.");
    } else if (!m_mate_separator_identified) {
      AR_REQUIRE(reads_1.size() == reads_2.size());
      // Mate separators are identified using the first block, in order to
      // reduce the need for locking in the post-processing step

      // Attempt to determine the mate separator character
      m_mate_separator = fastq::guess_mate_separator(reads_1, reads_2);
      m_mate_separator_identified = true;
    }
  }

  // Head must be checked after the first loop, to produce at least one chunk
  m_eof |= !m_head;
  chunk->eof = m_eof;
  chunk->mate_separator = m_mate_separator;
  chunk->first = m_first;
  m_first = false;

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));
  return chunks;
}

void
read_fastq::read_single_end(std::vector<fastq>& reads,
                            duplication_statistics& stats)
{
  for (; reads.size() < INPUT_READS && m_head && !m_eof; m_head--) {
    if (read_record(m_reader, reads)) {
      stats.process(reads.back());
    } else {
      m_eof = true;
    }
  }
}

void
read_fastq::read_interleaved(std::vector<fastq>& reads_1,
                             std::vector<fastq>& reads_2)
{
  for (; reads_1.size() < INPUT_READS && m_head && !m_eof; m_head--) {
    if (read_record(m_reader, reads_1)) {
      m_duplication_1->process(reads_1.back());
    } else {
      m_eof = true;
      break;
    }

    if (read_record(m_reader, reads_2)) {
      m_duplication_2->process(reads_2.back());
    }
  }
}

void
read_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'post_process_fastq_stats'

post_process_fastq::stats_pair::stats_pair(double sample_rate,
                                           unsigned int seed)
  : stats_1(std::make_unique<fastq_statistics>(sample_rate, seed))
  , stats_2(std::make_unique<fastq_statistics>(sample_rate, seed))
{
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
  , m_timer(std::make_unique<progress_timer>(config.log_progress))
{
  AR_REQUIRE(m_statistics_1 && m_statistics_2);

  for (size_t i = 0; i < config.max_threads; ++i) {
    m_stats.emplace_back(config.report_sample_rate, prng_seed());
  }
}

// Ensure that progress_timer definition is available to unique_ptr destructor
post_process_fastq::~post_process_fastq() = default;

chunk_vec
post_process_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<fastq_chunk>(data);
  AR_REQUIRE(chunk);
  auto& reads_1 = chunk->reads_1;
  auto& reads_2 = chunk->reads_2;

  auto stats = m_stats.acquire();

  AR_REQUIRE((reads_1.size() == reads_2.size()) || reads_2.empty());
  if (reads_2.empty()) {
    for (auto& read_1 : reads_1) {
      read_1.post_process(m_encoding);
      stats->stats_1->process(read_1);
    }
  } else {
    auto it_1 = reads_1.begin();
    auto it_2 = reads_2.begin();
    for (; it_1 != reads_1.end(); ++it_1, ++it_2) {
      fastq::validate_paired_reads(*it_1, *it_2, chunk->mate_separator);

      it_1->post_process(m_encoding);
      stats->stats_1->process(*it_1);

      it_2->post_process(m_encoding);
      stats->stats_2->process(*it_2);
    }
  }

  m_stats.release(std::move(stats));

  {
    std::unique_lock<std::mutex> lock(m_timer_lock);
    m_timer->increment(reads_1.size() + reads_2.size());
  }

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

void
post_process_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_timer_lock);

  while (auto it = m_stats.try_acquire()) {
    *m_statistics_1 += *it->stats_1;
    *m_statistics_2 += *it->stats_2;
  }

  m_timer->finalize();
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'split_fastq'

split_fastq::split_fastq(const userconfig& config,
                         const output_file& file,
                         size_t next_step)
  : analytical_step(processing_order::ordered, "split_fastq")
  , m_next_step(next_step)
  , m_isal_stream(is_isal_streaming_enabled(file, config.compression_level))
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
split_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<output_chunk>(data);
  AR_REQUIRE(chunk);
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!m_eof);
  m_eof = chunk->eof;

  chunk_vec chunks;
  for (const auto& src : chunk->buffers) {
    for (size_t src_offset = 0; src_offset < src.size();) {
      const auto n =
        std::min(src.size() - src_offset, BGZF_BLOCK_SIZE - m_buffer.size());
      m_buffer.append(src.data() + src_offset, n);

      src_offset += n;

      if (m_buffer.size() == BGZF_BLOCK_SIZE) {
        auto block = std::make_unique<output_chunk>();

        if (m_isal_stream) {
          m_isal_crc32 =
            libdeflate_crc32(m_isal_crc32, m_buffer.data(), m_buffer.size());
        }

        block->uncompressed_size = m_buffer.size();
        block->buffers.emplace_back(std::move(m_buffer));

        chunks.emplace_back(m_next_step, std::move(block));

        m_buffer = buffer();
        m_buffer.reserve(BGZF_BLOCK_SIZE);
      }
    }
  }

  if (m_eof) {
    auto block = std::make_unique<output_chunk>();
    block->eof = true;

    if (m_isal_stream) {
      m_isal_crc32 =
        libdeflate_crc32(m_isal_crc32, m_buffer.data(), m_buffer.size());
    }

    block->crc32 = m_isal_crc32;
    block->uncompressed_size = m_buffer.size();
    block->buffers.emplace_back(std::move(m_buffer));

    chunks.emplace_back(m_next_step, std::move(block));
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_split_fastq'

namespace {

size_t
isal_deflate_block(buffer& input_buffer,
                   buffer& output_buffer,
                   const size_t output_offset,
                   const bool eof)
{
  isal_zstream stream{};
  isal_deflate_stateless_init(&stream);

  stream.flush = FULL_FLUSH;
  stream.end_of_stream = eof;

  stream.level = ISAL_COMPRESSION_LEVEL;
  stream.level_buf_size = ISAL_BUFFER_SIZE;
  buffer level_buffer{ stream.level_buf_size };
  stream.level_buf = level_buffer.data();

  stream.avail_in = input_buffer.size();
  stream.next_in = input_buffer.data();
  stream.next_out = output_buffer.data() + output_offset;
  stream.avail_out = output_buffer.size() - output_offset;

  const auto ec = isal_deflate_stateless(&stream);
  switch (ec) {
    case COMP_OK:
      break;
    case INVALID_FLUSH:
      throw gzip_error("isal_deflate_stateless: invalid flush");
    case ISAL_INVALID_LEVEL:
      throw gzip_error("isal_deflate_stateless: invalid level");
    case ISAL_INVALID_LEVEL_BUF:
      throw gzip_error("isal_deflate_stateless: invalid buffer size");
    default: {
      std::ostringstream os;
      os << "isal_deflate_stateless: unknown error " << ec;

      throw gzip_error(os.str());
    }
  }

  // The easily compressible input should fit in a single output block
  AR_REQUIRE(stream.avail_in == 0);

  return stream.total_out;
}

} // namespace

gzip_split_fastq::gzip_split_fastq(const userconfig& config,
                                   const output_file& file,
                                   size_t next_step)
  : analytical_step(processing_order::unordered, "gzip_split_fastq")
  , m_config(config)
  , m_isal_stream(is_isal_streaming_enabled(file, config.compression_level))
  , m_format(file.format)
  , m_next_step(next_step)
{
}

chunk_vec
gzip_split_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<output_chunk>(data);
  AR_REQUIRE(chunk);
  AR_REQUIRE(chunk->buffers.size() == 1);

  buffer& input_buffer = chunk->buffers.front();
  buffer output_buffer;

  if (m_isal_stream) {
    output_buffer.resize(estimate_capacity(input_buffer.size(), false));
    const auto output_size =
      isal_deflate_block(input_buffer, output_buffer, 0, chunk->eof);
    output_buffer.resize(output_size);
  } else {
    if (m_format == output_format::ubam || m_config.compression_level == 0) {
      output_buffer.reserve(estimate_capacity(input_buffer.size(), chunk->eof));
      output_buffer.append(BGZF_HEADER);
      output_buffer.append_u8(1); // BFINAL=1, BTYPE=00; see RFC1951
      output_buffer.append_u16(input_buffer.size());
      output_buffer.append_u16(~input_buffer.size());
      output_buffer.append(input_buffer);
    } else if (m_config.compression_level == ISAL_COMPRESSION_LEVEL) {
      output_buffer.reserve(estimate_capacity(input_buffer.size(), chunk->eof));
      output_buffer.append(BGZF_HEADER);
      output_buffer.resize(output_buffer.capacity());

      const auto output_size = isal_deflate_block(input_buffer,
                                                  output_buffer,
                                                  BGZF_HEADER.size(),
                                                  true);

      // Resize the buffer to the actually used size
      output_buffer.resize(output_size + BGZF_HEADER.size());
    } else {
      // Libdeflate compression levels 1 to 12 are mapped onto 2 to 13
      AR_REQUIRE(m_config.compression_level >= 2 &&
                 m_config.compression_level <= 12);
      auto* compressor =
        libdeflate_alloc_compressor(m_config.compression_level);
      const auto output_bound =
        libdeflate_deflate_compress_bound(compressor, input_buffer.size());

      output_buffer.reserve(estimate_capacity(output_bound, chunk->eof));
      output_buffer.append(BGZF_HEADER);
      output_buffer.resize(output_buffer.capacity());

      const auto output_size =
        libdeflate_deflate_compress(compressor,
                                    input_buffer.data(),
                                    input_buffer.size(),
                                    output_buffer.data() + BGZF_HEADER.size(),
                                    output_buffer.size() - BGZF_HEADER.size());
      libdeflate_free_compressor(compressor);
      // The easily compressible input should fit in a single output block
      AR_REQUIRE(output_size);

      // Resize the buffer to the actually used size
      output_buffer.resize(output_size + BGZF_HEADER.size());
    }

    const auto input_crc32 =
      libdeflate_crc32(0, input_buffer.data(), input_buffer.size());
    output_buffer.append_u32(input_crc32);         // checksum of data
    output_buffer.append_u32(input_buffer.size()); // size of data

    AR_REQUIRE(output_buffer.size() <= BGZF_MAX_BLOCK_SIZE);
    // Write the final block size; -1 to fit 65536 in 16 bit
    output_buffer.put_u16(16, output_buffer.size() - 1);

    if (chunk->eof) {
      output_buffer.append(BGZF_EOF);
    }
  }

  // Enable reuse of the analytical_chunks
  std::swap(input_buffer, output_buffer);

  chunk_vec chunks;
  chunks.emplace_back(m_next_step, std::move(chunk));

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

write_fastq::write_fastq(const userconfig& config, const output_file& file)
  // Allow disk IO and writing to STDOUT at the same time
  : analytical_step(processing_order::ordered_io, "write_fastq")
  , m_output(file.name)
  , m_isal_stream(is_isal_streaming_enabled(file, config.compression_level))
{
  if (m_isal_stream) {
    buffer level_buf(ISAL_BUFFER_SIZE);
    buffer output_buf(OUTPUT_BLOCK_SIZE);

    struct isal_zstream stream = {};
    struct isal_gzip_header header = {};

    isal_gzip_header_init(&header);
    isal_deflate_init(&stream);
    stream.avail_in = 0;
    stream.flush = NO_FLUSH;
    stream.level = ISAL_COMPRESSION_LEVEL;
    stream.level_buf = level_buf.data();
    stream.level_buf_size = level_buf.size();
    stream.gzip_flag = IGZIP_GZIP_NO_HDR;
    stream.next_out = output_buf.data();
    stream.avail_out = output_buf.size();

    const auto ret = isal_write_gzip_header(&stream, &header);
    AR_REQUIRE(ret == 0, "buffer was not large enough for gzip header");

    output_buf.resize(stream.total_out);

    m_output.write(output_buf);
  }
}

chunk_vec
write_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<output_chunk>(data);
  AR_REQUIRE(chunk);
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!m_eof);

  try {
    m_eof = chunk->eof;
    m_uncompressed_bytes += chunk->uncompressed_size;

    const auto mode = (m_eof && !m_isal_stream) ? flush::on : flush::off;

    m_output.write(chunk->buffers, mode);

    if (m_eof && m_isal_stream) {
      buffer trailer;
      trailer.append_u32(chunk->crc32);
      trailer.append_u32(m_uncompressed_bytes);

      m_output.write(trailer, flush::on);
    }
  } catch (const io_error&) {
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
  } catch (const io_error&) {
    std::ostringstream msg;
    msg << "Error closing FASTQ file " << shell_escape(m_output.filename());

    throw io_error(msg.str(), errno);
  }
}

} // namespace adapterremoval
