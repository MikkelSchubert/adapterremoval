// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "fastq_io.hpp"          // declarations
#include "bgzf.hpp"              // for BGZF_EOF, BGZF_MAX_BLOCK_SIZE, ...
#include "buffer.hpp"            // for buffer
#include "commontypes.hpp"       // for output_format
#include "debug.hpp"             // for AR_REQUIRE, AR_REQUIRE_SINGLE_THREAD
#include "deflate.hpp"           // for ISAL_COMPRESSION_LEVEL, crc32, ...
#include "errors.hpp"            // for io_error, gzip_error, fastq_error
#include "fastq.hpp"             // for fastq
#include "fastq_enc.hpp"         // for MATE_SEPARATOR
#include "linereader_joined.hpp" // joined_line_readers
#include "logging.hpp"           // for log
#include "managed_io.hpp"        // for flush
#include "output.hpp"            // for output_file
#include "progress.hpp"          // for progress_timer
#include "scheduler.hpp"         // for provides analytical_step, chunk_ptr, ...
#include "statistics.hpp"        // for fastq_statistics, fastq_stats_ptr, ...
#include "strutils.hpp"          // for shell_escape, string_vec, ends_with
#include "threading.hpp"         // for threadsafe_data
#include "userconfig.hpp"        // for userconfig
#include "utilities.hpp"         // for prng_seed
#include <algorithm>             // for max, min
#include <cerrno>                // for errno
#include <cstddef>               // for size_t
#include <memory>                // for unique_ptr, make_unique, ...
#include <mutex>                 // for mutex, recursive_lock, unique_lock
#include <sstream>               // for basic_ostream, basic_ostringstream, ...
#include <string>                // for string<<
#include <string_view>           // for string_view
#include <utility>               // for move, swap
#include <vector>                // for vector

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// Helper function for isa-l

namespace {

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
read_record(joined_line_readers& reader, fastq& record)
{
  // Line numbers change as we attempt to read the record, and potentially
  // points to the next record in the case of invalid qualities/nucleotides
  const auto start = reader.position();

  try {
    return record.read_unsafe(reader);
  } catch (const fastq_error& error) {
    std::ostringstream stream;
    stream << "Error reading FASTQ record from "
           << reader.filenames(start, reader.position()) << ": "
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

class fastq_vec_cache
{
public:
  fastq_vec_cache() = default;
  fastq_vec_cache(const fastq_vec_cache&) = delete;
  fastq_vec_cache(fastq_vec_cache&&) = delete;
  fastq_vec_cache& operator=(const fastq_vec_cache&) = delete;
  fastq_vec_cache& operator=(fastq_vec_cache&&) = delete;
  ~fastq_vec_cache() = default;

  /** Attempt to acquire a previously used vector of fastq objects */
  void acquire(std::vector<fastq>& chunk, size_t capacity)
  {
    if (capacity) {
      const std::unique_lock lock{ s_chunks_lock };
      if (s_chunks_cache.empty()) {
        s_chunks_created++;
      } else {
        chunk = std::move(s_chunks_cache.back());
        s_chunks_cache.pop_back();
      }
    }

    chunk.resize(capacity);
  }

  /** Release a vector of fastq objects for reuse */
  void release(std::vector<fastq>&& chunk)
  {
    if (chunk.capacity()) {
      const std::unique_lock lock{ s_chunks_lock };
      s_chunks_cache.emplace_back(std::move(chunk));
    }
  }

  /** Verify that there were no leaks or lost chunks */
  void finalize()
  {
    const std::unique_lock lock{ s_chunks_lock };
    AR_REQUIRE(s_chunks_cache.size() == s_chunks_created);
    s_chunks_cache.clear();
    s_chunks_created = 0;
  }

private:
  //! Lock used to control access to cache of previously allocated FASTQ chunks
  std::mutex s_chunks_lock{};
  //! Cache of FASTQ chunks that have previously been returned by `release`
  std::vector<std::vector<fastq>> s_chunks_cache{};
  //! #chunks created, to check that `release` and `acquire` are called 1-to-1
  size_t s_chunks_created = 0;
};

//! Global cache of FASTQ vectors, used to reduce the number of allocs
fastq_vec_cache s_fastq_vec_cache;

} // namespace

read_fastq::read_fastq(const userconfig& config,
                       const step_id next_step,
                       const read_fastq::file_type mode)
  : analytical_step(processing_order::ordered, "read_fastq", next_step)
  , m_reader(select_filenames(config, mode))
  , m_mode(mode)
  , m_head(config.head)
{
}

void
read_fastq::add_steps(scheduler& sch,
                      const userconfig& config,
                      step_id next_step)
{
  const auto add_step = [&sch, &config](auto next, auto type) {
    return sch.add<read_fastq>(config, next, type);
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
read_fastq::signal()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_mode != file_type::read_2);

  chunk_vec chunks;
  if (m_eof) {
    chunks.emplace_back(scheduler_signal::suspend);
    return chunks;
  }

  auto chunk = std::make_unique<fastq_chunk>();

  if (m_mode == file_type::read_1) {
    read_single_end(chunk->reads_1);
  } else if (m_mode == file_type::interleaved) {
    read_interleaved(chunk->reads_1, chunk->reads_2);
  } else {
    AR_FAIL("invalid file_type value");
  }

  // Head must be checked after the first loop, to produce at least one chunk
  m_eof |= !m_head;
  chunk->eof = m_eof;
  chunk->first = m_first;
  m_first = false;

  if (m_eof) {
    chunks.emplace_back(scheduler_signal::suspend);
  }

  chunks.emplace_back(std::move(chunk));
  return chunks;
}

chunk_vec
read_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<fastq_chunk>(data);
  AR_REQUIRE(chunk);
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_mode == file_type::read_2);

  // The mate 2 reader only needs to track whether EOF was reached
  AR_REQUIRE(!m_eof);
  m_eof |= chunk->eof;

  read_single_end(chunk->reads_2);
  if (chunk->reads_1.size() != chunk->reads_2.size()) {
    throw fastq_error("Found unequal number of mate 1 and mate 2 reads; "
                      "input files may be truncated. Please fix before "
                      "continuing.");
  }

  chunk_vec chunks;
  chunks.emplace_back(std::move(chunk));
  return chunks;
}

void
read_fastq::read_single_end(std::vector<fastq>& reads)
{
  AR_REQUIRE(reads.empty());

  const size_t to_read = std::min<size_t>(INPUT_READS, m_head);
  s_fastq_vec_cache.acquire(reads, to_read);

  size_t nread = 0;
  for (; nread < to_read; ++nread) {
    if (!read_record(m_reader, reads.at(nread))) {
      m_eof = true;
      break;
    }
  }

  m_head -= nread;
  reads.resize(nread);
}

void
read_fastq::read_interleaved(std::vector<fastq>& reads_1,
                             std::vector<fastq>& reads_2)
{
  AR_REQUIRE(reads_1.empty() && reads_2.empty());

  const size_t to_read = std::min<size_t>(INPUT_READS, m_head);
  s_fastq_vec_cache.acquire(reads_1, to_read);
  s_fastq_vec_cache.acquire(reads_2, to_read);

  size_t nread = 0;
  for (; nread < to_read; ++nread) {
    if (!read_record(m_reader, reads_1.at(nread))) {
      m_eof = true;
      break;
    } else if (!read_record(m_reader, reads_2.at(nread))) {
      throw fastq_error("Found uneven number of interleaved reads; input files "
                        " may be truncated. Please fix before continuing.");
    }
  }

  m_head -= nread;
  reads_1.resize(nread);
  reads_2.resize(nread);
}

void
read_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);

  s_fastq_vec_cache.finalize();
}

void
read_fastq::release(std::vector<fastq>&& chunk)
{
  s_fastq_vec_cache.release(std::move(chunk));
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'scan_fastq'

scan_fastq::scan_fastq(const userconfig& config,
                       const step_id next_step,
                       statistics& stats)
  : analytical_step(processing_order::ordered, "scan_fastq", next_step)
  , m_mate_separator(config.input_mate_separator)
  , m_mate_separator_identified(config.input_mate_separator)
  , m_pre_trim_poly_x(config.pre_trim_poly_x)
  , m_post_trim_poly_x(config.post_trim_poly_x)
  , m_duplication_enabled(config.report_duplication)
  , m_duplication_1(stats.duplication_1)
  , m_duplication_2(stats.duplication_2)
{
  AR_REQUIRE(m_duplication_1 && m_duplication_2);

  m_2_color_checked = (*m_pre_trim_poly_x.get_reader() != "auto") &&
                      (*m_post_trim_poly_x.get_reader() != "auto");
}

chunk_vec
scan_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<fastq_chunk>(data);
  AR_REQUIRE(chunk);

  if (m_duplication_enabled) {
    for (const auto& read : chunk->reads_1) {
      m_duplication_1->process(read);
    }

    for (const auto& read : chunk->reads_2) {
      m_duplication_2->process(read);
    }
  }

  if (!m_mate_separator_identified) {
    // Mate separators are identified using an already ordered block, in order
    // to reduce the need for locking in the post-processing step
    if (chunk->reads_2.empty()) {
      m_mate_separator = fastq::guess_mate_separator(chunk->reads_1);
    } else {
      m_mate_separator =
        fastq::guess_mate_separator(chunk->reads_1, chunk->reads_2);
    }

    m_mate_separator_identified = true;
  }

  if (!m_2_color_checked) {
    detect_two_color(chunk->reads_1);
    m_2_color_checked = true;
  }

  chunk->mate_separator = m_mate_separator;

  chunk_vec chunks;
  chunks.emplace_back(std::move(chunk));
  return chunks;
}

namespace {

void
update_poly_x_nucleotides(threadsafe_data<std::string>& sink,
                          const std::string& nucleotides)
{
  auto writer = sink.get_writer();
  if (*writer == "auto") {
    *writer = nucleotides;
  }
}

} // namespace

void
scan_fastq::detect_two_color(const std::vector<fastq>& reads)
{
  // Check first available read; input is assumed to be homogeneous
  std::string nucleotides;
  for (const auto& read : reads) {
    if (read.is_two_color()) {
      log::info() << "Automatic trimming of poly-G tails enabled";
      nucleotides = "G";
    }
    break;
  }

  update_poly_x_nucleotides(m_pre_trim_poly_x, nucleotides);
  update_poly_x_nucleotides(m_post_trim_poly_x, nucleotides);
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
                                       const step_id next_step,
                                       statistics& stats)
  : analytical_step(processing_order::unordered,
                    "post_process_fastq",
                    next_step)
  , m_statistics_1(stats.input_1)
  , m_statistics_2(stats.input_2)
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
    const std::unique_lock lock(m_timer_lock);
    m_timer->increment(reads_1.size() + reads_2.size());
  }

  chunk_vec chunks;
  chunks.emplace_back(std::move(chunk));

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
                         const step_id next_step)
  : analytical_step(processing_order::ordered, "split_fastq", next_step)
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
      const auto n = std::min(src.size() - src_offset,
                              BGZF_MAX_INPUT_SIZE - m_buffer.size());
      m_buffer.append(src.data() + src_offset, n);

      src_offset += n;

      if (m_buffer.size() == BGZF_MAX_INPUT_SIZE) {
        auto block = std::make_unique<output_chunk>();

        if (m_isal_stream) {
          m_isal_crc32 = calculate_crc32(m_buffer, m_isal_crc32);
        }

        block->uncompressed_size = m_buffer.size();
        block->buffers.emplace_back(std::move(m_buffer));

        chunks.emplace_back(std::move(block));

        m_buffer = buffer();
        m_buffer.reserve(BGZF_MAX_INPUT_SIZE);
      }
    }
  }

  if (m_eof) {
    auto block = std::make_unique<output_chunk>();
    block->eof = true;

    if (m_isal_stream) {
      m_isal_crc32 = calculate_crc32(m_buffer, m_isal_crc32);
    }

    block->crc32 = m_isal_crc32;
    block->uncompressed_size = m_buffer.size();
    block->buffers.emplace_back(std::move(m_buffer));

    chunks.emplace_back(std::move(block));
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'gzip_split_fastq'

gzip_split_fastq::gzip_split_fastq(const userconfig& config,
                                   const output_file& file,
                                   const step_id next_step)
  : analytical_step(processing_order::unordered, "gzip_split_fastq", next_step)
  , m_isal_stream(is_isal_streaming_enabled(file, config.compression_level))
  , m_format(file.format)
  , m_compression_level(
      m_format == output_format::ubam ? 0 : config.compression_level)
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
    output_buffer.resize(BGZF_MAX_BLOCK_SIZE);
    deflate_block_isal_1(input_buffer, output_buffer, 0, chunk->eof);
  } else {
    bgzf_block_compress(input_buffer, output_buffer, m_compression_level);

    if (chunk->eof) {
      // May cause a realloc, but only on the very last chunk, so doesn't matter
      output_buffer.append(BGZF_EOF);
    }
  }

  // Enable reuse of the analytical_chunks
  std::swap(input_buffer, output_buffer);

  chunk_vec chunks;
  chunks.emplace_back(std::move(chunk));

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'write_fastq'

write_fastq::write_fastq(const userconfig& config, const output_file& file)
  // Allow disk IO and writing to STDOUT at the same time
  : analytical_step(processing_order::ordered_io,
                    "write_fastq",
                    step_id::UNDEFINED)
  , m_output(file.name)
  , m_isal_stream(is_isal_streaming_enabled(file, config.compression_level))
{
  if (m_isal_stream) {
    auto header = gzip_stream_header_isal();
    m_output.write(header);
  }
}

chunk_vec
write_fastq::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<output_chunk>(data);
  AR_REQUIRE(chunk);
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(!m_eof);

  m_eof = chunk->eof;
  m_uncompressed_bytes += chunk->uncompressed_size;

  if (m_eof && m_isal_stream) {
    auto tail = gzip_stream_tail_isal(chunk->crc32, m_uncompressed_bytes);
    chunk->buffers.emplace_back(std::move(tail));
  }

  try {
    m_output.write(chunk->buffers);

    // Close file to trigger any exceptions due to fclose, etc. This also forces
    // creation of the file, if no non-empty buffers have been written.
    if (m_eof) {
      m_output.finalize();
    }
  } catch (const io_error&) {
    std::ostringstream msg;
    msg << "Error writing to FASTQ file " << log_escape(m_output.filename());

    throw io_error(msg.str(), errno);
  }

  return {};
}

void
write_fastq::finalize()
{
  AR_REQUIRE_SINGLE_THREAD(m_lock);
  AR_REQUIRE(m_eof);
}

} // namespace adapterremoval
