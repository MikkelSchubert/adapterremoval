// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "buffer.hpp"            // for buffer, buffer_vec
#include "commontypes.hpp"       // for fastq_vec
#include "fastq.hpp"             // for fastq
#include "fastq_enc.hpp"         // for fastq_encoding
#include "linereader_joined.hpp" // for joined_line_readers
#include "managed_io.hpp"        // for managed_writer
#include "progress.hpp"          // for progress_timer
#include "scheduler.hpp"         // for analytical_step, chunk_ptr, chunk_vec
#include "statistics.hpp"        // for fastq_stats_ptr
#include <cstddef>               // for size_t
#include <cstdint>               // for uint32_t, uint64_t
#include <limits>                // for numeric_limits
#include <memory>                // for unique_ptr
#include <mutex>                 // for mutex
#include <string>                // for string

namespace adapterremoval {

class analytical_chunk;
class userconfig;
struct output_file;

using chunk_ptr = std::unique_ptr<analytical_chunk>;

//! Number of reads to load every cycle
const size_t INPUT_READS = 1024;
//! Size of chunks of when performing block compression
const size_t BGZF_BLOCK_SIZE = 0xff00;
//! Maximum size of BGZF blocks including header and tail
const size_t BGZF_MAX_BLOCK_SIZE = 0x10000;
//! Size of blocks to generate before writing to output
const size_t OUTPUT_BLOCK_SIZE = 4LLU * 64 * 1024;

/**
 * Simple file reading step.
 *
 * Reads from the mate 1 and the mate 2 files, storing the reads in a
 * fastq_file_chunk. Once the EOF has been reached, a single empty chunk will
 * be returned, marked using the 'eof' property.
 */
class read_fastq : public analytical_step
{
public:
  enum class file_type;

  read_fastq(const userconfig& config, size_t next_step, file_type mode);
  ~read_fastq() override = default;

  /** Adds read_fastq step(s) to the scheduler depending on current input */
  static void add_steps(scheduler& sch,
                        const userconfig& config,
                        size_t next_step);

  /** Reads lines from the input file and saves them in an fastq_file_chunk. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  read_fastq(const read_fastq&) = delete;
  read_fastq(read_fastq&&) = delete;
  read_fastq& operator=(const read_fastq&) = delete;
  read_fastq& operator=(read_fastq&&) = delete;

private:
  /** Fills a chunk with SE reads; stops on EOF or after `head` reads. */
  void read_single_end(fastq_vec& reads);
  /** Fills a chunk of interleaved reads; stops on EOF or after `head` reads. */
  void read_interleaved(fastq_vec& reads_1, fastq_vec& reads_2);

  //! The underlying file reader for mate 1 (and possibly mate 2) reads
  joined_line_readers m_reader;

  const size_t m_next_step;

  //! The kind of data read by the reader (SE, PE, interleaved)
  file_type m_mode;
  //! Used to track whether the first block has been read
  bool m_first = true;
  //! Used to track whether an EOF block has been received.
  bool m_eof = false;
  //! Number of reads to process
  uint64_t m_head = std::numeric_limits<uint64_t>::max();
  //! Character used to join read-names with mate numbers, e.g. '/'
  char m_mate_separator;
  //! Indicates if the mate separator is known / has been attempted identififed
  bool m_mate_separator_identified;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock{};
};

/**
 * Class responsible for collecting statistics on input FASTQ. This is split
 * into a separate step to minimize the amount of time that IO is blocked by the
 * FASTQ reading step.
 */
class post_process_fastq : public analytical_step
{
public:
  /** Constructor. */
  post_process_fastq(const userconfig& config,
                     size_t next_step,
                     statistics& stats);

  ~post_process_fastq() override = default;

  /** Reads lines from the input file and saves them in an fastq_file_chunk. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  post_process_fastq(const post_process_fastq&) = delete;
  post_process_fastq(post_process_fastq&&) = delete;
  post_process_fastq& operator=(const post_process_fastq&) = delete;
  post_process_fastq& operator=(post_process_fastq&&) = delete;

private:
  /** Linked mate 1/2 statistics to synchronize sampling of reads */
  struct stats_pair
  {
    stats_pair(double sample_rate, uint32_t seed)
      : stats_1(sample_rate, seed)
      , stats_2(sample_rate, seed)
    {
    }

    fastq_statistics stats_1;
    fastq_statistics stats_2;
  };

  //! Per thread statistics collected from raw reads
  threadstate<stats_pair> m_stats{};
  //! Destination for statistics collected from raw mate 1 reads
  fastq_stats_ptr m_statistics_1;
  //! Destination for statistics collected from raw mate 2 reads
  fastq_stats_ptr m_statistics_2;
  //! The analytical step following this step
  const size_t m_next_step;
  //! Encoding used to parse FASTQ reads.
  const fastq_encoding m_encoding;

  //! Lock used to control access to progress timer
  std::mutex m_timer_lock{};
  //! Timer for displaying read progress.
  progress_timer m_timer;
};

/**
 * Splits input into chunks that can be GZipped in parallel.
 */
class split_fastq : public analytical_step
{
public:
  /** Constructor; 'next_step' sets the destination of compressed chunks. */
  split_fastq(const userconfig& config,
              const output_file& file,
              size_t next_step);

  ~split_fastq() override = default;

  chunk_vec process(chunk_ptr chunk) override;
  void finalize() override;

  split_fastq(const split_fastq&) = delete;
  split_fastq(split_fastq&&) = delete;
  split_fastq& operator=(const split_fastq&) = delete;
  split_fastq& operator=(split_fastq&&) = delete;

private:
  //! The analytical step following this step
  const size_t m_next_step;
  //! Buffer used to store partial blocks
  buffer m_buffer{};

  //! Set if streaming compression is carried out using isa-l
  const bool m_isal_stream;
  //! CRC32 calculated across all input data (if using isa-l)
  uint32_t m_isal_crc32 = 0;

  //! Used to track whether an EOF block has been received.
  bool m_eof = false;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock{};
};

/**
 * GZip compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class gzip_split_fastq : public analytical_step
{
public:
  /** Constructor; 'next_step' sets the destination of compressed chunks. */
  gzip_split_fastq(const userconfig& config,
                   const output_file& file,
                   size_t next_step);

  ~gzip_split_fastq() override = default;

  /** Compresses input lines, saving compressed chunks to chunk->buffers. */
  chunk_vec process(chunk_ptr chunk) override;

  gzip_split_fastq(const gzip_split_fastq&) = delete;
  gzip_split_fastq(gzip_split_fastq&&) = delete;
  gzip_split_fastq& operator=(const gzip_split_fastq&) = delete;
  gzip_split_fastq& operator=(gzip_split_fastq&&) = delete;

private:
  const userconfig& m_config;
  //! Set if streaming compression is carried out using isa-l
  const bool m_isal_stream;
  //! Format of file being compressed
  const output_format m_format;
  //! The analytical step following this step
  const size_t m_next_step;
};

/**
 * Simple file reading step.
 *
 * The 'process' function takes a fastq_file_chunk object and writes the lines
 * at the offset corresponding to the 'type' argument to the corresponding
 * output file. The list of lines is cleared upon writing.
 */
class write_fastq : public analytical_step
{
public:
  /** Opens the specified file for writing using a managed writer */
  write_fastq(const userconfig& config, const output_file& file);

  /** Writes the reads of the type specified in the constructor. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Flushes the output file and prints progress report (if enabled). */
  void finalize() override;

private:
  //! Lazily opened / automatically closed handle
  managed_writer m_output;

  //! Set if streaming compression is carried out using isa-l
  const bool m_isal_stream;
  //! Used to track the total number of (uncompressed) bytes written
  size_t m_uncompressed_bytes = 0;
  //! Used to track whether an EOF block has been received.
  bool m_eof = false;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock{};
};

} // namespace adapterremoval
