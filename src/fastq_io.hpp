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
#pragma once

#include "buffer.hpp"            // for buffer, buffer_vec
#include "commontypes.hpp"       // for fastq_vec
#include "fastq.hpp"             // for fastq
#include "fastq_enc.hpp"         // for fastq_encoding
#include "linereader_joined.hpp" // for joined_line_readers
#include "managed_writer.hpp"    // for managed_writer
#include "progress.hpp"          // for progress_timer
#include "scheduler.hpp"         // for analytical_step, chunk_ptr, chunk_vec
#include "statistics.hpp"        // for fastq_stats_ptr
#include <memory>                // for unique_ptr
#include <mutex>                 // for mutex
#include <stddef.h>              // for size_t
#include <stdint.h>              // for uint32_t, uint64_t
#include <string>                // for string

namespace adapterremoval {

class fastq_output_chunk;
class fastq_read_chunk;
class userconfig;

using output_chunk_ptr = std::unique_ptr<fastq_output_chunk>;
using read_chunk_ptr = std::unique_ptr<fastq_read_chunk>;

//! Rough number of nucleotides to read every cycle
const size_t INPUT_BLOCK_SIZE = 4 * 64 * 1024;
//! Size of chunks of when performing block compression
const size_t GZIP_BLOCK_SIZE = 64 * 1024;
//! Size of blocks to generate before writing to output
const size_t OUTPUT_BLOCK_SIZE = 4 * 64 * 1024;

/**
 * Container object for (demultiplexed) reads.
 */
class fastq_read_chunk : public analytical_chunk
{
public:
  /** Create chunk representing lines starting at line offset (1-based). */
  explicit fastq_read_chunk(bool eof_ = false);

  ~fastq_read_chunk() override = default;

  //! Indicates that EOF has been reached.
  bool eof;

  //! Total number of nucleotides in this chunk
  size_t nucleotides;

  //! Lines read from the mate 1 files
  fastq_vec reads_1;
  //! Lines read from the mate 2 files
  fastq_vec reads_2;
};

/**
 * Container object for processed reads.
 */
class fastq_output_chunk : public analytical_chunk
{
public:
  /** Constructor; does nothing. */
  explicit fastq_output_chunk(bool eof_ = false, uint32_t crc32 = 0);

  ~fastq_output_chunk() override = default;

  /** Add FASTQ read to output buffer. */
  void add(const fastq& read);

  //! Indicates that EOF has been reached.
  bool eof;
  //! CRC32 of (uncompressed) data; only set if eof is true
  uint32_t crc32;

  //! Total number of nucleotides in this chunk
  size_t nucleotides;

  //! Encoded FASTQ reads
  std::string reads;

  //! Buffers of (compressed) FASTQ reads
  buffer_vec buffers;

  //! Size of (uncompressed) data in buffers;
  size_t uncompressed_size;
};

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
  /**
   * Constructor.
   */
  read_fastq(const userconfig& config, size_t next_step);

  /** Reads lines from the input file and saves them in an fastq_file_chunk. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  //! Copy construction not supported
  read_fastq(const read_fastq&) = delete;
  //! Assignment not supported
  read_fastq& operator=(const read_fastq&) = delete;

private:
  /** Fills a chunk with SE reads; stops on EOF or after `head` reads. */
  bool read_single_end(fastq_vec& reads_1);
  /** Fills a chunk with PE reads; stops on EOF or after `head` reads. */
  bool read_paired_end(fastq_vec& reads_1, fastq_vec& reads_2);
  /** Attempt to identify the mate separator for the provided reads */
  char identify_mate_separators(const fastq_vec& reads_1,
                                const fastq_vec& reads_2) const;

  //! The underlying file reader for mate 1 (and possibly mate 2) reads
  joined_line_readers m_io_input_1_base;
  //! The underlying file reader for mate 2 read (if not interleaved)
  joined_line_readers m_io_input_2_base;
  //! The reader used to read mate 1 reads.
  joined_line_readers* m_io_input_1;
  //! The reader used to read mate 1 reads; may be equal to m_io_input_1.
  joined_line_readers* m_io_input_2;
  //! The analytical step following this step
  const size_t m_next_step;

  //! Character used to join read-names with mate numbers, e.g. '/'
  char m_mate_separator;
  //! True if input is single-end
  bool m_single_end;
  //! Used to track whether an EOF block has been received.
  bool m_eof;

  //! Timer for displaying read progress.
  progress_timer m_timer;

  //! Number of reads to process
  uint64_t m_head;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
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
  post_process_fastq(size_t next_step,
                     statistics& stats,
                     const fastq_encoding& encoding);

  /** Reads lines from the input file and saves them in an fastq_file_chunk. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  //! Copy construction not supported
  post_process_fastq(const post_process_fastq&) = delete;
  //! Assignment not supported
  post_process_fastq& operator=(const post_process_fastq&) = delete;

private:
  //! Statistics collected from raw mate 1 reads
  fastq_stats_ptr m_statistics_1;
  //! Statistics collected from raw mate 2 reads
  fastq_stats_ptr m_statistics_2;
  //! The analytical step following this step
  const size_t m_next_step;
  //! Encoding used to parse FASTQ reads.
  const fastq_encoding m_encoding;

  //! Used to track whether an EOF block has been received.
  bool m_eof;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
};

/**
 * GZip compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class gzip_fastq : public analytical_step
{
public:
  /** Constructor; 'next_step' sets the destination of compressed chunks. */
  gzip_fastq(const userconfig& config, size_t next_step);

  /** Compresses input lines, saving compressed chunks to chunk->buffers. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Checks that all input has been processed and frees stream. */
  void finalize() override;

  //! Copy construction not supported
  gzip_fastq(const gzip_fastq&) = delete;
  //! Assignment not supported
  gzip_fastq& operator=(const gzip_fastq&) = delete;

private:
  //! The analytical step following this step
  const size_t m_next_step;
  //! Used to track whether an EOF block has been received.
  bool m_eof;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
};

/**
 * Splits input into chunks that can be GZipped in parallel.
 */
class split_fastq : public analytical_step
{
public:
  /** Constructor; 'next_step' sets the destination of compressed chunks. */
  split_fastq(const userconfig& config,
              const std::string& filename,
              size_t next_step);

  chunk_vec process(chunk_ptr chunk) override;
  void finalize() override;

  //! Copy construction not supported
  split_fastq(const split_fastq&) = delete;
  //! Assignment not supported
  split_fastq& operator=(const split_fastq&) = delete;

private:
  //! The analytical step following this step
  const size_t m_next_step;
  //! Buffer used to store partial blocks
  buffer m_buffer;
  //! Offset in current buffer
  size_t m_offset;

  //! Set if compression is carried out using isa-l
  bool m_isal_enabled;
  //! CRC32 calculated across all input data (if using isa-l)
  uint32_t m_isal_crc32;

  //! Used to track whether an EOF block has been received.
  bool m_eof;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
};

/**
 * GZip compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class gzip_split_fastq : public analytical_step
{
public:
  /** Constructor; 'next_step' sets the destination of compressed chunks. */
  gzip_split_fastq(const userconfig& config,
                   const std::string& filename,
                   size_t next_step);

  /** Compresses input lines, saving compressed chunks to chunk->buffers. */
  chunk_vec process(chunk_ptr chunk) override;

  //! Copy construction not supported
  gzip_split_fastq(const gzip_split_fastq&) = delete;
  //! Assignment not supported
  gzip_split_fastq& operator=(const gzip_split_fastq&) = delete;

private:
  const userconfig& m_config;
  //! Set if compression is carried out using isa-l
  bool m_isal_enabled;
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
  /**
   * Opens the specified file for writing using a managed writer.
   *
   * The filename "/dev/stdout" is handled specially and does not count as IO
   * for the purpose of scheduling the pipeline.
   */
  write_fastq(const userconfig& config, const std::string& filename);

  /** Writes the reads of the type specified in the constructor. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Flushes the output file and prints progress report (if enabled). */
  void finalize() override;

private:
  //! Lazily opened / automatically closed handle
  managed_writer m_output;

  //! Set if compression is carried out using isa-l
  bool m_isal_enabled;
  //! Used to track the total number of (uncompressed) bytes written
  size_t m_uncompressed_bytes;
  //! Used to track whether an EOF block has been received.
  bool m_eof;

  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
};

} // namespace adapterremoval
