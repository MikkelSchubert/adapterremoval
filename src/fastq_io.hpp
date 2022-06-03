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
#pragma once

#include <memory>   // for unique_ptr
#include <mutex>    // for mutex
#include <stddef.h> // for size_t
#include <string>   // for string
#include <zlib.h>   // for z_stream

#include "commontypes.hpp"       // for fastq_vec, string_vec
#include "fastq_enc.hpp"         // for fastq_encoding
#include "linereader_joined.hpp" // for joined_line_readers
#include "managed_writer.hpp"    // for buffer_ptr, buffer_vec, managed_writer
#include "scheduler.hpp"         // for analytical_step, chunk_vec, analyti...
#include "statistics.hpp"        // for fastq_stats_ptr, fastq_statistics, ...
#include "timer.hpp"             // for progress_timer

class fastq;
class fastq_encoding;
class fastq_output_chunk;
class fastq_read_chunk;
class userconfig;

typedef std::unique_ptr<fastq_output_chunk> output_chunk_ptr;
typedef std::unique_ptr<fastq_read_chunk> read_chunk_ptr;

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
  fastq_read_chunk(bool eof_ = false);

  virtual ~fastq_read_chunk() override;

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
  fastq_output_chunk(bool eof_ = false);

  virtual ~fastq_output_chunk() override;

  /** Add FASTQ read to output buffer. */
  void add(const fastq& read);

  //! Indicates that EOF has been reached.
  bool eof;

  //! Total number of nucleotides in this chunk
  size_t nucleotides;

  //! Encoded FASTQ reads
  std::string reads;

  //! Buffers of (compressed) FASTQ reads
  buffer_vec buffers;
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
  void read_single_end(read_chunk_ptr& chunk);
  /** Fills a chunk with PE reads; stops on EOF or after `head` reads. */
  void read_paired_end(read_chunk_ptr& chunk);

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

  //! True if input is single-end
  bool m_single_end;
  //! Used to track whether an EOF block has been received.
  bool m_eof;

  //! Timer for displaying read progress.
  progress_timer m_timer;

  //! Number of reads to process
  unsigned m_head;

#ifdef DEBUG
  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
#endif
};

/**
 * Class responsible for validating FASTQ records and collecting statistics.
 * This is split into a separate step to minimize the amount of time that IO
 * is blocked by the FASTQ reading step.
 */
class post_process_fastq : public analytical_step
{
public:
  /** Constructor. */
  post_process_fastq(const userconfig& config,
                     size_t next_step,
                     statistics* stats = nullptr);

  /** Reads lines from the input file and saves them in an fastq_file_chunk. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  //! Copy construction not supported
  post_process_fastq(const post_process_fastq&) = delete;
  //! Assignment not supported
  post_process_fastq& operator=(const post_process_fastq&) = delete;

private:
  void process_single_end(fastq_vec& reads_1);
  void process_paired_end(fastq_vec& reads_1, fastq_vec& reads_2);

  //! Encoding used to parse FASTQ reads.
  const fastq_encoding m_encoding;
  //! Character used to join read-names with mate numbers, e.g. '/'
  char m_mate_separator;
  //! Statistics collected from raw mate 1 reads
  fastq_stats_ptr m_statistics_1;
  //! Statistics collected from raw mate 2 reads
  fastq_stats_ptr m_statistics_2;
  //! The analytical step following this step
  const size_t m_next_step;

  //! Used to track whether an EOF block has been received.
  bool m_eof;

#ifdef DEBUG
  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
#endif
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
  //! GZip stream object
  z_stream m_stream;
  //! Used to track whether an EOF block has been received.
  bool m_eof;

#ifdef DEBUG
  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
#endif
};

/**
 * Splits input into chunks that can be GZipped in parallel.
 */
class split_fastq : public analytical_step
{
public:
  /** Constructor; 'next_step' sets the destination of compressed chunks. */
  split_fastq(size_t next_step);

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
  buffer_ptr m_buffer;
  //! Offset in current buffer
  size_t m_offset;

  //! Used to track whether an EOF block has been received.
  bool m_eof;

#ifdef DEBUG
  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
#endif
};

/**
 * GZip compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class gzip_split_fastq : public analytical_step
{
public:
  /** Constructor; 'next_step' sets the destination of compressed chunks. */
  gzip_split_fastq(const userconfig& config, size_t next_step);

  /** Compresses input lines, saving compressed chunks to chunk->buffers. */
  chunk_vec process(chunk_ptr chunk) override;

  //! Copy construction not supported
  gzip_split_fastq(const gzip_split_fastq&) = delete;
  //! Assignment not supported
  gzip_split_fastq& operator=(const gzip_split_fastq&) = delete;

private:
  const userconfig& m_config;
  //! The analytical step following this step
  const size_t m_next_step;
  //! Already allocated buffers usable for compressed output
  threadstate<unsigned char[]> m_buffers;
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
  write_fastq(const std::string& filename);

  /** Writes the reads of the type specified in the constructor. */
  chunk_vec process(chunk_ptr chunk) override;

  /** Flushes the output file and prints progress report (if enabled). */
  void finalize() override;

private:
  //! Lazily opened / automatically closed handle
  managed_writer m_output;

  //! Used to track whether an EOF block has been received.
  bool m_eof;

#ifdef DEBUG
  //! Lock used to verify that the analytical_step is only run sequentially.
  std::mutex m_lock;
#endif
};
