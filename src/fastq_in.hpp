// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "buffer.hpp"     // for buffer
#include "managed_io.hpp" // for managed_reader
#include "scheduler.hpp"  // for analytical_step, chunk_ptr, chunk_vec
#include "threading.hpp"  // for reusable_objects
#include <atomic>         // for atomic_bool
#include <iosfwd>         // for ostream
#include <memory>         // for unique_ptr, shared_ptr
#include <string>         // for string
#include <vector>         // for vector

namespace adapterremoval {

class reader_chunk;
struct reader_state;

using reader_chunk_ptr = std::unique_ptr<reader_chunk>;
using reader_state_ptr = std::shared_ptr<reader_state>;
using reusable_buffers = reusable_objects<buffer>;
using reusable_fastqs = reusable_objects<std::vector<fastq>>;

using shared_filename = std::shared_ptr<const std::string>;

/** State shared between all reader steps */
struct reader_state
{
  /** Create shared state for mate 1, 2, or interleaved reads */
  reader_state(reusable_buffers buffers_, reusable_fastqs fastqs_);

  /** Default destructor defined out of line */
  ~reader_state();

  //! Cache of reusable buffers
  reusable_buffers buffers;
  //! Cache of reusable FASTQ vectors
  reusable_fastqs fastqs;

  /** Is the block reader still producing input chunks? */
  [[nodiscard]] bool is_reading() const { return m_reading; }

  /** Used by the block reader to indicate that it no longer produces chunks */
  void mark_reading_finished() { m_reading = false; }

  enum class parsing
  {
    //! Still parsing
    yes,
    //! Terminated parsing due to having parsed --head records
    reached_head,
    //! Terminated parsing due to reaching EOF
    reached_eof,
  };

  /** Is the FASTQ parser still producing FASTQ chunks? */
  [[nodiscard]] bool is_parsing() const { return m_parsing == parsing::yes; }

  /** Returns true if the FASTQ parser has generated `--head` records */
  [[nodiscard]] bool is_parsing_at_head() const
  {
    return m_parsing == parsing::reached_head;
  }

  /**
   * Used by the FASTQ parser to indicate that it no longer produces chunks.
   * This function may only be called once, and not with `yes`.
   */
  void mark_parsing(parsing value);

private:
  //! Indicates if the reader is still producing chunks from input files
  std::atomic_bool m_reading{ true };
  //! Indicates if FASTQ chunks are being produced, capped at `--head` reads
  std::atomic<parsing> m_parsing{ parsing::yes };
};

/** Chunk used to pass raw/compressed/decompressed data between reader steps */
class reader_chunk : public analytical_chunk
{
public:
  /** Creates reader chunk that automatically releases the associated buffer */
  explicit reader_chunk(reusable_buffers buffers_);

  /** Releases the buffer, if it has capacity */
  ~reader_chunk() override;

  /** The type of data carried by the chunk, if any */
  enum class type
  {
    //! Unknown data, assumed to be uncompressed text, BGZF, or plain gzip
    unknown,
    //! BGZF compressed data
    bgzf,
    //! Uncompressed data, assumed to be text
    text
  };

  //! Name of current source
  shared_filename filename{};
  //! Compressed data or text. MUST be acquired from the same `reusable_buffers`
  //! passed to the constructor of the chunk
  buffer data{};
  //! Type of the current content in `data` buffer
  type data_type = type::unknown;
  //! Last block from current source
  bool end_of_file = false;

private:
  //! Cache used to release the buffer assigned to this chunk
  reusable_buffers buffers;
};

/** Converts reader_chunk::type values to string for debug logging */
std::ostream&
operator<<(std::ostream& os, reader_chunk::type value);

/**
 * Reads fixed-size chunks of data from one or more input files until EOF; the
 * final chunk size will be in the range 0 to X - 1, where X is the default
 * capacity for the reusable buffers cache, with the `end_of_file` flag set.
 *
 * The reader returns blank chunks after EOF, with neither data capacity nor
 * a filename.
 */
class block_reader : public analytical_step
{
public:
  block_reader(const std::vector<std::string>& filenames,
               reader_state_ptr state,
               step_id next_step);
  ~block_reader() override;

  /** Reads a chunk of data from the current file and returns it in a buffer */
  [[nodiscard]] chunk_vec signal() override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  block_reader(const block_reader&) = delete;
  block_reader(block_reader&&) = delete;
  block_reader& operator=(const block_reader&) = delete;
  block_reader& operator=(block_reader&&) = delete;

protected:
  /**
   * Overridable helper function that opens the specified file; must return a
   * reader that indicates EOF by returning less than the number of bytes
   * requested
   */
  virtual std::unique_ptr<managed_reader> open_reader(std::string filename);

  /** Close the current reader, if any */
  void close_current_reader();

  //! State shared between this step and other reader steps
  const reader_state_ptr m_state;

  //! Currently open input file
  std::unique_ptr<managed_reader> m_reader{};
  //! Reversed list of remaining files to read, including the current file
  std::vector<shared_filename> m_filenames{};
};

} // namespace adapterremoval
