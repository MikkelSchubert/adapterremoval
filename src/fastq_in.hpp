// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "buffer.hpp"        // for buffer
#include "buffer_reader.hpp" // for buffer_reader
#include "debug.hpp"         // for AR_REQUIRE
#include "fastq.hpp"         // for fastq
#include "managed_io.hpp"    // for managed_reader
#include "scheduler.hpp"     // for analytical_step, chunk_ptr, chunk_vec
#include "threading.hpp"     // for reusable_objects
#include <array>             // for array
#include <atomic>            // for atomic_bool
#include <cstddef>           // for size_t
#include <cstdint>           // for uint8_t, uint64_t
#include <deque>             // for deque
#include <iosfwd>            // for ostream
#include <limits>            // for numeric_limits
#include <memory>            // for unique_ptr, shared_ptr
#include <string>            // for string
#include <vector>            // for vector

struct inflate_state;

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

/** Describes the content of data derived from FASTQ files */
enum class fastq_input_type
{
  //! File containing only mate 1 reads
  mate_1,
  //! File containing only mate 2 reads
  mate_2,
  //! File containing only interleaved mate 1 and 2 reads
  interleaved,
};

/** FASTQ chunk that additionally tracks the file source */
class fastq_parser_chunk : public fastq_chunk
{
public:
  /** Was these reads read from file1 or file2 */
  fastq_input_type type = fastq_input_type::mate_1;
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

/** Decompresses plain Gzip, forwarding BGZF to parallel decompression */
class serial_decompressor : public analytical_step
{
public:
  serial_decompressor(reader_state_ptr state, step_id next_step);
  ~serial_decompressor() override;

  /** Emit decompressed chunk; signalling upstream if needed */
  [[nodiscard]] chunk_vec signal() override;

  /**
   * Classifies raw (unknown) data chunks and returns up to `m_pending_chunks`
   * chunks. Chunks may be decompressed (if gzip), forwarded as is (if not
   * compressed), or forwarded for parallel decompression (if bgzf).
   */
  [[nodiscard]] chunk_vec process(chunk_ptr data) override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  serial_decompressor(const serial_decompressor&) = delete;
  serial_decompressor(serial_decompressor&&) = delete;
  serial_decompressor& operator=(const serial_decompressor&) = delete;
  serial_decompressor& operator=(serial_decompressor&&) = delete;

private:
  /** Emits up to `m_pending_chunks` reader chunks */
  chunk_vec emit();

  /** Attempts to identify the next chunk type, returning true if successful */
  bool identify_next_chunk();

  /** Returns the next chunk, or nullptr if not enough data */
  reader_chunk_ptr next_chunk();
  /** Returns chunk of uncompressed data, or nullptr if not enough data */
  reader_chunk_ptr next_text_chunk();
  /** Returns chunk of decompressed gzip data, or nullptr if not enough data */
  reader_chunk_ptr next_gzip_chunk();
  /** Returns a BGZF block in a chunk, or nullptr if not enough data */
  reader_chunk_ptr next_bgzf_chunk();

  /**
   * Copies a BGZF block into `m_chunk`, returning true if `m_chunk` contains a
   * complete chunk or if EOF was reached in the process. Error handling is left
   * to the BGZF decompression step in the case of EOF.
   */
  bool build_bgzf_block();

  /** Warns if a run of BGZF blocks is not terminated by an EOF block */
  void bgzf_require_eof_block(shared_filename filename);
  /** Warns that a BGZF block BSIZE is smaller than BGZF_EOF */
  void bgzf_warning_small_bsize(shared_filename filename, size_t bsize);
  /** Warns that a BGZF block has an ISIZE larger than BGZF_MAX_BLOCK_SIZE */
  void bgzf_warning_large_isize(shared_filename filename, size_t isize);

  //! State shared between this step and other reader steps
  const reader_state_ptr m_state;

  /** Specifies the type of data currently being read */
  enum class data_type
  {
    //! Data type is unknown and must be identified to proceed
    unknown,
    //! Data is not Gzip/BGZF and is therefore presumed to be text
    text,
    //! Data is a Gzip stream, but does not follow the BGZF specification
    gzip,
    //! Data is a BGZF compressed block
    bgzf,
    //! Either a gzip stream or a BGZF block has been read; may either be
    //! followed by another gzip stream/BGZF block, or garbage (anything else)
    after_gzip,
    //! Trailing garbage after gzip files
    garbage,
  };

  //! Type of data currently being processed
  data_type m_data_type = data_type::unknown;
  //! Size of the current BGZF block, if m_data_type == data_type::bgzf;
  size_t m_bgzf_bsize = 0;
  //! End of decompressed data in m_chunk, if m_data_type == data_type::gzip;
  size_t m_gzip_offset = 0;
  //! The current BGZF/uncompressed chunk being built, if any
  buffer m_chunk{};
  //! Reusable buffer for gathering BGZF/gzip headers
  buffer m_header{};

  //! Size of the last BGZF block, used to check for EOF block in BGZF streams
  size_t m_last_bgzf_isize = 0;
  //! Indicates filename for which the last missing EOF marker was logged
  shared_filename m_bgzf_warned_for_eof{};
  //! Indicates filename for which the BSIZE related warning was logged
  shared_filename m_bgzf_warned_for_bsize{};
  //! Indicates filename for which the ISIZE related warning was logged
  shared_filename m_bgzf_warned_for_isize{};

  /** Wrapper class that tracks consumption of a `reader_chunk` */
  struct cached_chunk
  {
    cached_chunk(reader_chunk_ptr&& chunk_)
      : chunk(std::move(chunk_))
    {
      AR_REQUIRE(chunk);
    }

    /** Returns pointer to remaining data */
    [[nodiscard]] const uint8_t* avail_data() const
    {
      return chunk->data.data() + offset;
    }

    /** Returns length of remaining data */
    [[nodiscard]] size_t avail() const { return chunk->data.size() - offset; }

    /** Acts as if the original reader_chunk_ptr was dereferenced */
    const reader_chunk* operator->() const { return chunk.get(); }

    //! The underlying chunk
    reader_chunk_ptr chunk;
    //! Offset into the chunk's data, if it has already been partially processed
    size_t offset = 0;
  };

  //! Unprocessed (possibly compressed) chunks
  std::deque<cached_chunk> m_cache{};
  // Inflater used for plain gzip files/chunks
  std::unique_ptr<inflate_state> m_gzip_stream{};

  //! The number of chunks requested via `signal()` but not yet returned
  size_t m_pending_chunks = 0;
};

/** Decompresses BGZF chunks, and forwards other chunk-types as-is */
class parallel_decompressor : public analytical_step
{
public:
  parallel_decompressor(reader_state_ptr state, step_id next_step);
  ~parallel_decompressor() override;

  /**
   * Decompresses BGZF chunks, forwarding text chunks as is; once parsing is
   * done, both BGZF and text chunks are forwarded as is. Unknown chunks are
   * not permitted.
   */
  [[nodiscard]] chunk_vec process(chunk_ptr data) override;

  parallel_decompressor(const parallel_decompressor&) = delete;
  parallel_decompressor(parallel_decompressor&&) = delete;
  parallel_decompressor& operator=(const parallel_decompressor&) = delete;
  parallel_decompressor& operator=(parallel_decompressor&&) = delete;

private:
  //! State shared between this step and other reader steps
  const reader_state_ptr m_state;
};

/** Parses FASTQ reads from decompressed buffers */
class fastq_parser : public analytical_step
{
public:
  fastq_parser(std::vector<std::string> filenames,
               fastq_input_type type,
               uint64_t head,
               reader_state_ptr state,
               step_id next_step);
  ~fastq_parser() override;

  /** Returns chunk or signals for more data */
  [[nodiscard]] chunk_vec signal() override;

  /** Reads lines from the input file and saves them in an fastq_parser_chunk */
  [[nodiscard]] chunk_vec process(chunk_ptr data) override;

  /** Finalizer; checks that all input has been processed. */
  void finalize() override;

  fastq_parser(const fastq_parser&) = delete;
  fastq_parser(fastq_parser&&) = delete;
  fastq_parser& operator=(const fastq_parser&) = delete;
  fastq_parser& operator=(fastq_parser&&) = delete;

private:
  using result = buffer_reader::result;

  /** Emits up to `m_pending_chunks` FASTQ chunks */
  chunk_vec emit();

  //! Attempt to parse a chunk of FASTQ reads
  chunk_ptr parse_chunk();
  /** Attempt to parse up to `to_parse` FASTQ records from a single source */
  result read_single_end(uint64_t to_parse);
  /** Attempt to parse up to `to_parse` FASTQ interleaved record pairs */
  result read_interleaved(uint64_t to_parse);
  //! Attempt to parse a FASTQ record from available buffers
  buffer_reader::result parse_fastq(fastq& dst);

  //! State shared between this step and other reader steps
  const reader_state_ptr m_state;
  //! The input being type parsed
  const fastq_input_type m_type;

  //! Line reader on sequence of decompressed blocks
  buffer_reader m_reader;
  //! Start position of current read
  size_t m_position = 0;

  //! Current chunks being built for mate 1 and mate 2 reads
  std::array<std::vector<fastq>, 2> m_chunks{};
  //! The mate for which the current fastq record is being parsed (0 or 1)
  size_t m_chunk_mate = 0;
  //! Current fastq record being parsed
  size_t m_chunk_offset = 0;
  //! Remaining number of reads to parse
  uint64_t m_head = std::numeric_limits<uint64_t>::max();

  enum class fastq_parse
  {
    header,
    sequence,
    separator,
    qualities,
  };

  //! Current stage of FASTQ parsing
  fastq_parse m_chunk_parse = fastq_parse::header;
  //! Reused buffer for reading separator line
  std::string m_separator{};

  //! The number of chunks requested via `signal()` but not yet returned
  size_t m_pending_chunks = 0;
};

} // namespace adapterremoval
