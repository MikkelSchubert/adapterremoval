// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "fastq_in.hpp"      // declarations
#include "bgzf.hpp"          // for BGZF_HEADER, BGZF_MAX_BLOCK_SIZE, ...
#include "buffer.hpp"        // for buffer
#include "buffer_reader.hpp" // for buffer_reader
#include "commontypes.hpp"   // for read_mate
#include "debug.hpp"         // for AR_FAIL, AR_REQUIRE
#include "deflate.hpp"       // for checked_isal_inflate, inflate_state
#include "errors.hpp"        // for fastq_error, gzip_error
#include "fastq.hpp"         // for fastq
#include "logging.hpp"       // for info, warn
#include "managed_io.hpp"    // for managed_reader
#include "scheduler.hpp"     // for analytical_chunk, analytical_step, ...
#include "strutils.hpp"      // for log_escape
#include "userconfig.hpp"    // for userconfig
#include "utilities.hpp"     // for dynamic_cast_unique
#include <algorithm>         // for min, reverse, swap
#include <cstddef>           // for size_t
#include <cstdint>           // for uint8_t
#include <deque>             // for deque
#include <isa-l.h>           // for inflate_state
#include <memory>            // for make_shared, make_unique
#include <sstream>           // for ostringstream

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// reader_state

reader_state::reader_state(const reusable_buffers& buffers_,
                           const reusable_fastqs& fastqs_)
  : buffers(buffers_)
  , fastqs(fastqs_)
{
}

reader_state::~reader_state() = default;

void
reader_state::mark_parsing(parsing value)
{
  auto initial = parsing::yes;
  AR_REQUIRE(value != initial, "~2fb83d9febd61f6f~");
  const auto exchanged = m_parsing.compare_exchange_strong(initial, value);
  AR_REQUIRE(exchanged, "~6e7c8cea6a4745bc~");
}

////////////////////////////////////////////////////////////////////////////////
// reader_chunk

reader_chunk::reader_chunk(const reusable_buffers& buffers_)
  : buffers(buffers_)
{
}

reader_chunk::~reader_chunk()
{
  buffers.release(std::move(data));
}

std::ostream&
operator<<(std::ostream& os, reader_chunk::type value)
{
  switch (value) {
    case reader_chunk::type::unknown:
      return os << "reader_chunk::type::unknown";
    case reader_chunk::type::bgzf:
      return os << "reader_chunk::type::bgzf";
    case reader_chunk::type::text:
      return os << "reader_chunk::type::text";
    default:                                 // GCOVR_EXCL_LINE
      AR_FAIL("invalid reader_chunk::type"); // GCOVR_EXCL_LINE
  }
}

////////////////////////////////////////////////////////////////////////////////
// block_reader

block_reader::block_reader(const std::vector<std::string>& filenames,
                           reader_state_ptr state,
                           step_id next_step)
  : analytical_step(processing_order::ordered_io, "block_reader", next_step)
  , m_state(std::move(state))
{
  AR_REQUIRE(m_state, "~1cab54b36e0edde6~");
  AR_REQUIRE(!filenames.empty(), "~b5bb8b2a98f75df5~");

  for (auto it = filenames.rbegin(); it != filenames.rend(); ++it) {
    m_filenames.emplace_back(std::make_shared<const std::string>(*it));
  }
}

block_reader::~block_reader() = default;

chunk_vec
block_reader::signal()
{
  chunk_vec chunks;
  if (!m_state->is_parsing()) {
    AR_REQUIRE(!m_state->is_reading() || m_state->is_parsing_at_head(),
               "~4a20e925ef1d96b4~");
    close_current_reader();
    m_state->mark_reading_finished();

    chunks.emplace_back(scheduler_signal::suspend);
    return chunks;
  }

  // The reader must return chunks even if there is no data left to read, since
  // these signal to the recipient that another chunk should be emitted
  auto chunk = std::make_unique<reader_chunk>(m_state->buffers);
  if (m_state->is_reading()) {
    AR_REQUIRE(!m_filenames.empty());
    auto filename = m_filenames.back();

    if (!m_reader) {
      log::info() << "Opening FASTQ file " << log_escape(*filename);
      m_reader = open_reader(*filename);
      AR_REQUIRE(m_reader);
    }

    chunk->data_type = reader_chunk::type::unknown;
    chunk->filename = std::move(filename);
    chunk->data = m_state->buffers.acquire();
    const auto nread = m_reader->read(chunk->data.data(), chunk->data.size());
    if (nread != chunk->data.size()) {
      close_current_reader();

      chunk->data.resize(nread);
      chunk->end_of_file = true;
    }
  }

  chunks.emplace_back(std::move(chunk));

  return chunks;
}

void
block_reader::finalize()
{
  AR_REQUIRE((!m_reader && m_filenames.empty()) ||
               m_state->is_parsing_at_head(),
             "~fd72e01d62e4c694~");

  close_current_reader();
}

std::unique_ptr<managed_reader>
block_reader::open_reader(std::string filename)
{
  return std::make_unique<managed_reader>(std::move(filename));
}

void
block_reader::close_current_reader()
{
  if (auto reader = std::move(m_reader)) {
    log::info() << "Closing FASTQ file " << log_escape(*m_filenames.back());
    m_filenames.pop_back();
    if (m_filenames.empty()) {
      m_state->mark_reading_finished();
    }

    reader->close(); // called outside of destructor to allow exceptions
  }
}

////////////////////////////////////////////////////////////////////////////////
// serial_decompressor

serial_decompressor::serial_decompressor(reader_state_ptr state,
                                         step_id next_step)
  : analytical_step(processing_order::ordered, "serial_decompressor", next_step)
  , m_state(std::move(state))
{
  AR_REQUIRE(m_state, "~79a956c1fd6f2ac5~");

  m_header.reserve(BGZF_HEADER.size());
  m_gzip_stream = std::make_unique<inflate_state>();
  isal_inflate_init(m_gzip_stream.get());
  m_gzip_stream->crc_flag = ISAL_GZIP;
}

serial_decompressor::~serial_decompressor() = default;

chunk_vec
serial_decompressor::signal()
{
  chunk_vec chunks;
  if (m_state->is_parsing()) {
    m_pending_chunks++;

    chunks = emit();
  } else {
    chunks.emplace_back(scheduler_signal::suspend);
  }

  return chunks;
}

chunk_vec
serial_decompressor::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<reader_chunk>(data);
  AR_REQUIRE(chunk && chunk->data_type == reader_chunk::type::unknown,
             "~47f0d52258dba7cb~");

  if (!m_state->is_parsing()) {
    return {};
  } else if (chunk->data.capacity()) {
    AR_REQUIRE(chunk->filename, "~ad8643ad5a07417a~");
    m_cache.emplace_back(std::move(chunk));
  }

  // Return chunk or request more data; one signal -> one chunk
  return emit();
}

void
serial_decompressor::finalize()
{
  AR_REQUIRE(m_chunk.capacity() == 0 || m_state->is_parsing_at_head(),
             "~ab0c8ee7277efef1~");
  AR_REQUIRE(m_cache.empty() || m_state->is_parsing_at_head(),
             "~b194f127cdcb5d3e~");

  // Release any acquired buffers to `reusable_objects`
  m_state->buffers.release(std::move(m_chunk));
  m_cache.clear();
}

chunk_vec
serial_decompressor::emit()
{
  chunk_vec chunks;
  while (m_pending_chunks && m_state->is_parsing()) {
    if (auto output = next_chunk()) {
      chunks.emplace_back(std::move(output));
      m_pending_chunks--;
    } else if (chunks.empty() && m_state->is_reading()) {
      chunks.emplace_back(signal_target());
      break;
    } else {
      break;
    }
  }

  return chunks;
}

bool
serial_decompressor::identify_next_chunk()
{
  m_header.resize(0);
  for (const auto& chunk : m_cache) {
    if (size_t avail = chunk.avail()) {
      avail = std::min(avail, BGZF_HEADER.size() - m_header.size());
      m_header.append(chunk.avail_data(), avail);
    }

    if (m_header.size() == BGZF_HEADER.size() || chunk->end_of_file) {
      if (const auto block_size = bgzf_block_size(m_header)) {
        // BGZF_EOF represents the minimum possible size for a valid BGZF block
        if (block_size >= BGZF_EOF.size()) {
          m_data_type = data_type::bgzf;
          m_bgzf_bsize = block_size;
        } else {
          bgzf_warning_small_bsize(chunk->filename, block_size);
          m_data_type = data_type::gzip;
          // The file/stream does not meet BGZF specifications, so it probably
          // doesn't make sense to also warn if it is missing an EOF marker
          m_last_bgzf_isize = 0;
        }
      } else if (m_header.size() >= 2 && m_header.at(0) == 0x1f &&
                 m_header.at(1) == 0x8b) {
        m_data_type = data_type::gzip;
      } else if (m_data_type == data_type::after_gzip && m_header.size() != 0) {
        // TODO: Add "strict" mode treating this as an error
        log::warn() << "Found trailing garbage at end of "
                    << log_escape(*m_cache.front()->filename);
        m_data_type = data_type::garbage;
      } else {
        m_data_type = data_type::text;
      }

      bgzf_require_eof_block(chunk->filename);

      return true;
    }
  }

  return false;
}

reader_chunk_ptr
serial_decompressor::next_chunk()
{
  if (m_cache.empty()) {
    return {};
  }

  switch (m_data_type) {
    case data_type::after_gzip:
    case data_type::unknown: {
      if (identify_next_chunk()) {
        return next_chunk();
      } else {
        return {};
      }
    }
    case data_type::garbage: {
      do {
        if (m_cache.front()->end_of_file) {
          // Signal to downstream that this file is done
          m_data_type = data_type::text;
          m_cache.front().offset = 0;
          m_cache.front().chunk->data.resize(0);

          return next_text_chunk();
        }

        m_cache.pop_front();
      } while (!m_cache.empty());

      return {};
    }
    case data_type::text:
      return next_text_chunk();
    case data_type::gzip:
      return next_gzip_chunk();
    case data_type::bgzf:
      return next_bgzf_chunk();
    default:                              // GCOVR_EXCL_LINE
      AR_FAIL("invalid data_type value"); // GCOVR_EXCL_LINE
  }
}

reader_chunk_ptr
serial_decompressor::next_text_chunk()
{
  AR_REQUIRE(m_data_type == data_type::text);
  AR_REQUIRE(!m_cache.empty());
  AR_REQUIRE(m_cache.front().offset == 0);

  auto chunk = std::move(m_cache.front().chunk);
  m_cache.pop_front();

  chunk->data_type = reader_chunk::type::text;
  if (chunk->end_of_file) {
    m_data_type = data_type::unknown;
  }

  return chunk;
}

reader_chunk_ptr
serial_decompressor::next_gzip_chunk()
{
  AR_REQUIRE(m_data_type == data_type::gzip);
  AR_REQUIRE(!m_cache.empty());
  AR_REQUIRE(m_chunk.capacity() == 0 ||
             m_chunk.size() == m_state->buffers.capacity());
  if (m_chunk.capacity() == 0) {
    m_chunk = m_state->buffers.acquire();
  }

  do {
    auto& chunk_in = m_cache.front();
    auto filename = chunk_in->filename;

    m_gzip_stream->avail_in = chunk_in.avail();
    m_gzip_stream->next_in = const_cast<uint8_t*>(chunk_in.avail_data());

    AR_REQUIRE(m_gzip_offset <= m_chunk.size());
    m_gzip_stream->next_out = m_chunk.data() + m_gzip_offset;
    m_gzip_stream->avail_out = m_chunk.size() - m_gzip_offset;

    const size_t avail_in_before = m_gzip_stream->avail_in;
    const size_t avail_out_before = m_gzip_stream->avail_out;
    checked_isal_inflate(m_gzip_stream.get(),
                         chunk_in->end_of_file,
                         *chunk_in->filename);
    AR_REQUIRE(avail_out_before != m_gzip_stream->avail_out ||
               avail_in_before != m_gzip_stream->avail_in ||
               m_gzip_stream->block_state ==
                 isal_block_state::ISAL_BLOCK_FINISH);

    m_gzip_offset += avail_out_before - m_gzip_stream->avail_out;
    chunk_in.offset += avail_in_before - m_gzip_stream->avail_in;

    bool end_of_block = false;
    if (m_gzip_stream->block_state == isal_block_state::ISAL_BLOCK_FINISH) {
      // Prepare for subsequent/adjacent streams; end this chunk for simplicity
      isal_inflate_reset(m_gzip_stream.get());
      m_data_type = data_type::after_gzip;
      end_of_block = true;
    }

    bool end_of_file = false;
    if (chunk_in.avail() == 0) {
      if (chunk_in->end_of_file) {
        m_data_type = data_type::unknown;
        end_of_file = true;
      }

      m_gzip_stream->next_in = nullptr;
      m_gzip_stream->avail_in = 0;
      m_cache.pop_front();
    }

    if (end_of_block || end_of_file || m_gzip_stream->avail_out == 0) {
      auto chunk_out = std::make_unique<reader_chunk>(m_state->buffers);
      chunk_out->filename = filename;
      chunk_out->end_of_file = end_of_file;
      chunk_out->data_type = reader_chunk::type::text;
      chunk_out->data = std::move(m_chunk);
      chunk_out->data.resize(m_gzip_offset);
      m_gzip_offset = 0;

      return chunk_out;
    }
  } while (!m_cache.empty());

  return {};
}

reader_chunk_ptr
serial_decompressor::next_bgzf_chunk()
{
  AR_REQUIRE(m_data_type == data_type::bgzf);
  AR_REQUIRE(m_bgzf_bsize >= BGZF_EOF.size());
  AR_REQUIRE(!m_cache.empty());
  if (m_chunk.capacity() == 0) {
    // acquire empty buffer with at least the default capacity
    m_chunk = m_state->buffers.acquire(0);
  }

  reader_chunk_ptr chunk;
  if (!build_bgzf_block()) {
    return chunk;
  } else if (m_chunk.size() != m_bgzf_bsize) {
    throw gzip_error("error while decompressing " +
                     log_escape(*m_cache.front()->filename) +
                     ": file is truncated");
  }

  // If a complete BGZF block has too much content, then warn and fall back to
  // plain gzip decompression. This is closer to previous behavior and allows us
  // to read some files that are invalid BGZF but valid gzip
  m_last_bgzf_isize = m_chunk.get_u32(m_chunk.size() - 4);
  if (m_last_bgzf_isize > BGZF_MAX_BLOCK_SIZE) {
    auto& in_chunk = m_cache.front();
    bgzf_warning_large_isize(in_chunk->filename, m_last_bgzf_isize);
    // The file/stream does not meet BGZF specifications, so it probably doesn't
    // make sense to also warn if it is missing an EOF marker
    m_last_bgzf_isize = 0;

    if (in_chunk.avail() == 0) {
      // the current chunk is reused to avoid an extra empty chunk at EOF
      m_state->buffers.release(std::move(in_chunk.chunk->data));
      in_chunk.chunk->data = std::move(m_chunk);
      in_chunk.offset = 0;
    } else {
      // end_of_file is not set, since it's still set on the next chunk
      auto ptr = std::make_unique<reader_chunk>(m_state->buffers);
      ptr->data = std::move(m_chunk);
      ptr->filename = in_chunk->filename;
      m_cache.emplace_front(std::move(ptr));
    }

    m_data_type = data_type::gzip;
    m_bgzf_bsize = 0;

    // assign to `chunk` to allow copy elision
    chunk = next_gzip_chunk();
    return chunk;
  }

  AR_REQUIRE(!m_cache.empty());
  chunk = std::make_unique<reader_chunk>(m_state->buffers);
  chunk->filename = m_cache.front()->filename;
  chunk->data = std::move(m_chunk);
  chunk->data_type = reader_chunk::type::bgzf;

  if (m_cache.front().avail() == 0) {
    chunk->end_of_file = m_cache.front()->end_of_file;
    m_cache.pop_front();
  }

  m_bgzf_bsize = 0;
  if (chunk->end_of_file) {
    m_data_type = data_type::unknown;
    bgzf_require_eof_block(chunk->filename);
  } else {
    m_data_type = data_type::after_gzip;
  }

  return chunk;
}

bool
serial_decompressor::build_bgzf_block()
{
  AR_REQUIRE(m_data_type == data_type::bgzf);

  bool eof = false;
  for (; !m_cache.empty(); m_cache.pop_front()) {
    auto& chunk = m_cache.front();
    const auto to_copy = std::min(chunk.avail(), m_bgzf_bsize - m_chunk.size());
    m_chunk.append(chunk.avail_data(), to_copy);

    eof = chunk->end_of_file;
    chunk.offset += to_copy;
    if (m_chunk.size() == m_bgzf_bsize || eof) {
      break;
    }
  }

  return eof || m_chunk.size() == m_bgzf_bsize;
}

void
serial_decompressor::bgzf_require_eof_block(shared_filename filename)
{
  if (m_last_bgzf_isize && m_data_type != data_type::bgzf) {
    // only warn once per input file
    if (m_bgzf_warned_for_eof != filename) {
      // TODO: Add "strict" mode treating this as an error
      log::warn() << "BGZF stream in " << log_escape(*filename)
                  << " is missing an EOF marker; file may be truncated";
      m_bgzf_warned_for_eof = filename;
    }

    m_last_bgzf_isize = 0;
  }
}

void
serial_decompressor::bgzf_warning_small_bsize(shared_filename filename,
                                              size_t bsize)
{
  AR_REQUIRE(bsize < BGZF_EOF.size());

  if (m_bgzf_warned_for_bsize != filename) {
    m_bgzf_warned_for_bsize = filename;
    log::warn() << "BGZF block with a BSIZE that is "
                << (BGZF_EOF.size() - bsize)
                << " bytes smaller than the minimum " << BGZF_EOF.size()
                << " found in " << log_escape(*filename)
                << "; will attempt to decompress as plain gzip";
  }
}

void
serial_decompressor::bgzf_warning_large_isize(shared_filename filename,
                                              size_t isize)
{
  AR_REQUIRE(isize > BGZF_MAX_BLOCK_SIZE);

  if (m_bgzf_warned_for_isize != filename) {
    m_bgzf_warned_for_isize = filename;
    log::warn() << "BGZF block with ISIZE " << (isize - BGZF_MAX_BLOCK_SIZE)
                << " bytes larger than the maximum " << BGZF_MAX_BLOCK_SIZE
                << " found in " << log_escape(*filename)
                << "; will attempt to decompress as plain gzip";
  }
}

////////////////////////////////////////////////////////////////////////////////
// parallel_decompressor

parallel_decompressor::parallel_decompressor(reader_state_ptr state,
                                             step_id next_step)
  : analytical_step(processing_order::unordered,
                    "parallel_decompressor",
                    next_step)
  , m_state(std::move(state))
{
  AR_REQUIRE(m_state, "~c4b40d0026ce83ec~");
  AR_REQUIRE(m_state->buffers.capacity() >= BGZF_MAX_BLOCK_SIZE,
             "~27c490c08fca8088~");
}

parallel_decompressor::~parallel_decompressor() = default;

chunk_vec
parallel_decompressor::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<reader_chunk>(data);
  AR_REQUIRE(chunk, "~8a5ad2119fe9e996~");

  AR_REQUIRE(chunk->data_type != reader_chunk::type::unknown,
             "~aa2642c2b9c3b1d6~");
  if (m_state->is_parsing() && chunk->data_type == reader_chunk::type::bgzf) {
    AR_REQUIRE(chunk->filename, "~245a6e94975b7266~");
    AR_REQUIRE(chunk->data.size() >= BGZF_EOF.size(), "~410943f18b07dde0~");

    buffer buf = m_state->buffers.acquire();
    try {
      bgzf_block_decompress(chunk->data, buf, 0, *chunk->filename);
    } catch (...) {
      m_state->buffers.release(std::move(buf));
      throw;
    }

    std::swap(chunk->data, buf);
    m_state->buffers.release(std::move(buf));
    chunk->data_type = reader_chunk::type::text;
  }

  // Chunks must always be forwarded, since unordered steps must produce one
  // output per input, or produce no output at all. Otherwise downstream steps
  // may wait for the missing chunk, if a later chunk had already arrived.
  chunk_vec chunks;
  chunks.emplace_back(std::move(chunk));
  return chunks;
}

////////////////////////////////////////////////////////////////////////////////
// fastq_parser

fastq_parser::fastq_parser(std::vector<std::string> filenames,
                           fastq_input_type type,
                           uint64_t head,
                           reader_state_ptr state,
                           step_id next_step)
  : analytical_step(processing_order::ordered, "fastq_parser", next_step)
  , m_state(std::move(state))
  , m_type(type)
  , m_reader(std::move(filenames))
  , m_chunk_mate(type == fastq_input_type::mate_2 ? 1 : 0)
  , m_head(head)
{
  AR_REQUIRE(m_state, "~5a42dfdc29f8abdb~");

  m_reader.set_release_function(
    [this](buffer&& buf) mutable { m_state->buffers.release(std::move(buf)); });
}

fastq_parser::~fastq_parser() = default;

chunk_vec
fastq_parser::signal()
{
  chunk_vec chunks;
  if (m_state->is_parsing()) {
    m_pending_chunks++;

    chunks = emit();
  } else {
    chunks.emplace_back(scheduler_signal::suspend);
  }

  return chunks;
}

chunk_vec
fastq_parser::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<reader_chunk>(data);
  AR_REQUIRE(chunk, "~8ed559b0c3086856~");

  // BGZF chunks that have not been decompressed may arrive after parsing ended
  if (!m_state->is_parsing()) {
    return {};
  }

  AR_REQUIRE(chunk->data.capacity() && chunk->filename);
  AR_REQUIRE(chunk->data_type == reader_chunk::type::text,
             "~a12ddd10369e1b5c~");
  m_reader.add(std::move(chunk->data), chunk->end_of_file);

  return emit();
}

void
fastq_parser::finalize()
{
  AR_REQUIRE(m_reader.eof() || m_state->is_parsing_at_head(),
             "~1939d44ea0df7a34~");
  // Release any pending buffers to `reusable_objects`
  m_reader.close();
}

chunk_vec
fastq_parser::emit()
{
  chunk_vec chunks;
  while (m_pending_chunks && m_state->is_parsing()) {
    if (auto output = parse_chunk()) {
      chunks.emplace_back(std::move(output));
      m_pending_chunks--;
    } else if (chunks.empty()) {
      chunks.emplace_back(signal_target());
      break;
    } else {
      break;
    }
  }

  return chunks;
}

chunk_ptr
fastq_parser::parse_chunk()
{
  const auto to_parse = std::min<uint64_t>(m_head, m_state->fastqs.capacity());
  fastq_parser::result end;
  try {
    if (m_type == fastq_input_type::interleaved) {
      end = read_interleaved(to_parse);
    } else {
      end = read_single_end(to_parse);
    }
  } catch (const fastq_error& error) {
    std::ostringstream stream;
    stream << "Error reading FASTQ record from "
           << m_reader.names(m_position, m_reader.position()) << ": "
           << error.what();

    throw fastq_error(stream.str());
  }

  std::unique_ptr<fastq_parser_chunk> chunk;
  if (end != buffer_reader::result::partial) {
    AR_REQUIRE(m_chunk_offset == to_parse || end == buffer_reader::result::eof);
    // Check done here to avoid adding "Error reading..." prefix
    if (m_type == fastq_input_type::interleaved && m_chunk_mate == 1) {
      throw fastq_error("Found uneven number of interleaved reads; input files "
                        "may be truncated. Please fix before continuing.");
    }

    chunk = std::make_unique<fastq_parser_chunk>();
    chunk->type = m_type;

    // The move must always be done, so acquired buffers can be released later
    chunk->reads_1 = std::move(m_chunks.at(0));
    chunk->reads_2 = std::move(m_chunks.at(1));

    m_head -= m_chunk_offset;
    m_chunk_offset = 0;

    if (end == buffer_reader::result::eof) {
      m_state->mark_parsing(reader_state::parsing::reached_eof);
      chunk->eof = true;
    } else if (m_head == 0) {
      m_state->mark_parsing(reader_state::parsing::reached_head);
      chunk->eof = true;
    }
  }

  return chunk;
}

fastq_parser::result
fastq_parser::read_single_end(const uint64_t to_parse)
{
  auto& dst = m_chunks.at(m_chunk_mate);
  if (dst.empty()) {
    dst = m_state->fastqs.acquire(to_parse);
  }

  auto end = fastq_parser::result::complete;
  while (m_chunk_offset < to_parse) {
    end = parse_fastq(dst.at(m_chunk_offset));

    if (end == buffer_reader::result::complete) {
      m_position = m_reader.position();
      m_chunk_offset++;
    } else {
      break;
    }
  }

  if (end != buffer_reader::result::partial) {
    dst.resize(std::min(dst.size(), m_chunk_offset));
  }

  return end;
}

fastq_parser::result
fastq_parser::read_interleaved(const uint64_t to_parse)
{
  if (m_chunks.front().empty()) {
    m_chunks.at(0) = m_state->fastqs.acquire(to_parse);
    m_chunks.at(1) = m_state->fastqs.acquire(to_parse);
  }

  AR_REQUIRE(m_chunk_mate <= 1);
  auto end = fastq_parser::result::complete;
  while (m_chunk_offset < to_parse) {
    end = parse_fastq(m_chunks.at(m_chunk_mate).at(m_chunk_offset));

    if (end == buffer_reader::result::complete) {
      m_position = m_reader.position();
      m_chunk_offset += m_chunk_mate; // inc after reading mate 2
      m_chunk_mate ^= 1;
    } else {
      break;
    }
  }

  if (end != buffer_reader::result::partial) {
    m_chunks.at(0).resize(std::min(m_chunks.at(0).size(), m_chunk_offset));
    m_chunks.at(1).resize(std::min(m_chunks.at(1).size(), m_chunk_offset));
  }

  return end;
}

buffer_reader::result
fastq_parser::parse_fastq(fastq& dst)
{
  while (true) {
    switch (m_chunk_parse) {
      case fastq_parse::header: {
        const auto result = m_reader.getline(dst.m_header);
        if (result == buffer_reader::result::complete) {
          if (dst.m_header.empty()) {
            continue; // skip empty lines between records
          } else if (dst.m_header.size() < 2 || dst.m_header.at(0) != '@') {
            throw fastq_error("Malformed or empty FASTQ header");
          }

          m_chunk_parse = fastq_parse::sequence;
          break;
        } else {
          return result;
        }
      }
      case fastq_parse::sequence: {
        switch (const auto result = m_reader.getline(dst.m_sequence)) {
          case buffer_reader::result::complete:
            if (dst.m_sequence.empty()) {
              throw fastq_error("sequence is empty");
            }

            m_chunk_parse = fastq_parse::separator;
            break;
          case buffer_reader::result::eof:
            throw fastq_error("partial FASTQ record; cut off after header");
          case buffer_reader::result::partial:
          default:
            return result;
        }
        break;
      }
      case fastq_parse::separator: {
        switch (const auto result = m_reader.getline(m_separator)) {
          case buffer_reader::result::complete:
            if (m_separator.empty() || m_separator.front() != '+') {
              throw fastq_error("FASTQ record lacks separator character (+)");
            }

            m_chunk_parse = fastq_parse::qualities;
            break;
          case buffer_reader::result::eof:
            throw fastq_error("partial FASTQ record; cut off after sequence");
          case buffer_reader::result::partial:
          default:
            return result;
        }
        break;
      }
      case fastq_parse::qualities: {
        switch (const auto result = m_reader.getline(dst.m_qualities)) {
          case buffer_reader::result::complete:
            if (dst.m_qualities.length() != dst.m_sequence.length()) {
              throw fastq_error("sequence/quality lengths do not match");
            }

            m_chunk_parse = fastq_parse::header;
            return buffer_reader::result::complete;
          case buffer_reader::result::eof:
            throw fastq_error("partial FASTQ record; cut off after separator");
          case buffer_reader::result::partial:
          default:
            return result;
        }
      }
      default:                                // GCOVR_EXCL_LINE
        AR_FAIL("invalid fastq_parse value"); // GCOVR_EXCL_LINE
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// fastq_source

fastq_source::fastq_source(reader_state_ptr state_1,
                           reader_state_ptr state_2,
                           const step_id next_step)
  : analytical_step(processing_order::ordered, "fastq_source", next_step)
  , m_state_1(std::move(state_1))
  , m_state_2(std::move(state_2))
  , m_capacity(m_state_1 ? m_state_1->fastqs.capacity() : 0)
{
  AR_REQUIRE(m_state_1, "~b80bef5cfb52f770~");
  AR_REQUIRE(m_state_1 != m_state_2, "~a202e7872d3196ac~");
  AR_REQUIRE(!m_state_2 ||
               (m_state_1->fastqs.capacity() == m_state_2->fastqs.capacity()),
             "~b4c19603b9864dc5~");
}

fastq_source::~fastq_source() = default;

chunk_vec
fastq_source::process(chunk_ptr data)
{
  auto chunk = dynamic_cast_unique<fastq_parser_chunk>(data);
  AR_REQUIRE(chunk, "~57885f33c3bd81e8~");
  AR_REQUIRE(m_emitted_chunks < m_requested_chunks, "~e99a9768075c5401~");
  AR_REQUIRE(!m_eof_emitted, "~0c59fee205f2822d~");

  if (chunk->type != fastq_input_type::mate_2) {
    AR_REQUIRE(chunk->type == fastq_input_type::interleaved ||
                 chunk->reads_2.capacity() == 0,
               "~896425048fcc8b67~");
    AR_REQUIRE((chunk->reads_1.size() == m_capacity) ||
                 (chunk->eof && chunk->reads_1.size() < m_capacity),
               "~b9ed5c418b66d2d1~");
  }

  if (chunk->type != fastq_input_type::mate_1) {
    AR_REQUIRE(chunk->type == fastq_input_type::interleaved ||
                 chunk->reads_1.capacity() == 0,
               "~65c9baffc77c7168~");
    AR_REQUIRE((chunk->reads_2.size() == m_capacity) ||
                 (chunk->eof && chunk->reads_2.size() < m_capacity),
               "~226f53009d225eff~");
  }

  if (m_state_2) {
    AR_REQUIRE(chunk->type != fastq_input_type::interleaved,
               "~fc08efb2ce0f1222~");

    if (m_cache.empty() || chunk->type == m_cache.front()->type) {
      AR_REQUIRE(m_cache.size() < m_requested_chunks - m_emitted_chunks,
                 "~dc5ae4cd07004e15~");
      m_cache.emplace_back(std::move(chunk));
    } else {
      if (chunk->type == fastq_input_type::mate_1) {
        chunk->reads_2 = std::move(m_cache.front()->reads_2);
      } else {
        chunk->reads_1 = std::move(m_cache.front()->reads_1);
      }

      if (chunk->reads_1.size() != chunk->reads_2.size() ||
          chunk->eof != m_cache.front()->eof) {
        throw fastq_error("Found unequal number of mate 1 and mate 2 reads; "
                          "input files may be truncated. Please fix before "
                          "continuing.");
      }

      m_cache.pop_front();
    }
  } else if (chunk->type == fastq_input_type::interleaved) {
    AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size(),
               "~31eb2d253b3dd6cd~");
  } else if (chunk->type == fastq_input_type::mate_2) {
    AR_FAIL("~c62527b26be6b7dc~");
  }

  chunk_vec chunks;
  if (chunk) {
    m_emitted_chunks++;
    m_eof_emitted = chunk->eof;
    chunk->first = (m_emitted_chunks == 1);
    chunks.emplace_back(std::move(chunk));
  }

  if (m_suspended) {
    auto can_signal = can_signal_upstream();
    if (can_signal.first || can_signal.second) {
      chunks.emplace_back(scheduler_signal::resume);
      m_suspended = false;
    }
  }

  return chunks;
}

chunk_vec
fastq_source::signal()
{
  AR_REQUIRE(!m_suspended, "~c22ba037ab9d2c3c~");

  chunk_vec chunks;
  auto can_signal = can_signal_upstream();
  if (can_signal.first || can_signal.second) {
    m_requested_chunks++;

    if (can_signal.first) {
      chunks.emplace_back(signal_target(read_mate::_1));
    }

    if (can_signal.second) {
      chunks.emplace_back(signal_target(read_mate::_2));
    }
  }

  if (chunks.empty()) {
    chunks.emplace_back(scheduler_signal::suspend);
    m_suspended = true;
  }

  return chunks;
}

std::pair<bool, bool>
fastq_source::can_signal_upstream() const
{
  std::pair<bool, bool> can_signal;
  // Prevent a slow reader from causing build-up of reads for either mate
  // TODO: limit should be configurable, e.g. 2 * threads, with min value?
  if (m_requested_chunks < m_emitted_chunks + MAX_OUTSTANDING_CHUNKS) {
    // Each signal sent upstream should result in one FASTQ chunk until EOF, but
    // additional signals are allowed as a guard against desynchronization
    can_signal.first = m_state_1->is_parsing();
    can_signal.second = m_state_2 && m_state_2->is_parsing();
  }

  return can_signal;
}

void
fastq_source::finalize()
{
  AR_REQUIRE(!m_state_1->is_parsing(), "~2054659a9668ea62~");
  AR_REQUIRE(!m_state_2 || !m_state_2->is_parsing(), "~197fbe47828be02b~");
  AR_REQUIRE(m_cache.empty(), "~20fc12a3b1803c8c~");
  AR_REQUIRE(m_eof_emitted, "~2341150d380155b0~");
}

////////////////////////////////////////////////////////////////////////////////

void
add_input_steps(scheduler& sch,
                const userconfig& config,
                const reusable_buffers& buffers,
                const reusable_fastqs& fastqs,
                step_id next_step)
{
  const auto& files_1 = config.input_files_1;
  const auto& files_2 = config.input_files_2;
  const bool two_inputs = config.paired_ended_mode && !config.interleaved_input;

  reader_state_ptr state_1 = std::make_shared<reader_state>(buffers, fastqs);
  reader_state_ptr state_2 = two_inputs
                               ? std::make_shared<reader_state>(buffers, fastqs)
                               : reader_state_ptr{};

  next_step = sch.add<fastq_source>(state_1, state_2, next_step);
  sch.set_source_step(next_step);

  const auto file_1_type = config.interleaved_input
                             ? fastq_input_type::interleaved
                             : fastq_input_type::mate_1;
  step_id parse_step = sch.add<fastq_parser>(files_1,
                                             file_1_type,
                                             config.head,
                                             state_1,
                                             next_step);

  auto bgzf_step = sch.add<parallel_decompressor>(state_1, parse_step);
  auto gzip_step = sch.add<serial_decompressor>(state_1, bgzf_step);
  auto read_step = sch.add<block_reader>(files_1, state_1, gzip_step);
  sch.set_signal_target(parse_step, gzip_step);
  sch.set_signal_target(gzip_step, read_step);
  sch.set_signal_target(next_step, parse_step, read_mate::_1);

  if (two_inputs) {
    parse_step = sch.add<fastq_parser>(files_2,
                                       fastq_input_type::mate_2,
                                       config.head,
                                       state_2,
                                       next_step);

    bgzf_step = sch.add<parallel_decompressor>(state_2, parse_step);
    gzip_step = sch.add<serial_decompressor>(state_2, bgzf_step);
    read_step = sch.add<block_reader>(files_2, state_2, gzip_step);
    sch.set_signal_target(parse_step, gzip_step);
    sch.set_signal_target(gzip_step, read_step);
    sch.set_signal_target(next_step, parse_step, read_mate::_2);
  }
}

} // namespace adapterremoval
