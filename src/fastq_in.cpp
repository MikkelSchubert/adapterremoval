// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "fastq_in.hpp"   // declarations
#include "bgzf.hpp"       // for BGZF_HEADER, BGZF_MAX_BLOCK_SIZE, ...
#include "buffer.hpp"     // for buffer
#include "debug.hpp"      // for AR_FAIL, AR_REQUIRE
#include "deflate.hpp"    // for checked_isal_inflate, inflate_state
#include "errors.hpp"     // for fastq_error, gzip_error
#include "logging.hpp"    // for info, warn
#include "managed_io.hpp" // for managed_reader
#include "scheduler.hpp"  // for analytical_chunk, analytical_step, ...
#include "strutils.hpp"   // for log_escape
#include "utilities.hpp"  // for dynamic_cast_unique
#include <algorithm>      // for min
#include <cstddef>        // for size_t
#include <cstdint>        // for uint8_t
#include <deque>          // for deque
#include <isa-l.h>        // for inflate_state
#include <memory>         // for make_shared, make_unique

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// reader_state

reader_state::reader_state(reusable_buffers buffers_, reusable_fastqs fastqs_)
  : buffers(std::move(buffers_))
  , fastqs(std::move(fastqs_))
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

reader_chunk::reader_chunk(reusable_buffers buffers_)
  : buffers(std::move(buffers_))
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

} // namespace adapterremoval
