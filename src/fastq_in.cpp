// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "fastq_in.hpp"   // declarations
#include "buffer.hpp"     // for buffer
#include "debug.hpp"      // for AR_FAIL, AR_REQUIRE
#include "logging.hpp"    // for log
#include "managed_io.hpp" // for managed_reader
#include "scheduler.hpp"  // for analytical_chunk, analytical_step, ...
#include "strutils.hpp"   // for log_escape
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

} // namespace adapterremoval
