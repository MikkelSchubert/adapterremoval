// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "../scheduler_test.hpp" // for step_signal_to_chunk
#include "buffer.hpp"            // for buffer
#include "errors.hpp"            // for assert_failed
#include "fastq_in.hpp"          // for block_reader, reader_chunk_ptr, ...
#include "fastq_in_test.hpp"     // for DEFAULT_BUFFER_CAPACITY, ...
#include "logging.hpp"           // for log_capture
#include "managed_io.hpp"        // for managed_reader
#include "scheduler.hpp"         // for step_id
#include "strutils.hpp"          // for ends_with
#include "testing.hpp"           // for GENERATE, REQUIRE, REQUIRE_FALSE, ..
#include "testutils.hpp"         // for random_data_buffer, temporary_file
#include <filesystem>            // for temp_directory_path
#include <memory>                // for make_unique, unique_ptr
#include <random>                // for random_device
#include <sstream>               // for ostringstream
#include <string>                // for string
#include <system_error>          // for error_code
#include <tuple>                 // for ignore
#include <utility>               // for move, pair
#include <vector>                // for vector

namespace adapterremoval {

namespace {

/** Block reader that reads content from temporary files created on demand */
class mock_reader : public block_reader
{
public:
  mock_reader(const std::vector<std::string>& filenames,
              std::vector<buffer> buffers,
              reader_state_ptr state,
              step_id next_step)
    : block_reader(filenames, std::move(state), next_step)
    , m_buffers(std::move(buffers))
    , m_buffer(m_buffers.begin())
  {
    AR_REQUIRE(m_filenames.size() == m_buffers.size());
  }

  buffer total_input()
  {
    buffer copy;
    for (const auto& it : m_buffers) {
      copy.append(it);
    }

    return copy;
  }

protected:
  std::unique_ptr<managed_reader> open_reader(std::string filename) override
  {
    AR_REQUIRE(!m_filenames.empty());
    AR_REQUIRE(*m_filenames.back() == filename);
    AR_REQUIRE(m_buffer != m_buffers.end());

    FILE* handle = temporary_file(*m_buffer++);
    AR_REQUIRE(handle);

    return std::make_unique<managed_reader>(handle);
  }

  const std::vector<buffer> m_buffers{};
  std::vector<buffer>::const_iterator m_buffer{};
};

auto
make_reader(const std::vector<std::string>& filenames,
            std::vector<buffer> buffers)
{
  auto state = make_state();
  return std::pair(std::make_unique<mock_reader>(filenames,
                                                 std::move(buffers),
                                                 state,
                                                 DEFAULT_NEXT_STEP),
                   std::move(state));
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// block_reader

TEST_CASE("block_reader initialization", "[block_reader]")
{
  auto state = make_state();
  block_reader step{ { "filename" }, state, DEFAULT_NEXT_STEP };
  REQUIRE(step.next_step() == DEFAULT_NEXT_STEP);
  REQUIRE(step.name() == "block_reader");
  REQUIRE(step.ordering() == processing_order::ordered_io);
}

TEST_CASE("block_reader requires non-empty filenames", "[block_reader]")
{
  auto state = make_state();
  REQUIRE_THROWS_CONTAINS(mock_reader({}, {}, state, DEFAULT_NEXT_STEP),
                          assert_failed,
                          "~b5bb8b2a98f75df5~");
}

TEST_CASE("block_reader asserts on nullptr state", "[block_reader]")
{
  std::vector<std::string> filenames{ "dummy" };
  std::vector<buffer> buffers;
  buffers.emplace_back();

  REQUIRE_THROWS_CONTAINS(
    mock_reader(filenames, std::move(buffers), nullptr, DEFAULT_NEXT_STEP),
    assert_failed,
    "~1cab54b36e0edde6~");
}

TEST_CASE("block_reader asserts on process", "[block_reader]")
{
  std::vector<std::string> filenames{ "dummy" };
  std::vector<buffer> buffers;
  buffers.resize(1);
  auto [reader, state] = make_reader(filenames, std::move(buffers));

  REQUIRE_THROWS_CONTAINS(reader->process(nullptr),
                          assert_failed,
                          "process() called on block_reader");
}

#define PARAMETERIZE_FILE_SIZES                                                \
  GENERATE(0,                                                                  \
           1,                                                                  \
           DEFAULT_BUFFER_CAPACITY - 1,                                        \
           DEFAULT_BUFFER_CAPACITY,                                            \
           DEFAULT_BUFFER_CAPACITY + 1,                                        \
           (2 * DEFAULT_BUFFER_CAPACITY) - 1,                                  \
           (2 * DEFAULT_BUFFER_CAPACITY),                                      \
           (2 * DEFAULT_BUFFER_CAPACITY) + 1)

TEST_CASE("block_reader across files and blocks", "[block_reader]")
{
  const size_t file_count = GENERATE(1, 2);
  std::vector<std::string> filenames;
  std::vector<buffer> buffers;
  std::vector<size_t> sizes;

  sizes.push_back(PARAMETERIZE_FILE_SIZES);
  if (file_count > 1) {
    sizes.push_back(PARAMETERIZE_FILE_SIZES);
  }

  for (size_t i = 0; i < file_count; ++i) {
    const auto file_size = sizes.at(i);
    filenames.emplace_back("filename" + stringify(i));
    buffers.emplace_back(random_data_buffer(file_size));
  }

  auto [reader, state] = make_reader(filenames, std::move(buffers));
  buffer total_output;

  reader_chunk_ptr chunk;
  for (size_t i = 0; i < file_count; ++i) {
    const auto& filename = filenames.at(i);
    const auto filesize = sizes.at(i);

    size_t nread = 0;
    size_t nchunks = 0;
    while (nread < filesize) {
      chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
      REQUIRE(chunk->filename);
      REQUIRE(*chunk->filename == filename);
      REQUIRE(chunk->data_type == reader_chunk::type::unknown);
      REQUIRE(chunk->data.capacity() != 0);

      if (chunk->end_of_file) {
        REQUIRE(chunk->data.size() == filesize % DEFAULT_BUFFER_CAPACITY);
      } else {
        REQUIRE(chunk->data.size() == DEFAULT_BUFFER_CAPACITY);
      }

      nchunks++;
      nread += chunk->data.size();
      total_output.append(chunk->data);
    }

    if (filesize % DEFAULT_BUFFER_CAPACITY == 0) { // this includes empty files
      chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
      REQUIRE(chunk->filename);
      REQUIRE(*chunk->filename == filename);
      REQUIRE(chunk->data_type == reader_chunk::type::unknown);
      REQUIRE(chunk->data.capacity() != 0);
      REQUIRE(chunk->data.size() == 0);
      nchunks++;
    }

    REQUIRE(nchunks == (filesize / DEFAULT_BUFFER_CAPACITY) + 1);
    REQUIRE(nread == filesize);
    REQUIRE(chunk->end_of_file);
    const bool last_chunk = i + 1 == file_count;
    REQUIRE(state->is_reading() != last_chunk);
  }

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(!chunk->filename);
  REQUIRE(chunk->data_type == reader_chunk::type::unknown);
  REQUIRE(chunk->data.capacity() == 0);
  REQUIRE_FALSE(chunk->end_of_file);
  REQUIRE_FALSE(state->is_reading());

  chunk.reset();
  REQUIRE_NOTHROW(state->buffers.finalize());
  REQUIRE(total_output == reader->total_input());
}

TEST_CASE("block_reader returns blank chunks after EOF", "[block_reader]")
{
  auto [reader, state] = make_reader({ "dummy" }, std::vector<buffer>(1));

  auto chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(chunk->filename);
  REQUIRE(chunk->data_type == reader_chunk::type::unknown);
  REQUIRE(chunk->data.capacity() != 0);
  REQUIRE(chunk->data.size() == 0);
  REQUIRE(chunk->end_of_file);
  REQUIRE_FALSE(state->is_reading());

  for (size_t i = 0; i < 3; ++i) {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE(!chunk->filename);
    REQUIRE(chunk->data_type == reader_chunk::type::unknown);
    REQUIRE(chunk->data.capacity() == 0);
    REQUIRE_FALSE(chunk->end_of_file);
    REQUIRE_FALSE(state->is_reading());
  }

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("block_reader stops reading when parsing is done", "[block_reader]")
{
  std::vector<buffer> buffers;
  buffers.emplace_back(
    random_data_buffer((DEFAULT_BUFFER_CAPACITY * 3) - 1024));
  auto [reader, state] = make_reader({ "file1" }, std::move(buffers));

  log::log_capture cap;
  reader_chunk_ptr chunk;
  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(*chunk->filename == "file1");
  REQUIRE(chunk->data.size() == DEFAULT_BUFFER_CAPACITY);
  REQUIRE_FALSE(chunk->end_of_file);
  REQUIRE(state->is_reading());

  state->mark_parsing(reader_state::parsing::reached_head);

  if (GENERATE(true, false)) {
    const auto signal = STEP_SIGNAL_TO(scheduler_signal, *reader);
    REQUIRE(signal == scheduler_signal::suspend);
    REQUIRE(cap.str() == "[INFO] Opening FASTQ file 'file1'\n"
                         "[INFO] Closing FASTQ file 'file1'\n");
  }

  REQUIRE_NOTHROW(reader->finalize());
  REQUIRE(cap.str() == "[INFO] Opening FASTQ file 'file1'\n"
                       "[INFO] Closing FASTQ file 'file1'\n");

  chunk.reset();
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("block_reader handles signal after parser reaches EOF",
          "[block_reader]")
{
  std::vector<buffer> buffers;
  buffers.emplace_back(
    random_data_buffer(GENERATE(0, DEFAULT_BUFFER_CAPACITY - 1)));
  auto [reader, state] = make_reader({ "file1" }, std::move(buffers));

  std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(!state->is_reading());

  state->mark_parsing(reader_state::parsing::reached_eof);
  auto signal = STEP_SIGNAL_TO(scheduler_signal, *reader);
  REQUIRE(signal == scheduler_signal::suspend);
}

TEST_CASE("block_reader throws if parser reaches EOF first", "[block_reader]")
{
  std::vector<buffer> buffers;
  buffers.emplace_back(random_data_buffer(DEFAULT_BUFFER_CAPACITY * 2));
  auto [reader, state] = make_reader({ "file1" }, std::move(buffers));

  if (GENERATE(true, false)) {
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE(state->is_reading());
  }

  state->mark_parsing(reader_state::parsing::reached_eof);
  REQUIRE_THROWS_CONTAINS(reader->signal(),
                          assert_failed,
                          "~4a20e925ef1d96b4~");
}

TEST_CASE("block_reader logs opening/closing files", "[block_reader]")
{
  std::vector<std::string> filenames{ "file\r1", "file 2", "file3" };
  std::vector<buffer> buffers;
  buffers.resize(3);
  buffers.at(0) = random_data_buffer((DEFAULT_BUFFER_CAPACITY * 2) - 1);
  buffers.at(2) = random_data_buffer((DEFAULT_BUFFER_CAPACITY * 2) - 1);
  auto [reader, state] = make_reader(filenames, std::move(buffers));

  log::log_capture cap;

  reader_chunk_ptr chunk;
  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(*chunk->filename == "file\r1");
  REQUIRE(cap.str() == "[INFO] Opening FASTQ file 'file\\r1'\n");
  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(*chunk->filename == "file\r1");
  REQUIRE(ends_with(cap.str(),
                    "[INFO] Opening FASTQ file 'file\\r1'\n"
                    "[INFO] Closing FASTQ file 'file\\r1'\n"));

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(*chunk->filename == "file 2");
  REQUIRE(ends_with(cap.str(),
                    "[INFO] Closing FASTQ file 'file\\r1'\n"
                    "[INFO] Opening FASTQ file 'file 2'\n"
                    "[INFO] Closing FASTQ file 'file 2'\n"));

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(*chunk->filename == "file3");
  REQUIRE(ends_with(cap.str(),
                    "[INFO] Closing FASTQ file 'file 2'\n"
                    "[INFO] Opening FASTQ file 'file3'\n"));
  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(*chunk->filename == "file3");
  REQUIRE(ends_with(cap.str(),
                    "[INFO] Opening FASTQ file 'file3'\n"
                    "[INFO] Closing FASTQ file 'file3'\n"));

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  REQUIRE(!chunk->filename);
  REQUIRE(ends_with(cap.str(),
                    "[INFO] Opening FASTQ file 'file3'\n"
                    "[INFO] Closing FASTQ file 'file3'\n"));

  REQUIRE_FALSE(state->is_reading());
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("block_reader finalize requires EOF or parsing is done",
          "[block_reader]")
{
  std::vector<std::string> filenames{ "file1", "file2" };
  std::vector<buffer> buffers;
  buffers.emplace_back(random_data_buffer((DEFAULT_BUFFER_CAPACITY * 2) - 1));
  buffers.emplace_back(random_data_buffer(DEFAULT_BUFFER_CAPACITY - 1));
  auto [reader, state] = make_reader(filenames, std::move(buffers));

  SECTION("throws if neither EOF nor parsing is done")
  {
    REQUIRE_THROWS_CONTAINS(reader->finalize(),
                            assert_failed,
                            "~fd72e01d62e4c694~");
    return;
  }

  SECTION("throws until EOF")
  {
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE_THROWS_CONTAINS(reader->finalize(),
                            assert_failed,
                            "~fd72e01d62e4c694~");
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE_THROWS_CONTAINS(reader->finalize(),
                            assert_failed,
                            "~fd72e01d62e4c694~");
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE_NOTHROW(reader->finalize());
    REQUIRE_NOTHROW(state->buffers.finalize());
  }

  SECTION("parsing reached head before opening first file")
  {
    log::log_capture cap;
    state->mark_parsing(reader_state::parsing::reached_head);
    auto signal = STEP_SIGNAL_TO(scheduler_signal, *reader);

    REQUIRE(signal == scheduler_signal::suspend);
    REQUIRE_NOTHROW(reader->finalize());
    REQUIRE(cap.str() == ""); // NOLINT(readability-container-size-empty)
    REQUIRE_NOTHROW(state->buffers.finalize());
  }

  SECTION("parsing done inside first file")
  {
    log::log_capture cap;
    auto chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE(*chunk->filename == "file1");
    REQUIRE_FALSE(chunk->end_of_file);
    chunk.reset();

    state->mark_parsing(reader_state::parsing::reached_head);
    auto signal = STEP_SIGNAL_TO(scheduler_signal, *reader);
    REQUIRE(signal == scheduler_signal::suspend);

    REQUIRE_NOTHROW(reader->finalize());
    REQUIRE(cap.str() == "[INFO] Opening FASTQ file 'file1'\n"
                         "[INFO] Closing FASTQ file 'file1'\n");
    REQUIRE_NOTHROW(state->buffers.finalize());
  }

  SECTION("parsing done between files")
  {
    log::log_capture cap;
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    auto chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE(*chunk->filename == "file1");
    REQUIRE(chunk->end_of_file);
    chunk.reset();

    state->mark_parsing(reader_state::parsing::reached_head);
    auto signal = STEP_SIGNAL_TO(scheduler_signal, *reader);
    REQUIRE(signal == scheduler_signal::suspend);

    REQUIRE_NOTHROW(reader->finalize());
    REQUIRE(cap.str() == "[INFO] Opening FASTQ file 'file1'\n"
                         "[INFO] Closing FASTQ file 'file1'\n");
    REQUIRE_NOTHROW(state->buffers.finalize());
  }

  SECTION("parsing done after all files")
  {
    log::log_capture cap;
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    auto chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
    REQUIRE(*chunk->filename == "file2");
    REQUIRE(chunk->end_of_file);
    chunk.reset();

    state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));
    auto signal = STEP_SIGNAL_TO(scheduler_signal, *reader);
    REQUIRE(signal == scheduler_signal::suspend);

    REQUIRE_NOTHROW(reader->finalize());
    REQUIRE(cap.str() == "[INFO] Opening FASTQ file 'file1'\n"
                         "[INFO] Closing FASTQ file 'file1'\n"
                         "[INFO] Opening FASTQ file 'file2'\n"
                         "[INFO] Closing FASTQ file 'file2'\n");
    REQUIRE_NOTHROW(state->buffers.finalize());
  }

  REQUIRE_FALSE(state->is_reading());
}

TEST_CASE("block_reader keeps returning suspend after parsing is done",
          "[block_reader]")
{
  std::vector<std::string> filenames{ "file1" };
  std::vector<buffer> buffers;
  buffers.emplace_back(random_data_buffer((DEFAULT_BUFFER_CAPACITY * 2) - 1));
  auto [reader, state] = make_reader(filenames, std::move(buffers));

  if (GENERATE(true, false)) {
    std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, *reader);
  }

  REQUIRE(state->is_reading());
  state->mark_parsing(reader_state::parsing::reached_head);
  for (size_t i = 0; i < 3; ++i) {
    const auto signal = STEP_SIGNAL_TO(scheduler_signal, *reader);
    REQUIRE(signal == scheduler_signal::suspend);
  }

  REQUIRE_NOTHROW(reader->finalize());
  REQUIRE_NOTHROW(state->buffers.finalize());
}

namespace {

/** Returns a random, temporary filename */
std::string
temp_filename()
{
  // Catch::rng() is intentionally not used here, to avoid repeat filenames
  std::random_device rd;
  std::ostringstream os;
  for (size_t i = 0; i < 4; i++) {
    os << rd() << "_";
  }

  return (std::filesystem::temp_directory_path() / os.str()).string();
}

/** Reader class that opens files as write only */
class bad_reader : public block_reader
{
public:
  using block_reader::block_reader;

  ~bad_reader() override
  {
    if (!m_created.empty()) {
      m_reader.reset();
      // use error_code version to prevent remove from throwing
      std::error_code ec;
      // must be done after the reader has closed the handle
      std::filesystem::remove(m_created, ec);
    }
  }

protected:
  std::unique_ptr<managed_reader> open_reader(std::string filename) override
  {
    AR_REQUIRE(m_created.empty());
    FILE* handle = std::fopen(filename.c_str(), "wbx"); // Unreadable handle
    AR_REQUIRE(handle);
    m_created = filename;

    return std::make_unique<managed_reader>(handle);
  }

  std::string m_created{};
};

} // namespace

TEST_CASE("block_reader throws on missing file", "[block_reader]")
{
  auto state = make_state();
  auto filename = temp_filename();
  block_reader reader{ { filename }, state, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(reader.signal(),
                          io_error,
                          "failed to open " + log_escape(filename));
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("block_reader throws on IO errors", "[block_reader]")
{
  auto state = make_state();
  auto filename = temp_filename();

  bad_reader reader{ { filename }, state, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(reader.signal(),
                          io_error,
                          "error reading '<unknown file>': ");
  REQUIRE_NOTHROW(state->buffers.finalize());
}

} // namespace adapterremoval
