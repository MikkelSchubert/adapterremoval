// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "../scheduler_test.hpp" // for  STEP_PROCESS_TO, ...
#include "bgzf.hpp"              // for BGZF_EOF, BGZF_HEADER, ...
#include "commontypes.hpp"       // for read_mate
#include "deflate.hpp"           // for calculate_crc32, deflate_block_isal_1
#include "errors.hpp"            // for assert_failed
#include "fastq_in.hpp"          // for block_reader, reader_chunk_ptr, ...
#include "fastq_in_test.hpp"     // for DEFAULT_BUFFER_SIZE, chunk_eof, ...
#include "logging.hpp"           // for log_capture
#include "scheduler.hpp"         // for chunk_ptr, fastq_chunk, ...
#include "testing.hpp"           // for GENERATE, REQUIRE, REQUIRE_FALSE, ..
#include "testutils.hpp"         // for random_data_buffer, random_text_buffer
#include "utilities.hpp"         // for dynamic_cast_unique
#include <array>                 // for array
#include <memory>                // for make_shared, make_unique, shared_ptr
#include <sstream>               // for ostringstream
#include <string>                // for string

namespace adapterremoval {

namespace {

constexpr size_t DEFAULT_COMPRESSION_LEVEL = 4;

reader_chunk_ptr
make_unknown_chunk(reader_state_ptr& state,
                   chunk_eof eof,
                   std::shared_ptr<const std::string> filename = {})
{
  return make_reader_chunk(state, eof, reader_chunk::type::unknown, filename);
}

void
append_output_to(const reader_chunk_ptr& current, buffer& output)
{
  REQUIRE(current.get());
  if (current->data_type == reader_chunk::type::bgzf) {
    buffer dst;
    bgzf_block_decompress(current->data, dst);
    output.append(dst);
  } else {
    output.append(current->data);
  }
}

} // namespace

TEST_CASE("serial_decompressor initialization", "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  CHECK(step.name() == "serial_decompressor");
  CHECK(step.next_step() == DEFAULT_NEXT_STEP);
  CHECK(step.ordering() == processing_order::ordered);
}

TEST_CASE("serial_decompressor asserts on nullptr state",
          "[serial_decompressor]")
{
  REQUIRE_THROWS_CONTAINS(serial_decompressor(nullptr, DEFAULT_NEXT_STEP),
                          assert_failed,
                          "~79a956c1fd6f2ac5~");
}

TEST_CASE("serial_decompressor requires non-nullptr chunk",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(step.process(nullptr),
                          assert_failed,
                          "~47f0d52258dba7cb~");
}

TEST_CASE("serial_decompressor requires non-nullptr filename",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  auto chunk = make_unknown_chunk(state, chunk_eof::yes);
  chunk->filename.reset();

  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~ad8643ad5a07417a~");
}

TEST_CASE("serial_decompressor::signal sends signals upstream initially",
          "[serial_decompressor]")
{
  const auto target = make_step_id_for_testing(7913);
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  set_signal_target_for_testing(step,
                                make_step_id_for_testing(7913),
                                read_mate::_1);

  const auto signal = STEP_SIGNAL_TO(step_id, step);
  REQUIRE(signal == target);
}

TEST_CASE("serial_decompressor::signal returns chunk while possible",
          "[serial_decompressor]")
{
  const auto target = make_step_id_for_testing(7913);
  auto state = make_state(64);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  set_signal_target_for_testing(step,
                                make_step_id_for_testing(7913),
                                read_mate::_1);

  auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  auto text = random_text_buffer(1024);
  auto chunk = make_unknown_chunk(state, eof);
  gzip_and_append_data_to_buffer(chunk->data, text);

  buffer output;
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  do {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    output.append(chunk->data);
  } while (output.size() < 1024);

  REQUIRE(chunk->end_of_file == (eof == chunk_eof::yes));
  auto signal = STEP_SIGNAL_TO(step_id, step);
  REQUIRE(signal == target);
  REQUIRE(output == text);
}

TEST_CASE("serial_decompressor::signal returns suspend if parsing is done",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                               reader_state::parsing::reached_eof));
  const auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);
}

TEST_CASE("serial_decompressor::signal signals upstream if waiting for chunks",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  const auto signal_target = make_step_id_for_testing(7913);
  set_signal_target_for_testing(step, signal_target, read_mate::_1);

  const auto id = STEP_SIGNAL_TO(step_id, step);
  REQUIRE(id == signal_target);
}

TEMPLATE_TEST_CASE("serial_decompressor::process throws on non-reader_chunk",
                   "[serial_decompressor]",
                   fastq_chunk,
                   output_chunk)
{
  auto state = make_state();
  serial_decompressor proc{ state, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(proc.process(std::make_unique<TestType>()),
                          assert_failed,
                          "~47f0d52258dba7cb~");
}

TEST_CASE("serial_decompressor::process requires unknown data-type",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  auto chunk = std::make_unique<reader_chunk>(state->buffers);
  chunk->data = state->buffers.acquire(0);
  chunk->data_type =
    GENERATE(reader_chunk::type::bgzf, reader_chunk::type::text);

  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~47f0d52258dba7cb~");
}

TEST_CASE("serial_decompressor::process returns nothing if parsing is done",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  auto chunk = std::make_unique<reader_chunk>(state->buffers);

  state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                               reader_state::parsing::reached_eof));
  auto chunks = step.process(std::move(chunk));
  // NOLINTNEXTLINE(readability-container-size-empty)
  REQUIRE(chunks == chunk_vec{});
}

TEST_CASE("serial_decompressor::process drops chunks after done parsing",
          "[serial_decompressor]")
{
  auto state = make_state(1024);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  auto filename = std::make_shared<std::string>("dummy");

  auto chunk = make_unknown_chunk(state, chunk_eof::yes, filename);
  chunk->data.append(random_text_buffer(1024));
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));
  std::ignore = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);

  state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                               reader_state::parsing::reached_eof));
  for (size_t i = 0; i < 3; ++i) {
    chunk = make_unknown_chunk(state, chunk_eof::no, filename);
    chunk->data.resize(chunk->data.capacity());
    STEP_PROCESS_TO_NOTHING(step, std::move(chunk));
    auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
    REQUIRE(signal == scheduler_signal::suspend);
  }

  log::log_capture cap;
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("serial_decompressor skips trailing garbage after gzip",
          "[serial_decompressor]")
{
  constexpr size_t chunk_size = 1024;

  log::log_capture cap;
  auto state = make_state(chunk_size);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  // trailing garbage in the final block with gzip data
  const bool has_padding = GENERATE(true, false);
  // has full-sized garbage chunk after final gzip data
  const size_t has_n_chunks = GENERATE(0, 1, 2);
  // some small amount of garbage after the end
  const bool has_tiny_tail = GENERATE(true, false);
  //! gzip compresesd vs bgzf compressed
  const bool gzip_compressed = GENERATE(true, false);

  std::vector<reader_chunk_ptr> chunks;

  auto eof = (has_n_chunks || has_tiny_tail) ? chunk_eof::no : chunk_eof::yes;
  auto filename = std::make_shared<std::string>("first");
  buffer text = random_text_buffer(chunk_size * 2);

  std::string eof_warning;
  auto chunk = make_unknown_chunk(state, eof, filename);
  if (gzip_compressed) {
    gzip_and_append_data_to_buffer(chunk->data, text);
  } else {
    bgzf_block_compress(text, chunk->data, DEFAULT_COMPRESSION_LEVEL);
    if (GENERATE(true, false)) {
      chunk->data.append(BGZF_EOF);
    } else {
      eof_warning = "[WARNING] BGZF stream in 'first' is missing an EOF "
                    "marker; file may be truncated\n";
    }
  }

  if (has_padding) {
    size_t offset = chunk->data.size();
    REQUIRE(offset + 1 < chunk_size);
    chunk->data.append(random_data_buffer(chunk_size - offset));
    // Ensure that the header can never match by chance
    chunk->data.put_u8(offset, chunk->data.at(offset) | 0x80);
  }

  chunks.emplace_back(std::move(chunk));

  for (size_t i = 1; i <= has_n_chunks; ++i) {
    eof = (has_tiny_tail || i != has_n_chunks) ? chunk_eof::no : chunk_eof::yes;
    chunk = make_unknown_chunk(state, eof, filename);
    chunk->data.append(random_data_buffer(1024));
    // Ensure that the header can never match by chance
    chunk->data.put_u8(0, chunk->data.at(0) | 0x80);

    chunks.emplace_back(std::move(chunk));
  }

  // final (trailing) block
  if (has_tiny_tail) {
    size_t junk_size = GENERATE(0, 1, 2);
    if (!has_padding && !has_n_chunks) {
      junk_size++;
    }

    chunk = make_unknown_chunk(state, chunk_eof::yes, filename);
    chunk->data.append(random_data_buffer(junk_size));

    chunks.emplace_back(std::move(chunk));
  }

  filename = std::make_shared<std::string>("second");
  chunk = make_unknown_chunk(state, chunk_eof::yes, filename);
  chunk->data.append(random_text_buffer(chunk_size / 2));
  text.append(chunk->data);

  if (GENERATE(true, false)) {
    auto gzip = gzip_data_into_buffer(chunk->data);
    chunk->data.resize(0);
    chunk->data.append(gzip);
  }

  chunks.emplace_back(std::move(chunk));

  buffer total_output;
  bool last_chunk_eof = true;
  std::shared_ptr<const std::string> last_chunk_filename;

  for (auto& it : chunks) {
    STEP_PROCESS_TO_NOTHING(step, std::move(it));
  }

  do {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    append_output_to(chunk, total_output);

    REQUIRE((chunk->data_type == reader_chunk::type::text ||
             (*chunk->filename == "first" && !gzip_compressed)));
    REQUIRE((chunk->filename == last_chunk_filename || last_chunk_eof));
    last_chunk_filename = chunk->filename;
    last_chunk_eof = chunk->end_of_file;
  } while (!chunk->end_of_file || *chunk->filename != "second");

  CHECK(chunk->end_of_file);
  CHECK(total_output == text);
  if (has_padding || has_n_chunks || has_tiny_tail) {
    CHECK(cap.str() ==
          "[WARNING] Found trailing garbage at end of 'first'\n" + eof_warning);
  } else {
    // NOLINTNEXTLINE(readability-container-size-empty)
    REQUIRE(cap.str() == eof_warning);
  }
}

TEST_CASE("serial_decompressor handles gzip headers spanning chunks",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto chunk_1 = make_unknown_chunk(state, chunk_eof::no);
  auto chunk_2 = make_unknown_chunk(state, chunk_eof::yes);
  buffer gzip_1_text = random_text_buffer(128);
  chunk_1->data.append(gzip_data_into_buffer(gzip_1_text));

  buffer gzip_2_text = random_text_buffer(128);
  buffer gzip_2 = gzip_data_into_buffer(gzip_2_text);

  const int overlap = GENERATE(range(1, 17));

  SECTION("gzip 1 overlaps chunk 2")
  {
    const size_t new_size = chunk_1->data.size() - overlap;

    chunk_2->data.append(chunk_1->data.data() + new_size, overlap);
    chunk_2->data.append(gzip_2);
    chunk_1->data.resize(new_size);

    STEP_PROCESS_TO_NOTHING(step, std::move(chunk_1));
    STEP_PROCESS_TO_NOTHING(step, std::move(chunk_2));

    auto chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    REQUIRE(chunk->data_type == reader_chunk::type::text);
    REQUIRE_FALSE(chunk->end_of_file);
    REQUIRE(chunk->data == gzip_1_text);

    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    REQUIRE(chunk->data_type == reader_chunk::type::text);
    REQUIRE(chunk->end_of_file);
    REQUIRE(chunk->data == gzip_2_text);
  }

  SECTION("gzip 2 overlaps chunk 1")
  {
    chunk_1->data.append(gzip_2.data(), overlap);
    chunk_2->data.append(gzip_2.data() + overlap, gzip_2.size() - overlap);

    STEP_PROCESS_TO_NOTHING(step, std::move(chunk_1));
    STEP_PROCESS_TO_NOTHING(step, std::move(chunk_2));

    auto chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    REQUIRE(chunk->data_type == reader_chunk::type::text);
    REQUIRE_FALSE(chunk->end_of_file);
    REQUIRE(chunk->data == gzip_1_text);

    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    REQUIRE(chunk->data_type == reader_chunk::type::text);
    REQUIRE(chunk->end_of_file);
    REQUIRE(chunk->data == gzip_2_text);
  }
}

TEST_CASE("serial_decompressor handles bgzip headers spanning chunks",
          "[serial_decompressor]")
{
  log::log_capture cap; // silence EOF warnings
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  buffer bgzip_1 = bgzf_block_into_buffer(random_text_buffer(128));
  buffer bgzip_2 = bgzf_block_into_buffer(random_text_buffer(128));

  auto chunk_1 = make_unknown_chunk(state, chunk_eof::no);
  auto chunk_2 = make_unknown_chunk(state, chunk_eof::yes);

  const int overlap = GENERATE(range(1, 20)); // BGZF_EOF.size() + 2

  SECTION("bgzf 1 overlaps chunk 2")
  {
    chunk_1->data.append(bgzip_1.data(), bgzip_1.size() - overlap);
    chunk_2->data.append(bgzip_1.data() + bgzip_1.size() - overlap, overlap);
    chunk_2->data.append(bgzip_2);
    CHECK(true); // mandatory assert
  }

  SECTION("bgzf 2 overlaps chunk 1")
  {
    chunk_1->data.append(bgzip_1);
    chunk_1->data.append(bgzip_2.data(), overlap);
    chunk_2->data.append(bgzip_2.data() + overlap, bgzip_2.size() - overlap);
    CHECK(true); // mandatory assert
  }

  STEP_PROCESS_TO_NOTHING(step, std::move(chunk_1));
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk_2));

  chunk_1 = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
  REQUIRE(chunk_1->data_type == reader_chunk::type::bgzf);
  REQUIRE_FALSE(chunk_1->end_of_file);
  REQUIRE(chunk_1->data == bgzip_1);

  chunk_2 = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
  REQUIRE(chunk_2->data_type == reader_chunk::type::bgzf);
  REQUIRE(chunk_2->end_of_file);
  REQUIRE(chunk_2->data == bgzip_2);
}

TEST_CASE("serial_decompressor allows empty chunks after final gzip",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto text = random_text_buffer(1024);
  auto chunk = make_unknown_chunk(state, chunk_eof::no);
  chunk->data.append(gzip_data_into_buffer(text));
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));
  STEP_PROCESS_TO_NOTHING(step, make_unknown_chunk(state, chunk_eof::yes));

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
  REQUIRE_FALSE(chunk->end_of_file);
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data == text);

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
  REQUIRE(chunk->end_of_file);
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.capacity() != 0);
  REQUIRE(chunk->data.size() == 0);
}

TEST_CASE("serial_decompressor allows empty chunks after final bgzf",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto data = bgzf_block_into_buffer(random_text_buffer(1024));
  auto chunk = make_unknown_chunk(state, chunk_eof::no);
  chunk->data.append(data);
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));
  STEP_PROCESS_TO_NOTHING(step, make_unknown_chunk(state, chunk_eof::yes));

  log::log_capture cap;

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
  REQUIRE_FALSE(chunk->end_of_file);
  REQUIRE(chunk->data_type == reader_chunk::type::bgzf);
  REQUIRE(chunk->data == data);

  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
  REQUIRE(chunk->end_of_file);
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.capacity() != 0);
  REQUIRE(chunk->data.size() == 0);

  REQUIRE(cap.str() ==
          "[WARNING] BGZF stream in 'dummy filename' is missing an EOF "
          "marker; file may be truncated\n");
}

TEST_CASE("serial_decompressor handles files containing a mix of gzip/bgzf",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  auto chunk = make_unknown_chunk(state, chunk_eof::yes);

  enum class format
  {
    gzip,
    gzip_stream,
    bgzf,
  };

  std::vector<format> blocks;
  blocks.push_back(GENERATE(format::gzip, format::gzip_stream, format::bgzf));
  blocks.push_back(GENERATE(format::gzip, format::gzip_stream, format::bgzf));
  blocks.push_back(GENERATE(format::gzip, format::gzip_stream, format::bgzf));

  buffer input;
  for (auto fmt : blocks) {
    buffer text = random_text_buffer(1024);
    input.append(text);

    switch (fmt) {
      case format::gzip:
        gzip_and_append_data_to_buffer(chunk->data, text);
        break;
      case format::gzip_stream: {
        std::string_view view{ reinterpret_cast<const char*>(text.data()),
                               text.size() };
        chunk->data.append(gzip_stream_into_buffer({
          view.substr(0, 256),
          view.substr(256, 456),
          view.substr(712),
        }));
        break;
      }
      case format::bgzf:
        chunk->data.append(bgzf_block_into_buffer(text));
        break;
      default:
        FAIL("invalid format");
    }
  }

  log::log_capture cap; // silence EOF warnings

  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  buffer output;
  do {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    append_output_to(chunk, output);
  } while (!chunk->end_of_file);

  REQUIRE(output == input);
}

TEST_CASE("serial_decompressor returns full blocks when decompressing gzip",
          "[serial_decompressor]")
{
  auto state = make_state(1024);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto input = random_text_buffer(1280);
  auto gzip = gzip_data_into_buffer(input);
  const auto offset = gzip.size() / 2;

  STEP_SIGNAL_TO(step_id, step);
  auto chunk = make_unknown_chunk(state, chunk_eof::no);
  chunk->data.append(gzip.data(), offset);
  std::ignore = STEP_PROCESS_TO(step_id, step, std::move(chunk));

  chunk = make_unknown_chunk(state, chunk_eof::no);
  chunk->data.append(gzip.data() + offset, gzip.size() - offset);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));

  buffer expected;
  expected.append(input.data(), 1024LLU);
  REQUIRE(chunk->data == expected);

  // reader currently stops at block boundaries
  chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
  expected.resize(0);
  expected.append(input.data() + 1024, input.size() - 1024);
  REQUIRE(chunk->data == expected);
}

TEST_CASE("serial_decompressor returns multiple pending chunks",
          "[serial_decompressor]")
{
  auto state = make_state(512);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  const size_t pending = GENERATE(range(1, 5));
  for (size_t i = 0; i < pending; ++i) {
    std::ignore = STEP_SIGNAL_TO(step_id, step);
  }

  auto input = random_text_buffer((512 * 4) - 1);
  auto gzip = gzip_data_into_buffer(input);
  auto chunk = make_unknown_chunk(state, chunk_eof::yes);
  chunk->data.append(gzip);

  buffer total_output;
  auto chunks = step.process(std::move(chunk));
  REQUIRE(chunks.size() == pending);
  for (auto& variant : chunks) {
    auto* ptr = std::get_if<chunk_ptr>(&variant);
    REQUIRE(ptr);
    chunk = dynamic_cast_unique<reader_chunk>(*ptr);
    REQUIRE(chunk.get());
    total_output.append(chunk->data);
  }

  REQUIRE(total_output.size() <= input.size());
  buffer input_prefix;
  input_prefix.append(input.data(), total_output.size());
  REQUIRE(total_output == input_prefix);
}

namespace {

std::pair<reader_chunk_ptr, buffer>
create_undersized_bgzf_chunk(reader_state_ptr state,
                             chunk_eof eof,
                             shared_filename filename,
                             size_t undersized_by)
{
  REQUIRE(undersized_by <= 27);
  auto chunk = make_unknown_chunk(state, eof, filename);
  auto text = random_text_buffer(1024);
  bgzf_block_compress(text, chunk->data, DEFAULT_COMPRESSION_LEVEL);
  REQUIRE(chunk->data.get_u16(16) + 1LLU == chunk->data.size());
  chunk->data.put_u16(16, BGZF_EOF.size() - 1 - undersized_by);

  return { std::move(chunk), std::move(text) };
}

} // namespace

TEST_CASE("serial_decompressor attempts gzip if BSIZE is too small",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto [chunk, text] =
    create_undersized_bgzf_chunk(state,
                                 GENERATE(chunk_eof::yes, chunk_eof::no),
                                 {},
                                 GENERATE(range(1, 2, 27)));

  log::log_capture cap;
  std::ignore = STEP_SIGNAL_TO(step_id, step);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == text.size());
  REQUIRE(chunk->data == text);
}

TEST_CASE("serial_decompressor shows single BSIZE warning per file",
          "[serial_decompressor]")
{
  auto state = make_state(0x20000);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  log::log_capture cap;

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto filename = std::make_shared<std::string>("dummy filename");
  auto bgzf = create_undersized_bgzf_chunk(state, chunk_eof::no, filename, 1);
  auto chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(bgzf.first));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == bgzf.second.size());
  REQUIRE(chunk->data == bgzf.second);

  // possible regular chunk in between
  if (GENERATE(true, false)) {
    auto data = bgzf_block_into_buffer(random_text_buffer(14057));
    chunk = make_unknown_chunk(state, chunk_eof::no);
    chunk->data.append(data);

    std::ignore = STEP_SIGNAL_TO(step_id, step);
    chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));
    REQUIRE(chunk->data_type == reader_chunk::type::bgzf);
    REQUIRE(chunk->data.size() == data.size());
    REQUIRE(chunk->data == data);
  }

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  bgzf = create_undersized_bgzf_chunk(state, chunk_eof::yes, filename, 27);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(bgzf.first));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == bgzf.second.size());
  REQUIRE(chunk->data == bgzf.second);

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  filename = std::make_shared<std::string>("another filename");
  bgzf = create_undersized_bgzf_chunk(state, chunk_eof::yes, filename, 3);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(bgzf.first));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == bgzf.second.size());
  REQUIRE(chunk->data == bgzf.second);

  // only one warning per filename
  std::ostringstream os;
  os << "[WARNING] BGZF block with a BSIZE that is 1 bytes smaller than the "
        "minimum 28 found in 'dummy filename'; will attempt to decompress as "
        "plain gzip\n"
        "[WARNING] BGZF block with a BSIZE that is 3 bytes smaller than the "
        "minimum 28 found in 'another filename'; will attempt to decompress as "
        "plain gzip\n";

  REQUIRE(cap.str() == os.str());
}

namespace {

std::pair<reader_chunk_ptr, buffer>
create_oversized_bgzf_chunk(reader_state_ptr state,
                            chunk_eof eof,
                            shared_filename filename,
                            size_t oversized_by)
{
  REQUIRE(oversized_by > 0);
  auto chunk = make_unknown_chunk(state, eof, filename);
  chunk->data.append(BGZF_HEADER);
  chunk->data.resize(chunk->data.capacity());

  auto text = random_text_buffer(BGZF_MAX_BLOCK_SIZE + oversized_by);
  deflate_block_isal_1(text, chunk->data, BGZF_HEADER.size(), true);
  chunk->data.append_u32(calculate_crc32(text)); // CRC32
  chunk->data.append_u32(text.size());           // ISIZE

  AR_REQUIRE(chunk->data.size() <= BGZF_MAX_BLOCK_SIZE);
  // Write the final block size; -1 to fit 1..65536 in 16 bit
  chunk->data.put_u16(16, chunk->data.size() - 1);

  return std::pair{ std::move(chunk), std::move(text) };
}

} // namespace

TEST_CASE("serial_decompressor attempts gzip if ISIZE is too large",
          "[serial_decompressor]")
{
  auto state = make_state(0x20000);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  const auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  const auto oversized_by = GENERATE(1, 204);

  auto [chunk, text] =
    create_oversized_bgzf_chunk(state, eof, {}, oversized_by);

  log::log_capture cap;
  std::ignore = STEP_SIGNAL_TO(step_id, step);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == text.size());
  REQUIRE(chunk->data == text);
}

TEST_CASE("serial_decompressor shows single ISIZE warning per file",
          "[serial_decompressor]")
{
  auto state = make_state(0x20000);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  log::log_capture cap;

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto filename = std::make_shared<std::string>("dummy filename");
  auto bgzf = create_oversized_bgzf_chunk(state, chunk_eof::no, filename, 1);
  auto chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(bgzf.first));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == bgzf.second.size());
  REQUIRE(chunk->data == bgzf.second);

  // possible regular chunk in between
  if (GENERATE(true, false)) {
    auto data = bgzf_block_into_buffer(random_text_buffer(14057));
    chunk = make_unknown_chunk(state, chunk_eof::no);
    chunk->data.append(data);

    std::ignore = STEP_SIGNAL_TO(step_id, step);
    chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));
    REQUIRE(chunk->data_type == reader_chunk::type::bgzf);
    REQUIRE(chunk->data.size() == data.size());
    REQUIRE(chunk->data == data);
  }

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  bgzf = create_oversized_bgzf_chunk(state, chunk_eof::yes, filename, 1234);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(bgzf.first));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == bgzf.second.size());
  REQUIRE(chunk->data == bgzf.second);

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  filename = std::make_shared<std::string>("another filename");
  bgzf = create_oversized_bgzf_chunk(state, chunk_eof::yes, filename, 512);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(bgzf.first));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.size() == bgzf.second.size());
  REQUIRE(chunk->data == bgzf.second);

  REQUIRE(cap.str() == "[WARNING] BGZF block with ISIZE 1 bytes larger than "
                       "the maximum 65536 found in 'dummy filename'; will "
                       "attempt to decompress as plain gzip\n"
                       "[WARNING] BGZF block with ISIZE 512 bytes larger than "
                       "the maximum 65536 found in 'another filename'; will "
                       "attempt to decompress as plain gzip\n");
}

TEST_CASE("serial_decompressor only attempts gzip for whole BGZF chunks",
          "[serial_decompressor]")
{
  auto state = make_state(0x20000);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto text = random_text_buffer(BGZF_MAX_BLOCK_SIZE + 1);
  auto chunk = make_unknown_chunk(state, chunk_eof::yes);
  chunk->data.append(BGZF_HEADER);
  chunk->data.resize(chunk->data.capacity());
  deflate_block_isal_1(text, chunk->data, BGZF_HEADER.size(), true);
  chunk->data.append_u32(calculate_crc32(text)); // CRC32
  chunk->data.append_u32(text.size());           // ISIZE

  AR_REQUIRE(chunk->data.size() <= BGZF_MAX_BLOCK_SIZE);
  // Write the final block size; -1 to fit 1..65536 in 16 bit
  chunk->data.put_u16(16, chunk->data.size() - 1);

  log::log_capture cap;
  std::ignore = STEP_SIGNAL_TO(step_id, step);

  SECTION("warns if chunk is whole")
  {
    REQUIRE_NOTHROW(step.process(std::move(chunk)));

    std::ostringstream os;
    os << "[WARNING] BGZF block with ISIZE 1 bytes larger than the maximum "
          "65536 found in 'dummy filename'; will attempt to decompress as "
          "plain gzip\n";

    REQUIRE(cap.str() == os.str());
  }

  SECTION("does not warn if chunk is partial")
  {
    chunk->data.put_u16(16, chunk->data.size());
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            gzip_error,
                            "file is truncated");

    // NOLINTNEXTLINE(readability-container-size-empty)
    REQUIRE(cap.str() == "");
  }
}

TEST_CASE("serial_decompressor throws gzip_error on trailing invalid gzip",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  auto text = random_text_buffer(1024);
  auto chunk = make_unknown_chunk(state, eof);
  chunk->data.append(gzip_data_into_buffer(text));
  chunk->data.append_u8(0x1f); // Gzip header
  chunk->data.append_u8(0x8b);
  chunk->data.append_u8(0xFF); // compression method
  chunk->data.append(random_data_buffer(1024));

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data == text);
  REQUIRE_FALSE(chunk->end_of_file);

  REQUIRE_THROWS_AS(step.signal(), gzip_error);
}

TEST_CASE("serial_decompressor throws gzip_error on truncated gzip",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto chunk = make_unknown_chunk(state, chunk_eof::yes);
  chunk->data.append(gzip_data_into_buffer(random_text_buffer(1024)));
  chunk->data.resize(chunk->data.size() - 32);

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          gzip_error,
                          "file is truncated");
}

TEST_CASE("serial_decompressor::finalize throws if a chunk is being built",
          "[serial_decompressor]")
{
  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto text = random_text_buffer(1024);
  auto chunk = make_unknown_chunk(state, chunk_eof::no);
  if (GENERATE(true, false)) {
    gzip_and_append_data_to_buffer(chunk->data, text);
  } else {
    bgzf_block_compress(text, chunk->data, DEFAULT_COMPRESSION_LEVEL);
  }

  // incomplete stream
  REQUIRE(chunk->data.size() > 32);
  chunk->data.resize(chunk->data.size() - 32);

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  std::ignore = STEP_PROCESS_TO(step_id, step, std::move(chunk));

  // should have no effect
  if (GENERATE(true, false)) {
    state->mark_reading_finished();
  }

  SECTION("still parsing")
  {
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~ab0c8ee7277efef1~");
  }

  SECTION("parser reached EOF first")
  {
    state->mark_parsing(reader_state::parsing::reached_eof);
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~ab0c8ee7277efef1~");
  }

  SECTION("done parsing")
  {
    state->mark_parsing(reader_state::parsing::reached_head);
    REQUIRE_NOTHROW(step.finalize());

    log::log_capture cap;
    REQUIRE_NOTHROW(state->buffers.finalize()); // everything should be released
  }
}

TEST_CASE("serial_decompressor::finalize throws if not all data has been used",
          "[serial_decompressor]")
{
  constexpr size_t buffer_size = 1024;
  auto state = make_state(buffer_size);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto chunk = make_unknown_chunk(state, chunk_eof::no);
  if (GENERATE(true, false)) {
    auto text = random_text_buffer(2 * buffer_size);
    gzip_and_append_data_to_buffer(chunk->data, text);
  } else {
    chunk->data.append(bgzf_block_into_buffer(random_text_buffer(buffer_size)));
    chunk->data.append(bgzf_block_into_buffer(random_text_buffer(buffer_size)));
  }

  STEP_SIGNAL_TO(step_id, step);
  std::ignore = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));

  // should have no effect
  if (GENERATE(true, false)) {
    state->mark_reading_finished();
  }

  SECTION("still parsing")
  {
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~b194f127cdcb5d3e~");
  }

  SECTION("parser reached EOF first")
  {
    state->mark_parsing(reader_state::parsing::reached_eof);
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~b194f127cdcb5d3e~");
  }

  SECTION("done parsing")
  {
    state->mark_parsing(reader_state::parsing::reached_head);
    REQUIRE_NOTHROW(step.finalize());

    log::log_capture cap;
    REQUIRE_NOTHROW(state->buffers.finalize()); // everything should be released
  }
}

TEST_CASE("serial_decompressor handles mixed file input",
          "[serial_decompressor]")
{
  buffer total_input;
  auto make_chunk = [&total_input](reader_state_ptr& state, int type) {
    auto chunk = make_unknown_chunk(state, chunk_eof::yes);
    if (type == 0) {
      // empty file
    } else if (type == 5 || type == 6) {
      // treated as any other file
      chunk->data.append(BGZF_EOF);
    } else {
      buffer text = random_text_buffer(1024);
      total_input.append(text);
      if (type == 1 || type == 2) {
        gzip_and_append_data_to_buffer(chunk->data, text);
      } else if (type == 3 || type == 4) {
        bgzf_block_compress(text, chunk->data, DEFAULT_COMPRESSION_LEVEL);
      } else if (type == 7) {
        chunk->data.append(text);
      } else {
        FAIL("invalid file_type");
      }
    }

    if (type == 2 || type == 4 || type == 6) {
      chunk->data.append_u8(0x80); // must not match gzip magic
      chunk->data.append(random_data_buffer(10));
    }

    return chunk;
  };

  const std::array file_types = { GENERATE(range(0, 8)),
                                  GENERATE(range(0, 8)),
                                  GENERATE(range(0, 8)) };

  auto state = make_state();
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  for (const auto type : file_types) {
    STEP_PROCESS_TO_NOTHING(step, make_chunk(state, type));
  }

  log::log_capture cap;
  buffer total_output;
  for (const auto type : file_types) {
    auto chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    append_output_to(chunk, total_output);

    // trailing garbage results in empty (EOF) chunk
    if (type == 2 || type == 4 || type == 6) {
      REQUIRE_FALSE(chunk->end_of_file);
      chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
      REQUIRE(chunk->data.size() == 0);
    }

    REQUIRE(chunk->end_of_file);
  }

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  REQUIRE(total_input == total_output);
}

TEST_CASE("serial_decompressor handles mixed gzip/BGZF blocks",
          "[serial_decompressor]")
{
  auto make_block = [](buffer& total_input, int type) {
    buffer block;
    if (type == 0) {
      // treated as any other file
      block.append(BGZF_EOF);
    } else {
      buffer text = random_text_buffer(1024);
      total_input.append(text);
      if (type == 1) {
        gzip_and_append_data_to_buffer(block, text);
      } else if (type == 2) {
        bgzf_block_compress(text, block, DEFAULT_COMPRESSION_LEVEL);
      } else {
        FAIL("invalid file_type");
      }
    }

    return block;
  };

  const std::array block_types = { GENERATE(range(0, 3)),
                                   GENERATE(range(0, 3)),
                                   GENERATE(range(0, 3)) };

  buffer total_input;
  auto state = make_state();
  auto chunk = make_unknown_chunk(state, chunk_eof::yes);
  for (auto type : block_types) {
    chunk->data.append(make_block(total_input, type));
  }

  if (GENERATE(true, false)) {
    chunk->data.append_u8(0x80); // must not match gzip magic
    chunk->data.append(random_data_buffer(10));
  }

  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  log::log_capture cap;
  buffer total_output;

  do {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    append_output_to(chunk, total_output);
  } while (!chunk->end_of_file);

  REQUIRE(total_input == total_output);
}

TEST_CASE("serial_decompressor forwards text blocks as is",
          "[serial_decompressor]")
{
  auto state = make_state(1024);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  auto filename = std::make_shared<std::string>("some filename");

  buffer total_output;
  for (const size_t size : std::array{ 123, 0, 1024, 2048, 1025 }) {
    auto text = random_text_buffer(size);
    auto eof = size == 1025 ? chunk_eof::yes : chunk_eof::no;
    auto chunk = make_unknown_chunk(state, eof, filename);
    chunk->data.append(text);

    std::ignore = STEP_SIGNAL_TO(step_id, step);
    chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));
    REQUIRE(chunk->data_type == reader_chunk::type::text);
    REQUIRE(chunk->data == text);
    REQUIRE(chunk->end_of_file == (size == 1025));
    REQUIRE(chunk->filename == filename);
  }

  // EOF should reset filetype after multiple blocks
  buffer bgzf;
  bgzf.append(BGZF_EOF);
  auto chunk = make_unknown_chunk(state, chunk_eof::yes);
  chunk->data.append(bgzf);
  std::ignore = STEP_SIGNAL_TO(step_id, step);
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, step, std::move(chunk));
  REQUIRE(chunk->data_type == reader_chunk::type::bgzf);
  REQUIRE(chunk->data == bgzf);
  REQUIRE(chunk->end_of_file);
}

TEST_CASE("serial_decompressor warns if BGZF stream lacks EOF marker",
          "[serial_decompressor]")
{
  auto state = make_state(1024);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };
  auto text = random_text_buffer(1234);
  auto bgzf = bgzf_block_into_buffer(text);

  std::string message;
  if (GENERATE(true, false)) {
    bgzf.append(BGZF_EOF);
  } else {
    message = "[WARNING] BGZF stream in 'dummy filename' is missing an EOF "
              "marker; file may be truncated\n";
  }

  switch (GENERATE(0, 1, 2)) {
    case 0: // nothing before EOF
      break;
    case 1: // trailing garbage
      bgzf.append_u8(0x80);
      bgzf.append(random_data_buffer(123));
      message.insert(
        0,
        "[WARNING] Found trailing garbage at end of 'dummy filename'\n");
      break;
    case 2: { // trailing gzip file
      auto additional_text = random_text_buffer(123);
      bgzf.append(gzip_data_into_buffer(additional_text));
      text.append(additional_text);
      break;
    }
    default:
      FAIL("not implemented");
  }

  auto chunk = make_unknown_chunk(state, chunk_eof::yes);
  chunk->data.append(bgzf);
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  log::log_capture cap;
  buffer total_output;

  do {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    append_output_to(chunk, total_output);
  } while (!chunk->end_of_file);

  REQUIRE(cap.str() == message);
}

TEST_CASE("serial_decompressor shows single missing-EOF warning per file",
          "[serial_decompressor]")
{
  auto state = make_state(1024);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto text = random_text_buffer(1234);
  auto bgzf = bgzf_block_into_buffer(text);
  auto additional_text = random_text_buffer(123);
  bgzf.append(gzip_data_into_buffer(additional_text));
  text.append(additional_text);
  additional_text = random_text_buffer(123);
  bgzf.append(bgzf_block_into_buffer(additional_text));
  text.append(additional_text);

  auto filename = std::make_shared<std::string>("filename1");
  auto chunk = make_unknown_chunk(state, chunk_eof::yes, filename);
  chunk->data.append(bgzf);
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  filename = std::make_shared<std::string>("filename2");
  chunk = make_unknown_chunk(state, chunk_eof::yes, filename);
  additional_text = random_text_buffer(123);
  chunk->data.append(bgzf_block_into_buffer(additional_text));
  text.append(additional_text);
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  log::log_capture cap;
  buffer total_output;

  do {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    append_output_to(chunk, total_output);
  } while (!chunk->end_of_file || *chunk->filename != "filename2");

  REQUIRE(cap.str() == "[WARNING] BGZF stream in 'filename1' is missing an EOF "
                       "marker; file may be truncated\n"
                       "[WARNING] BGZF stream in 'filename2' is missing an EOF "
                       "marker; file may be truncated\n");
}

TEST_CASE("serial_decompressor skips EOF warning after ISIZE warning",
          "[serial_decompressor]")
{
  auto state = make_state(1024);
  serial_decompressor step{ state, DEFAULT_NEXT_STEP };

  auto text = random_text_buffer(1234);
  auto bgzf = bgzf_block_into_buffer(text);
  auto additional_text = random_text_buffer(123);
  bgzf.append(gzip_data_into_buffer(additional_text));
  text.append(additional_text);
  if (GENERATE(true, false)) {
    additional_text = random_text_buffer(123);
    bgzf.append(bgzf_block_into_buffer(additional_text));
    text.append(additional_text);
  }

  auto filename = std::make_shared<std::string>("filename1");
  auto chunk = make_unknown_chunk(state, chunk_eof::yes, filename);
  chunk->data.append(bgzf);
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  filename = std::make_shared<std::string>("filename2");
  chunk = make_unknown_chunk(state, chunk_eof::yes, filename);
  additional_text = random_text_buffer(123);
  chunk->data.append(bgzf_block_into_buffer(additional_text));
  text.append(additional_text);
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));

  log::log_capture cap;
  buffer total_output;

  do {
    chunk = STEP_SIGNAL_TO_CHUNK(reader_chunk, step);
    append_output_to(chunk, total_output);
  } while (!chunk->end_of_file || *chunk->filename != "filename2");

  REQUIRE(cap.str() == "[WARNING] BGZF stream in 'filename1' is missing an EOF "
                       "marker; file may be truncated\n"
                       "[WARNING] BGZF stream in 'filename2' is missing an EOF "
                       "marker; file may be truncated\n");
}

} // namespace adapterremoval
