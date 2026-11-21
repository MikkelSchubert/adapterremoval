// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "../scheduler_test.hpp" // for step_process_to_chunk
#include "bgzf.hpp"              // for BGZF_EOF, BGZF_MAX_BLOCK_SIZE, ...
#include "errors.hpp"            // for assert_failed, gzip_error
#include "fastq_in.hpp"          // for parallel_decompressor, reader_chunk, ...
#include "fastq_in_test.hpp"     // for DEFAULT_NEXT_STEP, make_state, ...
#include "logging.hpp"           // for level, log_capture, set_level
#include "scheduler.hpp"         // for chunk_vec, fastq_chunk, ...
#include "testing.hpp"           // for TEST_CASE, GENERATE, REQUIRE, ..
#include "testutils.hpp"         // for bgzf_block_into_buffer, ...
#include <cstddef>               // for size_t
#include <memory>                // for make_unique, shared_ptr
#include <string>                // for string

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////

namespace {

reader_chunk_ptr
make_compressed_chunk(std::shared_ptr<reader_state>& state,
                      chunk_eof eof,
                      std::shared_ptr<std::string> filename = {})
{
  return make_reader_chunk(state, eof, reader_chunk::type::bgzf, filename);
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// parallel_decompressor

TEST_CASE("parallel_decompressor throws on nullptr state",
          "[parallel_decompressor]")
{
  std::shared_ptr<reader_state> state;
  REQUIRE_THROWS_CONTAINS(parallel_decompressor(state, DEFAULT_NEXT_STEP),
                          assert_failed,
                          "~c4b40d0026ce83ec~");
}

TEST_CASE("parallel_decompressor sets analytical_step properties",
          "[parallel_decompressor]")
{
  auto state = make_state();
  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };

  REQUIRE(proc.next_step() == DEFAULT_NEXT_STEP);
  REQUIRE(proc.name() == "parallel_decompressor");
  REQUIRE(proc.ordering() == processing_order::unordered);
}

TEST_CASE("parallel_decompressor throws on signal", "[parallel_decompressor]")
{
  auto state = make_state();
  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(proc.signal(),
                          assert_failed,
                          "signal() called on parallel_decompressor");
}

TEST_CASE("parallel_decompressor throws on process(nullptr)",
          "[parallel_decompressor]")
{
  auto state = make_state();
  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(proc.process({}),
                          assert_failed,
                          "~8a5ad2119fe9e996~");
}

TEMPLATE_TEST_CASE("parallel_decompressor throws on process(non-reader_chunk)",
                   "[parallel_decompressor]",
                   fastq_chunk,
                   output_chunk)
{
  auto state = make_state();
  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(proc.process(std::make_unique<TestType>()),
                          assert_failed,
                          "~8a5ad2119fe9e996~");
}

TEST_CASE("parallel_decompressor throws on unknown chunk",
          "[parallel_decompressor]")
{
  auto state = make_state();

  auto chunk = make_compressed_chunk(state, chunk_eof::no);
  chunk->data_type = reader_chunk::type::unknown;

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  REQUIRE_THROWS_CONTAINS(proc.process(std::move(chunk)),
                          assert_failed,
                          "~aa2642c2b9c3b1d6~");

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor forwards uncompressed reader_chunks",
          "[parallel_decompressor]")
{
  // sizes must not trigger min/max-BGZF-size checks
  const auto size =
    GENERATE(0, 1, 27, 28, 1024, BGZF_MAX_BLOCK_SIZE, BGZF_MAX_BLOCK_SIZE + 1);
  const auto text = random_text_buffer(size);
  const auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  const auto filename = std::make_shared<std::string>("some filename");

  auto state = make_state();
  auto chunk = make_compressed_chunk(state, eof, filename);
  chunk->data_type = reader_chunk::type::text;
  chunk->data.append(text);
  const auto* data_ptr = chunk->data.data();

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));
  REQUIRE(chunk->filename == filename);
  REQUIRE(chunk->end_of_file == (eof == chunk_eof::yes));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.data() == data_ptr);
  REQUIRE(chunk->data == text);

  chunk.reset();

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor forwards chunks if parsing is done",
          "[parallel_decompressor]")
{
  auto state = make_state();
  state->mark_parsing(GENERATE(reader_state::parsing::reached_eof,
                               reader_state::parsing::reached_head));

  buffer data = random_text_buffer(1024);
  auto data_type = reader_chunk::type::text;
  auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  auto filename = std::make_shared<std::string>("some filename");

  if (GENERATE(true, false)) {
    data_type = reader_chunk::type::bgzf;
    data = bgzf_block_into_buffer(data);
  }

  auto chunk = make_compressed_chunk(state, eof, filename);
  chunk->data_type = data_type;
  chunk->data.append(data);

  const auto* chunk_inner_ptr = chunk.get();

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));
  REQUIRE(chunk.get() == chunk_inner_ptr);
  REQUIRE(chunk->data_type == data_type);
  REQUIRE(chunk->data == data);
  REQUIRE(chunk->end_of_file == (eof == chunk_eof::yes));
  REQUIRE(chunk->filename == filename);
  chunk.reset();

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor decompresses BGZF chunks",
          "[parallel_decompressor]")
{
  const auto text = random_text_buffer(1024);
  const auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  const auto filename = std::make_shared<std::string>("some filename");

  auto state = make_state();
  auto chunk = make_compressed_chunk(state, eof, filename);
  chunk->data.append(bgzf_block_into_buffer(text));

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));
  REQUIRE(chunk->filename == filename);
  REQUIRE(chunk->end_of_file == (eof == chunk_eof::yes));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data == text);
  chunk.reset();

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor acquires chunks", "[parallel_decompressor]")
{
  auto state = make_state();
  auto text = random_text_buffer(1024);
  auto chunk = make_compressed_chunk(state, chunk_eof::no);
  chunk->data.append(bgzf_block_into_buffer(text));

  auto prealloc = state->buffers.acquire();
  const auto* prealloc_ptr = prealloc.data();
  state->buffers.release(std::move(prealloc));

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));
  REQUIRE(chunk->data.data() == prealloc_ptr);
  chunk.reset();

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor requires minimum buffer capacity for chunks",
          "[parallel_decompressor]")
{
  const size_t buffer_capacity =
    GENERATE(1, 2, BGZF_MAX_BLOCK_SIZE - 2, BGZF_MAX_BLOCK_SIZE - 1);
  auto state = make_state(buffer_capacity);

  REQUIRE_THROWS_CONTAINS(parallel_decompressor(state, DEFAULT_NEXT_STEP),
                          assert_failed,
                          "~27c490c08fca8088~");
}

TEST_CASE("parallel_decompressor uses default capacity for acquired chunks",
          "[parallel_decompressor]")
{
  const size_t buffer_capacity =
    BGZF_MAX_BLOCK_SIZE + GENERATE(0, 1, 1024, 4096);
  auto state = make_state(buffer_capacity);
  auto text = random_text_buffer(1024);
  auto chunk = make_compressed_chunk(state, chunk_eof::no);
  chunk->data.append(bgzf_block_into_buffer(text));

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));
  REQUIRE(chunk->data.capacity() == buffer_capacity);
  chunk.reset();

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor releases chunks", "[parallel_decompressor]")
{
  auto state = make_state();
  auto text = random_text_buffer(1024);
  auto chunk = make_compressed_chunk(state, chunk_eof::no);
  chunk->data.append(bgzf_block_into_buffer(text));
  const auto* data_ptr = chunk->data.data();

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));

  auto released = state->buffers.acquire();
  REQUIRE(released.data() == data_ptr);
  state->buffers.release(std::move(released));
  chunk.reset();

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor requires filename (any value)",
          "[parallel_decompressor]")
{
  auto state = make_state();
  auto chunk = make_compressed_chunk(state, chunk_eof::no);
  chunk->filename.reset();

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  REQUIRE_THROWS_CONTAINS(proc.process(std::move(chunk)),
                          assert_failed,
                          "~245a6e94975b7266~");

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor includes filename in decompression errors",
          "[parallel_decompressor]")
{
  auto state = make_state();
  auto chunk = make_compressed_chunk(state, chunk_eof::no);
  chunk->filename =
    std::make_shared<const std::string>("b8ec7a5a-ca0c-\r-943a-f8a4c2388aac");

  if (GENERATE(true, false)) {
    chunk->data.append(random_data_buffer(1024)); // format error
  } else {
    chunk->data.append(bgzf_block_into_buffer(random_text_buffer(1024)));
    chunk->data.put_u8(32, chunk->data.at(32) ^ 0xFF); // deflate error
  }

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  REQUIRE_THROWS_CONTAINS(proc.process(std::move(chunk)),
                          gzip_error,
                          "error while decompressing "
                          "'b8ec7a5a-ca0c-\\r-943a-f8a4c2388aac'");

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor reuses chunk objects",
          "[parallel_decompressor]")
{
  auto state = make_state();
  if (GENERATE(true, false)) {
    state->mark_parsing(GENERATE(reader_state::parsing::reached_eof,
                                 reader_state::parsing::reached_head));
  }

  auto chunk = make_compressed_chunk(state, chunk_eof::no);
  chunk->data_type =
    GENERATE(reader_chunk::type::text, reader_chunk::type::bgzf);
  chunk->data.append(BGZF_EOF);
  const auto* chunk_inner_ptr = chunk.get();

  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };
  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));
  REQUIRE(chunk.get() == chunk_inner_ptr);
}

TEST_CASE("parallel_decompressor decompresses BGZF EOF",
          "[parallel_decompressor]")
{
  auto state = make_state();
  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };

  auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  auto filename = std::make_shared<std::string>("some filename");
  auto chunk = make_compressed_chunk(state, eof, filename);
  chunk->data.append(BGZF_EOF);

  chunk = STEP_PROCESS_TO_CHUNK(reader_chunk, proc, std::move(chunk));
  REQUIRE(chunk->data_type == reader_chunk::type::text);
  REQUIRE(chunk->data.capacity() != 0);
  REQUIRE(chunk->data.size() == 0);
  REQUIRE(chunk->filename == filename);
  REQUIRE(chunk->end_of_file == (eof == chunk_eof::yes));
  chunk.reset();

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("parallel_decompressor rejects BGZF buffers smaller than BGZF EOF",
          "[parallel_decompressor]")
{
  auto state = make_state();
  parallel_decompressor proc{ state, DEFAULT_NEXT_STEP };

  auto eof = GENERATE(chunk_eof::yes, chunk_eof::no);
  auto filename = std::make_shared<std::string>("some filename");
  auto chunk = make_compressed_chunk(state, eof, filename);
  chunk->data.append(BGZF_EOF);
  chunk->data.resize(BGZF_EOF.size() - GENERATE(1, 2, 27, 28));

  REQUIRE_THROWS_CONTAINS(proc.process(std::move(chunk)),
                          assert_failed,
                          "~410943f18b07dde0~");

  REQUIRE_NOTHROW(state->buffers.finalize());
}

} // namespace adapterremoval
