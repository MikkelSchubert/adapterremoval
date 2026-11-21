// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "fastq_in_test.hpp" // declarations
#include "debug.hpp"         // for AR_REQUIRE
#include "errors.hpp"        // for assert_failed
#include "fastq_in.hpp"      // for reader_state
#include "testing.hpp"       // for make_step_id_for_testing
#include <cstddef>           // for size_t

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<reader_state>
make_state(const size_t buffer_capacity, const size_t fastq_capacity)
{
  reusable_buffers buffers{ "buffers", buffer_capacity };
  reusable_fastqs fastqs{ "fastqs", fastq_capacity };

  return std::make_shared<reader_state>(buffers, fastqs);
}

std::shared_ptr<reader_state>
make_state(const std::shared_ptr<reader_state>& source)
{
  AR_REQUIRE(source);
  return std::make_shared<reader_state>(source->buffers, source->fastqs);
}

reader_chunk_ptr
make_reader_chunk(const reader_state_ptr& state,
                  chunk_eof eof,
                  reader_chunk::type type,
                  std::shared_ptr<const std::string> filename)
{
  if (!filename) {
    filename = std::make_shared<std::string>("dummy filename");
  }

  auto chunk = std::make_unique<reader_chunk>(state->buffers);
  chunk->data_type = type;
  chunk->filename = std::move(filename);
  chunk->end_of_file = eof == chunk_eof::yes;
  chunk->data = state->buffers.acquire(0);

  return chunk;
}

////////////////////////////////////////////////////////////////////////////////
// reader_state

TEST_CASE("reader_state does not permit setting parsing::yes", "[reader_state]")
{
  reusable_buffers buffers{ "buffers", DEFAULT_BUFFER_CAPACITY };
  reusable_fastqs fastqs{ "fastqs", DEFAULT_FASTQ_CAPACITY };
  reader_state state{ buffers, fastqs };

  const auto initial = GENERATE(reader_state::parsing::yes,
                                reader_state::parsing::reached_head,
                                reader_state::parsing::reached_eof);
  if (initial != reader_state::parsing::yes) {
    REQUIRE_NOTHROW(state.mark_parsing(initial));
  }

  REQUIRE_THROWS_CONTAINS(state.mark_parsing(reader_state::parsing::yes),
                          assert_failed,
                          "~2fb83d9febd61f6f~");
}

TEST_CASE("reader_state does not permit setting multiple parsing states",
          "[reader_state]")
{
  reusable_buffers buffers{ "buffers", DEFAULT_BUFFER_CAPACITY };
  reusable_fastqs fastqs{ "fastqs", DEFAULT_FASTQ_CAPACITY };
  reader_state state{ buffers, fastqs };

  const auto first = GENERATE(reader_state::parsing::reached_head,
                              reader_state::parsing::reached_eof);
  const auto second = GENERATE(reader_state::parsing::reached_head,
                               reader_state::parsing::reached_eof);

  REQUIRE_NOTHROW(state.mark_parsing(first));
  REQUIRE_THROWS_CONTAINS(state.mark_parsing(second),
                          assert_failed,
                          "~6e7c8cea6a4745bc~");
  REQUIRE_THROWS_CONTAINS(state.mark_parsing(reader_state::parsing::yes),
                          assert_failed,
                          "~2fb83d9febd61f6f~");
}

////////////////////////////////////////////////////////////////////////////////
// reader_chunk

TEST_CASE("reader_chunk::type to string", "[reader_chunk]")
{
  using ::Catch::fallbackStringifier;
  CHECK(fallbackStringifier(reader_chunk::type::unknown) ==
        "reader_chunk::type::unknown");
  CHECK(fallbackStringifier(reader_chunk::type::text) ==
        "reader_chunk::type::text");
  CHECK(fallbackStringifier(reader_chunk::type::bgzf) ==
        "reader_chunk::type::bgzf");
}

} // namespace adapterremoval
