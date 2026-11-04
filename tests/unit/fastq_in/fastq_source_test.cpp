// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "../scheduler_test.hpp" // for STEP_PROCESS_TO_CHUNK, ...
#include "commontypes.hpp"       // for read_mate
#include "errors.hpp"            // for assert_failed, fastq_error
#include "fastq_in.hpp"          // for MAX_OUTSTANDING_CHUNKS, ...
#include "fastq_in_test.hpp"     // for DEFAULT_BUFFER_CAPACITY, ...
#include "scheduler.hpp"         // for chunk_ptr, chunk_variant, ...
#include "strutils.hpp"          // for stringify
#include "testing.hpp"           // for GENERATE, REQUIRE, REQUIRE_FALSE, ...
#include "testutils.hpp"         // for random_fastq_records
#include <memory>                // for make_unique, unique_ptr
#include <tuple>                 // for ignore
#include <utility>               // for move, pair, swap
#include <vector>                // for vector

namespace adapterremoval {

namespace {

//! Default signal target for mate 1 upstream/source
const step_id DEFAULT_SIGNAL_TARGET_1 = make_step_id_for_testing(0xDEAD);
//! Default signal target for mate 2 upstream/source
const step_id DEFAULT_SIGNAL_TARGET_2 = make_step_id_for_testing(0xBEEF);

std::unique_ptr<fastq_parser_chunk>
make_fastq_parser_chunk(const reader_state_ptr& state,
                        fastq_input_type type,
                        chunk_eof eof = chunk_eof::no)
{
  AR_REQUIRE(state);
  auto chunk = std::make_unique<fastq_parser_chunk>();
  chunk->eof = (eof == chunk_eof::yes);
  chunk->type = type;

  if (type != fastq_input_type::mate_2) {
    chunk->reads_1 = random_fastq_records(state->fastqs.capacity(), "/1");
  }

  if (type != fastq_input_type::mate_1) {
    chunk->reads_2 = random_fastq_records(state->fastqs.capacity(), "/2");
  }

  return chunk;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// fastq_source

TEST_CASE("fastq_source initialization", "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  REQUIRE(step.next_step() == DEFAULT_NEXT_STEP);
  REQUIRE(step.name() == "fastq_source");
  REQUIRE(step.ordering() == processing_order::ordered);
}

TEST_CASE("fastq_source throws on nullptr state 1", "[fastq_source]")
{
  REQUIRE_THROWS_CONTAINS(fastq_source(nullptr, nullptr, DEFAULT_NEXT_STEP),
                          assert_failed,
                          "~b80bef5cfb52f770~");

  auto state_2 = make_state();
  REQUIRE_THROWS_CONTAINS(fastq_source(nullptr, state_2, DEFAULT_NEXT_STEP),
                          assert_failed,
                          "~b80bef5cfb52f770~");
}

TEST_CASE("fastq_source permits nullptr state 2", "[fastq_source]")
{
  auto state_1 = make_state();
  REQUIRE_NOTHROW(fastq_source(state_1, nullptr, DEFAULT_NEXT_STEP));
}

TEST_CASE("fastq_source throws if state 1 is state 2", "[fastq_source]")
{
  auto state = make_state();
  REQUIRE_THROWS_CONTAINS(fastq_source(state, state, DEFAULT_NEXT_STEP),
                          assert_failed,
                          "~a202e7872d3196ac~");
}

TEST_CASE("fastq_source throws if fastq capacities differ", "[fastq_source]")
{
  SECTION("different objects, same capacity is OK")
  {
    auto state_1 = make_state();
    auto state_2 = make_state();
    REQUIRE_NOTHROW(fastq_source{ state_1, state_2, DEFAULT_NEXT_STEP });
  }

  SECTION("different objects, different capacity is rejected")
  {
    auto state_1 = make_state(DEFAULT_BUFFER_CAPACITY, 1234);
    auto state_2 = make_state();
    REQUIRE_THROWS_CONTAINS(fastq_source(state_1, state_2, DEFAULT_NEXT_STEP),
                            assert_failed,
                            "~b4c19603b9864dc5~");
  }
}

TEST_CASE("fastq_source throws on nullptr SE chunk", "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(step.process(nullptr),
                          assert_failed,
                          "~57885f33c3bd81e8~");
}

TEST_CASE("fastq_source throws on nullptr PE chunk", "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_CONTAINS(step.process(nullptr),
                          assert_failed,
                          "~57885f33c3bd81e8~");
}

TEMPLATE_TEST_CASE("fastq_source throws on process(non-fastq_parser_chunk)",
                   "[fastq_source]",
                   fastq_chunk,
                   output_chunk)
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  std::unique_ptr<fastq_source> step;
  if (GENERATE(true, false)) {
    step = std::make_unique<fastq_source>(state_1, nullptr, DEFAULT_NEXT_STEP);
  } else {
    step = std::make_unique<fastq_source>(state_1, state_2, DEFAULT_NEXT_STEP);
  }

  REQUIRE_THROWS_CONTAINS(step->process(std::make_unique<TestType>()),
                          assert_failed,
                          "~57885f33c3bd81e8~");
}

TEST_CASE("fastq_source throws on process(reader_chunk)", "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  std::unique_ptr<fastq_source> step;
  if (GENERATE(true, false)) {
    step = std::make_unique<fastq_source>(state_1, nullptr, DEFAULT_NEXT_STEP);
  } else {
    step = std::make_unique<fastq_source>(state_1, state_2, DEFAULT_NEXT_STEP);
  }

  auto chunk = std::make_unique<fastq_chunk>();
  REQUIRE_THROWS_CONTAINS(step->process(std::move(chunk)),
                          assert_failed,
                          "~57885f33c3bd81e8~");
}

TEST_CASE("fastq_source::signal returns up to N signals in SE mode",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_1, read_mate::_1);

  for (size_t i = 0; i < MAX_OUTSTANDING_CHUNKS; ++i) {
    const auto signal = STEP_SIGNAL_TO(step_id, step);
    REQUIRE(signal == DEFAULT_SIGNAL_TARGET_1);
  }

  const auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);

  REQUIRE_THROWS_CONTAINS(step.signal(), assert_failed, "~c22ba037ab9d2c3c~");
}

TEST_CASE("fastq_source::signal resumes after receiving SE chunks",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_1, read_mate::_1);

  for (size_t i = 0; i < MAX_OUTSTANDING_CHUNKS; ++i) {
    std::ignore = STEP_SIGNAL_TO(step_id, step);
  }

  auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);

  auto chunk = make_fastq_parser_chunk(state, fastq_input_type::mate_1);
  auto output = step.process(std::move(chunk));
  REQUIRE(output.size() == 2);
  REQUIRE(output.at(1) == chunk_variant{ scheduler_signal::resume });

  const size_t finished_chunks = GENERATE(range(0, 4));
  for (size_t i = 0; i < finished_chunks; ++i) {
    chunk = make_fastq_parser_chunk(state, fastq_input_type::mate_1);
    std::ignore =
      STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(chunk));
  }

  for (size_t i = 0; i < finished_chunks + 1; ++i) {
    std::ignore = STEP_SIGNAL_TO(step_id, step);
  }

  signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);
}

TEST_CASE("fastq_source::signal returns up to N signals for in PE mode",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_1, read_mate::_1);
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_2, read_mate::_2);

  for (size_t i = 0; i < MAX_OUTSTANDING_CHUNKS; ++i) {
    const auto signals = step.signal();
    REQUIRE(signals.size() == 2);
    REQUIRE(signals.at(0) == chunk_variant{ DEFAULT_SIGNAL_TARGET_1 });
    REQUIRE(signals.at(1) == chunk_variant{ DEFAULT_SIGNAL_TARGET_2 });
  }

  const auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);

  REQUIRE_THROWS_CONTAINS(step.signal(), assert_failed, "~c22ba037ab9d2c3c~");
}

TEST_CASE("fastq_source::signal resumes after receiving chunks",
          "[fastq_source]")
{
  auto state_1 = make_state(DEFAULT_BUFFER_CAPACITY, 1);
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_1, read_mate::_1);
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_2, read_mate::_2);

  size_t i = 0;
  for (; i < MAX_OUTSTANDING_CHUNKS; ++i) {
    const auto signals = step.signal();
    REQUIRE(signals.size() == 2);
    REQUIRE(signals.at(0) == chunk_variant{ DEFAULT_SIGNAL_TARGET_1 });
    REQUIRE(signals.at(1) == chunk_variant{ DEFAULT_SIGNAL_TARGET_2 });
  }

  auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);

  auto chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
  // NOLINTNEXTLINE(readability-container-size-empty)
  REQUIRE(step.process(std::move(chunk)) == chunk_vec{});

  chunk = make_fastq_parser_chunk(state_2, fastq_input_type::mate_2);
  auto output = step.process(std::move(chunk));
  REQUIRE(output.size() == 2);
  REQUIRE(output.at(1) == chunk_variant{ scheduler_signal::resume });

  const size_t finished_chunks = GENERATE(range(0, 4));
  for (i = 0; i < finished_chunks; ++i) {
    chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
    // NOLINTNEXTLINE(readability-container-size-empty)
    REQUIRE(step.process(std::move(chunk)) == chunk_vec{});

    chunk = make_fastq_parser_chunk(state_2, fastq_input_type::mate_2);
    std::ignore =
      STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(chunk));
  }

  for (i = 0; i < finished_chunks + 1; ++i) {
    const auto signals = step.signal();
    REQUIRE(signals.size() == 2);
    REQUIRE(signals.at(0) == chunk_variant{ DEFAULT_SIGNAL_TARGET_1 });
    REQUIRE(signals.at(1) == chunk_variant{ DEFAULT_SIGNAL_TARGET_2 });
  }

  signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);
}

TEST_CASE("fastq_source::signal keeps signalling other target if one is done",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_1, read_mate::_1);
  set_signal_target_for_testing(step, DEFAULT_SIGNAL_TARGET_2, read_mate::_2);

  size_t i = 0;
  for (; i < MAX_OUTSTANDING_CHUNKS / 2; ++i) {
    const auto signals = step.signal();
    REQUIRE(signals.size() == 2);
    REQUIRE(signals.at(0) == chunk_variant{ DEFAULT_SIGNAL_TARGET_1 });
    REQUIRE(signals.at(1) == chunk_variant{ DEFAULT_SIGNAL_TARGET_2 });
  }

  step_id expected_target{ step_id::UNDEFINED };
  if (GENERATE(true, false)) {
    expected_target = DEFAULT_SIGNAL_TARGET_1;
    state_2->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                   reader_state::parsing::reached_eof));
  } else {
    expected_target = DEFAULT_SIGNAL_TARGET_2;
    state_1->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                   reader_state::parsing::reached_eof));
  }

  for (; i < MAX_OUTSTANDING_CHUNKS; ++i) {
    const auto signal = STEP_SIGNAL_TO(step_id, step);
    REQUIRE(signal == expected_target);
  }

  const auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(signal == scheduler_signal::suspend);

  REQUIRE_THROWS_CONTAINS(step.signal(), assert_failed, "~c22ba037ab9d2c3c~");
}

TEST_CASE("fastq_source::process caches and emits unpaired PE chunks in order",
          "[fastq_source]")
{
  constexpr auto MAX_CHUNKS = MAX_OUTSTANDING_CHUNKS;
  auto state_1 = make_state(DEFAULT_BUFFER_CAPACITY, 1);
  auto state_2 = make_state(state_1);

  // generate fastq chunks in random order, respecting chunk ordering
  auto& rng = Catch::rng();
  std::vector<chunk_ptr> chunks;
  for (size_t i = 0, j = 0; i < MAX_CHUNKS || j < MAX_CHUNKS;) {
    fastq_chunk_ptr chunk;
    if (i != MAX_CHUNKS && (j == MAX_CHUNKS || rng() & 1)) {
      chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
      chunk->reads_1.at(0) = fastq("@read:" + stringify(i) + "/1", "A", "1");
      chunk->eof = (i + 1 == MAX_CHUNKS);
      i++;
    } else {
      chunk = make_fastq_parser_chunk(state_2, fastq_input_type::mate_2);
      chunk->reads_2.at(0) = fastq("@read:" + stringify(j) + "/2", "C", "2");
      chunk->eof = (j + 1 == MAX_CHUNKS);
      j++;
    }

    chunks.emplace_back(std::move(chunk));
  }

  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  for (size_t i = 0; i < MAX_CHUNKS; ++i) {
    std::ignore = step.signal();
  }

  std::vector<fastq_chunk_ptr> outputs;
  for (auto& chunk : chunks) {
    auto output = step.process(std::move(chunk));
    for (auto& output_chunk : output) {
      auto* output_chunk_ptr = std::get_if<chunk_ptr>(&output_chunk);
      REQUIRE(output_chunk_ptr);
      auto fastq_chunk_ptr =
        dynamic_cast_unique<fastq_chunk>(*output_chunk_ptr);
      REQUIRE(fastq_chunk_ptr.get());
      outputs.emplace_back(std::move(fastq_chunk_ptr));
    }
  }

  REQUIRE(outputs.size() == MAX_CHUNKS);
  for (size_t i = 0; i < MAX_CHUNKS; ++i) {
    auto& chunk = *outputs.at(i);
    REQUIRE(
      chunk.reads_1 ==
      std::vector<fastq>{ fastq{ "@read:" + stringify(i) + "/1", "A", "1" } });
    REQUIRE(
      chunk.reads_2 ==
      std::vector<fastq>{ fastq{ "@read:" + stringify(i) + "/2", "C", "2" } });
    REQUIRE(chunk.eof == (i + 1 == MAX_CHUNKS));
  }
}

TEST_CASE("fastq_source::process emits SE chunks immediately", "[fastq_source]")
{
  constexpr size_t MAX_CHUNKS = 10;
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  for (size_t i = 0; i < MAX_CHUNKS; ++i) {
    std::ignore = STEP_SIGNAL_TO(step_id, step);
    auto chunk = make_fastq_parser_chunk(state, fastq_input_type::mate_1);
    chunk->eof = (i + 1 == MAX_CHUNKS);
    const auto reads = chunk->reads_1;
    auto out_chunk = STEP_PROCESS_TO_CHUNK(fastq_chunk, step, std::move(chunk));
    REQUIRE(out_chunk->reads_1 == reads);
    REQUIRE(out_chunk->eof == (i + 1 == MAX_CHUNKS));
  }
}

TEST_CASE("fastq_source::process emits interleaved chunks immediately",
          "[fastq_source]")
{
  constexpr size_t MAX_CHUNKS = 10;
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  for (size_t i = 0; i < MAX_CHUNKS; ++i) {
    std::ignore = STEP_SIGNAL_TO(step_id, step);
    auto chunk = make_fastq_parser_chunk(state, fastq_input_type::interleaved);
    chunk->eof = (i + 1 == MAX_CHUNKS);
    const auto reads_1 = chunk->reads_1;
    const auto reads_2 = chunk->reads_2;
    auto out_chunk = STEP_PROCESS_TO_CHUNK(fastq_chunk, step, std::move(chunk));
    REQUIRE(out_chunk->reads_1 == reads_1);
    REQUIRE(out_chunk->reads_2 == reads_2);
    REQUIRE(out_chunk->eof == (i + 1 == MAX_CHUNKS));
  }
}

TEST_CASE("fastq_source::process requires that interleaved read counts match",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  std::ignore = STEP_SIGNAL_TO(step_id, step);
  // other checks catch non-EOF blocks
  auto chunk = make_fastq_parser_chunk(state,
                                       fastq_input_type::interleaved,
                                       chunk_eof::yes);

  if (GENERATE(true, false)) {
    chunk->reads_1.pop_back();
  } else {
    chunk->reads_2.pop_back();
  }

  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~31eb2d253b3dd6cd~");
}

TEST_CASE("fastq_source::process requires that mate 1/2 read counts match",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  using pair = std::pair<size_t, size_t>;
  const auto reads =
    GENERATE(pair{ 0, 1 },
             pair{ 1, 0 },
             pair{ DEFAULT_FASTQ_CAPACITY, DEFAULT_FASTQ_CAPACITY - 1 },
             pair{ DEFAULT_FASTQ_CAPACITY - 1, DEFAULT_FASTQ_CAPACITY });

  auto chunk_1 =
    make_fastq_parser_chunk(state_1, fastq_input_type::mate_1, chunk_eof::yes);
  chunk_1->reads_1.resize(reads.first);
  auto chunk_2 =
    make_fastq_parser_chunk(state_2, fastq_input_type::mate_2, chunk_eof::yes);
  chunk_2->reads_2.resize(reads.second);

  if (GENERATE(true, false)) {
    std::swap(chunk_1, chunk_2);
  }

  std::ignore = step.signal();
  std::ignore = step.process(std::move(chunk_1));
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk_2)),
                          fastq_error,
                          "Found unequal number of mate 1 and mate 2 reads");
}

TEST_CASE("fastq_source::process rejects too large SE chunks", "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };

  const size_t reads =
    DEFAULT_FASTQ_CAPACITY + GENERATE(1, 2, DEFAULT_FASTQ_CAPACITY);

  auto chunk = make_fastq_parser_chunk(state, fastq_input_type::mate_1);
  chunk->reads_1.resize(reads);
  chunk->eof = GENERATE(true, false);

  std::ignore = step.signal();
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~b9ed5c418b66d2d1~");
}

TEST_CASE("fastq_source::process rejects too large interleaved chunks",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };

  const size_t reads =
    DEFAULT_FASTQ_CAPACITY + GENERATE(1, 2, DEFAULT_FASTQ_CAPACITY);

  auto chunk = make_fastq_parser_chunk(state, fastq_input_type::interleaved);
  chunk->reads_1.resize(reads);
  chunk->reads_2.resize(reads);
  chunk->eof = GENERATE(true, false);

  std::ignore = step.signal();
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~b9ed5c418b66d2d1~");
}

TEST_CASE("fastq_source::process rejects too large PE chunks", "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  std::ignore = step.signal();

  const auto eof = GENERATE(true, false);
  const size_t reads =
    DEFAULT_FASTQ_CAPACITY + GENERATE(1, 2, DEFAULT_FASTQ_CAPACITY);

  SECTION("mate 1")
  {
    auto chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
    chunk->reads_1.resize(reads);
    chunk->eof = eof;

    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~b9ed5c418b66d2d1~");
  }

  SECTION("mate 2")
  {
    auto chunk = make_fastq_parser_chunk(state_2, fastq_input_type::mate_2);
    chunk->reads_2.resize(reads);
    chunk->eof = eof;

    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~226f53009d225eff~");
  }
}

TEST_CASE("fastq_source::process requires fixed-size PE chunks before EOF",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  std::ignore = step.signal();

  const int delta = GENERATE(-DEFAULT_FASTQ_CAPACITY, -32, -1, 1, 32);
  const size_t reads = static_cast<int>(DEFAULT_FASTQ_CAPACITY) + delta;

  SECTION("mate 1")
  {
    auto chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
    chunk->reads_1.resize(reads);
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~b9ed5c418b66d2d1~");
  }

  SECTION("mate 2")
  {
    auto chunk = make_fastq_parser_chunk(state_2, fastq_input_type::mate_2);
    chunk->reads_2.resize(reads);
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~226f53009d225eff~");
  }
}

TEST_CASE("fastq_source::process requires that mate 1/2 EOF matches",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  auto chunk_1 = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
  auto chunk_2 = make_fastq_parser_chunk(state_2, fastq_input_type::mate_2);

  chunk_1->eof = GENERATE(true, false);
  chunk_2->eof = !chunk_1->eof;

  if (GENERATE(true, false)) {
    std::swap(chunk_1, chunk_2);
  }

  std::ignore = step.signal();
  std::ignore = step.process(std::move(chunk_1));
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk_2)),
                          fastq_error,
                          "Found unequal number of mate 1 and mate 2 reads");
}

TEST_CASE("fastq_source::process rejects SE chunks if none were requested",
          "[fastq_source]")
{
  const auto mate =
    GENERATE(fastq_input_type::mate_1, fastq_input_type::interleaved);

  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  auto chunk = make_fastq_parser_chunk(state, mate);

  SECTION("signalled")
  {
    std::ignore = step.signal();
    REQUIRE_NOTHROW(step.process(std::move(chunk)));
  }

  SECTION("not signalled")
  {
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~e99a9768075c5401~");
  }
}

TEST_CASE("fastq_source::process rejects PE chunks if none were requested",
          "[fastq_source]")
{
  const auto mate =
    GENERATE(fastq_input_type::mate_1, fastq_input_type::mate_2);
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  auto chunk = make_fastq_parser_chunk(state_1, mate);

  SECTION("signalled")
  {
    std::ignore = step.signal();
    REQUIRE_NOTHROW(step.process(std::move(chunk)));
  }

  SECTION("not signalled")
  {
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~e99a9768075c5401~");
  }
}

TEST_CASE("fastq_source::process does not allow chunk cache to grow forever",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  const auto mate =
    GENERATE(fastq_input_type::mate_1, fastq_input_type::mate_2);

  std::ignore = step.signal();
  REQUIRE_NOTHROW(step.process(make_fastq_parser_chunk(state_1, mate)));
  REQUIRE_THROWS_CONTAINS(step.process(make_fastq_parser_chunk(state_1, mate)),
                          assert_failed,
                          "~dc5ae4cd07004e15~");
}

TEST_CASE("fastq_source::process rejects mate 2 chunks if SE mode",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  auto chunk = make_fastq_parser_chunk(state, fastq_input_type::mate_2);

  std::ignore = step.signal();
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~c62527b26be6b7dc~");
}

TEST_CASE("fastq_source::process rejects interleaved chunks if PE mode",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  auto chunk = make_fastq_parser_chunk(state_1, fastq_input_type::interleaved);
  std::ignore = step.signal();

  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~fc08efb2ce0f1222~");
}

TEST_CASE("fastq_source::process rejects unexpected reads capacity",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  std::ignore = step.signal();

  auto chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
  chunk->reads_1.reserve(state_1->fastqs.capacity());
  chunk->reads_2.reserve(state_2->fastqs.capacity());

  SECTION("mate 1")
  {
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~896425048fcc8b67~");
  }

  SECTION("mate 2")
  {
    chunk->type = fastq_input_type::mate_2;
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~65c9baffc77c7168~");
  }
}

TEST_CASE("fastq_source::process rejects reads after EOF", "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, nullptr, DEFAULT_NEXT_STEP };

  const auto initial_eof = GENERATE(chunk_eof::no, chunk_eof::yes);
  const auto type =
    GENERATE(fastq_input_type::mate_1, fastq_input_type::interleaved);
  std::ignore = step.signal();
  std::ignore = step.process(make_fastq_parser_chunk(state, type, initial_eof));

  std::ignore = step.signal();
  auto eof = GENERATE(chunk_eof::no, chunk_eof::yes);
  auto chunk = make_fastq_parser_chunk(state, type, eof);

  if (initial_eof == chunk_eof::yes) {
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~0c59fee205f2822d~");
  } else {
    REQUIRE_NOTHROW(step.process(std::move(chunk)));
  }
}

TEST_CASE("fastq_source::process rejects reads after EOF in PE mode",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  const auto initial_eof = GENERATE(chunk_eof::no, chunk_eof::yes);
  std::ignore = step.signal();
  std::ignore = step.process(
    make_fastq_parser_chunk(state_1, fastq_input_type::mate_1, initial_eof));
  std::ignore = step.process(
    make_fastq_parser_chunk(state_1, fastq_input_type::mate_2, initial_eof));

  std::ignore = step.signal();
  auto eof = GENERATE(chunk_eof::no, chunk_eof::yes);
  auto mate = GENERATE(fastq_input_type::mate_1, fastq_input_type::mate_2);
  auto chunk = make_fastq_parser_chunk(state_1, mate, eof);

  if (initial_eof == chunk_eof::yes) {
    REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                            assert_failed,
                            "~0c59fee205f2822d~");
  } else {
    REQUIRE_NOTHROW(step.process(std::move(chunk)));
  }
}

TEST_CASE("fastq_source::finalize checks that parsing is done",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };

  auto type = GENERATE(fastq_input_type::mate_1, fastq_input_type::interleaved);
  std::ignore = step.signal();
  std::ignore =
    step.process(make_fastq_parser_chunk(state, type, chunk_eof::yes));

  SECTION("not done parsing")
  {
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~2054659a9668ea62~");
  }

  SECTION("done parsing")
  {
    state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));
    REQUIRE_NOTHROW(step.finalize());
  }
}

TEST_CASE("fastq_source::finalize checks that all parsing is done in PE mode",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  std::ignore = step.signal();
  std::ignore = step.process(
    make_fastq_parser_chunk(state_1, fastq_input_type::mate_1, chunk_eof::yes));
  std::ignore = step.process(
    make_fastq_parser_chunk(state_1, fastq_input_type::mate_2, chunk_eof::yes));

  SECTION("not done parsing")
  {
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~2054659a9668ea62~");
  }

  SECTION("only source 1 is done parsing")
  {
    state_1->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                   reader_state::parsing::reached_eof));
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~197fbe47828be02b~");
  }

  SECTION("only source 2 is done parsing")
  {
    state_2->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                   reader_state::parsing::reached_eof));
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~2054659a9668ea62~");
  }

  SECTION("done parsing")
  {
    // state mismatch is allowed, so that head can handle truncated inputs
    state_1->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                   reader_state::parsing::reached_eof));
    state_2->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                   reader_state::parsing::reached_eof));
    REQUIRE_NOTHROW(step.finalize());
  }
}

TEST_CASE("fastq_source::finalize checks that there are no cached chunks",
          "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };
  auto chunk = make_fastq_parser_chunk(
    state_1,
    GENERATE(fastq_input_type::mate_1, fastq_input_type::mate_2));

  std::ignore = step.signal();
  // NOLINTNEXTLINE(readability-container-size-empty)
  REQUIRE(step.process(std::move(chunk)) == chunk_vec{});

  state_1->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));
  state_2->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));

  REQUIRE_THROWS_CONTAINS(step.finalize(), assert_failed, "~20fc12a3b1803c8c~");
}

TEST_CASE("fastq_source::finalize checks that EOF was reached",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, nullptr, DEFAULT_NEXT_STEP };

  auto eof = GENERATE(true, false);
  auto type = GENERATE(fastq_input_type::mate_1, fastq_input_type::interleaved);
  auto chunk = make_fastq_parser_chunk(state, type);
  chunk->eof = eof;

  std::ignore = step.signal();
  std::ignore = step.process(std::move(chunk));

  state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                               reader_state::parsing::reached_eof));

  if (eof) {
    REQUIRE_NOTHROW(step.finalize());
  } else {
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~2341150d380155b0~");
  }
}

TEST_CASE("fastq_source::finalize checks that EOF was reached for PE reads",
          "[fastq_source]")
{
  const auto eof = GENERATE(chunk_eof::no, chunk_eof::yes);

  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  std::ignore = step.signal();
  std::ignore = step.process(
    make_fastq_parser_chunk(state_1, fastq_input_type::mate_1, eof));
  std::ignore = step.process(
    make_fastq_parser_chunk(state_2, fastq_input_type::mate_2, eof));

  state_1->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));
  state_2->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));

  if (eof == chunk_eof::yes) {
    REQUIRE_NOTHROW(step.finalize());
  } else {
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~2341150d380155b0~");
  }
}

TEST_CASE("fastq_source marks first chunk for SE/interleaved data",
          "[fastq_source]")
{
  auto state = make_state();
  fastq_source step{ state, {}, DEFAULT_NEXT_STEP };
  const auto mate =
    GENERATE(fastq_input_type::mate_1, fastq_input_type::interleaved);

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto in_chunk = make_fastq_parser_chunk(state, mate);
  auto output = STEP_PROCESS_TO_CHUNK(fastq_chunk, step, std::move(in_chunk));
  REQUIRE(output->first);

  for (size_t i = 0; i < 10; ++i) {
    std::ignore = STEP_SIGNAL_TO(step_id, step);
    in_chunk = make_fastq_parser_chunk(state, mate);
    output = STEP_PROCESS_TO_CHUNK(fastq_chunk, step, std::move(in_chunk));
    REQUIRE_FALSE(output->first);
  }
}

TEST_CASE("fastq_source marks first chunk for PE data", "[fastq_source]")
{
  auto state_1 = make_state();
  auto state_2 = make_state(state_1);
  fastq_source step{ state_1, state_2, DEFAULT_NEXT_STEP };

  std::ignore = step.signal();
  auto in_chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
  STEP_PROCESS_TO_NOTHING(step, std::move(in_chunk));
  in_chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_2);
  auto output = STEP_PROCESS_TO_CHUNK(fastq_chunk, step, std::move(in_chunk));
  REQUIRE(output->first);

  for (size_t i = 0; i < 10; ++i) {
    std::ignore = step.signal();
    in_chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_1);
    STEP_PROCESS_TO_NOTHING(step, std::move(in_chunk));
    in_chunk = make_fastq_parser_chunk(state_1, fastq_input_type::mate_2);
    output = STEP_PROCESS_TO_CHUNK(fastq_chunk, step, std::move(in_chunk));
    REQUIRE_FALSE(output->first);
  }
}

} // namespace adapterremoval
