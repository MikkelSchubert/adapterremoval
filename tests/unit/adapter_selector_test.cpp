// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_database.hpp" // for adapter_database
#include "adapter_selector.hpp" // for adapter_chunk, adapter_finalizer, ...
#include "catch.hpp"
#include "errors.hpp"        // for assertion_failed
#include "fastq.hpp"         // for fastq
#include "logging.hpp"       // for log_capture
#include "scheduler.hpp"     // for fastq_chunk
#include "sequence_sets.hpp" // for sample_set
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include "threading.hpp"     // for threadsafe_data
#include "utilities.hpp"     // for dynamic_cast_unique
#include <cstdint>
#include <memory> // for make_shared, make_unique

namespace adapterremoval {

namespace {

// Parameterize tests over supported SIMD instruction sets
#define PARAMETERIZE_IS GENERATE(from_range(simd::supported()))

constexpr size_t DEFAULT_NEXT_STEP = 12345;

const double DEFAULT_MISMATCH_THRESHOLD = 1.0 / 6.0;

using hits_vec = std::vector<adapter_detection_stats::hits>;

using Catch::Matchers::Contains;

}; // namespace

////////////////////////////////////////////////////////////////////////////////
// adapter_preselector

TEST_CASE("adapter_preselector requires fastq block")
{
  adapter_preselector step{ DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_AS(step.process(chunk_ptr{}), assert_failed);
  REQUIRE_THROWS_AS(step.process(std::make_unique<output_chunk>()),
                    assert_failed);
  REQUIRE_THROWS_AS(step.process(std::make_unique<adapter_chunk>()),
                    assert_failed);
}

TEST_CASE("process single, empty block")
{
  const bool eof = GENERATE(true, false);
  auto in_chunk = std::make_unique<fastq_chunk>();

  in_chunk->eof = eof;

  adapter_preselector step{ DEFAULT_NEXT_STEP };
  auto out = step.process(std::move(in_chunk));

  REQUIRE(out.size() == 1);
  REQUIRE(out.at(0).first == DEFAULT_NEXT_STEP);

  auto out_chunk = dynamic_cast_unique<adapter_chunk>(out.at(0).second);
  REQUIRE(out_chunk.get());
  REQUIRE(out_chunk->data->reads_1.empty());
  REQUIRE(out_chunk->data->eof == eof);
  REQUIRE(out_chunk->last_adapter_selection_block == eof);
}

TEST_CASE("process single block")
{
  const bool eof = GENERATE(true, false);
  auto in_chunk = std::make_unique<fastq_chunk>();
  in_chunk->reads_1.resize(100, fastq("read", "ACGTTGCTGATCAACTGGTACATATCAAG"));
  in_chunk->eof = eof;

  adapter_preselector step{ DEFAULT_NEXT_STEP };
  auto out = step.process(std::move(in_chunk));
  REQUIRE(out.size() == 1);
  REQUIRE(out.at(0).first == DEFAULT_NEXT_STEP);

  auto out_chunk = dynamic_cast_unique<adapter_chunk>(out.at(0).second);
  REQUIRE(out_chunk.get());
  REQUIRE(out_chunk->data->eof == eof);
  REQUIRE(out_chunk->data->reads_1.size() == 100);
  REQUIRE(out_chunk->last_adapter_selection_block == eof);
}

TEST_CASE("adapter_preselector checks that selection is done")
{
  adapter_preselector step{ DEFAULT_NEXT_STEP };

  SECTION("didn't finish")
  {
    REQUIRE_THROWS_AS(step.finalize(), assert_failed);
  }

  SECTION("reached eof")
  {
    auto chunk = std::make_unique<fastq_chunk>();
    chunk->eof = true;
    step.process(std::move(chunk));
    REQUIRE_NOTHROW(step.finalize());
  }

  SECTION("processed enough reads")
  {
    fastq tmpl{ "read", "ACGTTGCTGATCAACTGGTACATATCAAG" };
    auto in_chunk = std::make_unique<fastq_chunk>();
    in_chunk->reads_1.resize(1000, tmpl);
    for (size_t i = 0; i < 99; ++i) {
      auto out = step.process(std::move(in_chunk));
      auto out_chunk = dynamic_cast_unique<adapter_chunk>(out.at(0).second);
      REQUIRE(out_chunk.get());
      REQUIRE(!out_chunk->last_adapter_selection_block);
      REQUIRE(out_chunk->data.get());
      REQUIRE(out_chunk->data->reads() == 1000);
      in_chunk = std::move(out_chunk->data);
    }

    auto out = step.process(std::move(in_chunk));
    auto out_chunk = dynamic_cast_unique<adapter_chunk>(out.at(0).second);
    REQUIRE(out_chunk.get());
    REQUIRE(out_chunk->last_adapter_selection_block);

    REQUIRE_NOTHROW(step.finalize());
  }
}

TEST_CASE("adapter_preselector returns fastq chunks after selecting")
{
  adapter_preselector step{ DEFAULT_NEXT_STEP };
  auto in_chunk = std::make_unique<fastq_chunk>();
  const auto* in_chunk_ptr = in_chunk.get();
  in_chunk->reads_1.resize(1000,
                           fastq("read", "ACGTTGCTGATCAACTGGTACATATCAAG"));
  for (size_t i = 0; i < 100; ++i) {
    auto out = step.process(std::move(in_chunk));
    auto out_chunk = dynamic_cast_unique<adapter_chunk>(out.at(0).second);
    REQUIRE(out_chunk.get());
    REQUIRE(out_chunk->data.get() == in_chunk_ptr);
    in_chunk = std::move(out_chunk->data);
  }

  auto out = step.process(std::move(in_chunk));
  auto out_chunk = dynamic_cast_unique<fastq_chunk>(out.at(0).second);
  REQUIRE(out_chunk.get());
  REQUIRE(out_chunk.get() == in_chunk_ptr);

  REQUIRE_NOTHROW(step.finalize());
}

////////////////////////////////////////////////////////////////////////////////
// adapter_selector

TEST_CASE("adapter selector requires adapter or fastq block")
{
  adapter_database database;
  database.add_known();

  adapter_selector step{ database,
                         PARAMETERIZE_IS,
                         DEFAULT_MISMATCH_THRESHOLD,
                         DEFAULT_NEXT_STEP };

  REQUIRE_THROWS_AS(step.process(chunk_ptr{}), assert_failed);
  REQUIRE_THROWS_AS(step.process(std::make_unique<output_chunk>()),
                    assert_failed);
}

TEST_CASE("adapter selector returns fastq chunks as is")
{
  adapter_database database;
  database.add_known();
  adapter_selector step{ database,
                         PARAMETERIZE_IS,
                         DEFAULT_MISMATCH_THRESHOLD,
                         DEFAULT_NEXT_STEP };

  auto in_chunk = std::make_unique<fastq_chunk>();
  const auto* in_chunk_ptr = in_chunk.get();
  in_chunk->reads_1.resize(10, fastq("read", "ACGTTGCTGATCAACTGGTACATATCAAG"));

  auto out = step.process(std::move(in_chunk));
  auto out_chunk = dynamic_cast_unique<fastq_chunk>(out.at(0).second);
  REQUIRE(out_chunk.get() == in_chunk_ptr);
}

TEST_CASE("adapter selector on SE data")
{
  adapter_set adapters{
    { "GTTATTTA", "ACGTGTTA" },
    { "ACGGACGT", "GGCAGTTA" },
  };
  adapter_database database;
  database.add(adapters);

  adapter_selector step{ database,
                         PARAMETERIZE_IS,
                         DEFAULT_MISMATCH_THRESHOLD,
                         DEFAULT_NEXT_STEP };

  auto chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();
  chunk->data->reads_1 = {
    fastq{ "read1", "AGTGTTATTTAA" },
    fastq{ "read1", "ACGGACGTTTA" },
    fastq{ "read1", "ACGGACGTTTA" },
  };

  auto out = step.process(std::move(chunk));
  REQUIRE(out.size() == 1);
  auto out_chunk = dynamic_cast_unique<adapter_chunk>(out.at(0).second);
  REQUIRE(out_chunk.get());

  REQUIRE(out_chunk->adapters.mate_1() ==
          hits_vec{ { 2, 16 }, {}, {}, { 1, 8 } });
  REQUIRE(out_chunk->adapters.mate_2() == hits_vec{});
}

TEST_CASE("adapter selector on PE data")
{
  adapter_set adapters{
    { "GTTATTTA", "ACGTGTTA" },
    { "ACGGACGT", "GGCAGTTA" },
  };
  adapter_database database;
  database.add(adapters);

  adapter_selector step{ database,
                         PARAMETERIZE_IS,
                         DEFAULT_MISMATCH_THRESHOLD,
                         DEFAULT_NEXT_STEP };

  auto chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();
  chunk->data->reads_1 = {
    fastq{ "read1", "AGTGTTATTTAA" },
    fastq{ "read1", "ACGGACGTTTA" },
    fastq{ "read1", "ACGGACGTTTA" },
  };
  chunk->data->reads_2 = {
    fastq{ "read2", "ACTACGTGTAA" },
    fastq{ "read2", "ACTACGTGTTT" },
    fastq{ "read2", "ACTACGTGTTT" },
  };

  auto out = step.process(std::move(chunk));
  REQUIRE(out.size() == 1);
  auto out_chunk = dynamic_cast_unique<adapter_chunk>(out.at(0).second);
  REQUIRE(out_chunk.get());

  REQUIRE(out_chunk->adapters.mate_1() ==
          hits_vec{ { 2, 16 }, {}, {}, { 1, 8 } });
  REQUIRE(out_chunk->adapters.mate_2() == hits_vec{ {}, { 3, 24, 3 }, {}, {} });
}

////////////////////////////////////////////////////////////////////////////////
// adapter_finalizer

namespace {

chunk_vec
test_finalizer_on_chunks(threadsafe_data<sample_set> samples,
                         std::vector<chunk_ptr>& chunks)
{

  REQUIRE(samples.get_reader()->adapters() == adapter_set{});

  adapter_set adapters{
    { "GTTATTTA", "ACGTGTTA" },
    { "ACGGACGT", "GGCAGTTA" },
  };
  adapter_database database;
  database.add(adapters);

  adapter_finalizer step(database, samples, DEFAULT_NEXT_STEP);
  REQUIRE(samples.get_reader()->adapters() == adapter_set{});

  chunk_vec results;
  for (auto&& chunk : chunks) {
    for (auto&& it : step.process(std::move(chunk))) {
      results.emplace_back(std::move(it));
    }
  }

  return results;
}

/** Helper function to allow safe checks against chunk identity */
template<typename T>
uintptr_t
ptr_to_id(const std::unique_ptr<T>& data)
{
  return reinterpret_cast<uintptr_t>(data.get());
}

} // namespace

TEST_CASE("post adapter-selection chunks are passed as is")
{
  adapter_database database;
  threadsafe_data<sample_set> samples;
  adapter_finalizer step(database, samples, DEFAULT_NEXT_STEP);
  const std::vector<fastq> reads{
    fastq{ "read1", "AGTGTTATTTAA" },
    fastq{ "read1", "ACGGACGTTTA" },
    fastq{ "read1", "ACGGACGTTTA" },
  };

  auto data = std::make_unique<fastq_chunk>();
  data->reads_1 = reads;

  const auto before = reinterpret_cast<uintptr_t>(data.get());
  auto chunks = step.process(std::move(data));

  REQUIRE(chunks.size() == 1);
  data = dynamic_cast_unique<fastq_chunk>(chunks.at(0).second);
  const auto after = reinterpret_cast<uintptr_t>(data.get());

  REQUIRE(before == after);
  REQUIRE(data->reads_1 == reads);
}

TEST_CASE("no notifications before final chunk")
{
  log::log_capture cap;
  std::vector<chunk_ptr> chunks;
  chunk_vec results;

  threadsafe_data<sample_set> samples;
  auto chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();
  chunk->adapters = adapter_detection_stats{ 10,
                                             { {}, {}, {}, { 7, 0 } },
                                             { {}, { 3, 0 }, {}, {} } };
  chunks.emplace_back(std::move(chunk));
  results = test_finalizer_on_chunks(samples, chunks);

  REQUIRE_THAT(cap.str(), !Contains("adapter"));
}

TEST_CASE("user is notified when adapters are selected for SE data")
{
  log::log_capture cap;
  std::vector<chunk_ptr> chunks;
  chunk_vec results;

  threadsafe_data<sample_set> samples;
  auto chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();
  chunk->last_adapter_selection_block = true;
  const auto ptr_before = ptr_to_id(chunk->data);

  SECTION("not enough matches")
  {
    chunk->adapters = adapter_detection_stats{ 10, { {}, {}, {}, { 7, 0 } } };
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);
    REQUIRE_THAT(
      cap.str(),
      Contains("Could not identify adapter sequences automatically"));
    REQUIRE(samples.get_reader()->adapters() == adapter_set{ { "", "" } });
  }

  SECTION("sufficient hits for --adapter1")
  {
    chunk->adapters = adapter_detection_stats{ 10, { {}, {}, {}, { 10, 0 } } };
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);

    REQUIRE_THAT(cap.str(),
                 Contains("Using 'User adapters #1' read 1 adapter for "
                          "--adapter1: GTTATTTA"));
    REQUIRE_THAT(cap.str(), !Contains("--adapter2"));
    REQUIRE(samples.get_reader()->adapters() ==
            adapter_set{ { "GTTATTTA", "" } });
  }

  REQUIRE(results.size() == 1);
  REQUIRE(results.at(0).first == DEFAULT_NEXT_STEP);

  const auto ptr_after =
    reinterpret_cast<uintptr_t>(results.at(0).second.get());
  REQUIRE(ptr_before == ptr_after);
}

TEST_CASE("user is notified when adapters are selected for PE data")
{
  log::log_capture cap;
  std::vector<chunk_ptr> chunks;
  chunk_vec results;

  threadsafe_data<sample_set> samples;
  auto chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();
  chunk->last_adapter_selection_block = true;
  const auto ptr_before = ptr_to_id(chunk->data);

  SECTION("no reads, no messages")
  {
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);

    REQUIRE_THAT(cap.str(), !Contains("adapter"));
    REQUIRE(samples.get_reader()->adapters() == adapter_set{ { "", "" } });
  }

  SECTION("not enough matches")
  {
    chunk->adapters = adapter_detection_stats{ 10,
                                               { {}, {}, {}, { 7, 0 } },
                                               { {}, { 3, 0 }, {}, {} } };
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);
    REQUIRE_THAT(
      cap.str(),
      Contains("Could not identify adapter sequences automatically"));
    REQUIRE(samples.get_reader()->adapters() == adapter_set{ { "", "" } });
  }

  SECTION("sufficient hits for --adapter1")
  {
    chunk->adapters = adapter_detection_stats{ 10,
                                               { {}, {}, {}, { 10, 0 } },
                                               { {}, { 9, 0 }, {}, {} } };
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);

    REQUIRE_THAT(cap.str(),
                 Contains("Using 'User adapters #1' read 1 adapter for "
                          "--adapter1: GTTATTTA"));
    REQUIRE_THAT(cap.str(), Contains("Could not identify --adapter2 sequence"));
    REQUIRE(samples.get_reader()->adapters() ==
            adapter_set{ { "GTTATTTA", "" } });
  }

  SECTION("sufficient hits for --adapter2")
  {
    chunk->adapters = adapter_detection_stats{ 10,
                                               { {}, {}, {}, { 9, 0 } },
                                               { {}, { 10, 0 }, {}, {} } };
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);

    REQUIRE_THAT(cap.str(), Contains("Could not identify --adapter1 sequence"));
    REQUIRE_THAT(cap.str(),
                 Contains("Using 'User adapters #1' read 2 adapter for "
                          "--adapter2: ACGTGTTA"));
    REQUIRE(samples.get_reader()->adapters() ==
            adapter_set{ { "", "ACGTGTTA" } });
  }

  SECTION("sufficient hits for both")
  {
    chunk->adapters = adapter_detection_stats{ 10,
                                               { {}, {}, {}, { 10, 0 } },
                                               { {}, { 10, 0 }, {}, {} } };
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);

    REQUIRE_THAT(cap.str(),
                 Contains("Using 'User adapters #1' read 1 adapter for "
                          "--adapter1: GTTATTTA"));
    REQUIRE_THAT(cap.str(),
                 Contains("Using 'User adapters #1' read 2 adapter for "
                          "--adapter2: ACGTGTTA"));
    REQUIRE(samples.get_reader()->adapters() ==
            adapter_set{ { "GTTATTTA", "ACGTGTTA" } });
  }

  SECTION("sufficient hits for both, in wrong orientation")
  {
    chunk->adapters = adapter_detection_stats{ 10,
                                               { {}, { 10, 0 }, {}, {} },
                                               { {}, {}, {}, { 10, 0 } } };
    chunks.emplace_back(std::move(chunk));
    results = test_finalizer_on_chunks(samples, chunks);

    REQUIRE_THAT(cap.str(),
                 Contains("Using 'User adapters #1' read 2 adapter for "
                          "--adapter1: ACGTGTTA"));
    REQUIRE_THAT(cap.str(),
                 Contains("Using 'User adapters #1' read 1 adapter for "
                          "--adapter2: GTTATTTA"));
    REQUIRE_THAT(cap.str(),
                 Contains("The orientation of identified adapters does not "
                          "match expectations"));
    REQUIRE(samples.get_reader()->adapters() ==
            adapter_set{ { "ACGTGTTA", "GTTATTTA" } });
  }

  REQUIRE(results.size() == 1);
  REQUIRE(results.at(0).first == DEFAULT_NEXT_STEP);
  REQUIRE(ptr_before == ptr_to_id(results.at(0).second));
}

TEST_CASE("stats are merged across chunks")
{
  log::log_capture _cap;
  std::vector<chunk_ptr> chunks;

  auto chunk_1 = std::make_unique<adapter_chunk>();
  chunk_1->data = std::make_unique<fastq_chunk>();
  chunk_1->last_adapter_selection_block = false;
  chunk_1->adapters = adapter_detection_stats{ 10,
                                               { {}, {}, {}, { 7, 0 } },
                                               { {}, { 3, 0 }, {}, {} } };

  auto chunk_2 = std::make_unique<adapter_chunk>();
  chunk_2 = std::make_unique<adapter_chunk>();
  chunk_2->data = std::make_unique<fastq_chunk>();
  chunk_2->last_adapter_selection_block = true;
  chunk_2->adapters = adapter_detection_stats{ 10,
                                               { {}, {}, {}, { 4, 0 } },
                                               { {}, { 8, 0 }, {}, {} } };

  SECTION("before final chunk")
  {
    threadsafe_data<sample_set> samples;
    chunks.emplace_back(std::move(chunk_1));
    test_finalizer_on_chunks(samples, chunks);

    REQUIRE(samples.get_reader()->adapters() == adapter_set{});
  }

  SECTION("sufficient hits")
  {
    threadsafe_data<sample_set> samples;
    chunks.emplace_back(std::move(chunk_1));
    chunks.emplace_back(std::move(chunk_2));
    test_finalizer_on_chunks(samples, chunks);
    REQUIRE(samples.get_reader()->adapters() ==
            adapter_set{ { "GTTATTTA", "ACGTGTTA" } });
  }
}

TEST_CASE("final chunk must be final adapter_chunk")
{
  adapter_database database;
  database.add_known();

  threadsafe_data<sample_set> samples;
  adapter_finalizer step(database, samples, DEFAULT_NEXT_STEP);

  auto chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();
  chunk->last_adapter_selection_block = true;

  step.process(std::move(chunk));

  chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();
  chunk->last_adapter_selection_block = GENERATE(true, false);

  REQUIRE_THROWS_AS(step.process(std::move(chunk)), assert_failed);
}

TEST_CASE("finalize checks that last chunk was observed")
{
  adapter_database database;
  database.add_known();

  threadsafe_data<sample_set> samples;
  adapter_finalizer step(database, samples, DEFAULT_NEXT_STEP);

  SECTION("before final chunk #1")
  {
    REQUIRE_THROWS_AS(step.finalize(), assert_failed);
  }

  auto chunk = std::make_unique<adapter_chunk>();
  chunk->data = std::make_unique<fastq_chunk>();

  SECTION("before final chunk #2")
  {
    step.process(std::move(chunk));
    REQUIRE_THROWS_AS(step.finalize(), assert_failed);
  }

  SECTION("after final chunk")
  {
    chunk->last_adapter_selection_block = true;
    step.process(std::move(chunk));
    REQUIRE_NOTHROW(step.finalize());
  }
}

} // namespace adapterremoval