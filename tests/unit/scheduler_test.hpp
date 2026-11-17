// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp" // for read_mate
#include "scheduler.hpp"   // for analytical_step, chunk_ptr, chunk_vec
#include "testing.hpp"     // for REQUIRE
#include "utilities.hpp"   // for dynamic_cast_unique
#include <cstddef>         // for size_t
#include <functional>      // for function
#include <memory>          // for unique_ptr
#include <string_view>     // for string_view
#include <utility>         // for move

namespace adapterremoval {

// GCOVR_EXCL_BR_START

/** Helper that checks and returns a single expected value from a chunk_vec */
template<typename T>
T
chunk_vec_to(chunk_vec&& chunks,
             std::function<void(bool, std::string_view)> check)
{
  check(!chunks.empty(), "no chunks returned");
  check(chunks.size() == 1, "more than 1 chunk returned");
  auto* variant = std::get_if<T>(&chunks.at(0));
  check(variant != nullptr, "chunk is wrong variant");
  return std::move(*variant);
}

/** Helper that checks/casts/returns a single expected chunk from a chunk_vec */
template<typename T>
std::unique_ptr<T>
chunk_vec_to_chunk(chunk_vec&& chunks,
                   std::function<void(bool, std::string_view)> check)
{
  auto base_chunk = chunk_vec_to<chunk_ptr>(std::move(chunks), check);
  check(base_chunk.get(), "chunk is NULL");
  auto chunk = dynamic_cast_unique<T>(base_chunk);
  check(chunk.get(), "chunk could not be cast to target");

  return chunk;
}

// GCOVR_EXCL_BR_STOP

/** Helper macro that ensures that the actual source location is reported */
#define STEP_CHECK_WITH_LOCATION                                               \
  [](bool check_, std::string_view message_) {                                 \
    if (!check_) {                                                             \
      FAIL(message_);                                                          \
    }                                                                          \
  }

/** Process a chunk, for which no output is expected */
#define STEP_PROCESS_TO_NOTHING(step, chunk)                                   \
  REQUIRE((step).process(chunk).empty())

/** Process a single chunk and return the expected single chunk */
#define STEP_PROCESS_TO_CHUNK(T, step, chunk)                                  \
  chunk_vec_to_chunk<T>((step).process(chunk), STEP_CHECK_WITH_LOCATION)

/** Process a single chunk and return a value of the expected type */
#define STEP_PROCESS_TO(T, step, chunk)                                        \
  chunk_vec_to<T>((step).process(chunk), STEP_CHECK_WITH_LOCATION)

/** Send a signal for which no output is expected */
#define STEP_SIGNAL_TO_NOTHING(step) REQUIRE((step).signal().empty())

/** Signal a step and return the expected single chunk */
#define STEP_SIGNAL_TO_CHUNK(T, step)                                          \
  chunk_vec_to_chunk<T>((step).signal(), STEP_CHECK_WITH_LOCATION)

/** Signal a step and return the expected single signal */
#define STEP_SIGNAL_TO(T, step)                                                \
  chunk_vec_to<T>((step).signal(), STEP_CHECK_WITH_LOCATION)

/** Creates a step_id for testing */
step_id
make_step_id_for_testing(size_t id);

/** Sets the signal target of an analytical step for testing */
void
set_signal_target_for_testing(analytical_step& step,
                              step_id target,
                              read_mate mate);

} // namespace adapterremoval
