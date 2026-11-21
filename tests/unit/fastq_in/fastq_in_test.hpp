// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "fastq_in.hpp"       // for reader_state
#include "scheduler.hpp"      // for step_id
#include "scheduler_test.hpp" // for make_step_id_for_testing
#include <cstddef>            // for size_t
#include <memory>             // for shared_ptr

namespace adapterremoval {

constexpr size_t DEFAULT_BUFFER_CAPACITY = 0x10000;
constexpr size_t DEFAULT_FASTQ_CAPACITY = 0x80;
const step_id DEFAULT_NEXT_STEP = make_step_id_for_testing(12345);

/** Returns reader state with the specified default capacities */
std::shared_ptr<reader_state>
make_state(size_t buffer_capacity = DEFAULT_BUFFER_CAPACITY,
           size_t fastq_capacity = DEFAULT_FASTQ_CAPACITY);

/** Return reader state that shares reusable objects with `source` */
std::shared_ptr<reader_state>
make_state(const std::shared_ptr<reader_state>& source);

/** Indicates if a chunk is is at the end of file; enum for readability */
enum class chunk_eof
{
  no,
  yes
};

/** Returns a reader-chunk with an acquired, 0-sized buffer */
[[nodiscard]] reader_chunk_ptr
make_reader_chunk(const reader_state_ptr& state,
                  chunk_eof eof,
                  reader_chunk::type type = reader_chunk::type::unknown,
                  std::shared_ptr<const std::string> filename = {});

} // namespace adapterremoval
