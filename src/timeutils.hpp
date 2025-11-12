// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <chrono> // for system_clock
#include <string> // for string

namespace adapterremoval {

/** Returns the timepoint for start of execution  */
const std::chrono::system_clock::time_point&
start_time();

/**
 * Returns a timestamp in the specified format, see
 * https://en.cppreference.com/w/cpp/io/manip/put_time
 */
std::string
format_time(const std::chrono::system_clock::time_point& now,
            const char* format,
            bool milliseconds = false);

} // namespace adapterremoval
