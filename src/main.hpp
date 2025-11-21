// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#ifdef __FAST_MATH__
#error "AdapterRemoval cannot be compiled with -ffast-math"
#endif

#include <string_view> // for string_view

namespace adapterremoval {

class userconfig;

// See main_adapter_rm.cpp
int
remove_adapter_sequences(const userconfig& config);
// See main_adapter_id.cpp
int
generate_reports(const userconfig& config);
// See main_benchmark.cpp
int
benchmark(const userconfig& config);

/** Terminate on internal error; prints  */
[[noreturn]] void
terminate(std::string_view message);

} // namespace adapterremoval
