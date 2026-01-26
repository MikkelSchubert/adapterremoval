// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <string>      // for string
#include <string_view> // for string_view

namespace adapterremoval::program {

/** Returns program name */
std::string_view
name();

/** Returns version string without revision information */
std::string_view
short_version();

/** Returns complete version string, possibly including git tags */
std::string_view
long_version();

/** Returns program name plus short version string */
std::string_view
short_name();

/** Returns program name plus long version string */
std::string_view
long_name();

} // namespace adapterremoval::program
