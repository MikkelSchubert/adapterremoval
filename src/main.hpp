// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <string_view> // for string_view

#ifdef __FAST_MATH__
#error "AdapterRemoval cannot be compiled with -ffast-math"
#endif

namespace adapterremoval {

const std::string_view NAME = "AdapterRemoval";
const std::string_view VERSION = "3.0.0-alpha3";
const std::string_view FULL_NAME = "AdapterRemoval v3.0.0-alpha3";

} // namespace adapterremoval
