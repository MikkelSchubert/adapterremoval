// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <string_view> // for string_view

namespace adapterremoval {

class statistics;
class userconfig;

bool
write_json_report(const userconfig& config,
                  const statistics& stats,
                  std::string_view filename);

bool
write_html_report(const userconfig& config,
                  const statistics& stats,
                  std::string_view filename);

void
print_terminal_preamble(const userconfig& config);

void
print_terminal_postamble(const userconfig& config, bool any_errors = false);

} // namespace adapterremoval
