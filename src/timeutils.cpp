// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "timeutils.hpp" // declarations
#include "debug.hpp"     // for AR_REQUIRE
#include <cctype>        // for isprint, isalnum, tolower, toupper
#include <chrono>        // for system_clock
#include <cmath>         // for log10, pow, round
#include <ctime>         // for localtime_r, tm
#include <iomanip>       // for operator<<, setprecision
#include <sstream>       // for ostringstream, operator<<, basic_ostream, bas...

namespace adapterremoval {

static const auto s_start_timestamp = std::chrono::system_clock::now();

const std::chrono::system_clock::time_point&
start_time()
{
  return s_start_timestamp;
}

std::string
format_time(const std::chrono::system_clock::time_point& now,
            const char* format,
            const bool milliseconds)
{
  AR_REQUIRE(format);
  using namespace std::chrono;

  const auto in_time_t = system_clock::to_time_t(now);

  tm in_localtime{};
  std::ostringstream ss;

#if defined(_WIN32)
  if (!localtime_s(&in_localtime, &in_time_t)) {
#else
  if (localtime_r(&in_time_t, &in_localtime)) {
#endif
    ss << std::put_time(&in_localtime, format);

    if (milliseconds) {
      const auto ms =
        duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
      ss << '.' << std::setfill('0') << std::setw(3) << (ms.count() % 1000);
    }

    return ss.str();
  } else {
    return "<invalid time>";
  }
}

} // namespace adapterremoval
