// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "progress.hpp"
#include "debug.hpp"    // for AR_FAIL
#include "logging.hpp"  // for info, log_stream, cerr
#include "strutils.hpp" // for format_rough_number, format_thousand_sep
#include <chrono>       // for microseconds
#include <iomanip>      // for operator<<, setfill, setw
#include <sstream>      // for operator<<, basic_ostream, ostringstream
#include <string>       // for string, operator<<, basic_string, char_traits
#include <vector>       // for vector

namespace adapterremoval {

namespace {
//! Print progress report every N read
const size_t REPORT_EVERY_NTH_READ = 1e6;
//! Print an updated progress report every S seconds
const size_t REPORT_EVERY_NTH_LOOP = 10;
//! Increment the spinner every time-unit
const auto SPIN_EVERY = std::chrono::microseconds(100000);

std::string
format_time(double seconds)
{
  std::ostringstream stream;
  stream.precision(1);
  stream << std::setfill('0');

  if (seconds > 60 * 60) {
    stream << static_cast<size_t>(seconds) / (60LU * 60) << ":";
    stream << std::setw(2);
  }

  if (seconds > 60) {
    stream << (static_cast<size_t>(seconds) % (60LU * 60)) / 60 << ":";
    stream << std::setw(4) << std::setfill('0');
  }

  stream << std::fixed << (static_cast<size_t>(seconds * 100) % 6000) / 100.0
         << "s";

  return stream.str();
}

std::string
format_progress(size_t reads,
                double seconds,
                double rate,
                bool finalize = false)
{
  std::ostringstream ss;

  if (finalize) {
    ss << "Processed " << format_thousand_sep(reads) << " reads in "
       << format_time(seconds) << "; " << format_rough_number(rate)
       << " reads/s on average";
  } else {
    ss << "Processed " << format_rough_number(reads) << " reads in "
       << format_time(seconds) << "; " << format_rough_number(rate)
       << " reads/s";
  }

  return ss.str();
}

} // namespace

progress_timer::progress_timer(progress_type type)
  : m_type(type)
{
  start();
}

progress_timer::~progress_timer()
{
  stop();
}

void
progress_timer::increment(size_t inc)
{
  switch (m_type) {
    case progress_type::none: {
      m_total += inc;
      break;
    }

    case progress_type::simple: {
      m_total += inc;
      m_current += inc;

      if (m_current >= REPORT_EVERY_NTH_READ) {
        const double current_time = m_timer.duration();
        const double rate = m_current / (current_time - m_last_time);
        log::info() << format_progress(m_total, m_timer.duration(), rate);

        m_current = 0;
        m_last_time = current_time;
      }

      break;
    }

    case progress_type::spinner: {
      std::unique_lock<std::mutex> lock(m_lock);

      m_total += inc;
      break;
    }

    default:
      AR_FAIL("invalid progress report type");
  }
}

void
progress_timer::finalize()
{
  stop();

  const double rate = m_total / m_timer.duration();
  log::info() << format_progress(m_total, m_timer.duration(), rate, true);
}

void
progress_timer::start()
{
  if (m_type == progress_type::spinner && !m_spinner.joinable()) {
    m_spinning = true;
    m_spinner = std::thread(&progress_timer::loop, this);
  }
}

void
progress_timer::loop()
{
  double last_time = 0;
  size_t last_reads = 0;

  const std::vector<std::string> symbols = { "⠙", "⠸", "⢰", "⣠",
                                             "⣄", "⡆", "⠇", "⠋" };
  auto symbol_it = symbols.begin();
  std::string last_message;

  for (size_t loop = 0;; loop++) {
    std::this_thread::sleep_for(SPIN_EVERY);

    double current_time = 0;
    size_t current_reads = 0;

    {
      std::unique_lock<std::mutex> lock(m_lock);
      if (m_spinning) {
        current_time = m_timer.duration();
        current_reads = m_total;
      } else {
        break;
      }
    }

    if (loop % REPORT_EVERY_NTH_LOOP == 0) {
      const double rate =
        (current_reads - last_reads) / (current_time - last_time);

      last_message = format_progress(current_reads, current_time, rate);
      last_reads = current_reads;
      last_time = current_time;
    }

    log::cerr().transient() << *symbol_it++ << " " << last_message;

    if (symbol_it == symbols.end()) {
      symbol_it = symbols.begin();
    }
  }
}

void
progress_timer::stop()
{
  if (m_spinner.joinable()) {
    {
      std::unique_lock<std::mutex> lock(m_lock);
      m_spinning = false;
    }

    m_spinner.join();
  }
}

} // namespace adapterremoval
