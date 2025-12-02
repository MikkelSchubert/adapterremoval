// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "linereader_joined.hpp" // declarations
#include "debug.hpp"             // for AR_REQUIRE
#include "strutils.hpp"          // for log_escape
#include <limits>                // for numeric_limits
#include <memory>                // for make_unique
#include <sstream>               // for ostringstream
#include <vector>                // for vector

namespace adapterremoval {

joined_line_readers::joined_line_readers(std::vector<std::string> filenames)
{
  AR_REQUIRE(!filenames.empty());
  for (auto&& filename : filenames) {
    m_filenames.emplace_back(std::move(filename),
                             std::numeric_limits<size_t>::max());
  }
}

bool
joined_line_readers::getline(std::string& dst)
{
  while (true) {
    if (m_reader && m_reader->getline(dst)) {
      m_current_line++;
      return true;
    } else if (!open_next_file()) {
      return false;
    }
  }
}

std::string
joined_line_readers::filenames(size_t start, size_t end) const
{
  AR_REQUIRE(1 <= start && start <= end && end <= m_current_line);

  string_vec parts;
  while (start <= end) {
    bool found_file = false;
    size_t current_offset = 0;
    for (const auto& it : m_filenames) {
      if (start <= it.second) {
        size_t file_start = start - current_offset;
        size_t file_end = std::min(end, it.second) - current_offset;
        // For parsing errors we expect that start < end, but for I/O errors
        // the line number indicates the line at which the read failed
        if (file_start < file_end) {
          file_end--;
        }

        std::ostringstream os;
        if (file_start != file_end) {
          os << log_escape(it.first) << " at lines " << file_start << "-"
             << file_end;
        } else {
          os << log_escape(it.first) << " at line " << file_start;
        }

        parts.emplace_back(os.str());
        start = std::min(end, it.second) + 1;
        found_file = true;
        break;
      }

      current_offset = it.second;
    }

    AR_REQUIRE(found_file);
  }

  return join_text(parts, ", ", ", and ");
}

bool
joined_line_readers::open_next_file()
{
  if (m_reader) {
    // m_current_line is the the first line in the next file at this point
    m_filenames.at(m_current_file).second = m_current_line - 1;
    m_current_file++;
  }

  if (m_current_file >= m_filenames.size()) {
    m_reader.reset();
    return false;
  }

  m_reader =
    std::make_unique<line_reader>(m_filenames.at(m_current_file).first);

  return true;
}

} // namespace adapterremoval
