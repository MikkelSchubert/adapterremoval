// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "linereader_joined.hpp" // declarations
#include "debug.hpp"             // for AR_REQUIRE
#include "strutils.hpp"          // for log_escape
#include <algorithm>             // for min
#include <cstddef>               // for size_t
#include <limits>                // for numeric_limits
#include <memory>                // for make_unique
#include <sstream>               // for ostringstream
#include <string>                // string
#include <utility>               // for move
#include <vector>                // for vector

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// joined_filenames

joined_filenames::joined_filenames(std::vector<std::string> filenames)
{
  AR_REQUIRE(!filenames.empty());
  for (auto&& filename : filenames) {
    m_filenames.emplace_back(std::move(filename),
                             std::numeric_limits<size_t>::max());
  }
}

void
joined_filenames::inc_position(size_t n)
{
  AR_REQUIRE(m_current_file < m_filenames.size());
  m_position += n;
}

std::string
joined_filenames::filenames(size_t start, size_t end) const
{
  AR_REQUIRE(start <= end && end <= m_position);
  if (start == m_position && remaining_filenames() == 0) {
    return "end of input files";
  }

  size_t current_offset = 0;
  std::vector<std::string> parts;
  for (const auto& [filename, current_end] : m_filenames) {
    if (end < current_offset) {
      break;
    } else if (start >= current_offset && start < current_end) {
      const size_t file_start = start - current_offset;
      const size_t file_end = std::min(current_end, end) - current_offset;
      start += file_end - file_start;

      std::ostringstream os;
      // Reading failed at line, or parsing failed after reading 1 line in file
      if (file_start == file_end || file_start + 1 == file_end) {
        os << log_escape(filename) << " at line " << file_start + 1;
      } else {
        os << log_escape(filename) << " at lines " << file_start + 1 << "-"
           << file_end;
      }

      parts.emplace_back(os.str());
    } else if (start == current_offset && current_offset == current_end) {
      parts.emplace_back("empty file " + log_escape(filename));
    }

    current_offset = current_end;
  }

  return join_text(parts, ", ", ", and ");
}

const std::string&
joined_filenames::current_filename() const
{
  AR_REQUIRE(m_current_file < m_filenames.size());
  return m_filenames.at(m_current_file).first;
}

void
joined_filenames::next_file()
{
  AR_REQUIRE(m_current_file < m_filenames.size());
  m_filenames.at(m_current_file).second = m_position;
  m_current_file++;
}

////////////////////////////////////////////////////////////////////////////////
// joined_line_readers

joined_line_readers::joined_line_readers(std::vector<std::string> filenames)
  : m_filenames(std::move(filenames))
{
}

bool
joined_line_readers::getline(std::string& dst)
{
  while (true) {
    if (m_reader && m_reader->getline(dst)) {
      m_filenames.inc_position();
      return true;
    } else if (!open_next_file()) {
      return false;
    }
  }
}

size_t
joined_line_readers::position() const noexcept
{
  return m_filenames.position();
}

std::string
joined_line_readers::filenames(size_t start, size_t end) const
{
  return m_filenames.filenames(start, end);
}

bool
joined_line_readers::open_next_file()
{
  if (m_reader) {
    m_reader.reset();
    m_filenames.next_file();
  }

  if (m_filenames.remaining_filenames()) {
    m_reader = std::make_unique<line_reader>(m_filenames.current_filename());
    return true;
  }

  return false;
}

} // namespace adapterremoval
