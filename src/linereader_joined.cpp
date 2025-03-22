// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "linereader_joined.hpp"
#include <memory> // for make_unique
#include <vector> // for vector

namespace adapterremoval {

joined_line_readers::joined_line_readers(const string_vec& filenames)
  : m_filenames(filenames.rbegin(), filenames.rend())
{
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

bool
joined_line_readers::open_next_file()
{
  if (m_filenames.empty()) {
    m_reader.reset();
    return false;
  }

  m_reader = std::make_unique<line_reader>(m_filenames.back());
  m_filename = m_filenames.back();
  m_current_line = 1;

  m_filenames.pop_back();

  return true;
}

} // namespace adapterremoval
