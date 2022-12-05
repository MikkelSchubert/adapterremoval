/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
\*************************************************************************/
#include <memory> // for make_unique
#include <vector> // for vector

#include "linereader_joined.hpp"

namespace adapterremoval {

joined_line_readers::joined_line_readers(const string_vec& filenames)
  : m_filenames(filenames.rbegin(), filenames.rend())
  , m_reader()
  , m_filename()
  , m_current_line(0)
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

const std::string&
joined_line_readers::filename() const
{
  return m_filename;
}

size_t
joined_line_readers::linenumber() const
{
  return m_current_line;
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
