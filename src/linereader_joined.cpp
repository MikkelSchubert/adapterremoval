/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>

#include "debug.hpp"
#include "linereader_joined.hpp"
#include "threads.hpp"

joined_line_readers::joined_line_readers(const string_vec& filenames)
  : m_filenames(filenames.rbegin(), filenames.rend())
  , m_reader()
  , m_current_line(1)
{}

joined_line_readers::~joined_line_readers() {}

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
    return false;
  }

  auto filename = m_filenames.back();

  {
    print_locker lock;
    std::cerr << "Opening FASTQ file '" << filename
              << "', line numbers start at " << m_current_line << std::endl;
  }

  m_reader.reset(new line_reader(filename));
  m_filenames.pop_back();

  return true;
}
