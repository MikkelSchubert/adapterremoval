/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Stinus Lindgreen - stinus@binf.ku.dk            *
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
#pragma once

#include <fstream>
#include <string>
#include <vector>

class HTMLTmplHead
{
public:
  HTMLTmplHead();
  ~HTMLTmplHead();

  void set_name(const std::string& value);
  void set_version(const std::string& value);

  void write(std::ofstream& out);

private:
  bool m_written;
  std::string m_name;
  bool m_name_is_set;
  std::string m_version;
  bool m_version_is_set;
};

class HTMLTmplBody
{
public:
  HTMLTmplBody();
  ~HTMLTmplBody();

  void write(std::ofstream& out);

private:
  bool m_written;
};
