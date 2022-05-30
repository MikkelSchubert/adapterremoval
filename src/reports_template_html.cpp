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
#include "reports_template_html.hpp"
#include "debug.hpp"

size_t g_html_id = 1;

html_head::html_head()
  : m_written()
  , m_name()
  , m_name_is_set()
  , m_version()
  , m_version_is_set()
{
  //
}

html_head::~html_head()
{
  AR_REQUIRE(m_written, "template html_head was not written");
}

html_head&
html_head::set_name(const std::string& value)
{
  m_name = value;
  m_name_is_set = true;
  return *this;
}

html_head&
html_head::set_version(const std::string& value)
{
  m_version = value;
  m_version_is_set = true;
  return *this;
}

void
html_head::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_head already written");
  AR_REQUIRE(m_name_is_set, "html_head::name not set");
  AR_REQUIRE(m_version_is_set, "html_head::version not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "<!DOCTYPE html>\n";
  out << "<html lang='en'>\n";
  out << "\n";
  out << "<head>\n";
  out << "    <meta charset='utf-8'>\n";
  out << "    <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n";
  out << "    <title>" << m_name << " " << m_version << "</title>\n";
  out << "    <link rel='stylesheet' href='https://unpkg.com/purecss@2.1.0/build/pure-min.css'\n";
  out << "        integrity='sha384-yHIFVG6ClnONEA5yB5DJXfW2/KC173DIQrYoZMEtBvGzmf0PKiGyNEqe9N6BNDBH' crossorigin='anonymous'>\n";
  out << "    <style type='text/css'>\n";
  out << "\n";
  out << "    </style>\n";
  out << "</head>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_body::html_body()
  : m_written()
{
  //
}

html_body::~html_body()
{
  AR_REQUIRE(m_written, "template html_body was not written");
}

void
html_body::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_body already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "<body>\n";
  out << "    <div id='layout'>\n";
  out << "\n";
  out << "    </div>\n";
  out << "</body>\n";
  out << "\n";
  out << "</html>";
  // clang-format on
  m_written = true;
}
