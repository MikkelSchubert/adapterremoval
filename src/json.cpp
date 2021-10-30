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
#include <cmath>
#include <iomanip>
#include <sstream>

#include "json.hpp"

std::string
_escape(const std::string& value)
{
  std::stringstream stream;

  stream << '"';

  for (char c : value) {
    switch (c) {
      case '"':
        stream << "\\\"";
        break;
      case '\\':
        stream << "\\\\";
        break;
      case '\b':
        stream << "\\b";
        break;
      case '\f':
        stream << "\\f";
        break;
      case '\n':
        stream << "\\n";
        break;
      case '\r':
        stream << "\\r";
        break;
      case '\t':
        stream << "\\t";
        break;
      default: {
        if (c <= 0x1F) {
          stream << "\\u" << std::hex << static_cast<int>(c);
        } else {
          stream << c;
        }
      }
    }
  }

  stream << '"';

  return stream.str();
}

json_section::json_section(json_section&& other) noexcept
  : m_parent(nullptr)
{
  std::swap(m_parent, other.m_parent);
}

json_section::json_section(json_writer* parent)
  : m_parent(parent)
{
  m_parent->_write("{");
  m_parent->m_values = false;
  m_parent->m_indent++;
}

json_section::json_section(json_writer* parent, const std::string& key)
  : m_parent(parent)
{
  m_parent->_write(key, "{");
  m_parent->m_values = false;
  m_parent->m_indent++;
}

json_section::~json_section()
{
  if (m_parent) {
    m_parent->end('}');
  }
}

json_writer::json_writer(std::ostream& stream)
  : m_stream(stream)
  , m_indent(1)
  , m_values(false)
{
  m_stream << "{";
}

json_writer::~json_writer()
{
  while (m_indent) {
    end('}');
  }
}

json_section
json_writer::start()
{
  return json_section(this);
}

json_section
json_writer::start(const std::string& key)
{
  return json_section(this, key);
}

void
json_writer::end(char c)
{
  AR_DEBUG_ASSERT(m_indent);
  if (m_values) {
    m_stream << "\n" << std::setw(2 * m_indent - 1) << c;
  } else {
    m_stream << c;
  }

  m_values = true;
  m_indent--;
}

void
json_writer::start_list(const std::string& key)
{
  _write(key, "[");
  m_values = false;
  m_indent++;
}

void
json_writer::end_list()
{
  end(']');
}

void
json_writer::write(const std::string& key, const std::string& value)
{
  _write(key, _escape(value));
}

void
json_writer::write(const std::string& key,
                   const std::vector<std::string>& value)
{
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < value.size(); ++i) {
    if (i) {
      ss << ", ";
    }

    ss << _escape(value.at(i));
  }

  ss << "]";

  _write(key, ss.str());
}

void
json_writer::write_bool(const std::string& key, const bool value)
{
  _write(key, value ? "true" : "false");
}

void
json_writer::write_int(const std::string& key, const int64_t value)
{
  _write(key, std::to_string(value));
}

void
json_writer::write_float(const std::string& key, const double value)
{
  if (std::isnan(value)) {
    _write(key, "null");
  } else {
    std::ostringstream out;
    out.precision(3);
    out << std::fixed << value;

    _write(key, out.str());
  }
}

void
json_writer::write_null(const std::string& key)
{
  _write(key, "null");
}

void
json_writer::_write(const std::string& value)
{
  AR_DEBUG_ASSERT(m_indent);
  if (m_values) {
    m_stream << ",";
  }

  m_stream << "\n" << std::setw(2 * m_indent + value.size()) << value;
  m_values = true;
}

void
json_writer::_write(const std::string& key, const std::string& value)
{
  AR_DEBUG_ASSERT(m_indent);
  _write(_escape(key));

  m_stream << ": " << value;
}
