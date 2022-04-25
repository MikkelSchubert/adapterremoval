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
#include <algorithm> // for find
#include <cmath>     // for isnan
#include <iomanip>   // for operator<<, setw
#include <sstream>   // for operator<<, basic_ostream, stringstream, ostream
#include <utility>   // for swap

#include "debug.hpp"     // for AR_DEBUG_ASSERT
#include "json.hpp"      // for definitions
#include "utilities.hpp" // for make_shared

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

std::string
format_f64(double value)
{
  if (std::isnan(value)) {
    return "null";
  }

  std::ostringstream out;
  out.precision(3);
  out << std::fixed << value;

  return out.str();
}

////////////////////////////////////////////////////////////////////////////////

json_value::~json_value() {}

std::string
json_value::to_string() const
{
  std::stringstream ss;

  write(ss);

  return ss.str();
}

////////////////////////////////////////////////////////////////////////////////

json_token::json_token(const std::string& value)
  : m_value(value)
{
}

json_ptr
json_token::from_str(const std::string& value)
{
  return make_shared<json_token>(_escape(value));
}

json_ptr
json_token::from_str_vec(const string_vec& values)
{
  string_vec escaped;
  for (const auto& it : values) {
    escaped.push_back(_escape(it));
  }

  return json_token::from_raw_vec(escaped);
}

json_ptr
json_token::from_i64(const int64_t value)
{
  return make_shared<json_token>(std::to_string(value));
}

json_ptr
json_token::from_i64_vec(const counts& values)
{
  string_vec svalues;
  for (size_t i = 0; i < values.size(); ++i) {
    svalues.push_back(std::to_string(values.get(i)));
  }

  return json_token::from_raw_vec(svalues);
}

json_ptr
json_token::from_f64(const double value)
{
  return make_shared<json_token>(format_f64(value));
}

json_ptr
json_token::from_f64_vec(const rates& values)
{
  string_vec svalues;
  for (size_t i = 0; i < values.size(); ++i) {
    svalues.push_back(format_f64(values.get(i)));
  }

  return json_token::from_raw_vec(svalues);
}

json_ptr
json_token::from_boolean(const bool value)
{
  return make_shared<json_token>(value ? "true" : "false");
}

json_ptr
json_token::from_null()
{
  return make_shared<json_token>("null");
}

std::string
json_token::to_string() const
{
  return m_value;
}

void
json_token::write(std::ostream& out, size_t /* indent */) const
{
  out << m_value;
}

json_ptr
json_token::from_raw_vec(const string_vec& values)
{
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < values.size(); ++i) {
    if (i) {
      ss << ", ";
    }

    ss << values.at(i);
  }
  ss << "]";

  return make_shared<json_token>(ss.str());
}

////////////////////////////////////////////////////////////////////////////////

json_list::json_list()
  : m_values()
{
}

void
json_list::write(std::ostream& out, size_t indent_) const
{
  const auto indent = std::string(indent_, ' ');

  if (m_values.empty()) {
    out << "[]";
  } else {
    out << "[\n";

    for (size_t i = 0; i < m_values.size(); ++i) {
      out << indent << "  ";
      m_values.at(i)->write(out, indent_ + 2);

      if (i < m_values.size() - 1) {
        out << ",";
      }

      out << "\n";
    }

    out << indent << "]";
  }
}

json_dict_ptr
json_list::dict()
{
  auto ptr = make_shared<json_dict>();
  m_values.push_back(ptr);

  return ptr;
}

////////////////////////////////////////////////////////////////////////////////
json_dict::json_dict()
  : m_keys()
  , m_values()
{
}

void
json_dict::write(std::ostream& out, size_t indent_) const
{
  AR_DEBUG_ASSERT(m_keys.size() == m_values.size());
  const auto indent = std::string(indent_, ' ');

  if (m_keys.empty()) {
    out << "{}";
  } else {
    out << "{\n";

    for (size_t i = 0; i < m_keys.size(); ++i) {
      const auto it = m_values.find(m_keys.at(i));
      AR_DEBUG_ASSERT(it != m_values.end());

      out << indent << "  " << _escape(it->first) << ": ";
      it->second->write(out, indent_ + 2);
      if (i < m_keys.size() - 1) {
        out << ",";
      }

      out << "\n";
    }

    out << indent << "}";
  }
}

json_dict_ptr
json_dict::dict(const std::string& key)
{
  auto ptr = make_shared<json_dict>();
  _set(key, ptr);

  return ptr;
}

json_list_ptr
json_dict::list(const std::string& key)
{
  auto ptr = make_shared<json_list>();
  _set(key, ptr);

  return ptr;
}

void
json_dict::str(const std::string& key, const std::string& value)
{
  _set(key, json_token::from_str(value));
}

void
json_dict::str_vec(const std::string& key, const string_vec& values)
{
  _set(key, json_token::from_str_vec(values));
}

void
json_dict::i64(const std::string& key, const int64_t value)
{
  _set(key, json_token::from_i64(value));
}

void
json_dict::i64_vec(const std::string& key, const counts& values)
{
  _set(key, json_token::from_i64_vec(values));
}

void
json_dict::f64(const std::string& key, const double value)
{
  _set(key, json_token::from_f64(value));
}

void
json_dict::f64_vec(const std::string& key, const rates& values)
{
  _set(key, json_token::from_f64_vec(values));
}

void
json_dict::boolean(const std::string& key, const bool value)
{
  _set(key, json_token::from_boolean(value));
}

void
json_dict::null(const std::string& key)
{
  _set(key, json_token::from_null());
}

void
json_dict::_set(const std::string& key, const json_ptr& ptr)
{
  const auto it = std::find(m_keys.begin(), m_keys.end(), key);
  if (it == m_keys.end()) {
    m_keys.push_back(key);
  }

  m_values[key] = ptr;
}
