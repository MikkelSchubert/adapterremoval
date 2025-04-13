// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "read_group.hpp" // declarations
#include "debug.hpp"      // for AR_REQUIRE
#include "strutils.hpp"   // for string_vec, indent_lines
#include <stdexcept>      // for invalid_argument
#include <string>         // for string
#include <string_view>    // for string_view

namespace adapterremoval {

read_group::read_group()
  : m_header{ "@RG\tID:1" }
  , m_id{ "1" }
{
}

read_group::read_group(std::string_view value)
  : m_header{ "@RG" }
{
  using invalid = std::invalid_argument;

  // It's not unreasonable to except users to try to specify a full @RG line
  if (starts_with(value, "@RG\t")) {
    value = value.substr(4);
  }

  if (!value.empty()) {
    for (const auto& field : split_text(value, '\t')) {
      for (const auto c : field) {
        if (c < ' ' || c > '~') {
          throw invalid("only characters in the range ' ' to '~' are allowed");
        }
      }

      if (field.size() < 4) {
        throw invalid("tags must be at least 4 characters long");
      } else if (!is_ascii_letter(field.at(0))) {
        throw invalid("first character in tag name must be a letter");
      } else if (!is_ascii_letter_or_digit(field.at(1))) {
        throw invalid(
          "second character in tag name must be a letter or number");
      } else if (field.at(2) != ':') {
        throw invalid("third character in tag must be a colon");
      } else if (starts_with(field, "ID:")) {
        if (!m_id.empty()) {
          throw invalid("multiple ID tags found");
        }

        m_id = field.substr(3);
      }
    }
  }

  // Generic read-group key if the user did not supply one
  if (m_id.empty()) {
    m_id = "1";
    m_header += "\tID:1";
  }

  if (!value.empty()) {
    m_header += "\t";
    m_header += value;
  }
}

void
read_group::set_id(std::string_view id)
{
  AR_REQUIRE(!id.empty());
  update_tag("ID", id);
  m_id = id;
}

bool
read_group::operator==(const read_group& other) const
{
  return m_header == other.m_header && m_id == other.m_id;
}

void
read_group::update_tag(std::string_view key, std::string_view value)
{
  AR_REQUIRE(!key.empty());
  std::string cache;
  cache.push_back('\t');
  cache.append(key);
  cache.push_back(':');

  auto index = m_header.find(cache);
  if (index != std::string::npos) {
    if (value.empty()) {
      cache = m_header.substr(0, index);
    } else {
      cache = m_header.substr(0, index + cache.size());
      cache.append(value);
    }

    index = m_header.find('\t', index + 1);
    if (index != std::string::npos) {
      cache.append(m_header.substr(index));
    }

    m_header = cache;
  } else if (!value.empty()) {
    m_header.append(cache);
    m_header.append(value);
  }
}

std::ostream&
operator<<(std::ostream& os, const read_group& value)
{
  return os << "read_group{id=" << log_escape(value.id())
            << ", header=" << log_escape(value.header()) << "}";
}

} // namespace adapterremoval
