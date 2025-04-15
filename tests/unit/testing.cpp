// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "testing.hpp"  // declarations
#include "strutils.hpp" // for log_escape
#include <cctype>       // for isprint

namespace {

std::string
format_exception(std::string_view name, const std::exception& value)
{
  Catch::ReusableStringStream os;
  os << name << "{ " << adapterremoval::log_escape(value.what()) << " }";
  return os.str();
}

} // namespace

namespace Catch {

std::string
Catch::StringMaker<std::invalid_argument, void>::convert(
  const std::invalid_argument& value)
{
  return format_exception("invalid_argument", value);
}

std::string
StringMaker<std::string>::convert(const std::string& str)
{
  if (getCurrentContext().getConfig()->showInvisibles()) {
    return adapterremoval::log_escape(str);
  }

  for (const auto c : str) {
    // Escape any special characters that cannot be distinguished at a glance
    if (!std::isprint(c) || (std::isspace(c) && (c != ' ' && c != '\n'))) {
      return adapterremoval::log_escape(str);
    }
  }

  return '"' + str + '"';
}

} // namespace Catch
