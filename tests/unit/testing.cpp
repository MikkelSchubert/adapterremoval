// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "testing.hpp"  // declarations
#include "strutils.hpp" // for log_escape
#include <cctype>       // for isprint
#include <cstdint>      // for uint8_t
#include <exception>    // for exception
#include <stdexcept>    // for invalid_argument
#include <string>       // for string
#include <string_view>  // for string_view

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

// This definition replaces the one in `catch.hpp`, which has been commented out
std::string
StringMaker<std::string>::convert(const std::string& str)
{
  if (getCurrentContext().getConfig()->showInvisibles()) {
    return adapterremoval::log_escape(str);
  }

  for (const auto c : str) {
    // Escape any special characters that cannot be distinguished at a glance
    const auto uc = static_cast<uint8_t>(c);
    if (!std::isprint(uc) && c != '\n') {
      return adapterremoval::log_escape(str);
    }
  }

  return '"' + str + '"';
}

} // namespace Catch
