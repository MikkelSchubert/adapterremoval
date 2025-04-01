// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "testing.hpp"  // declarations
#include "strutils.hpp" // for log_escape

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

} // namespace Catch
