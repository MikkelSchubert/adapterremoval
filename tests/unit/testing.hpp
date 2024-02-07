/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include "catch.hpp"
#include "fastq.hpp"
#include "strutils.hpp"
#include <ostream>
#include <sstream>

namespace Catch {

template<>
struct StringMaker<adapterremoval::fastq::ntrimmed>
{
  static std::string convert(adapterremoval::fastq::ntrimmed const& value)
  {
    std::stringstream ss;
    ss << "ntrimmed{" << value.first << ", " << value.second << "}";
    return ss.str();
  }
};

template<>
struct StringMaker<adapterremoval::fastq>
{
  static std::string convert(adapterremoval::fastq const& value)
  {
    std::string s;
    value.into_string(s);
    return adapterremoval::shell_escape(s);
  }
};

} // namespace Catch

namespace adapterremoval {

inline std::ostream&
operator<<(std::ostream& stream, const adapterremoval::fastq& record)
{
  return stream << "'@" << record.header() << "\\n"
                << record.sequence() << "\\n+\\n"
                << record.qualities() << "\\n'";
}

} // namespace adapterremoval
