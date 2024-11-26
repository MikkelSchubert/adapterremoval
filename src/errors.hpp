/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include <ios>       // for ios_base::failure
#include <stdexcept> // for runtime_error
#include <string>    // for string

namespace adapterremoval {

std::string
format_io_error(const std::string& message, int error_number);

/** Exception explaining 'abort' calls when running unit-tests. */
class assert_failed : public std::logic_error
{
public:
  /** Creates exception with the specified error message. */
  explicit assert_failed(const std::string& what)
    : std::logic_error(what)
  {
  }
};

/** Represents errors during basic IO. */
class io_error : public std::ios_base::failure
{
public:
  /** Produces a combined error including a description of the error code */
  explicit io_error(const std::string& message, int error_number = 0)
    : std::ios_base::failure(message,
                             std::error_code(static_cast<int>(error_number),
                                             std::generic_category()))
  {
  }
};

/** Represents errors during GZip (de)compression. */
class gzip_error : public io_error
{
public:
  explicit gzip_error(const std::string& message)
    : io_error(message)
  {
  }
};

/** Exception raised for parsing and validation errors. */
class parsing_error : public std::runtime_error
{
public:
  explicit parsing_error(const std::string& message)
    : std::runtime_error(message)
  {
  }
};

/** Exception raised for FASTQ parsing and validation errors. */
class fastq_error : public parsing_error
{
public:
  explicit fastq_error(const std::string& message)
    : parsing_error(message)
  {
  }
};

} // namespace adapterremoval
