// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
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
