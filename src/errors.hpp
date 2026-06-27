// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <iosfwd>      // for ostream
#include <stdexcept>   // for runtime_error
#include <string>      // for string
#include <string_view> // for string_view

namespace adapterremoval {

std::string
format_io_error(std::string_view message, int error_number);

/** Exception explaining 'abort' calls when running unit-tests. */
class assert_failed : public std::logic_error
{
public:
  /** Creates exception with the specified error message. */
  explicit assert_failed(const std::string& what);

  friend std::ostream& operator<<(std::ostream& os, const assert_failed& value);
};

/** Base class for errors that should not represent bugs in the program */
class program_failure : public std::runtime_error
{
public:
  /** Produces an error without an associated error code */
  explicit program_failure(const std::string& message);

  /** Generic stream functions for program failure sub-classes */
  friend std::ostream& operator<<(std::ostream& os,
                                  const program_failure& value);

protected:
  //! The specific kind of error
  [[nodiscard]] virtual std::string_view kind() const = 0;
};

#define AR_DEFINE_FAILURE_FROM(subclass, name)                                 \
  class name : public subclass /** NOLINT(bugprone-macro-parentheses) */       \
  {                                                                            \
  public:                                                                      \
    explicit name(const std::string& message)                                  \
      : subclass(message)                                                      \
    {                                                                          \
    }                                                                          \
                                                                               \
  protected:                                                                   \
    [[nodiscard]] std::string_view kind() const override { return #name; }     \
  }

/** Represents errors during basic IO. */
class io_error : public program_failure
{
public:
  /** Produces an error without an associated error code */
  explicit io_error(const std::string& message);

  /** Produces a combined error including a description of the error code */
  explicit io_error(const std::string& message, int error_number);

protected:
  [[nodiscard]] std::string_view kind() const override { return "io_error"; }
};

/** Represents errors during GZip (de)compression. */
AR_DEFINE_FAILURE_FROM(io_error, gzip_error);

/** Exception raised for parsing and validation errors. */
AR_DEFINE_FAILURE_FROM(program_failure, parsing_error);

/** Exception raised for errors during serialization / encoding. */
AR_DEFINE_FAILURE_FROM(program_failure, serializing_error);

/** Exception raised for FASTQ parsing and validation errors. */
AR_DEFINE_FAILURE_FROM(parsing_error, fastq_error);

/** Exception raised to trigger a generic (testable) fatal error */
AR_DEFINE_FAILURE_FROM(program_failure, fatal_error);

} // namespace adapterremoval
