/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include "managed_io.hpp" // for managed_reader
#include <array>          // for array
#include <cstdio>         // for BUFSIZ
#include <ios>            // for ios_base, ios_base::failure
#include <isa-l.h>        // for ISAL_MAJOR_VERSION, ISAL_MINOR_VERSION
#include <memory>         // for shared_ptr, unique_ptr
#include <string>         // for string
#include <vector>         // for vector

struct inflate_state;
struct isal_gzip_header;

// Required due to 2.30 fixing issues when reading concatenated gzip files
static_assert((ISAL_MAJOR_VERSION > 2 ||
               (ISAL_MAJOR_VERSION == 2 && ISAL_MINOR_VERSION >= 30)),
              "isa-l v2.30+ required");

namespace adapterremoval {

//! Buffer used for compressed and uncompressed line data
using line_buffer = std::array<char, 10 * BUFSIZ>;

/** Base-class for line reading; used by receivers. */
class line_reader_base
{
public:
  /** Does nothing. */
  line_reader_base() = default;

  /** Closes the file, if still open. */
  virtual ~line_reader_base() = default;

  /** Reads a line into dst, returning false on EOF. */
  virtual bool getline(std::string& dst) = 0;

  //! Copy construction not supported
  line_reader_base(const line_reader_base&) = delete;
  //! Assignment not supported
  line_reader_base& operator=(const line_reader_base&) = delete;
};

/** Simple helper class for reading strings from memory */
class vec_reader : public line_reader_base
{
public:
  /** Constructs a reader from a set of lines (not copied) */
  vec_reader(const std::vector<std::string>& lines);

  /** Reads a line into dst, returning false on EOF. */
  bool getline(std::string& dst);

private:
  //! Reference to vector containing lines of text
  const std::vector<std::string>& m_lines;
  //! Current position in m_lines
  std::vector<std::string>::const_iterator m_it;
};

/**
 * Simple line reader.
 *
 * Currently reads
 *  - uncompressed files
 *  - gzip compressed files
 *
 * Errors are reported using either 'io_error' or 'gzip_error'.
 */
class line_reader : public line_reader_base
{
public:
  /** Creates a line handler from an existing handle */
  explicit line_reader(FILE* handle);
  /** Constructor; opens file and throws on errors. */
  explicit line_reader(std::string filename);

  /** Reads a line into dst, returning false on EOF. */
  bool getline(std::string& dst) override;

  //! Copy construction not supported
  line_reader(const line_reader&) = delete;
  //! Assignment not supported
  line_reader& operator=(const line_reader&) = delete;

private:
  //! Refills 'm_buffer' and sets 'm_buffer_ptr' and 'm_buffer_end'.
  void refill_buffers();

  //! Reader used to fill unparsed buffers
  managed_reader m_reader;
  /** Refills 'm_raw_buffer'; sets 'm_raw_buffer_ptr' and 'm_raw_buffer_end'. */
  void refill_raw_buffer(size_t avail_in = 0);
  /** Points 'm_buffer' and other points to corresponding 'm_raw_buffer's. */
  void refill_buffers_uncompressed();

  std::unique_ptr<inflate_state> m_gzip_stream;
  std::unique_ptr<isal_gzip_header> m_gzip_header;

  /** Returns true if the raw buffer contains gzip'd data. */
  bool is_raw_buffer_gzip() const;
  /** Initializes gzip stream and output buffers. */
  void initialize_buffers_gzip();
  /** Refills 'm_buffer' from compressed data; may refill raw buffers. */
  void refill_buffers_gzip();

  //! Pointer to buffer of decompressed data; may be equal to m_raw_buffer.
  std::shared_ptr<line_buffer> m_buffer;
  //! Pointer to current location in input buffer.
  line_buffer::iterator m_buffer_ptr;
  //! Pointer to end of current buffer.
  line_buffer::iterator m_buffer_end;

  //! Pointer to buffer of raw data.
  std::shared_ptr<line_buffer> m_raw_buffer;
  //! Pointer to end of current raw buffer.
  line_buffer::iterator m_raw_buffer_end;

  //! Indicates if a read across the EOF has been attempted.
  bool m_eof;
};

} // namespace adapterremoval
