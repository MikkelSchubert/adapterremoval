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

#include <cstdio>   // for FILE, fopen, fclose
#include <fstream>  // for ofstream
#include <memory>   // for unique_ptr
#include <stddef.h> // for size_t
#include <string>   // for string
#include <utility>  // for pair
#include <vector>   // for vector

#include "buffer.hpp" // for buffer_vec

namespace adapterremoval {

/**
 * Writer that manages open handles if open files exceeds ulimits.
 *
 * The file is lazily opened the first time a write is performed;
 * if the file cannot be opened due to the number of already open
 * files, the writer will close the least recently used handle and
 * retry.
 */
class managed_writer
{
public:
  explicit managed_writer(const std::string& filename);
  ~managed_writer();

  /**
   * Opens a file using fopen and returns the handle.
   *
   * If too many handles are used, this funtion will close writers until
   * the file can be succesfully opened.
   */
  static FILE* fopen(const std::string& filename, const char* mode);

  void write_buffer(const buffer& buf, bool flush);
  void write_buffers(const buffer_vec& buffers, bool flush);
  void write_string(const std::string& buffer, bool flush);

  void close();

  const std::string& filename() const;

  managed_writer(const managed_writer&) = delete;
  managed_writer& operator=(const managed_writer&) = delete;

private:
  /* Ensure that the writer is open, closing existing files if nessesary. */
  static void open_writer(managed_writer* ptr);
  /* Removes the writer from the list of open writers. */
  static void remove_writer(managed_writer* ptr);
  /* Sets the writer as the most recently used writer. */
  static void add_head_writer(managed_writer* ptr);
  /* Close the least recently used writer. */
  static void close_tail_writer();

  //! Destination filename; is created lazily.
  std::string m_filename;
  //! Lazily opened, managed handle; may be closed to free up handles.
  std::ofstream m_stream;
  //! Indicates if the file has been created
  bool m_created;

  //! Previous managed_writer; used more recently than this.
  managed_writer* m_prev;
  //! Next managed_writer; used before this.
  managed_writer* m_next;

  //! Most recently used managed_writer
  static managed_writer* s_head;
  //! Least recently used managed_writer
  static managed_writer* s_tail;
  //! Indicates if a performance warning has been printed
  static bool s_warning_printed;
};

} // namespace adapterremoval
