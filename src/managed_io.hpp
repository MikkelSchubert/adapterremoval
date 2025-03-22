// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "buffer.hpp" // for buffer_vec
#include <array>      // for array
#include <cstdio>     // for FILE
#include <string>     // for string

namespace adapterremoval {

/** Indicates if the writes should be flushed */
enum class flush
{
  on,
  off
};

/** Reader that closes unused writer handles if open files exceeds ulimits. */
class managed_reader
{
public:
  /** Take ownership of an existing file handle */
  explicit managed_reader(FILE* handle);
  /*
   * Opens a managed reader. If too many handles are used, this function will
   * close writers until the file can be successfully opened.
   */
  explicit managed_reader(std::string filename);

  /** Closes the handle if it has not been closed already */
  ~managed_reader();

  /** Closes the handle; no-op if the handle has been closed already */
  void close();

  /** Returns the filename of the (previously opened) file */
  const std::string& filename() const { return m_filename; }

  /** Reads `size` bytes into the destination buffer */
  size_t read(void* buffer, size_t size);

  managed_reader(const managed_reader&) = delete;
  managed_reader(managed_reader&& other) = delete;
  managed_reader& operator=(managed_reader&& other) = delete;
  managed_reader& operator=(const managed_reader&) = delete;

private:
  //! Source filename
  std::string m_filename{};
  //! File handle or nullptr if the file has been closed
  FILE* m_file = nullptr;

  friend class io_manager;
};

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
  /** Create lazy writer; does not open filename immediately */
  explicit managed_writer(std::string filename);
  /** Checks that the file handle has been closed */
  ~managed_writer();

  /** Write buffers, opening/creating the file as needed */
  void write(const buffer& buf, flush mode = flush::off);
  void write(const buffer_vec& buffers, flush mode = flush::off);
  void write(const std::string& buffer, flush mode = flush::off);

  /** Closes the handle; no-op if the handle has been closed already */
  void close();

  /** Returns the filename of the (previously opened) file */
  const std::string& filename() const { return m_filename; }

  managed_writer(const managed_writer&) = delete;
  managed_writer(managed_writer&& other) = delete;
  managed_writer& operator=(managed_writer&& other) = delete;
  managed_writer& operator=(const managed_writer&) = delete;

private:
  //! Destination filename; is created lazily.
  std::string m_filename{};
  //! Indicates if the file has been created/truncated
  bool m_created = false;
  //! Lazily opened, managed handle; may be closed to free up handles.
  FILE* m_file = nullptr;
  //! Indicates if the handle is a stream and can't be closed
  bool m_stream = false;

  //! Managed writer used more recently than this writer.
  managed_writer* m_prev = nullptr;
  //! Managed writer used prior to this writer.
  managed_writer* m_next = nullptr;

  friend class io_manager;
  friend class writer_lock;
};

} // namespace adapterremoval
