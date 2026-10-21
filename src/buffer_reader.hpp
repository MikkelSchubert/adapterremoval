// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "buffer.hpp"            // for buffer
#include "linereader_joined.hpp" // for joined_filenames
#include <cstddef>               // for size_t
#include <deque>                 // for deque
#include <functional>            // for function
#include <iosfwd>                // for ostream
#include <string>                // for string
#include <utility>               // for move, pair
#include <vector>                // for vector

namespace adapterremoval {

/**
 * Joined line-reader class that delegates reading of the file input to another
 * source; this source feeds buffers of raw text to the reader, which enables
 * calling the `getline` function until all complete lines have been read or EOF
 * is reached.
 *
 * Unlike the `line_reader` and `joined_line_reader` class, this reader does not
 * decompress compressed input. This task must be carried out prior to feeding
 * buffers to the `buffer_reader`.
 */
class buffer_reader
{
public:
  //! Enum indicating result of attempted getline call
  enum class result
  {
    //! Complete line read into buffer
    complete,
    //! Partial line (or nothing) was read; add more buffers to continue
    partial,
    //! No more lines to be read, across all files; dst is cleared
    eof,
  };

  /** Creates buffer-reader over multiple files in the specified order */
  explicit buffer_reader(std::vector<std::string> filenames);

  /** Releases remaining buffers */
  ~buffer_reader();

  /**
   * Add buffer to reader with an EOF flag indicating that the buffer represents
   * the tail end of the current file. It is not permitted to add buffers past
   * the EOF of the final, expected file specified in the constructor.
   */
  void add(buffer&& buf, bool eof = false);

  /**
   * Reads a (partial) line into `dst`. Returns `complete` if a full line has
   * been read, `partial` if only a partial line has been read, or `eof` if
   * nothing has been read into `dst` and no further data is available.
   *
   * `dst` is cleared when starting a new line and appended to when resuming a
   * partially read line. The caller should therefore pass the same `dst` until
   * `getline` returns `complete` or `eof`. `dst` is also cleared if `eof`.
   */
  result getline(std::string& dst);

  /** Returns the number of remaining buffers */
  [[nodiscard]] size_t buffer_count() const noexcept
  {
    return m_buffers.size();
  }

  /**
   * Returns true once the final buffer has been consumed; this may or may not
   * be preceded by `getline()` returning `eof`, depending on trailing newlines
   */
  [[nodiscard]] bool eof() const noexcept
  {
    return m_filenames.remaining_filenames() == 0;
  }

  /** Global number representing the current position across all read buffers */
  [[nodiscard]] size_t position() const noexcept
  {
    return m_filenames.position();
  }

  /** List names covering range of positions; must not exceed `position()` */
  [[nodiscard]] std::string names(size_t start, size_t end) const
  {
    return m_filenames.filenames(start, end);
  }

  /**
   * Set function used to release consumed buffers. This may for example be used
   * to return these to a global buffer.
   **/
  void set_release_function(std::function<void(buffer&&)> func)
  {
    m_release = std::move(func);
  }

  buffer_reader(const buffer_reader&) = delete;
  buffer_reader(buffer_reader&&) = delete;
  buffer_reader& operator=(const buffer_reader&) = delete;
  buffer_reader& operator=(buffer_reader&&) = delete;

private:
  /** Ready the next buffer, if any */
  void next_buffer();

  //! Buffers with EOF flags that have not yet been consumed
  std::deque<std::pair<buffer, bool>> m_buffers{};
  //! Current position in current buffer (if any)
  const char* m_buffer_ptr = nullptr;
  //! End position of current buffer (if any)
  const char* m_buffer_end = nullptr;
  //! The current size of the destination buffer
  size_t m_read_len = std::string::npos;

  //! Tracks line numbers across one or more buffers
  joined_filenames m_filenames;
  //! The remaining number of buffers with EOFs set
  size_t m_remaining_eof = 0;

  //! Optional function to release buffer, e.g. back to shared cache
  std::function<void(buffer&&)> m_release{};

  //! Used as memory location for empty buffers that may not have backing memory
  const char m_eof_char{};
};

std::ostream&
operator<<(std::ostream&, buffer_reader::result);

} // namespace adapterremoval
