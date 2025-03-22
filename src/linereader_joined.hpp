// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2017 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp" // for string_vec
#include "linereader.hpp"  // for line_reader_base
#include <cstddef>         // for size_t
#include <memory>          // for unique_ptr
#include <string>          // for string

namespace adapterremoval {

/**
 * Multi-file line-reader
 *
 * Wrapper around line_reader that automatically reads through one or more
 * files, returning the content as a contiguous stream of lines. No assumptions
 * are made about the format of the individual files.
 */
class joined_line_readers : public line_reader_base
{
public:
  /** Creates line-reader over multiple files in the specified order. */
  explicit joined_line_readers(const string_vec& filenames);

  /**
   * Reads a line from the currently open file; if EOF is encountered, the
   * currently open file is closed and the next file is opened. Returns true
   * if a line was successfully read, or false if no files remain.
   */
  bool getline(std::string& dst) override;

  /** Currently open file; empty if no file is open. */
  const std::string& filename() const { return m_filename; }

  /** Line number in the current file (1-based); 0 if no file is open. */
  size_t linenumber() const { return m_current_line; }

  joined_line_readers(const joined_line_readers&) = delete;
  joined_line_readers(joined_line_readers&&) = delete;
  joined_line_readers& operator=(const joined_line_readers&) = delete;
  joined_line_readers& operator=(joined_line_readers&&) = delete;

private:
  /**
   * Open the next file, removes it from the queue, and returns true; returns
   * false if no files remain to be processed.
   */
  bool open_next_file();

  //! Files left to read; stored in reverse order.
  string_vec m_filenames{};
  //! Currently open file, if any.
  std::unique_ptr<line_reader> m_reader{};
  //! The currently open file
  std::string m_filename{};
  //! Current line across all files.
  size_t m_current_line = 0;
};

} // namespace adapterremoval
