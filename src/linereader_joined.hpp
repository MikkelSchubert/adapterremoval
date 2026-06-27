// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2017 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "linereader.hpp" // for line_reader_base
#include <cstddef>        // for size_t
#include <memory>         // for unique_ptr
#include <string>         // for string
#include <vector>         // for vector

namespace adapterremoval {

/**
 * Class for tracking line-numbers across filenames and generating line ranges
 * that may span these. This is used for generating error messages for
 * multi-line records that may span file boundaries
 */
class joined_filenames
{
public:
  /** Creates list of filenames with unknown linecounts */
  explicit joined_filenames(std::vector<std::string> filenames);

  /** Increment the position in the current file */
  void inc_position(size_t n = 1);

  /** Returns the position across all files. Should not be used in logs */
  [[nodiscard]] size_t position() const noexcept { return m_position; }

  /** List filenames covered by an inclusive range of positions */
  [[nodiscard]] std::string filenames(size_t start, size_t end) const;

  /** Returns the remaining filenames */
  [[nodiscard]] size_t remaining_filenames() const noexcept
  {
    return m_filenames.size() - m_current_file;
  }

  /** Returns the current filename */
  [[nodiscard]] const std::string& current_filename() const;

  /** Set the line count for the currently active file and switch file */
  void next_file();

private:
  //! Files (to be) read and the number of lines read, if any
  std::vector<std::pair<std::string, size_t>> m_filenames{};

  //! Current active file
  size_t m_current_file = 0;
  //! Current position across all files
  size_t m_position = 0;
};

/**
 * Multi-file line-reader
 *
 * Wrapper around line_reader that automatically reads through one or more
 * files, returning the content as a contiguous stream of lines. No
 * assumptions are made about the format of the individual files.
 */
class joined_line_readers : public line_reader_base
{
public:
  /** Creates line-reader over multiple files in the specified order. */
  explicit joined_line_readers(std::vector<std::string> filenames);

  ~joined_line_readers() override = default;

  /**
   * Reads a line from the currently open file; if EOF is encountered, the
   * currently open file is closed and the next file is opened. Returns true
   * if a line was successfully read, or false if no files remain.
   */
  bool getline(std::string& dst) override;

  /** Global line number representing the current position across all files. */
  [[nodiscard]] size_t position() const noexcept;

  /** List filenames covering range of lines; must not exceed `position()` */
  [[nodiscard]] std::string filenames(size_t start, size_t end) const;

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

  //! Currently open file, if any.
  std::unique_ptr<line_reader> m_reader{};
  //! The list of filenames to process
  joined_filenames m_filenames;
};

} // namespace adapterremoval
