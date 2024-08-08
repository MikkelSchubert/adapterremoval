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

#include "commontypes.hpp" // for string_vec
#include <cstddef>         // for size_t
#include <functional>      // for function

namespace adapterremoval {

class fastq;
class buffer;
enum class output_format;

enum class fastq_flags
{
  //! SE read
  se,
  //! SE read that failed QC
  se_fail,
  //! Mate 1 read
  pe_1,
  //! Mate 1 read that failed QC
  pe_1_fail,
  //! Mate 1 read
  pe_2,
  //! Mate 1 read that failed QC
  pe_2_fail,

};

/** Contains SAM/BAM read-group information */
class read_group
{
public:
  read_group();

  /**
   * Parses a read-group string in the form "ID:1\tSM:sample" (optionally
   * including a leading "@RG\t"). Throws std::invalid_argument if the value
   * is invalid.
   */
  explicit read_group(std::string_view value);

  /** Returns the read-group ID for use in per-read 'RG' tags */
  [[nodiscard]] std::string_view id() const { return m_id; }

  /** Returns the full @RG header, not including a trailing new-line */
  [[nodiscard]] std::string_view header() const { return m_header; }

private:
  //! The full read_group header, including leading `@RG\t`
  std::string m_header{};
  //! Value mapping reads (via `RG:Z:${ID}`) to the @RG header
  std::string m_id{};
};

class fastq_serializer
{
public:
  static void header(buffer& buf, const string_vec& args, const read_group& rg);
  static void record(buffer& buf,
                     const fastq& record,
                     fastq_flags flags,
                     char mate_separator,
                     const read_group& rg);
};

class sam_serializer
{
public:
  static void header(buffer& buf, const string_vec& args, const read_group& rg);
  static void record(buffer& buf,
                     const fastq& record,
                     fastq_flags flags,
                     char mate_separator,
                     const read_group& rg);
};

class bam_serializer
{
public:
  static void header(buffer& buf, const string_vec& args, const read_group& rg);
  static void record(buffer& buf,
                     const fastq& record,
                     fastq_flags flags,
                     char mate_separator,
                     const read_group& rg);
};

class serializer
{
public:
  explicit serializer(output_format format);

  /** Set the read group; used to write headers/record for SAM/BAM */
  void set_read_group(const read_group& rg) { m_read_group = rg; }

  /** Set the mate separator; used to trim mate information for SAM/BAM */
  void set_mate_separator(char value) { m_mate_separator = value; }

  void header(buffer& buf, const string_vec& args) const;
  void record(buffer& buf, const fastq& record, fastq_flags flags) const;

private:
  //! Read group information to be added to SAM/BAM records
  read_group m_read_group{};
  //! Mate separator in processed reads
  char m_mate_separator = '\0';
  //! Function used to serialize file headers for processed reads
  std::function<decltype(fastq_serializer::header)> m_header{};
  //! Function used to serialize processed reads
  std::function<decltype(fastq_serializer::record)> m_record{};
};

} // namespace adapterremoval
