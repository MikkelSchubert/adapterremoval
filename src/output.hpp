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

#include "commontypes.hpp" // for string_vec, read_type, merge_strategy
#include <array>           // for array
#include <cstddef>         // for size_t
#include <memory>          // for unique_ptr
#include <string>          // for string
#include <string_view>     // for string_view
#include <vector>          // for vector

namespace adapterremoval {

class fastq;
class buffer;
class scheduler;
class userconfig;
class analytical_chunk;
enum class fastq_flags;

using chunk_ptr = std::unique_ptr<analytical_chunk>;
using chunk_ptr_vec = std::vector<chunk_ptr>;
using chunk_pair = std::pair<size_t, chunk_ptr>;
using chunk_vec = std::vector<chunk_pair>;

//! Path used to indicate that a file is not needed
const std::string DEV_NULL = "/dev/null";

struct output_file
{
  std::string name{};
  output_format format = output_format::fastq_gzip;

  bool operator==(const output_file& other) const
  {
    return name == other.name && format == other.format;
  }
};

//! Implemented in main_adapter_rm.cpp
size_t
add_write_step(scheduler& sch,
               const userconfig& config,
               const output_file& file);

/** Per sample output filenames / steps  */
class sample_output_files
{
public:
  sample_output_files();

  /** Sets the output file for a given read type. */
  void set_file(read_type rtype, output_file file);

  /** Unique output filenames, indexed using `offset` */
  [[nodiscard]] const std::vector<output_file>& files() const
  {
    return m_files;
  }

  /** Unique pipeline steps, indexed using `offset` */
  [[nodiscard]] const std::vector<size_t>& pipeline_steps() const
  {
    return m_pipeline_steps;
  }

  /** Returns the offset to the pipeline step/filename for a given read type. */
  [[nodiscard]] size_t offset(read_type value) const;

  /** Returns the format for a given offset */
  [[nodiscard]] output_format format(size_t offset) const
  {
    return m_files.at(offset).format;
  }

  /** Returns the filename for a given offset */
  [[nodiscard]] const std::string& filename(size_t offset) const
  {
    return m_files.at(offset).name;
  }

  //! Constant used to represent disabled output files/steps.
  static const size_t disabled;

private:
  friend class output_files;

  //! Unique output files. Multiple read types may be mapped to a file
  std::vector<output_file> m_files{};
  //! Unique pipeline steps IDs. Multiple read types may be mapped to a step
  std::vector<size_t> m_pipeline_steps{};
  //! Mapping of read types to filenames/steps.
  std::array<size_t, static_cast<size_t>(read_type::max)> m_offsets{};
};

/** Class used to organize filenames of output files. */
class output_files
{
public:
  output_files() = default;

  static bool parse_extension(const std::string& filename, output_format& sink);
  static std::string_view file_extension(output_format format);

  /** Constant indicating that a step has been disabled. */
  static const size_t disabled;

  //! JSON file containing settings / statistics
  std::string settings_json{};
  //! HTML file containing settings / statistics / plots
  std::string settings_html{};

  //! Filename for unidentified mate 1 reads (demultiplexing)
  output_file unidentified_1 = output_file{ DEV_NULL };
  //! Pipeline step responsible for compresssing/writing unidentified 1 reads
  size_t unidentified_1_step = disabled;
  //! Filename for unidentified mate 1 reads (demultiplexing)
  output_file unidentified_2 = output_file{ DEV_NULL };
  //! Pipeline step responsible for compresssing/writing unidentified 2 reads
  size_t unidentified_2_step = disabled;

  using sample_output_vec = std::vector<sample_output_files>;

  void add_sample(sample_output_files&& sample) { m_samples.push_back(sample); }

  /** Adds write steps for each  */
  void add_write_steps(scheduler& sch, const userconfig& config);

  [[nodiscard]] const sample_output_files& get_sample(size_t idx) const
  {
    return m_samples.at(idx);
  }

  [[nodiscard]] const sample_output_vec& samples() const { return m_samples; }

private:
  sample_output_vec m_samples{};
};

/** Helper class used to generate per file-type chunks for processed reads . */
class processed_reads
{
public:
  explicit processed_reads(const sample_output_files& map, bool first);

  /** Adds a read of the given type to be processed */
  void add(const fastq& read, read_type type, fastq_flags flags);

  /** Returns a chunk for each generated type of processed reads. */
  chunk_vec finalize(bool eof);

private:
  const sample_output_files& m_map;

  //! A set output chunks being created; typically fewer than read_type::max.
  chunk_ptr_vec m_chunks{};
};

/** Map of samples to downstream FASTQ processing/writing steps. */
class post_demux_steps
{
public:
  explicit post_demux_steps(const output_files& output);

  //! Step used to write unidentified mate 1 reads
  const size_t unidentified_1;
  //! Format of unidentified 1 reads
  const output_format unidentified_1_format;
  //! Step used to write unidentified mate 2 reads; may be disabled
  const size_t unidentified_2;
  //! Format of unidentified 1 reads
  const output_format unidentified_2_format;

  /* Processing step for each sample. */
  std::vector<size_t> samples{};

  /** Constant indicating that a step has been disabled. */
  static const size_t disabled;
};

/** Helper class used to generate per file-type chunks for processed reads . */
class demultiplexed_reads
{
public:
  explicit demultiplexed_reads(const post_demux_steps& steps);

  void add_unidentified_1(const fastq& read, fastq_flags flags);
  void add_unidentified_2(const fastq& read, fastq_flags flags);
  void add_read_1(fastq&& read, size_t sample);
  void add_read_2(fastq&& read, size_t sample);

  /** Returns a chunk for each generated type of processed reads. */
  chunk_vec flush(bool eof);

private:
  const post_demux_steps& m_steps;

  //! Cache of demultiplex reads; used to reduce the number of output chunks
  //! generated from each processed chunk, which would otherwise increase
  //! linearly with the number of barcodes.
  chunk_ptr_vec m_samples{};
  //! Cache of unidentified mate 1 reads
  chunk_ptr m_unidentified_1{};
  //! Format of unidentified 1 reads
  const output_format m_unidentified_1_format;
  //! Cache of unidentified mate 2 reads
  chunk_ptr m_unidentified_2{};
  //! Format of unidentified 1 reads
  const output_format m_unidentified_2_format;
};

} // namespace adapterremoval
