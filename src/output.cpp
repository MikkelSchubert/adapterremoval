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
#include "output.hpp"    // declarations
#include "buffer.hpp"    // for buffer
#include "debug.hpp"     // for AR_REQUIRE, AR_FAIL
#include "fastq.hpp"     // for fastq
#include "fastq_io.hpp"  // for write_fastq, split_fastq, gzip_split_fastq
#include "scheduler.hpp" // for analytical_chunk
#include "strutils.hpp"  // for ends_with
#include <limits>        // for numeric_limits

namespace adapterremoval {

class userconfig;

////////////////////////////////////////////////////////////////////////////////
// Implementations for `sample_output_files`

const size_t sample_output_files::disabled = std::numeric_limits<size_t>::max();

size_t
add_write_step(scheduler& sch,
               const userconfig& config,
               const output_file& file)
{
  AR_REQUIRE(file.name != DEV_NULL);
  size_t step_id = sch.add<write_fastq>(config, file);

  switch (file.format) {
    case output_format::fastq:
      break;
    case output_format::fastq_gzip: {
      step_id = sch.add<gzip_split_fastq>(config, file, step_id);
      step_id = sch.add<split_fastq>(config, file, step_id);
      break;
    }
    default:
      AR_FAIL("invalid output format");
  }

  return step_id;
}

sample_output_files::sample_output_files()
{
  m_offsets.fill(sample_output_files::disabled);
}

void
sample_output_files::set_file(const read_type rtype, output_file file)
{
  const auto index = static_cast<size_t>(rtype);
  AR_REQUIRE(m_offsets.at(index) == sample_output_files::disabled);

  // If the file type isn't being saved, then there is no need to process the
  // reads. This saves time especially when output compression is enabled.
  if (file.name != DEV_NULL) {
    // FIXME: This assumes that filesystem is case sensitive
    auto it = m_files.begin();
    for (; it != m_files.end(); ++it) {
      if (file.name == it->name) {
        break;
      }
    }

    if (it == m_files.end()) {
      m_files.push_back(std::move(file));

      m_offsets.at(index) = m_files.size() - 1;
    } else {
      AR_REQUIRE(file == *it);
      const auto existing_index = it - m_files.begin();
      m_offsets.at(index) = existing_index;
    }
  }
}

size_t
sample_output_files::offset(read_type value) const
{
  return m_offsets.at(static_cast<size_t>(value));
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `output_files`

const size_t output_files::disabled = static_cast<size_t>(-1);

bool
output_files::parse_extension(const std::string& filename, output_format& sink)
{
  const std::string value = "." + to_lower(filename);
  if (ends_with(value, ".fq.gz") || ends_with(value, ".fastq.gz")) {
    sink = output_format::fastq_gzip;
    return true;
  } else if (ends_with(value, ".fq") || ends_with(value, ".fastq")) {
    sink = output_format::fastq;
    return true;
  }

  return false;
}

std::string_view
output_files::file_extension(const output_format format)
{
  switch (format) {
    case output_format::fastq:
      return ".fastq";
      break;
    case output_format::fastq_gzip:
      return ".fastq.gz";
      break;
    default:
      AR_FAIL("invalid output format");
  }
}

void
output_files::add_write_steps(scheduler& sch, const userconfig& config)
{
  for (auto& sample : m_samples) {
    AR_REQUIRE(sample.m_pipeline_steps.empty());
    for (const auto& file : sample.m_files) {
      sample.m_pipeline_steps.push_back(add_write_step(sch, config, file));
    }
  }

  if (unidentified_1.name != DEV_NULL) {
    unidentified_1_step = add_write_step(sch, config, unidentified_1);
  }

  if (unidentified_1.name != unidentified_2.name &&
      unidentified_2.name != DEV_NULL) {
    unidentified_2_step = add_write_step(sch, config, unidentified_2);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `processed_reads`

processed_reads::processed_reads(const sample_output_files& map, const bool eof)
  : m_map(map)
{
  const auto pipeline_steps = map.pipeline_steps().size();
  AR_REQUIRE(map.files().size() == pipeline_steps);

  for (size_t i = 0; i < pipeline_steps; ++i) {
    m_chunks.emplace_back(std::make_unique<analytical_chunk>(eof));
  }
}

void
processed_reads::add(const fastq& read, const read_type type)
{
  const size_t offset = m_map.offset(type);
  if (offset != sample_output_files::disabled) {
    m_chunks.at(offset)->add(read);
  }
}

chunk_vec
processed_reads::finalize()
{
  chunk_vec chunks;

  for (size_t i = 0; i < m_chunks.size(); ++i) {
    chunks.emplace_back(m_map.pipeline_steps().at(i),
                        std::move(m_chunks.at(i)));
  }

  return chunks;
}

} // namespace adapterremoval
