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
#include "output.hpp"     // declarations
#include "buffer.hpp"     // for buffer
#include "debug.hpp"      // for AR_REQUIRE, AR_FAIL
#include "fastq.hpp"      // for fastq
#include "fastq_io.hpp"   // for write_fastq, split_fastq, gzip_split_fastq
#include "scheduler.hpp"  // for analytical_chunk
#include "serializer.hpp" // for fastq_serialzier
#include "strutils.hpp"   // for ends_with
#include <limits>         // for numeric_limits

namespace adapterremoval {

class userconfig;

namespace {

buffer&
get_buffer(chunk_ptr& chunk)
{
  AR_REQUIRE(chunk);
  if (chunk->buffers.empty()) {
    chunk->buffers.emplace_back();
  }

  return chunk->buffers.back();
}

void
serialize_header(chunk_ptr& chunk, const output_format format)
{
  fastq_serializer::header(get_buffer(chunk), format);
}

void
serialize_record(chunk_ptr& chunk,
                 const fastq& record,
                 const output_format format,
                 const fastq_flags flags)
{
  fastq_serializer::record(get_buffer(chunk), record, format, flags);
  chunk->nucleotides += record.length();
}

} // namespace

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
    case output_format::sam:
      break;
    case output_format::fastq_gzip:
    case output_format::sam_gzip: {
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
  } else if (ends_with(value, ".fq") || ends_with(value, ".fastq")) {
    sink = output_format::fastq;
  } else if (ends_with(value, ".sam.gz")) {
    sink = output_format::sam_gzip;
  } else if (ends_with(value, ".sam")) {
    sink = output_format::sam;
  } else {
    return false;
  }

  return true;
}

std::string_view
output_files::file_extension(const output_format format)
{
  switch (format) {
    case output_format::fastq:
      return ".fastq";
    case output_format::fastq_gzip:
      return ".fastq.gz";
    case output_format::sam:
      return ".sam";
    case output_format::sam_gzip:
      return ".sam.gz";
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

  if (unidentified_1.name == unidentified_2.name) {
    unidentified_2_step = unidentified_1_step;
  } else if (unidentified_2.name != DEV_NULL) {
    unidentified_2_step = add_write_step(sch, config, unidentified_2);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `processed_reads`

processed_reads::processed_reads(const sample_output_files& map, bool first)
  : m_map(map)
{
  const auto pipeline_steps = map.pipeline_steps().size();
  AR_REQUIRE(map.files().size() == pipeline_steps);

  for (size_t i = 0; i < pipeline_steps; ++i) {
    m_chunks.emplace_back(std::make_unique<analytical_chunk>());
    if (first) {
      serialize_header(m_chunks.back(), m_map.format(i));
    }
  }
}

void
processed_reads::add(const fastq& read,
                     const read_type type,
                     const fastq_flags flags)
{
  const size_t offset = m_map.offset(type);
  if (offset != sample_output_files::disabled) {
    serialize_record(m_chunks.at(offset), read, m_map.format(offset), flags);
  }
}

chunk_vec
processed_reads::finalize(bool eof)
{
  chunk_vec chunks;

  for (size_t i = 0; i < m_chunks.size(); ++i) {
    const auto next_step = m_map.pipeline_steps().at(i);
    auto& chunk = m_chunks.at(i);
    chunk->eof = eof;

    chunks.emplace_back(next_step, std::move(chunk));
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `post_demux_steps`

const size_t post_demux_steps::disabled = output_files::disabled;

post_demux_steps::post_demux_steps(const output_files& output)
  : unidentified_1(output.unidentified_1_step)
  , unidentified_1_format(output.unidentified_1.format)
  , unidentified_2(output.unidentified_2_step)
  , unidentified_2_format(output.unidentified_2.format)
{
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `demultiplexed_reads`

namespace {

template<typename T>
void
flush_chunk(chunk_vec& output, std::unique_ptr<T>& ptr, size_t step, bool eof)
{
  if (eof || ptr->nucleotides >= INPUT_BLOCK_SIZE) {
    ptr->eof = eof;
    output.push_back(chunk_pair(step, std::move(ptr)));
    ptr = std::make_unique<T>();
  }
}

} // namespace

demultiplexed_reads::demultiplexed_reads(const post_demux_steps& steps)
  : m_steps(steps)
  , m_unidentified_1_format(steps.unidentified_1_format)
  , m_unidentified_2_format(steps.unidentified_2_format)
{
  if (m_steps.unidentified_1 != post_demux_steps::disabled) {
    m_unidentified_1 = std::make_unique<analytical_chunk>();
    serialize_header(m_unidentified_1, m_unidentified_1_format);
  }

  if (m_steps.unidentified_2 != post_demux_steps::disabled &&
      m_steps.unidentified_1 != m_steps.unidentified_2) {
    m_unidentified_2 = std::make_unique<analytical_chunk>();
    serialize_header(m_unidentified_2, m_unidentified_2_format);
  }

  for (const auto next_step : m_steps.samples) {
    AR_REQUIRE(next_step != post_demux_steps::disabled);

    m_samples.push_back(std::make_unique<analytical_chunk>());
  }
}

void
demultiplexed_reads::add_unidentified_1(const fastq& read,
                                        const fastq_flags flags)
{
  if (m_unidentified_1) {
    serialize_record(m_unidentified_1, read, m_unidentified_1_format, flags);
  }
}

void
demultiplexed_reads::add_unidentified_2(const fastq& read,
                                        const fastq_flags flags)
{
  if (m_unidentified_2) {
    serialize_record(m_unidentified_2, read, m_unidentified_1_format, flags);
  } else if (m_steps.unidentified_1 == m_steps.unidentified_2) {
    if (m_unidentified_1) {
      serialize_record(m_unidentified_1, read, m_unidentified_2_format, flags);
    }
  }
}

void
demultiplexed_reads::add_read_1(fastq&& read, size_t sample)
{
  m_samples.at(sample)->reads_1.push_back(std::move(read));
}

void
demultiplexed_reads::add_read_2(fastq&& read, size_t sample)
{
  m_samples.at(sample)->reads_2.push_back(std::move(read));
}

chunk_vec
demultiplexed_reads::flush(bool eof)
{
  chunk_vec output;

  if (m_unidentified_1) {
    flush_chunk(output, m_unidentified_1, m_steps.unidentified_1, eof);
  }

  if (m_unidentified_2) {
    flush_chunk(output, m_unidentified_2, m_steps.unidentified_2, eof);
  }

  for (size_t nth = 0; nth < m_samples.size(); ++nth) {
    flush_chunk(output, m_samples.at(nth), m_steps.samples.at(nth), eof);
  }

  return output;
}

} // namespace adapterremoval
