// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "output.hpp"     // declarations
#include "buffer.hpp"     // for buffer
#include "debug.hpp"      // for AR_REQUIRE, AR_FAIL
#include "fastq.hpp"      // for fastq
#include "fastq_io.hpp"   // for write_fastq, split_fastq, gzip_split_fastq
#include "scheduler.hpp"  // for analytical_chunk
#include "serializer.hpp" // for read_meta
#include "strutils.hpp"   // for ends_with
#include <limits>         // for numeric_limits

namespace adapterremoval {

class userconfig;

namespace {

buffer&
get_buffer(output_chunk_ptr& chunk)
{
  AR_REQUIRE(chunk);
  if (chunk->buffers.empty()) {
    chunk->buffers.emplace_back();
  }

  return chunk->buffers.back();
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
    case output_format::sam_gzip:
    case output_format::bam:
    case output_format::ubam: {
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
sample_output_files::set_file(const read_file rtype, output_file file)
{
  const auto index = static_cast<size_t>(rtype);
  AR_REQUIRE(m_offsets.at(index) == sample_output_files::disabled);

  // If the file type isn't being saved, then there is no need to process the
  // reads. This saves time especially when output compression is enabled.
  if (file.name != DEV_NULL) {
    // FIXME: This assumes that filesystem is case sensitive
    auto it = m_output.begin();
    for (; it != m_output.end(); ++it) {
      if (file.name == it->file.name) {
        break;
      }
    }

    if (it == m_output.end()) {
      m_output.push_back({ std::move(file), disabled });
      m_offsets.at(index) = m_output.size() - 1;
    } else {
      AR_REQUIRE(file == it->file);
      const auto existing_index = it - m_output.begin();
      m_offsets.at(index) = existing_index;
    }
  }
}

void
sample_output_files::set_step(read_file rtype, size_t step)
{
  const auto index = static_cast<size_t>(rtype);
  const auto offset = m_offsets.at(index);
  AR_REQUIRE(offset != disabled);

  AR_REQUIRE(m_output.at(offset).step == disabled);
  m_output.at(offset).step = step;
}

size_t
sample_output_files::offset(read_file value) const
{
  return m_offsets.at(static_cast<size_t>(value));
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `output_files`

const size_t output_files::disabled = static_cast<size_t>(-1);

bool
output_files::parse_format(std::string_view filename, output_format& sink)
{
  const std::string value = to_lower(std::string(filename));
  // ubam is tested separately, to avoid supporting non-standard file extensions
  if (value == "ubam") {
    sink = output_format::ubam;
    return true;
  }

  return parse_extension("." + value, sink);
}

bool
output_files::parse_extension(std::string_view filename, output_format& sink)
{
  const std::string value = to_lower(std::string(filename));
  if (ends_with(value, ".fq.gz") || ends_with(value, ".fastq.gz")) {
    sink = output_format::fastq_gzip;
  } else if (ends_with(value, ".fq") || ends_with(value, ".fastq")) {
    sink = output_format::fastq;
  } else if (ends_with(value, ".sam.gz")) {
    sink = output_format::sam_gzip;
  } else if (ends_with(value, ".sam")) {
    sink = output_format::sam;
  } else if (ends_with(value, ".bam")) {
    sink = output_format::bam;
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
    case output_format::bam:
    case output_format::ubam:
      return ".bam";
    default:
      AR_FAIL("invalid output format");
  }
}

void
output_files::add_write_steps(scheduler& sch, const userconfig& config)
{
  AR_REQUIRE(unidentified_1_step == disabled &&
             unidentified_2_step == disabled);

  for (auto& sample : m_samples) {
    for (auto& it : sample.m_output) {
      AR_REQUIRE(it.step == disabled);
      it.step = add_write_step(sch, config, it.file);
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

processed_reads::processed_reads(const sample_output_files& map)
  : m_map(map)
{
  for (size_t i = 0; i < map.size(); ++i) {
    m_chunks.emplace_back(std::make_unique<output_chunk>());
    m_serializers.emplace_back(map.format(i));
  }
}

processed_reads::~processed_reads() = default;

void
processed_reads::set_sample(const sample& value)
{
  for (auto& it : m_serializers) {
    it.set_sample(value);
  }
}

void
processed_reads::set_mate_separator(char value)
{
  m_mate_separator = value;
  for (auto& it : m_serializers) {
    it.set_mate_separator(value);
  }
}

void
processed_reads::set_demultiplexing_only(bool value)
{
  for (auto& it : m_serializers) {
    it.set_demultiplexing_only(value);
  }
}

void
processed_reads::write_headers(const string_vec& args)
{
  for (size_t i = 0; i < m_chunks.size(); ++i) {
    m_serializers.at(i).header(get_buffer(m_chunks.at(i)), args);
  }
}

void
processed_reads::add(const fastq& read, read_type flags, size_t barcode)
{
  add(read, read_meta(flags).barcode(barcode));
}

void
processed_reads::add(const fastq& read, const read_meta& meta)
{
  const size_t offset = m_map.offset(meta.get_file());
  if (offset != sample_output_files::disabled) {
    auto& buffer = get_buffer(m_chunks.at(offset));
    m_serializers.at(offset).record(buffer, read, meta);
    m_chunks.at(offset)->reads++;
  }
}

chunk_vec
processed_reads::finalize(bool eof)
{
  chunk_vec chunks;

  for (size_t i = 0; i < m_chunks.size(); ++i) {
    auto& chunk = m_chunks.at(i);
    chunk->eof = eof;

    chunks.emplace_back(m_map.step(i), std::move(chunk));
  }

  return chunks;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for `post_demux_steps`

const size_t post_demux_steps::disabled = output_files::disabled;

///////////////////////////////////////////////////////////////////////////////
// Implementations for `demultiplexed_reads`

namespace {

void
flush_chunk(chunk_vec& output,
            fastq_chunk_ptr& ptr,
            size_t step,
            const bool eof,
            const char mate_separator)
{
  if (eof || ptr->reads() >= INPUT_READS) {
    ptr->eof = eof;
    ptr->mate_separator = mate_separator;
    output.emplace_back(step, std::move(ptr));
    ptr = std::make_unique<fastq_chunk>();
  }
}

} // namespace

demultiplexed_reads::demultiplexed_reads(const post_demux_steps& steps)
  : m_steps(steps)
{
  // Position 0 is used for unidentified reads
  m_cache.push_back(std::make_unique<fastq_chunk>());

  for (const auto next_step : m_steps.samples) {
    AR_REQUIRE(next_step != post_demux_steps::disabled);

    m_cache.push_back(std::make_unique<fastq_chunk>());
  }

  // The first chunks should have headers written depending on format
  for (auto& it : m_cache) {
    it->first = true;
  }
}

void
demultiplexed_reads::add_unidentified_1(fastq&& read)
{
  m_cache.at(0)->reads_1.push_back(std::move(read));
}

void
demultiplexed_reads::add_unidentified_2(fastq&& read)
{
  m_cache.at(0)->reads_2.push_back(std::move(read));
}

void
demultiplexed_reads::add_read_1(fastq&& read, size_t sample, size_t barcode)
{
  auto& chunk = *m_cache.at(sample + 1);
  chunk.reads_1.push_back(std::move(read));
  chunk.barcodes.push_back(barcode);
}

void
demultiplexed_reads::add_read_2(fastq&& read, size_t sample)
{
  m_cache.at(sample + 1)->reads_2.push_back(std::move(read));
}

chunk_vec
demultiplexed_reads::flush(bool eof, char mate_separator)
{
  chunk_vec output;

  // Unidentified reads; these are treated as a pseudo-sample for simplicity
  flush_chunk(output, m_cache.at(0), m_steps.unidentified, eof, mate_separator);

  for (size_t i = 1; i < m_cache.size(); ++i) {
    flush_chunk(output,
                m_cache.at(i),
                m_steps.samples.at(i - 1),
                eof,
                mate_separator);
  }

  return output;
}

} // namespace adapterremoval
