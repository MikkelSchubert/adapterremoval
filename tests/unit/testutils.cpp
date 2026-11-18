// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "testutils.hpp" // declarations
#include "bgzf.hpp"      // for bgzf_block_compress
#include "buffer.hpp"    // for buffer
#include "deflate.hpp"   // for calculate_crc32, deflate_block_isal_1, ...
#include "fastq.hpp"     // for fastq
#include "strutils.hpp"  // for stringify
#include "testing.hpp"   // for REQUIRE, RandomIntegerGenerator, ...
#include <cstdint>       // for uint8_t
#include <cstdio>        // for FILE
#include <libdeflate.h>  // for libdeflate_alloc_compressor, ...
#include <limits>        // for numeric_limits
#include <memory>        // for unique_ptr
#include <string>        // for string
#include <string_view>   // for string_view
#include <vector>        // for vector

namespace adapterremoval {

namespace {

using IntRng = Catch::Generators::RandomIntegerGenerator<size_t>;

/** Return FASTQ read with 150 bp of random DNA and 150 random quality scores */
fastq
random_fastq_impl(IntRng& rng, const size_t nth, const std::string_view mate)
{
  constexpr size_t RECORD_SIZE = 150;

  auto name = "@read" + stringify(nth);
  name.append(mate);

  std::string sequence;
  for (size_t j = 0; j < RECORD_SIZE; ++j, rng.next()) {
    sequence.push_back("ACGT"[rng.get()]);
  }

  std::string qualities;
  for (size_t j = 0; j < RECORD_SIZE; ++j, rng.next()) {
    qualities.push_back("!;@I"[rng.get()]);
  }

  return fastq{ std::move(name), std::move(sequence), std::move(qualities) };
}

buffer
random_fastq_buffer_impl(const size_t max_records, const size_t max_buffer_size)
{
  // Generates values in the range [0; 3]
  IntRng rng{ 0, 3 };

  buffer buf;
  for (size_t i = 1; buf.size() < max_buffer_size && i <= max_records; ++i) {
    auto record = random_fastq_impl(rng, i, {});
    buf.append(record.header());
    buf.append("\n");
    buf.append(record.sequence());
    buf.append("\n+\n");
    buf.append(record.qualities());
    buf.append("\n");
  }

  if (buf.size() > max_buffer_size) {
    buf.resize(max_buffer_size);
  }

  return buf;
}

void
gzip_and_append_data_to_buffer(buffer& dst,
                               const uint8_t* data,
                               const size_t len)
{
  using compressor_ptr = std::unique_ptr<libdeflate_compressor,
                                         decltype(&libdeflate_free_compressor)>;
  thread_local compressor_ptr compressor{ libdeflate_alloc_compressor(4),
                                          libdeflate_free_compressor };

  const size_t offset = dst.size();
  const size_t bound = libdeflate_gzip_compress_bound(compressor.get(), len);
  dst.resize(offset + bound);
  auto size = libdeflate_gzip_compress(compressor.get(),
                                       data,
                                       len,
                                       dst.data() + offset,
                                       dst.size() - offset);
  REQUIRE(size > 0);
  dst.resize(size + offset);
}

} // namespace

buffer
random_data_buffer(size_t size)
{
  buffer buf;
  buf.reserve(((size + 3) / 4) * 4);

  auto& rng = Catch::rng();
  for (size_t i = 0; i < size; i += 4) {
    buf.append_u32(rng());
  }

  buf.resize(size);
  return buf;
}

buffer
random_text_buffer(size_t size)
{
  return random_fastq_buffer_impl(std::numeric_limits<size_t>::max(), size);
}

buffer
random_fastq_buffer(size_t n)
{
  return random_fastq_buffer_impl(n, std::numeric_limits<size_t>::max());
}

std::vector<fastq>
random_fastq_records(size_t n, std::string_view mate)
{
  // Generates values in the range [0; 3]
  IntRng rng{ 0, 3 };

  std::vector<fastq> records;
  records.reserve(n);
  for (size_t i = 1; i <= n; ++i) {
    records.emplace_back(random_fastq_impl(rng, i, mate));
  }

  return records;
}

FILE*
temporary_file(const buffer& buf)
{
  FILE* handle = std::tmpfile();
  REQUIRE(handle);

  if (buf.size()) {
    const auto n = std::fwrite(buf.data(), 1, buf.size(), handle);
    REQUIRE(n == buf.size());

    REQUIRE(std::fseek(handle, 0, SEEK_SET) == 0);
  }

  return handle;
}

buffer
gzip_data_into_buffer(const buffer& data)
{
  buffer buf;
  gzip_and_append_data_to_buffer(buf, data.data(), data.size());
  return buf;
}

buffer
gzip_data_into_buffer(std::string_view data)
{
  buffer buf;
  gzip_and_append_data_to_buffer(buf,
                                 reinterpret_cast<const uint8_t*>(data.data()),
                                 data.size());
  return buf;
}

buffer
bgzf_block_into_buffer(const buffer& data)
{
  buffer buf;
  bgzf_block_compress(data, buf, 4);
  return buf;
}

void
gzip_and_append_data_to_buffer(buffer& dst, const buffer& data)
{
  gzip_and_append_data_to_buffer(dst, data.data(), data.size());
}

void
gzip_and_append_data_to_buffer(buffer& dst, std::string_view data)
{
  gzip_and_append_data_to_buffer(dst,
                                 reinterpret_cast<const uint8_t*>(data.data()),
                                 data.size());
}

buffer
gzip_stream_into_buffer(std::initializer_list<std::string_view> blocks)
{
  size_t size = 0;
  uint32_t crc32 = 0;
  buffer buf = gzip_stream_header_isal();

  for (const auto* it = blocks.begin(); it != blocks.end(); ++it) {
    buffer src;
    src.append(*it);
    buffer dst{ it->size() + 64 };
    deflate_block_isal_1(src, dst, 0, it + 1 == blocks.end());
    buf.append(dst);

    size += it->size();
    crc32 = calculate_crc32(src, crc32);
  }

  buf.append(gzip_stream_tail_isal(crc32, size));

  return buf;
}

} // namespace adapterremoval
