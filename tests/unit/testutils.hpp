// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef>          // for size_t
#include <cstdio>           // for FILE
#include <initializer_list> // for initializer_list
#include <string_view>      // for string_view
#include <vector>           // for vector

namespace adapterremoval {

class buffer;
class fastq;

/** Returns buffer filled with `size` bytes of random data. Not thread-safe */
buffer
random_data_buffer(size_t size);

/** Returns buffer with `size` bytes of FASTQ like text. Not thread-safe */
buffer
random_text_buffer(size_t size);

/** Returns buffer containing `n` FASTQ records. Not thread-safe */
buffer
random_fastq_buffer(size_t n);

/** Returns vector containing `n` FASTQ records. Not thread-safe */
std::vector<fastq>
random_fastq_records(size_t n, std::string_view mate = {});

/** Returns handle to temporary file with the content set to `buf` */
FILE*
temporary_file(const buffer& buf);

/** Gzip compresses the provided data and returns the stream in a buffer */
buffer
gzip_data_into_buffer(const buffer& data);

/** Gzip compresses the provided data and returns the stream in a buffer */
buffer
gzip_data_into_buffer(std::string_view data);

/**
 * Creates single bgzf block from `data` and returns it in a buffer. The size of
 * `data` must be no larger than BGZF_MAX_INPUT_SIZE.
 */
buffer
bgzf_block_into_buffer(const buffer& data);

/** Gzip compresses the provided data and append it to a buffer */
void
gzip_and_append_data_to_buffer(buffer& dst, const buffer& data);

/** Gzip compresses the provided data and append it to a buffer */
void
gzip_and_append_data_to_buffer(buffer& dst, std::string_view data);

/** Create a gzip stream with N data blocks and return it in a buffer */
buffer
gzip_stream_into_buffer(std::initializer_list<std::string_view> blocks);

} // namespace adapterremoval
