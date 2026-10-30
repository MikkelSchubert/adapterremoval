// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "buffer.hpp"       // for buffer
#include "deflate.hpp"      // for calculate_crc32, deflate_block_isal_1, ...
#include "errors.hpp"       // for io_error
#include "linereader.hpp"   // declarations
#include "logging.hpp"      // for log_capture
#include "strutils.hpp"     // for split_lines, string_vec
#include "testing.hpp"      // for TEST_CASE, REQUIRE, ...
#include <cstdint>          // for uint32_t
#include <cstdio>           // for fclose, fwrite, rewind, tmpfile, FILE, fread
#include <initializer_list> // for initializer_list
#include <isa-l.h>          // for ISAL_GZIP, inflate_state, isal_inflate_ini
#include <libdeflate.h>     // for libdeflate_alloc_compressor, ...
#include <memory>           // for make_unique, unique_ptr
#include <string>           // for string
#include <string_view>      // for string_view, ""sv
#include <tuple>            // for tuple_size_v
#include <vector>           // for vector

using namespace std::string_view_literals;

namespace adapterremoval {

namespace {

template<typename T>
buffer
gzip_data(const T& data)
{
  using compressor_ptr = std::unique_ptr<libdeflate_compressor,
                                         decltype(&libdeflate_free_compressor)>;
  thread_local compressor_ptr compressor{ libdeflate_alloc_compressor(4),
                                          libdeflate_free_compressor };

  buffer buf{ data.size() + 64 };
  auto size = libdeflate_gzip_compress(compressor.get(),
                                       data.data(),
                                       data.size(),
                                       buf.data(),
                                       buf.size());
  REQUIRE(size > 0);
  buf.resize(size);

  return buf;
}

buffer
gzip_stream_strings(std::initializer_list<std::string_view> strings)
{
  size_t size = 0;
  uint32_t crc32 = 0;
  buffer buf = gzip_stream_header_isal();

  for (const auto* it = strings.begin(); it != strings.end(); ++it) {
    buffer src;
    src.append(*it);
    buffer dst{ it->size() + 64 };
    deflate_block_isal_1(src, dst, 0, it + 1 == strings.end());
    buf.append(dst);

    size += it->size();
    crc32 = calculate_crc32(src, crc32);
  }

  buf.append(gzip_stream_tail_isal(crc32, size));

  return buf;
}

FILE*
temporary_file(const buffer& buf)
{
  FILE* handle = std::tmpfile();
  REQUIRE(handle);

  const auto n = std::fwrite(buf.data(), 1, buf.size(), handle);
  REQUIRE(n == buf.size());

  REQUIRE(std::fseek(handle, 0, SEEK_SET) == 0);

  return handle;
}

string_vec
read_lines(line_reader& reader)
{
  string_vec lines;
  std::string buffer;

  while (reader.getline(buffer)) {
    lines.push_back(buffer);
  }

  return lines;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// Helper functions

TEST_CASE("temporary_file/gzip_data returns gzipped file")
{
  const size_t data_size = 1024;
  const auto data = random_text_buffer(data_size);

  const auto gzip = gzip_data(data);
  const auto gzip_size = gzip.size();
  FILE* handle = temporary_file(gzip);

  buffer src{ gzip_size + 1 }; // +1 to catch extra data being included
  const auto n = std::fread(src.data(), 1, src.size(), handle);
  REQUIRE(std::fclose(handle) == 0);
  REQUIRE(n == gzip_size);
  src.resize(n);

  buffer dst{ data_size };
  auto state = std::make_unique<inflate_state>();
  isal_inflate_init(state.get());
  state->crc_flag = ISAL_GZIP;
  state->next_in = src.data();
  state->avail_in = src.size();
  state->next_out = dst.data();
  state->avail_out = dst.size();
  checked_isal_inflate_stateless(state.get(), data_size);

  REQUIRE(dst == data);
}

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("line_reader throws on null")
{
  REQUIRE_THROWS_AS(line_reader(nullptr), assert_failed);
}

TEST_CASE("line_reader throws on missing file")
{
  // This file hopefully does not exist ...
  REQUIRE_THROWS_WITH(line_reader("eb8d324f-bae2-4484-894a-883b6abacc6b"),
                      Catch::Contains("failed to open file:"));
}

TEST_CASE("line_reader on empty file returns no lines")
{
  line_reader reader(temporary_file(buffer{}));
  REQUIRE(read_lines(reader).empty());
}

TEST_CASE("line_reader on empty line returns single, empty line")
{
  buffer input;
  input.append(GENERATE("\n", "\r\n"));
  if (GENERATE(true, false)) {
    input = gzip_data(input);
  }

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == string_vec({ "" }));
}

TEST_CASE("line_reader does not require trailing newline")
{
  buffer input = "HnSw4n WJAsds32emFfQdwO6Z\n"
                 "tiuSGNsadNZ\n"
                 "pAehX8GBlmyOPR"_buffer;
  input.append(GENERATE("", "\n"));

  if (GENERATE(true, false)) {
    input = gzip_data(input);
  }

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4n WJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR" });
}

TEST_CASE("line_reader does not require trailing newline (CRLF)")
{
  buffer input = "HnSw4n WJAsds32emFfQdwO6Z\r\n"
                 "tiuSGNsadNZ\r\n"
                 "pAehX8GBlmyOPR"_buffer;
  input.append(GENERATE("", "\r\n"));

  if (GENERATE(true, false)) {
    input = gzip_data(input);
  }

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4n WJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR" });
}

TEST_CASE("line_reader ignores trailing newline, but keep empty lines")
{
  buffer input = "HnSw4nWJAsds32emFfQdwO6Z\n"
                 "tiuSGNsadNZ\n"
                 "pAehX8GBlmyOPR\n"
                 "\n"
                 "\n"_buffer;

  if (GENERATE(true, false)) {
    input = gzip_data(input);
  }

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR",
                                            "",
                                            "" });
}

TEST_CASE("line_reader ignores trailing newline, but keep empty lines (CRLF)")
{
  buffer input = "HnSw4nWJAsds32emFfQdwO6Z\r\n"
                 "tiuSGNsadNZ\r\n"
                 "pAehX8GBlmyOPR\r\n"
                 "\r\n"
                 "\r\n"_buffer;

  if (GENERATE(true, false)) {
    input = gzip_data(input);
  }

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR",
                                            "",
                                            "" });
}

TEST_CASE("line_reader handles mixed LF/CRLF")
{
  buffer input = "HnSw4nWJAsds32emFfQdwO6Z\n"
                 "tiuSGNsadNZ\r\n"
                 "pAehX8GBlmyOPR\n"
                 "\n"
                 "\r\n"_buffer;

  if (GENERATE(true, false)) {
    input = gzip_data(input);
  }

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR",
                                            "",
                                            "" });
}

TEST_CASE("line_reader handles gzipped file in blocks")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ\n",
                             "pAehX8GBlmyOPR\n" };

  buffer buf;
  for (const auto& it : input) {
    buf.append(gzip_data(it));
  }

  line_reader reader(temporary_file(buf));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR" });
}

TEST_CASE("line_reader handles streamed gzipped file")
{
  buffer buf = gzip_stream_strings({
    "HnSw4nWJAsds32emFfQdwO6Z\n",
    "tiuSGNsadNZ\n",
    "pAehX8GBlmyOPR\n",
  });

  line_reader reader(temporary_file(buf));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR" });
}

TEST_CASE("line_reader handles gzipped file in blocks, with partial lines")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ",
                             "pAehX8GBlmyOPR\n" };

  buffer buf;
  for (const auto& it : input) {
    buf.append(gzip_data(it));
  }

  line_reader reader(temporary_file(buf));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZpAehX8GBlmyOPR" });
}

TEST_CASE("line_reader handles streamed gzipped file with partial lines")
{
  buffer buf = gzip_stream_strings({
    "HnSw4nWJAsds32emFfQdwO6Z\n",
    "tiuSGNsadNZ",
    "pAehX8GBlmyOPR\n",
  });

  line_reader reader(temporary_file(buf));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZpAehX8GBlmyOPR" });
}

TEST_CASE("line_reader handles gzipped file with empty blocks")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "",
                             "tiuSGNsadNZ\n",
                             "pAehX8GBlmyOPR\n" };

  buffer buf;
  for (const auto& it : input) {
    buf.append(gzip_data(it));
  }

  line_reader reader(temporary_file(buf));
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZ",
                                            "pAehX8GBlmyOPR" });
}

TEST_CASE("line_reader handles gzipped file with trailing junk")
{
  buffer buf;
  buf.append(gzip_data("HnSw4nWJAsds32emFfQdwO6Z\n"_buffer));
  buf.append("trailing junk");
  FILE* handle = temporary_file(buf);

  const log::log_capture ss;
  line_reader reader(handle);
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z" });
  REQUIRE_THAT(ss.str(),
               Catch::Matchers::Contains(
                 "[WARNING] Ignoring trailing garbage at the end of"));
}

TEST_CASE("line_reader handles multiple gzipped file with trailing junk")
{
  buffer buf;
  buf.append(gzip_data("HnSw4nWJAsds32emFfQdwO6Z\n"sv));
  buf.append(gzip_data("tiuSGNsadNZpAehX8GBlmyOPR\n"sv));
  buf.append("trailing junk");
  FILE* handle = temporary_file(buf);

  const log::log_capture ss;
  line_reader reader(handle);
  REQUIRE(read_lines(reader) == string_vec{ "HnSw4nWJAsds32emFfQdwO6Z",
                                            "tiuSGNsadNZpAehX8GBlmyOPR" });
  REQUIRE_THAT(ss.str(),
               Catch::Matchers::Contains(
                 "[WARNING] Ignoring trailing garbage at the end of"));
}

TEST_CASE("line_reader handles empty gzipped file with trailing junk")
{
  buffer buf;
  buf.append(gzip_data(buffer{}));
  buf.append("trailing junk");
  FILE* handle = temporary_file(buf);

  const log::log_capture ss;
  line_reader reader(handle);
  REQUIRE(read_lines(reader).empty());
  REQUIRE_THAT(ss.str(),
               Catch::Matchers::Contains(
                 "[WARNING] Ignoring trailing garbage at the end of"));
}

TEST_CASE("line reader always clears buffer")
{
  buffer input = "HnSw4nWJAsds32emFfQdwO6Z\n"
                 "tiuSGNsadNZ\n"
                 "pAehX8GBlmyOPR"_buffer;
  input.append(GENERATE("", "\n"));

  if (GENERATE(true, false)) {
    input = gzip_data(input);
  }

  line_reader reader(temporary_file(std::move(input)));

  std::string buffer;
  while (reader.getline(buffer)) {
    REQUIRE(!buffer.empty());
  }

  REQUIRE(buffer.empty());
  REQUIRE_FALSE(reader.getline(buffer));
  REQUIRE(buffer.empty());
}

TEST_CASE("line_reader handles blocks ends/starts straddling buffers")
{
  const auto buffer_size = std::tuple_size_v<line_reader::line_buffer>;

  const buffer header = gzip_stream_header_isal();
  const buffer input_1 = random_text_buffer(0xFFFF);

  buffer gzip;
  gzip.append(header);
  gzip.resize(buffer_size);
  deflate_block_uncompressed(input_1, gzip, header.size());
  gzip.append(gzip_stream_tail_isal(calculate_crc32(input_1), input_1.size()));

  const int shift = GENERATE(range(-8, 8));
  // remaining size - gzip header - deflate header - gzip tail - shift
  const size_t input_2_size =
    buffer_size - gzip.size() - header.size() - 5 - 8 - shift;
  const buffer input_2 = random_text_buffer(input_2_size);

  gzip.append(header);
  size_t offset = gzip.size();
  gzip.resize(buffer_size);
  deflate_block_uncompressed(input_2, gzip, offset);
  gzip.append(gzip_stream_tail_isal(calculate_crc32(input_2), input_2.size()));

  REQUIRE(gzip.size() == buffer_size - shift);

  // input_3 will start just after or partially overlapping the first raw buffer
  buffer input_3 = random_text_buffer(1024);
  gzip.append(gzip_data(input_3));

  line_reader reader(temporary_file(gzip));

  buffer input;
  input.append(input_1);
  input.append(input_2);
  input.append(input_3);
  auto expected_lines =
    split_lines(std::string_view{ reinterpret_cast<const char*>(input.data()),
                                  input.size() });
  // split_lines will include an empty string if the final line ends with '\n'
  if (expected_lines.back().empty()) {
    expected_lines.pop_back();
  }

  const auto lines = read_lines(reader);
  REQUIRE(lines.size() == expected_lines.size());
  REQUIRE(lines == expected_lines);
}

////////////////////////////////////////////////////////////////////////////////
// Error handling

TEST_CASE("line_reader throws on truncated files")
{
  buffer data = gzip_data(random_text_buffer(1024));

  data.resize(data.size() - GENERATE(1, 5, 9, 10));
  line_reader reader{ temporary_file(data) };

  REQUIRE_THROWS_CONTAINS(read_lines(reader), gzip_error, "file is truncated");
}

TEST_CASE("line_reader throws on truncated files (non-final)")
{
  buffer data;
  data.reserve(0x20000);
  data.append(gzip_stream_header_isal());
  const size_t header_size = data.size();
  data.resize(data.capacity());

  const buffer text = random_data_buffer(1024);
  deflate_block_isal_1(text, data, header_size, false);

  data.resize(data.size() - GENERATE(0, 1, 5, 9, 10));
  line_reader reader{ temporary_file(data) };

  REQUIRE_THROWS_CONTAINS(read_lines(reader), gzip_error, "file is truncated");
}

TEST_CASE("line_reader throws on invalid gzip header")
{
  buffer data = gzip_data(random_data_buffer(1024));
  data.put_u8(2, 0xFF); // CM

  line_reader reader{ temporary_file(data) };

  REQUIRE_THROWS_CONTAINS(read_lines(reader),
                          gzip_error,
                          "not a valid gzip file");
}

TEST_CASE("line_reader throws on invalid deflate header")
{
  buffer data = gzip_data(random_data_buffer(1024));
  data.put_u8(10, 0xFF);

  line_reader reader{ temporary_file(data) };

  REQUIRE_THROWS_CONTAINS(read_lines(reader),
                          gzip_error,
                          "invalid deflate block");
}

TEST_CASE("line_reader throws on corrupt files")
{
  buffer data = gzip_data(random_data_buffer(1024));
  data.put_u8(data.size() - 32, data.at(data.size() - 32) ^ 0xFF);

  line_reader reader{ temporary_file(data) };

  REQUIRE_THROWS_AS(read_lines(reader), gzip_error);
}

TEST_CASE("line_reader throws on wrong isize/crc32")
{
  buffer data = gzip_data(random_data_buffer(1024));
  const int offset = data.size() - GENERATE(4, 8); // isize (-4) or crc32 (-8)
  data.put_u32(offset, data.get_u32(offset) + GENERATE(-1, 1));

  line_reader reader{ temporary_file(data) };

  // isa-l treats both CRC32/ISIZE mismatches as checksum errors
  REQUIRE_THROWS_CONTAINS(read_lines(reader),
                          gzip_error,
                          "incorrect content checksum or ISIZE");
}

} // namespace adapterremoval
