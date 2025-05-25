// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"     // for io_error
#include "linereader.hpp" // declarations
#include "logging.hpp"    // for log_capture
#include "testing.hpp"    // for TEST_CASE, REQUIRE, ...
#include <array>          // for array
#include <cstdio>         // for fclose, fwrite, rewind, tmpfile, FILE, fread
#include <libdeflate.h>   // for libdeflate_alloc_compressor, libdeflate_free...
#include <string>         // for string, basic_string, operator==
#include <vector>         // for operator==, vector

namespace adapterremoval {

using string_vec = std::vector<std::string>;

FILE*
temporary_file(const string_vec& lines)
{
  FILE* handle = std::tmpfile();

  for (const auto& line : lines) {
    const auto n = std::fwrite(line.data(), sizeof(char), line.size(), handle);
    REQUIRE(n == line.size());
  }

  std::rewind(handle);

  return handle;
}

FILE*
temporary_gzip_file(const string_vec& blocks)
{
  FILE* handle = std::tmpfile();

  auto* compressor = libdeflate_alloc_compressor(6);
  for (const auto& block : blocks) {
    std::array<char, 1024> buffer{};
    auto size = libdeflate_gzip_compress(compressor,
                                         block.data(),
                                         block.size(),
                                         buffer.data(),
                                         buffer.size());

    const auto n = std::fwrite(buffer.data(), sizeof(char), size, handle);
    REQUIRE(n == size);
  }

  libdeflate_free_compressor(compressor);

  std::rewind(handle);

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

  REQUIRE(buffer.empty());

  return lines;
}

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

////////////////////////////////////////////////////////////////////////////////
// Newlines (LF and CRLF)

TEST_CASE("temporary_gzip_file returns gzipped file")
{
  FILE* handle = temporary_gzip_file({ "foo" });
  std::array<char, 10> buffer{};
  auto n = std::fread(buffer.data(), sizeof(char), buffer.size(), handle);

  REQUIRE_FALSE(std::fclose(handle));
  REQUIRE(n == buffer.size());
  REQUIRE(buffer.at(0) == '\x1f');
  REQUIRE(buffer.at(1) == '\x8b');
}

TEST_CASE("line_reader on empty file returns no lines")
{
  line_reader reader(temporary_file({}));
  REQUIRE(read_lines(reader).empty());
}

TEST_CASE("line_reader on empty line returns single, empty line")
{
  line_reader reader(temporary_file({ "\n" }));
  REQUIRE(read_lines(reader) == string_vec({ "" }));
}

TEST_CASE("line_reader does not require trailing newline")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ\n",
                             "pAehX8GBlmyOPR" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZ",
                                "pAehX8GBlmyOPR" };

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader ignores trailing newline")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ\n",
                             "pAehX8GBlmyOPR\n" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZ",
                                "pAehX8GBlmyOPR" };

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader ignores trailing newline (CRLF)")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\r\n",
                             "tiuSGNsadNZ\r\n",
                             "pAehX8GBlmyOPR\r\n" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZ",
                                "pAehX8GBlmyOPR" };

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader ignores trailing newline, but keep empty lines")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ\n",
                             "pAehX8GBlmyOPR\n",
                             "\n",
                             "\n" };
  const string_vec expected = {
    "HnSw4nWJAsds32emFfQdwO6Z", "tiuSGNsadNZ", "pAehX8GBlmyOPR", "", "",
  };

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader ignores trailing newline, but keep empty lines (CRLF)")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\r\n",
                             "tiuSGNsadNZ\r\n",
                             "pAehX8GBlmyOPR\r\n",
                             "\r\n",
                             "\r\n" };
  const string_vec expected = {
    "HnSw4nWJAsds32emFfQdwO6Z", "tiuSGNsadNZ", "pAehX8GBlmyOPR", "", "",
  };

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader handles mixed LF/CRLF")
{
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ\r\n",
                             "pAehX8GBlmyOPR\n",
                             "\n",
                             "\r\n" };
  const string_vec expected = {
    "HnSw4nWJAsds32emFfQdwO6Z", "tiuSGNsadNZ", "pAehX8GBlmyOPR", "", "",
  };

  line_reader reader(temporary_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader handles gzipped file")
{
  // One long string/single block:
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n"
                             "tiuSGNsadNZ\n"
                             "pAehX8GBlmyOPR\n" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZ",
                                "pAehX8GBlmyOPR" };

  line_reader reader(temporary_gzip_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader handles gzipped file in blocks")
{
  // Multiple strings/blocks:
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ\n",
                             "pAehX8GBlmyOPR\n" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZ",
                                "pAehX8GBlmyOPR" };

  line_reader reader(temporary_gzip_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader handles gzipped file in blocks, with partial lines")
{
  // Multiple strings/blocks:
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZ",
                             "pAehX8GBlmyOPR\n" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZpAehX8GBlmyOPR" };

  line_reader reader(temporary_gzip_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader handles gzipped file with empty blocks")
{
  // Multiple strings/blocks:
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "",
                             "tiuSGNsadNZ\n",
                             "pAehX8GBlmyOPR\n" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZ",
                                "pAehX8GBlmyOPR" };

  line_reader reader(temporary_gzip_file(input));
  REQUIRE(read_lines(reader) == expected);
}

TEST_CASE("line_reader handles gzipped file with trailing junk")
{
  // Multiple strings/blocks:
  const string_vec input = { "HnSw4nWJAsds32emFfQdwO6Z\n",
                             "tiuSGNsadNZpAehX8GBlmyOPR\n" };
  const string_vec expected = { "HnSw4nWJAsds32emFfQdwO6Z",
                                "tiuSGNsadNZpAehX8GBlmyOPR" };

  FILE* handle = temporary_gzip_file(input);
  REQUIRE(std::fseek(handle, 0, SEEK_END) == 0);
  REQUIRE(std::fwrite("trailing junk", sizeof(char), 13, handle) == 13);
  std::rewind(handle);

  log::log_capture ss;
  line_reader reader(handle);
  REQUIRE(read_lines(reader) == expected);
  REQUIRE_THAT(ss.str(),
               Catch::Matchers::Contains(
                 "[WARNING] Ignoring trailing garbage at the end of"));
}

} // namespace adapterremoval
