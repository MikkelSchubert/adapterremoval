// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "buffer.hpp"        // for buffer
#include "buffer_reader.hpp" // for buffer_reader
#include "catch.hpp"
#include "errors.hpp"  // for assert failed
#include "testing.hpp" // for TEST_CASE, REQUIRE, CHECK, ...

namespace adapterremoval {

TEST_CASE("buffer_reader requires non-empty filename list")
{
  REQUIRE_THROWS_AS(buffer_reader({}), assert_failed);
}

TEST_CASE("buffer_reader returns partial from start")
{
  buffer_reader reader({ "reads.fastq" });
  REQUIRE(reader.buffer_count() == 0);
  REQUIRE_FALSE(reader.eof());
  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(reader.position() == 0);
  REQUIRE_FALSE(reader.eof());
  REQUIRE(line.empty());
}

TEST_CASE("buffer_reader with one filename")
{
  buffer_reader reader({ "reads.fastq" });
  REQUIRE(reader.buffer_count() == 0);
  REQUIRE_FALSE(reader.eof());
  reader.add(buffer(), true);
  REQUIRE(reader.buffer_count() == 1);
  REQUIRE_FALSE(reader.eof());

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(reader.buffer_count() == 0);
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans empty buffers")
{
  buffer_reader reader{ { "file1" } };
  reader.add(""_buffer);
  reader.add("line "_buffer);
  reader.add(""_buffer);
  reader.add("without "_buffer);
  reader.add(""_buffer);
  reader.add("newline"_buffer);
  reader.add(""_buffer);
  reader.add("\n"_buffer);
  reader.add(""_buffer);
  reader.add(""_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line without newline");
  REQUIRE_FALSE(reader.eof());

  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans empty buffers, with partial")
{
  buffer_reader reader{ { "file1" } };
  reader.add(""_buffer);
  reader.add("line "_buffer);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(line == "line ");
  REQUIRE_FALSE(reader.eof());

  reader.add(""_buffer);
  reader.add("without "_buffer);
  reader.add(""_buffer);
  reader.add("newline"_buffer);

  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(line == "line without newline");

  reader.add(""_buffer);
  reader.add("\n"_buffer);
  reader.add(""_buffer);
  reader.add(""_buffer, true);

  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line without newline");

  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans empty buffers for EOF")
{
  buffer_reader reader{ { "file1" } };
  reader.add("line"_buffer);
  reader.add(""_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line");
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans empty buffers for EOF, with partial")
{
  buffer_reader reader{ { "file1" } };
  reader.add("line"_buffer);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(line == "line");
  REQUIRE_FALSE(reader.eof());

  reader.add(""_buffer, true);
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line");
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans empty buffers ending on EOF")
{
  buffer_reader reader{ { "file1" } };
  reader.add("line\n"_buffer);
  reader.add(""_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans empty buffers ending on EOF, with partial")
{
  buffer_reader reader{ { "file1" } };
  reader.add("line\n"_buffer);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line");
  REQUIRE_FALSE(reader.eof());

  reader.add(""_buffer, true);
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader does not span EOF")
{
  buffer_reader reader{ { "file1", "file2" } };
  reader.add("line without newline"_buffer, true);
  reader.add("proper line\n"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line without newline");
  REQUIRE_FALSE(reader.eof());

  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "proper line");
  REQUIRE_FALSE(reader.eof());

  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans multiple files with newline")
{
  buffer_reader reader{ { "file1", "file2" } };
  reader.add("line 1\nline 2\n"_buffer, true);
  reader.add("line 3\nline 4\n"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 2");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 3");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 4");
  REQUIRE_FALSE(reader.eof());

  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans multiple files with newline, with partial")
{
  buffer_reader reader{ { "file1", "file2" } };
  reader.add("line 1\nline 2\n"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 2");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(line.empty());
  REQUIRE_FALSE(reader.eof());

  reader.add("line 3\nline 4\n"_buffer, true);

  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 3");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 4");
  REQUIRE_FALSE(reader.eof());

  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
  REQUIRE(reader.eof());
}

TEST_CASE("buffer_reader spans multiple files without newline")
{
  buffer_reader reader{ { "file1", "file2" } };
  reader.add("line 1\nline 2"_buffer, true);
  reader.add("line 3\nline 4"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 2");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 3");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 4");
  REQUIRE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
}

TEST_CASE("buffer_reader spans multiple files without newline, with partial")
{
  buffer_reader reader{ { "file1", "file2" } };
  reader.add("line 1\nline 2"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 2");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(line.empty());
  REQUIRE_FALSE(reader.eof());

  reader.add("line 3\nline 4"_buffer, true);

  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 3");
  REQUIRE_FALSE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 4");
  REQUIRE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
}

TEST_CASE("buffer_reader returns empty lines, except trailing")
{
  buffer_reader reader{ { "dummy" } };
  reader.add("1\n\n2\n\n\n3\n"_buffer, true);

  std::string line;
  REQUIRE(reader.position() == 0);
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 1);
  REQUIRE(line == "1");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 2);
  REQUIRE(line.empty());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 3);
  REQUIRE(line == "2");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 4);
  REQUIRE(line.empty());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 5);
  REQUIRE(line.empty());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 6);
  REQUIRE(line == "3");
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(reader.position() == 6); // unchanged
  REQUIRE(line.empty());
}

TEST_CASE("buffer_reader trims carriage returns")
{
  buffer_reader reader{ { "filename" } };
  reader.add("line 1\r\nline 2\r"_buffer);
  reader.add("\nline 3"_buffer);
  reader.add("\r\nline 4"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 2");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 3");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 4");
}

// This documents current behavior; this should perhaps be changed/fixed
TEST_CASE("buffer_reader does not trim trailing carriage return")
{
  buffer_reader reader{ { "filename" } };
  reader.add("line 1\r"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1\r");
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
}

TEST_CASE("buffer_reader does not accept buffers after last EOF")
{
  buffer_reader reader{ { "filename" } };
  reader.add(buffer(), true);

  // VARIANTS: Neither content nor eof flag should affect this check
  buffer buf;
  if (GENERATE(true, false)) {
    buf.append("content");
  }
  const auto eof = GENERATE(true, false);
  REQUIRE_THROWS_AS(reader.add(std::move(buf), eof), assert_failed);
}

TEST_CASE("buffer_reader asserts if destination is modified")
{
  buffer_reader reader{ { "filename" } };
  reader.add("part 1"_buffer);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(line == "part 1");

  reader.add("|part 2"_buffer, true);

  SECTION("reading again succeeds")
  {
    REQUIRE(reader.getline(line) == buffer_reader::result::complete);
    REQUIRE(line == "part 1|part 2");
  }

  SECTION("shorter dest fails")
  {
    line.pop_back();
    REQUIRE_THROWS_AS(reader.getline(line), assert_failed);
    REQUIRE(line == "part ");
  }

  SECTION("longer dest fails")
  {
    line.push_back('x');
    REQUIRE_THROWS_AS(reader.getline(line), assert_failed);
    REQUIRE(line == "part 1x");
  }
}

////////////////////////////////////////////////////////////////////////////////
// buffer_reader -- positions and filenames

TEST_CASE("buffer_reader position in empty file")
{
  buffer_reader reader{ { "dummy" } };
  reader.add(buffer{}, true);

  CHECK(reader.position() == 0);
  CHECK(reader.names(0, 0) == "'dummy' at line 1");

  std::string line;
  CHECK(reader.getline(line) == buffer_reader::result::eof);
  CHECK(reader.position() == 0);
  CHECK(reader.names(0, 0) == "end of input files");
}

TEST_CASE("buffer_reader position in file with trailing newline")
{
  buffer_reader reader{ { "dummy" } };
  reader.add("line 1\nline 2\nline 3\n"_buffer, true);

  std::string line;
  CHECK(reader.position() == 0);
  CHECK(reader.names(0, 0) == "'dummy' at line 1");
  CHECK(reader.getline(line) == buffer_reader::result::complete);
  CHECK(reader.position() == 1);
  CHECK(reader.names(1, 1) == "'dummy' at line 2");
  CHECK(reader.getline(line) == buffer_reader::result::complete);
  CHECK(reader.position() == 2);
  CHECK(reader.names(2, 2) == "'dummy' at line 3");
  CHECK(reader.getline(line) == buffer_reader::result::complete);
  CHECK(reader.position() == 3);
  CHECK(reader.names(3, 3) == "'dummy' at line 4");
  CHECK(reader.getline(line) == buffer_reader::result::eof);
  CHECK(reader.position() == 3); // unchanged
  CHECK(reader.names(3, 3) == "end of input files");
}

TEST_CASE("buffer_reader position in file without trailing newline")
{
  buffer_reader reader{ { "dummy" } };
  reader.add("line 1\nline 2\nline 3"_buffer, true);

  std::string line;
  CHECK(reader.position() == 0);
  CHECK(reader.names(0, 0) == "'dummy' at line 1");
  CHECK(reader.getline(line) == buffer_reader::result::complete);
  CHECK(reader.position() == 1);
  CHECK(reader.names(1, 1) == "'dummy' at line 2");
  CHECK(reader.getline(line) == buffer_reader::result::complete);
  CHECK(reader.position() == 2);
  CHECK(reader.names(2, 2) == "'dummy' at line 3");
  CHECK(reader.getline(line) == buffer_reader::result::complete); // eof
  CHECK(reader.position() == 3);
  CHECK(reader.names(3, 3) == "end of input files");
  CHECK(reader.getline(line) == buffer_reader::result::eof);
  CHECK(reader.position() == 3);                     // unchanged
  CHECK(reader.names(3, 3) == "end of input files"); // unchanged
}

TEST_CASE("buffer_reader position across multiple files")
{
  buffer_reader reader{ { "file1", "file2", "file3" } };
  reader.add("line 1\nline 2\n"_buffer, true);
  reader.add("line 3\nline "_buffer);
  reader.add("4\n"_buffer, true);
  reader.add("line"_buffer);
  reader.add(" 5\nline "_buffer);
  reader.add("6\n"_buffer, true);

  std::string line;
  REQUIRE(reader.position() == 0);
  REQUIRE(reader.names(0, 0) == "'file1' at line 1");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1");
  REQUIRE(reader.position() == 1);
  REQUIRE(reader.names(0, 0) == "'file1' at line 1");
  REQUIRE(reader.names(0, 1) == "'file1' at line 1");
  REQUIRE(reader.names(1, 1) == "'file1' at line 2");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 2");
  REQUIRE(reader.position() == 2);
  REQUIRE(reader.names(1, 1) == "'file1' at line 2");
  REQUIRE(reader.names(1, 2) == "'file1' at line 2");
  REQUIRE(reader.names(2, 2) == "'file1' at line 3"); // not yet eof
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 3");
  REQUIRE(reader.position() == 3);
  REQUIRE(reader.names(2, 2) == "'file2' at line 1"); // now past eof
  REQUIRE(reader.names(2, 3) == "'file2' at line 1");
  REQUIRE(reader.names(3, 3) == "'file2' at line 2");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 4");
  REQUIRE(reader.position() == 4);
  REQUIRE(reader.names(3, 3) == "'file2' at line 2");
  REQUIRE(reader.names(3, 4) == "'file2' at line 2");
  REQUIRE(reader.names(4, 4) == "'file2' at line 3");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 5");
  REQUIRE(reader.position() == 5);
  REQUIRE(reader.names(4, 4) == "'file3' at line 1");
  REQUIRE(reader.names(4, 5) == "'file3' at line 1");
  REQUIRE(reader.names(5, 5) == "'file3' at line 2");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 6");
  REQUIRE(reader.position() == 6);
  REQUIRE(reader.names(5, 5) == "'file3' at line 2");
  REQUIRE(reader.names(5, 6) == "'file3' at line 2");
  REQUIRE(reader.names(6, 6) == "'file3' at line 3");
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(reader.names(6, 6) == "end of input files");
}

TEST_CASE("buffer_reader position across empty files")
{
  buffer_reader reader{ { "file1", "file2", "file3" } };
  reader.add("line 1\n"_buffer, true);
  reader.add(""_buffer, true);
  reader.add("line 2"_buffer, true);

  std::string line;
  REQUIRE(reader.position() == 0);
  REQUIRE(reader.names(0, 0) == "'file1' at line 1");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 1");
  REQUIRE(reader.position() == 1);
  REQUIRE(reader.names(0, 0) == "'file1' at line 1");
  REQUIRE(reader.names(0, 1) == "'file1' at line 1");
  REQUIRE(reader.names(1, 1) == "'file1' at line 2");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "line 2");
  REQUIRE(reader.position() == 2);
  REQUIRE(reader.names(1, 1) == "empty file 'file2', and 'file3' at line 1");
  REQUIRE(reader.names(1, 2) == "empty file 'file2', and 'file3' at line 1");
  REQUIRE(reader.names(2, 2) == "end of input files");
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(reader.position() == 2);
}

////////////////////////////////////////////////////////////////////////////////
// buffer_reader::release

TEST_CASE("buffer_reader releases consumed on partial read")
{
  std::vector<buffer> cache{};
  const auto release = [&cache](buffer&& buf) {
    cache.emplace_back(std::move(buf));
  };

  auto buf = "line "_buffer;
  const auto* data = buf.data();

  buffer_reader reader{ { "unnamed" } };
  reader.set_release_function(release);
  reader.add(std::move(buf));

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(cache.size() == 1);
  REQUIRE(cache.back().data() == data);
}

TEST_CASE("buffer_reader releases consumed on complete read")
{
  std::vector<buffer> cache{};
  const auto release = [&cache](buffer&& buf) {
    cache.emplace_back(std::move(buf));
  };

  auto buf_1 = "line "_buffer;
  const auto* data_1 = buf_1.data();
  auto buf_2 = "1\nline 2"_buffer;

  buffer_reader reader{ { "unnamed" } };
  reader.set_release_function(release);
  reader.add(std::move(buf_1));
  reader.add(std::move(buf_2));

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(cache.size() == 1);
  REQUIRE(cache.back().data() == data_1);
}

TEST_CASE("buffer_reader does not release buffers with zero capacity")
{
  buffer_reader reader{ { "unnamed" } };
  reader.set_release_function([&](buffer&&) { FAIL(); });
  reader.add(buffer());
  reader.add(buffer(), true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
}

TEST_CASE("buffer_reader does release buffers with capacity")
{
  buffer buf;
  buf.reserve(8);
  const auto* data = buf.data();
  bool released = false;

  buffer_reader reader{ { "unnamed" } };
  reader.set_release_function([&](buffer&& b) {
    REQUIRE(!released);
    REQUIRE(b.data() == data);
    released = true;
  });

  // VARIANTS: The buffer should be released even if it is not added to the
  //           internal cache, i.e. in the case of non-EOF empty buffers
  reader.add(std::move(buf), GENERATE(true, false));

  std::string line;
  REQUIRE(reader.getline(line) != buffer_reader::result::complete);
  REQUIRE(released);
}

TEST_CASE("buffer_reader releases consumed buffers in order")
{
  std::vector<buffer> cache{};
  const auto release = [&cache](buffer&& buf) {
    cache.emplace_back(std::move(buf));
  };

  auto buf_1 = "line 1\nline "_buffer;
  const auto* data_1 = buf_1.data();
  auto buf_2 = "2\nline 3"_buffer;
  const auto* data_2 = buf_2.data();

  buffer_reader reader{ { "unnamed" } };
  reader.set_release_function(release);
  reader.add(std::move(buf_1));
  reader.add(std::move(buf_2));

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(cache.empty());
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(cache.size() == 1);
  REQUIRE(cache.back().data() == data_1);
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(cache.size() == 2);
  REQUIRE(cache.back().data() == data_2);
}

TEST_CASE("buffer_reader releases unconsumed buffers with capacity")
{
  std::vector<buffer> cache{};
  const auto release = [&cache](buffer&& buf) {
    cache.emplace_back(std::move(buf));
  };

  auto buf_1 = "line 1\nline "_buffer;
  const auto* data_1 = buf_1.data();
  auto buf_2 = "2\nline 3"_buffer;
  const auto* data_2 = buf_2.data();

  {
    buffer_reader reader{ { "unnamed" } };
    reader.set_release_function(release);
    reader.add(std::move(buf_1));
    reader.add(std::move(buf_2));

    // EOF blocks are kept, even if empty, but should not be released
    if (GENERATE(true, false)) {
      reader.add(buffer{}, true);
    }
  }

  REQUIRE(cache.size() == 2);
  REQUIRE(cache.at(0).data() == data_1);
  REQUIRE(cache.at(1).data() == data_2);
}

TEST_CASE("buffer_reader releases unconsumed, empty buffers with capacity")
{
  std::vector<buffer> cache{};
  const auto release = [&cache](buffer&& buf) {
    cache.emplace_back(std::move(buf));
  };

  auto buf_1 = "line 1\nline 2"_buffer;
  const auto* data_1 = buf_1.data();
  buffer buf_2;
  buf_2.reserve(17);
  const auto* data_2 = buf_2.data();

  {
    buffer_reader reader{ { "unnamed" } };
    reader.set_release_function(release);
    reader.add(std::move(buf_1));
    reader.add(std::move(buf_2), true);
  }

  REQUIRE(cache.size() == 2);
  REQUIRE(cache.at(0).data() == data_1);
  REQUIRE(cache.at(1).data() == data_2);
}

TEST_CASE("buffer_reader remembers current line if release throws")
{
  bool threw = false;
  buffer_reader reader{ { "file1", "file2" } };
  reader.add("line 1"_buffer, true);
  reader.add("line 2\n"_buffer, true);

  reader.set_release_function([&threw](buffer&&) {
    if (!threw) {
      threw = true;
      throw fastq_error("some error");
    }
  });

  std::string line;
  REQUIRE(reader.position() == 0);
  REQUIRE_THROWS_AS(reader.getline(line), fastq_error);
  REQUIRE(line == "line 1");
  REQUIRE(reader.position() == 0);
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 1);
  REQUIRE(line == "line 1");
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(reader.position() == 2);
  REQUIRE(line == "line 2");
}

////////////////////////////////////////////////////////////////////////////////
// buffer_reader::close

TEST_CASE("buffer_reader close on empty reader results in EOF")
{
  buffer_reader reader{ { "file1" } };

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE_FALSE(reader.eof());
  reader.close();
  REQUIRE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
}

TEST_CASE("buffer_reader close on non-empty reader results in EOF")
{
  buffer_reader reader{ { "file1" } };
  reader.add("line 1\nline "_buffer);
  reader.add("2\nline 3"_buffer);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE_FALSE(reader.eof());
  REQUIRE(line == "line 1");
  reader.close();
  REQUIRE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
}

TEST_CASE("buffer_reader closing at EOF is no-op")
{
  buffer_reader reader{ { "file1" } };
  reader.add("test"_buffer, true);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::complete);
  REQUIRE(line == "test");
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(reader.eof());
  reader.close();
  REQUIRE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
}

TEST_CASE("buffer_reader adding to closed reader is not allowed")
{
  buffer_reader reader{ { "file1" } };
  REQUIRE_FALSE(reader.eof());
  reader.close();
  REQUIRE(reader.eof());

  REQUIRE_THROWS_AS(reader.add("test"_buffer, false), assert_failed);
  REQUIRE_THROWS_AS(reader.add("test"_buffer, true), assert_failed);
}

TEST_CASE("buffer_reader closing terminates partial reads")
{
  buffer_reader reader{ { "file1" } };
  reader.add("test"_buffer, false);

  std::string line;
  REQUIRE(reader.getline(line) == buffer_reader::result::partial);
  REQUIRE(line == "test");
  REQUIRE_FALSE(reader.eof());
  reader.close();
  REQUIRE(reader.eof());
  REQUIRE(reader.getline(line) == buffer_reader::result::eof);
  REQUIRE(line.empty());
}

TEST_CASE("buffer_reader close releases remaining buffers with capacity")
{
  buffer_reader reader{ { "file1" } };
  reader.add("line 1\nline "_buffer);
  reader.add("2\nline 3"_buffer);
  reader.add(buffer{}, GENERATE(true, false));

  // Partially or entirely unconsumed buffers
  if (GENERATE(true, false)) {
    std::string line;
    REQUIRE(reader.getline(line) == buffer_reader::result::complete);
    REQUIRE(line == "line 1");
  }

  std::vector<buffer> cache;
  reader.set_release_function(
    [&cache](buffer&& buf) { cache.emplace_back(std::move(buf)); });
  REQUIRE_NOTHROW(reader.close());
  REQUIRE(cache.size() == 2);
  REQUIRE(cache.at(0) == "line 1\nline "_buffer);
  REQUIRE(cache.at(1) == "2\nline 3"_buffer);
}

TEST_CASE("buffer_reader close releases empty buffers with capacity")
{
  buffer added_buf;
  added_buf.reserve(17);
  const auto* added_data = added_buf.data();

  buffer_reader reader{ { "file1" } };
  reader.add(std::move(added_buf), true);

  std::vector<buffer> cache;
  reader.set_release_function(
    [&cache](buffer&& buf) { cache.emplace_back(std::move(buf)); });
  REQUIRE_NOTHROW(reader.close());
  REQUIRE(cache.size() == 1);
  REQUIRE(cache.at(0).data() == added_data);
}

TEST_CASE("buffer_reader close releases nothing if empty")
{
  buffer_reader reader{ { "file1" } };

  if (GENERATE(true, false)) {
    reader.add("line 1\nline "_buffer);
    reader.add("2\nline 3"_buffer, GENERATE(true, false));

    std::string line;
    REQUIRE(reader.getline(line) == buffer_reader::result::complete);
    REQUIRE(reader.getline(line) == buffer_reader::result::complete);
    REQUIRE(reader.getline(line) != buffer_reader::result::eof);
  }

  reader.set_release_function([](buffer&&) { FAIL("should not release"); });
  REQUIRE_NOTHROW(reader.close());
}

TEST_CASE("buffer_reader close handles exceptions in release function")
{
  size_t nreleased = 0;
  std::vector<buffer> cache;
  {
    buffer_reader reader{ { "file1" } };
    reader.add("line 1\n"_buffer);
    reader.add("line 2\n"_buffer);
    reader.add("line 3\n"_buffer, true);

    reader.set_release_function([&](buffer&& buf) {
      REQUIRE(buf.capacity());
      if (nreleased++ == 1) {
        throw fastq_error("some error");
      }

      cache.emplace_back(std::move(buf));
    });

    REQUIRE_FALSE(reader.eof());
    REQUIRE(reader.buffer_count() == 3);
    REQUIRE_THROWS_AS(reader.close(), fastq_error);
    REQUIRE(reader.buffer_count() == 2); // left for destructor
    REQUIRE(reader.eof());
    REQUIRE(cache.size() == 1);
    REQUIRE(cache.at(0) == "line 1\n"_buffer);
  }

  REQUIRE(cache.size() == 3);
  REQUIRE(cache.at(1) == "line 2\n"_buffer);
  REQUIRE(cache.at(2) == "line 3\n"_buffer);
}

TEST_CASE("buffer_reader close handles exceptions after move in release")
{
  size_t nreleased = 0;
  std::vector<buffer> cache;
  {
    buffer_reader reader{ { "file1" } };
    reader.add("line 1\n"_buffer);
    reader.add("line 2\n"_buffer);
    reader.add("line 3\n"_buffer, true);

    reader.set_release_function([&](buffer&& buf) {
      REQUIRE(buf.capacity());
      if (nreleased++ == 1) {
        buffer b = std::move(buf);
        throw fastq_error("some error");
      }

      cache.emplace_back(std::move(buf));
    });

    REQUIRE_FALSE(reader.eof());
    REQUIRE(reader.buffer_count() == 3);
    REQUIRE_THROWS_AS(reader.close(), fastq_error);
    REQUIRE(reader.buffer_count() == 2); // left for destructor
    REQUIRE(reader.eof());
    REQUIRE(cache.size() == 1);
    REQUIRE(cache.at(0) == "line 1\n"_buffer);
  }

  REQUIRE(cache.size() == 2);
  REQUIRE(cache.at(1) == "line 3\n"_buffer);
}

////////////////////////////////////////////////////////////////////////////////
// buffer_reader::result

TEST_CASE("buffer_reader::result to string")
{
  using Catch::fallbackStringifier;

  CHECK(fallbackStringifier(buffer_reader::result::complete) ==
        "buffer_reader::result::complete");
  CHECK(fallbackStringifier(buffer_reader::result::partial) ==
        "buffer_reader::result::partial");
  CHECK(fallbackStringifier(buffer_reader::result::eof) ==
        "buffer_reader::result::eof");
}

} // namespace adapterremoval
