// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"            // for assert_failed
#include "linereader_joined.hpp" // declarations
#include "testing.hpp"           // for TEST_CASE, REQUIRE, ...
#include <string>                // for string
#include <utility>               // for move
#include <vector>                // for vector

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// joined_filenames basics

TEST_CASE("joined_filenames requires non-empty list")
{
  REQUIRE_THROWS_AS(joined_filenames({}), assert_failed);
}

TEST_CASE("joined_filenames with one filename")
{
  joined_filenames f({ "reads.fastq" });

  CHECK(f.remaining_filenames() == 1);
  CHECK(f.current_filename() == "reads.fastq");

  f.next_file();

  CHECK(f.remaining_filenames() == 0);
  REQUIRE_THROWS_AS(f.current_filename(), assert_failed);
  REQUIRE_THROWS_AS(f.next_file(), assert_failed);
}

TEST_CASE("joined_filenames with multiple filenames")
{
  joined_filenames f({ "file1.fastq", "file2.fastq", "file3.fastq" });

  CHECK(f.remaining_filenames() == 3);
  CHECK(f.current_filename() == "file1.fastq");
  f.next_file();

  CHECK(f.remaining_filenames() == 2);
  CHECK(f.current_filename() == "file2.fastq");
  f.next_file();

  CHECK(f.remaining_filenames() == 1);
  CHECK(f.current_filename() == "file3.fastq");
  f.next_file();

  CHECK(f.remaining_filenames() == 0);
  REQUIRE_THROWS_AS(f.current_filename(), assert_failed);
  REQUIRE_THROWS_AS(f.next_file(), assert_failed);
}

TEST_CASE("joined_filenames line numbers increase globally")
{
  joined_filenames f({ "file1.fastq", "file2.fastq" });

  CHECK(f.remaining_filenames() == 2);
  CHECK(f.current_filename() == "file1.fastq");
  f.inc_position();
  CHECK(f.position() == 1);
  f.next_file();
  CHECK(f.position() == 1);

  CHECK(f.remaining_filenames() == 1);
  CHECK(f.current_filename() == "file2.fastq");
  f.inc_position(17);
  CHECK(f.position() == 18);
  f.next_file();

  CHECK(f.position() == 18);
  CHECK(f.remaining_filenames() == 0);
  REQUIRE_THROWS_AS(f.current_filename(), assert_failed);
}

////////////////////////////////////////////////////////////////////////////////
// joined_filenames filenames

TEST_CASE("joined_filenames start must be <= than end")
{
  joined_filenames jf{ { "file1" } };
  jf.inc_position(10);

  REQUIRE_NOTHROW(jf.filenames(0, 0));
  REQUIRE_NOTHROW(jf.filenames(0, 10));
  REQUIRE_NOTHROW(jf.filenames(9, 10));
  REQUIRE_THROWS_AS(jf.filenames(10, 9), assert_failed);
  REQUIRE_THROWS_AS(jf.filenames(10, 0), assert_failed);
}

TEST_CASE("joined_filenames end must be <= position")
{
  joined_filenames jf{ { "file1" } };
  REQUIRE_NOTHROW(jf.filenames(0, 0));
  REQUIRE_THROWS_AS(jf.filenames(0, 1), assert_failed);
  REQUIRE_THROWS_AS(jf.filenames(1, 2), assert_failed);

  jf.inc_position(10);
  REQUIRE_NOTHROW(jf.filenames(0, 10));
  REQUIRE_NOTHROW(jf.filenames(9, 10));
  REQUIRE_THROWS_AS(jf.filenames(10, 11), assert_failed);
}

TEST_CASE("joined_filenames end must be <= position, multiple files")
{
  joined_filenames jf{ { "file1", "file2" } };
  REQUIRE_NOTHROW(jf.filenames(0, 0));
  REQUIRE_THROWS_AS(jf.filenames(0, 1), assert_failed);
  REQUIRE_THROWS_AS(jf.filenames(1, 2), assert_failed);

  jf.inc_position(10);

  REQUIRE_NOTHROW(jf.filenames(0, 10));
  REQUIRE_NOTHROW(jf.filenames(9, 10));
  REQUIRE_THROWS_AS(jf.filenames(10, 11), assert_failed);

  jf.next_file();

  REQUIRE_NOTHROW(jf.filenames(0, 10));
  REQUIRE_NOTHROW(jf.filenames(9, 10));
  REQUIRE_THROWS_AS(jf.filenames(10, 11), assert_failed);

  jf.inc_position(1);

  REQUIRE_NOTHROW(jf.filenames(10, 11));
}

TEST_CASE("joined_filenames no closed files")
{
  std::vector<std::string> filenames{ "file1" };
  if (GENERATE(false, true)) {
    filenames.emplace_back("file2");
  }

  joined_filenames jf{ std::move(filenames) };
  jf.inc_position(1000000);

  CHECK(jf.filenames(0, 0) == "'file1' at line 1");
  CHECK(jf.filenames(0, 10) == "'file1' at lines 1-10");
  CHECK(jf.filenames(1, 10) == "'file1' at lines 2-10");
  CHECK(jf.filenames(10, 10) == "'file1' at line 11");
  CHECK(jf.filenames(10, 1000000) == "'file1' at lines 11-1000000");
}

TEST_CASE("joined_filenames single file")
{
  joined_filenames jf{ { "file1" } };
  jf.inc_position(1234);

  if (GENERATE(true, false)) {
    jf.next_file();
  }

  CHECK(jf.filenames(0, 0) == "'file1' at line 1");
  CHECK(jf.filenames(0, 1234) == "'file1' at lines 1-1234");
}

TEST_CASE("joined_filenames multiple files")
{
  joined_filenames jf{ { "file1", "file2" } };
  jf.inc_position(1000);
  jf.next_file();
  jf.inc_position(1234);

  if (GENERATE(true, false)) {
    jf.next_file();
  }

  CHECK(jf.filenames(0, 0) == "'file1' at line 1");
  CHECK(jf.filenames(0, 1) == "'file1' at line 1");
  CHECK(jf.filenames(0, 750) == "'file1' at lines 1-750");
  CHECK(jf.filenames(749, 1000) ==
        "'file1' at lines 750-1000, and 'file2' at line 1");
  CHECK(jf.filenames(1000, 1000) == "'file2' at line 1");
  CHECK(jf.filenames(749, 1234) ==
        "'file1' at lines 750-1000, and 'file2' at lines 1-234");
  CHECK(jf.filenames(1001, 1234) == "'file2' at lines 2-234");
}

TEST_CASE("joined_filenames multiple files, one empty")
{
  joined_filenames jf{ { "file1", "file2", "file3" } };
  jf.inc_position(1000);
  jf.next_file();
  jf.next_file();
  jf.inc_position(1234);

  if (GENERATE(true, false)) {
    jf.next_file();
  }

  CHECK(jf.filenames(0, 0) == "'file1' at line 1");
  CHECK(jf.filenames(0, 750) == "'file1' at lines 1-750");
  CHECK(jf.filenames(749, 1000) == "'file1' at lines 750-1000, empty file "
                                   "'file2', and 'file3' at line 1");
  CHECK(jf.filenames(1000, 1000) ==
        "empty file 'file2', and 'file3' at line 1");
  CHECK(jf.filenames(749, 1234) == "'file1' at lines 750-1000, empty file "
                                   "'file2', and 'file3' at lines 1-234");
  CHECK(jf.filenames(1001, 1234) == "'file3' at lines 2-234");
}

TEST_CASE("joined_filenames multiple files, one empty at front")
{
  joined_filenames jf{ { "file1", "file2" } };
  jf.next_file();
  jf.inc_position(1000);

  if (GENERATE(true, false)) {
    jf.next_file();
  }

  CHECK(jf.filenames(0, 0) == "empty file 'file1', and 'file2' at line 1");
  CHECK(jf.filenames(0, 1) == "empty file 'file1', and 'file2' at line 1");
  CHECK(jf.filenames(0, 750) == "empty file 'file1', and 'file2' at lines "
                                "1-750");
  CHECK(jf.filenames(1, 750) == "'file2' at lines 2-750");
}

TEST_CASE("joined_filenames multiple files, one empty at back")
{
  joined_filenames jf{ { "file1", "file2" } };
  jf.inc_position(1000);
  jf.next_file();
  jf.next_file();

  CHECK(jf.filenames(999, 1000) ==
        "'file1' at line 1000, and empty file 'file2'");
}

TEST_CASE("joined_filenames at end of files")
{
  std::vector<std::string> filenames{ "file1" };
  joined_filenames jf{ std::move(filenames) };
  jf.inc_position(100);

  CHECK(jf.filenames(99, 100) == "'file1' at line 100");
  jf.next_file();
  CHECK(jf.filenames(100, 100) == "end of input files");
}

TEST_CASE("joined_filenames at end of files with empty file")
{
  std::vector<std::string> filenames{ "file1", "file2" };
  joined_filenames jf{ std::move(filenames) };
  jf.inc_position(100);
  jf.next_file();
  jf.next_file();

  CHECK(jf.filenames(100, 100) == "end of input files");
}

TEST_CASE("joined_filenames at end of files with empty file only")
{
  std::vector<std::string> filenames{ "file1" };
  joined_filenames jf{ std::move(filenames) };
  jf.next_file();

  CHECK(jf.filenames(0, 0) == "end of input files");
}

} // namespace adapterremoval
