/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <limits>
#include <stdexcept>
#include <gtest/gtest.h>

#include "strutils.h"

namespace ar
{

///////////////////////////////////////////////////////////////////////////////
// Tests for 'indent_lines'

TEST(strutils_indent, empty_input_empty_output)
{
    ASSERT_EQ("", indent_lines(""));
}


TEST(strutils_indent, empty_line)
{
    ASSERT_EQ("\n", indent_lines("\n"));
}


TEST(strutils_indent, two_empty_lines)
{
    ASSERT_EQ("\n\n", indent_lines("\n\n"));
}


TEST(strutils_indent, three_empty_lines)
{
    ASSERT_EQ("\n\n\n", indent_lines("\n\n\n"));
}


TEST(strutils_indent, single_line_no_trailing_endline)
{
    ASSERT_EQ("    this is a test", indent_lines("this is a test"));
}


TEST(strutils_indent, single_line_trailing_endline)
{
    ASSERT_EQ("    this is a test\n", indent_lines("this is a test\n"));
}


TEST(strutils_indent, two_lines_no_trailing_endline)
{
    ASSERT_EQ("    this is a test\n    line #2",
              indent_lines("this is a test\nline #2"));
}


TEST(strutils_indent, two_lines_trailing_endline)
{
    ASSERT_EQ("    this is a test\n    line #2\n",
              indent_lines("this is a test\nline #2\n"));
}


TEST(strutils_indent, empty_line_middle)
{
    ASSERT_EQ("    this is a test\n\n    line #2\n",
              indent_lines("this is a test\n\nline #2\n"));
}


TEST(strutils_indent, empty_line_trailing)
{
    ASSERT_EQ("    this is a test\n    line #2\n\n",
              indent_lines("this is a test\nline #2\n\n"));
}


TEST(strutils_indent, empty_lines_trailing)
{
    ASSERT_EQ("    this is a test\n    line #2\n\n\n",
              indent_lines("this is a test\nline #2\n\n\n"));
}


///////////////////////////////////////////////////////////////////////////////
// Tests for 'columnize_'


///////////////////////////////////////////////////////////////////////////////
// Tests for 'indent_lines'

TEST(strutils_columnize, empty_lines)
{
    ASSERT_EQ("", columnize_text(""));
    ASSERT_EQ("", columnize_text("\n"));
    ASSERT_EQ("", columnize_text("\n\n"));
    ASSERT_EQ("", columnize_text("\n\n\n"));
}


TEST(strutils_columnize, text_within_empty_lines)
{
    ASSERT_EQ("foo bar", columnize_text("foo bar"));
    ASSERT_EQ("foo bar", columnize_text("foo\nbar"));
    ASSERT_EQ("foo bar", columnize_text("foo\nbar\n"));
    ASSERT_EQ("foo bar", columnize_text("\nfoo\nbar"));
    ASSERT_EQ("foo bar", columnize_text("foo\nbar\n\n"));
    ASSERT_EQ("foo bar", columnize_text("\nfoo\nbar\n"));
    ASSERT_EQ("foo bar", columnize_text("\n\nfoo\nbar\n"));
}


TEST(strutils_columnize, maximum_width)
{
    ASSERT_EQ("foo bar zood", columnize_text("foo bar\nzood", 12));
    ASSERT_EQ("foo bar\nzood", columnize_text("foo bar\nzood", 11));
    ASSERT_EQ("foo bar\nzood", columnize_text("foo bar\nzood", 7));
    ASSERT_EQ("foo\nbar\nzood", columnize_text("foo bar\nzood", 6));
    ASSERT_EQ("foo\nbar\nzood", columnize_text("foo bar\nzood", 3));
    ASSERT_EQ("foo\nbar\nzood", columnize_text("foo bar\nzood", 1));
    ASSERT_EQ("foo\nbar\nzood", columnize_text("foo bar\nzood", 0));
}


TEST(strutils_columnize, ljust)
{
    ASSERT_EQ("foo bar zood", columnize_text("foo bar\nzood", 12, 2));
    ASSERT_EQ("foo bar\n  zood", columnize_text("foo bar\nzood", 11, 2));
    ASSERT_EQ("foo bar\n  zood", columnize_text("foo bar\nzood", 7, 2));
    ASSERT_EQ("foo\n  bar\n  zood", columnize_text("foo bar\nzood", 6, 2));
    ASSERT_EQ("foo\n  bar\n  zood", columnize_text("foo bar\nzood", 3, 2));
    ASSERT_EQ("foo\n  bar\n  zood", columnize_text("foo bar\nzood", 1, 2));
    ASSERT_EQ("foo\n  bar\n  zood", columnize_text("foo bar\nzood", 0, 2));
}

} // namespace ar
