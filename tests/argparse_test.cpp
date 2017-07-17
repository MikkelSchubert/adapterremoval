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
#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "argparse.hpp"
#include "debug.hpp"

namespace ar
{

typedef std::unique_ptr<argparse::consumer_base> consumer_autoptr;


///////////////////////////////////////////////////////////////////////////////
// flag -- boolean

TEST(flag, defaults)
{
	consumer_autoptr ptr(new argparse::flag());
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(" ", ptr->metavar());
	ASSERT_EQ("", ptr->help());
	ASSERT_EQ("off", ptr->to_str());
}


TEST(flag, help)
{
	consumer_autoptr ptr(new argparse::flag(NULL, "help! help!"));
	ASSERT_EQ("help! help!", ptr->help());
}


TEST(flag, consumes_zero_arguments)
{
	string_vec arguments;
	arguments.push_back("--foo");
	consumer_autoptr ptr(new argparse::flag());
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(0, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
}


TEST(flag, consume_past_the_end)
{
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::flag());
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(0, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
}


TEST(flag, consume__with_sink__true)
{
	bool sink = true;
	consumer_autoptr ptr(new argparse::flag(&sink));
	ASSERT_EQ("on", ptr->to_str());
	const string_vec arguments;
	ASSERT_EQ(0, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_TRUE(sink);
	ASSERT_EQ("on", ptr->to_str());
}


TEST(flag, consume__with_sink__false)
{
	bool sink = false;
	consumer_autoptr ptr(new argparse::flag(&sink));
	ASSERT_EQ("off", ptr->to_str());
	const string_vec arguments;
	ASSERT_EQ(0, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_TRUE(sink);
	ASSERT_EQ("on", ptr->to_str());
}


///////////////////////////////////////////////////////////////////////////////
// any -- string

TEST(any, defaults)
{
	consumer_autoptr ptr(new argparse::any());
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("", ptr->metavar());
	ASSERT_EQ("", ptr->help());
	ASSERT_EQ("<not set>", ptr->to_str());
}


TEST(any, args_set)
{
	std::string sink = "kitchensink";
	consumer_autoptr ptr(new argparse::any(&sink, "a metavar", "help!"));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("a metavar", ptr->metavar());
	ASSERT_EQ("help!", ptr->help());
	ASSERT_EQ("kitchensink", ptr->to_str());
}


TEST(any, consumes_one_argument)
{
	string_vec arguments;
	arguments.push_back("foo");
	arguments.push_back("bar");
	consumer_autoptr ptr(new argparse::any());
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ("foo", ptr->to_str());
}


TEST(any, consume_past_the_end)
{
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::any());
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(static_cast<size_t>(-1), ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("<not set>", ptr->to_str());
}


TEST(any, consume__with_sink__empty)
{
	std::string sink;
	consumer_autoptr ptr(new argparse::any(&sink));
	ASSERT_EQ("<not set>", ptr->to_str());
	string_vec arguments;
	arguments.push_back("foo");
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ("foo", sink);
	ASSERT_EQ("foo", ptr->to_str());
}


TEST(any, consume__with_sink__preset)
{
	std::string sink = "kitchensink";
	consumer_autoptr ptr(new argparse::any(&sink));
	ASSERT_EQ("kitchensink", ptr->to_str());
	string_vec arguments;
	arguments.push_back("foo");
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ("foo", sink);
	ASSERT_EQ("foo", ptr->to_str());
}


///////////////////////////////////////////////////////////////////////////////
// many -- strings

TEST(many, defaults)
{
	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("", ptr->metavar());
	ASSERT_EQ("", ptr->help());
	ASSERT_EQ("<not set>", ptr->to_str());
}


TEST(many, args_set)
{
	string_vec sink;
	sink.push_back("kitchensink");

	consumer_autoptr ptr(new argparse::many(&sink, "a metavar", "help!"));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("a metavar", ptr->metavar());
	ASSERT_EQ("help!", ptr->help());
	ASSERT_EQ("kitchensink", ptr->to_str());
}


TEST(many, consumes_one_argument)
{
	string_vec arguments;
	arguments.push_back("foo");

	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ("foo", ptr->to_str());
}


TEST(many, consumes_two_arguments)
{
	string_vec arguments;
	arguments.push_back("foo");
	arguments.push_back("bar");

	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(2, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ("foo;bar", ptr->to_str());
}


TEST(many, consumes_until_next)
{
	string_vec arguments;
	arguments.push_back("foo");
	arguments.push_back("--zoo");
	arguments.push_back("bar");

	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ("foo", ptr->to_str());
}


TEST(many, consume_past_the_end)
{
	string_vec sink;
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::many(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(0, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ("<not set>", ptr->to_str());
}


TEST(many, consume__with_sink__empty)
{
	string_vec sink;
	string_vec expected;
	expected.push_back("foo");

	consumer_autoptr ptr(new argparse::many(&sink));
	ASSERT_EQ("<not set>", ptr->to_str());
	string_vec arguments;
	arguments.push_back("foo");
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ(expected, sink);
	ASSERT_EQ("foo", ptr->to_str());
}


TEST(many, consume__with_sink__preset)
{
	string_vec sink;
	sink.push_back("kitchensink");
	string_vec expected;
	expected.push_back("foo");

	consumer_autoptr ptr(new argparse::many(&sink));
	ASSERT_EQ("kitchensink", ptr->to_str());
	string_vec arguments;
	arguments.push_back("foo");
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ(expected, sink);
	ASSERT_EQ("foo", ptr->to_str());
}


///////////////////////////////////////////////////////////////////////////////
// knob -- unsigned

TEST(knob, defaults)
{
	unsigned sink = 0;
	consumer_autoptr ptr(new argparse::knob(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("", ptr->metavar());
	ASSERT_EQ("", ptr->help());
	ASSERT_EQ("0", ptr->to_str());
}


TEST(knob, defaults__sink_required)
{
	ASSERT_THROW(argparse::knob(NULL), assert_failed);
}


TEST(knob, args_set)
{
	unsigned sink = 7913;
	consumer_autoptr ptr(new argparse::knob(&sink, "a metavar", "help!"));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("a metavar", ptr->metavar());
	ASSERT_EQ("help!", ptr->help());
	ASSERT_EQ("7913", ptr->to_str());
}


TEST(knob, consumes_one_argument)
{
	unsigned sink = 0;
	string_vec arguments;
	arguments.push_back("47");
	arguments.push_back("bar");
	consumer_autoptr ptr(new argparse::knob(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ(47, sink);
	ASSERT_EQ("47", ptr->to_str());
}


TEST(knob, consume_past_the_end)
{
	unsigned sink = 13;
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(static_cast<size_t>(-1), ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(13, sink);
	ASSERT_EQ("13", ptr->to_str());
}


TEST(knob, reject_negative_values)
{
	unsigned sink = 13;
	string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	arguments.push_back("-47");
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(static_cast<size_t>(-1), ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(13, sink);
	ASSERT_EQ("13", ptr->to_str());
}


TEST(knob, accept_zero)
{
	unsigned sink = 13;
	string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	arguments.push_back("0");
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ(0, sink);
	ASSERT_EQ("0", ptr->to_str());
}


TEST(knob, upper_bound)
{
	unsigned sink = 13;
	string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	arguments.push_back("4294967295");
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ(4294967295, sink);
	ASSERT_EQ("4294967295", ptr->to_str());
}


TEST(knob, past_upper_bound)
{
    unsigned sink = 13;
    string_vec arguments;
    consumer_autoptr ptr(new argparse::knob(&sink));
    arguments.push_back("4294967296");
    ASSERT_FALSE(ptr->is_set());
    ASSERT_EQ(-1, ptr->consume(arguments.begin(), arguments.end()));
    ASSERT_FALSE(ptr->is_set());
}


TEST(knob, trailing_garbage)
{
    unsigned sink = 13;
    string_vec arguments;
    consumer_autoptr ptr(new argparse::knob(&sink));
    arguments.push_back("7913w");
    ASSERT_FALSE(ptr->is_set());
    ASSERT_EQ(-1, ptr->consume(arguments.begin(), arguments.end()));
    ASSERT_FALSE(ptr->is_set());
}


///////////////////////////////////////////////////////////////////////////////
// floaty_knob -- double

TEST(floaty_knob, defaults)
{
	double sink = 0;
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("", ptr->metavar());
	ASSERT_EQ("", ptr->help());
	ASSERT_EQ("0", ptr->to_str());
}


TEST(floaty_knob, default_nan)
{
	double sink = std::numeric_limits<double>::quiet_NaN();
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("", ptr->metavar());
	ASSERT_EQ("", ptr->help());
	ASSERT_EQ("<not set>", ptr->to_str());
}


TEST(floaty_knob, defaults__sink_required)
{
	ASSERT_THROW(argparse::floaty_knob(NULL), assert_failed);
}


TEST(floaty_knob, args_set)
{
	double sink = 3.142;
	consumer_autoptr ptr(new argparse::floaty_knob(&sink, "a metavar", "help!"));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ("a metavar", ptr->metavar());
	ASSERT_EQ("help!", ptr->help());
	ASSERT_EQ("3.142", ptr->to_str());
}


TEST(floaty_knob, consumes_one_argument)
{
	double sink = 47.0;
	string_vec arguments;
	arguments.push_back("-19.84");
	arguments.push_back("bar");
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(1, ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_TRUE(ptr->is_set());
	ASSERT_EQ(-19.84, sink);
	ASSERT_EQ("-19.84", ptr->to_str());
}


TEST(floaty_knob, consume_past_the_end)
{
	double sink = 13;
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(static_cast<size_t>(-1), ptr->consume(arguments.begin(), arguments.end()));
	ASSERT_FALSE(ptr->is_set());
	ASSERT_EQ(13, sink);
	ASSERT_EQ("13", ptr->to_str());
}


TEST(floaty_knob, trailing_garbage)
{
    double sink = 47.0;
    string_vec arguments;
    arguments.push_back("-19.84wat");
    consumer_autoptr ptr(new argparse::floaty_knob(&sink));
    ASSERT_FALSE(ptr->is_set());
    ASSERT_EQ(-1, ptr->consume(arguments.begin(), arguments.end()));
    ASSERT_FALSE(ptr->is_set());
}


///////////////////////////////////////////////////////////////////////////////
// parser

} // namespace ar
