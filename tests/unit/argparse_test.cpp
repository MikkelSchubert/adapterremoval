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

#include "testing.hpp"
#include "argparse.hpp"
#include "debug.hpp"

namespace ar
{

typedef std::unique_ptr<argparse::consumer_base> consumer_autoptr;


///////////////////////////////////////////////////////////////////////////////
// flag -- boolean

TEST_CASE("Flag defaults", "[argparse::flag]")
{
	consumer_autoptr ptr(new argparse::flag());
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == " ");
	REQUIRE(ptr->help() == "");
	REQUIRE(ptr->to_str() == "off");
}


TEST_CASE("Flag help", "[argparse::flag]")
{
	consumer_autoptr ptr(new argparse::flag(nullptr, "help! help!"));
	REQUIRE(ptr->help() == "help! help!");
}


TEST_CASE("Flag consumes no arguments", "[argparse::flag]")
{
	string_vec arguments;
	arguments.push_back("--foo");
	consumer_autoptr ptr(new argparse::flag());
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 0);
	CHECK(ptr->is_set());
}


TEST_CASE("Flag may be called on end of arguments", "[argparse::flag]")
{
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::flag());
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 0);
	CHECK(ptr->is_set());
}


TEST_CASE("Flag uses sink with true", "[argparse::flag]")
{
	bool sink = true;
	consumer_autoptr ptr(new argparse::flag(&sink));
	REQUIRE(ptr->to_str() == "on");
	const string_vec arguments;
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 0);
	CHECK(ptr->is_set());
	CHECK(sink);
	REQUIRE(ptr->to_str() == "on");
}


TEST_CASE("Flag uses sink with false", "[argparse::flag]")
{
	bool sink = false;
	consumer_autoptr ptr(new argparse::flag(&sink));
	REQUIRE(ptr->to_str() == "off");
	const string_vec arguments;
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 0);
	CHECK(ptr->is_set());
	CHECK(sink);
	REQUIRE(ptr->to_str() == "on");
}


///////////////////////////////////////////////////////////////////////////////
// any -- string

TEST_CASE("Any defaults", "[argparse::any]")
{
	consumer_autoptr ptr(new argparse::any());
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "");
	REQUIRE(ptr->help() == "");
	REQUIRE(ptr->to_str() == "<not set>");
}


TEST_CASE("Any value set", "[argparse::any]")
{
	std::string sink = "kitchensink";
	consumer_autoptr ptr(new argparse::any(&sink, "a metavar", "help!"));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "a metavar");
	REQUIRE(ptr->help() == "help!");
	REQUIRE(ptr->to_str() == "kitchensink");
}


TEST_CASE("Any consumes one argument", "[argparse::any]")
{
	string_vec arguments;
	arguments.push_back("foo");
	arguments.push_back("bar");
	consumer_autoptr ptr(new argparse::any());
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(ptr->to_str() == "foo");
}


TEST_CASE("Any cannot consume empty list of arguments", "[argparse::any]")
{
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::any());
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == static_cast<size_t>(-1));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->to_str() == "<not set>");
}


TEST_CASE("Any with empty sink", "[argparse::any]")
{
	std::string sink;
	consumer_autoptr ptr(new argparse::any(&sink));
	REQUIRE(ptr->to_str() == "<not set>");
	string_vec arguments;
	arguments.push_back("foo");
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == "foo");
	REQUIRE(ptr->to_str() == "foo");
}


TEST_CASE("Any with preset sink", "[argparse::any]")
{
	std::string sink = "kitchensink";
	consumer_autoptr ptr(new argparse::any(&sink));
	REQUIRE(ptr->to_str() == "kitchensink");
	string_vec arguments;
	arguments.push_back("foo");
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == "foo");
	REQUIRE(ptr->to_str() == "foo");
}


///////////////////////////////////////////////////////////////////////////////
// many -- strings

TEST_CASE("Many defaults", "[argparse::many]")
{
	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "");
	REQUIRE(ptr->help() == "");
	REQUIRE(ptr->to_str() == "<not set>");
}


TEST_CASE("Many with defaults", "[argparse::many]")
{
	string_vec sink;
	sink.push_back("kitchensink");

	consumer_autoptr ptr(new argparse::many(&sink, "a metavar", "help!"));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "a metavar");
	REQUIRE(ptr->help() == "help!");
	REQUIRE(ptr->to_str() == "kitchensink");
}


TEST_CASE("Many consumes one argument", "[argparse::many]")
{
	string_vec arguments;
	arguments.push_back("foo");

	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(ptr->to_str() == "foo");
}


TEST_CASE("Many consumes two arguments", "[argparse::many]")
{
	string_vec arguments;
	arguments.push_back("foo");
	arguments.push_back("bar");

	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 2);
	CHECK(ptr->is_set());
	REQUIRE(ptr->to_str() == "foo;bar");
}


TEST_CASE("Many consumes until next option", "[argparse::many]")
{
	string_vec arguments;
	arguments.push_back("foo");
	arguments.push_back("--zoo");
	arguments.push_back("bar");

	string_vec sink;
	consumer_autoptr ptr(new argparse::many(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(ptr->to_str() == "foo");
}


TEST_CASE("Many does not consume empty list of arguments", "[argparse::many]")
{
	string_vec sink;
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::many(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 0);
	CHECK(ptr->is_set());
	REQUIRE(ptr->to_str() == "<not set>");
}


TEST_CASE("Many with empty sink", "[argparse::many]")
{
	string_vec sink;
	string_vec expected;
	expected.push_back("foo");

	consumer_autoptr ptr(new argparse::many(&sink));
	REQUIRE(ptr->to_str() == "<not set>");
	string_vec arguments;
	arguments.push_back("foo");
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == expected);
	REQUIRE(ptr->to_str() == "foo");
}


TEST_CASE("many with preset sink", "[argparse::many]")
{
	string_vec sink;
	sink.push_back("kitchensink");
	string_vec expected;
	expected.push_back("foo");

	consumer_autoptr ptr(new argparse::many(&sink));
	REQUIRE(ptr->to_str() == "kitchensink");
	string_vec arguments;
	arguments.push_back("foo");
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == expected);
	REQUIRE(ptr->to_str() == "foo");
}


///////////////////////////////////////////////////////////////////////////////
// knob -- unsigned

TEST_CASE("Knob defaults", "[argparse::knob]")
{
	unsigned sink = 0;
	consumer_autoptr ptr(new argparse::knob(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "");
	REQUIRE(ptr->help() == "");
	REQUIRE(ptr->to_str() == "0");
}


TEST_CASE("Knob requires a sink", "[argparse::knob]")
{
	REQUIRE_THROWS_AS(argparse::knob(nullptr), assert_failed);
}


TEST_CASE("Knot values are set", "[argparse::knob]")
{
	unsigned sink = 7913;
	consumer_autoptr ptr(new argparse::knob(&sink, "a metavar", "help!"));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "a metavar");
	REQUIRE(ptr->help() == "help!");
	REQUIRE(ptr->to_str() == "7913");
}


TEST_CASE("Knob consumes one argument", "[argparse::knob]")
{
	unsigned sink = 0;
	string_vec arguments;
	arguments.push_back("47");
	arguments.push_back("bar");
	consumer_autoptr ptr(new argparse::knob(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == 47);
	REQUIRE(ptr->to_str() == "47");
}


TEST_CASE("Knob does not consume empty list of arguments", "[argparse::knob]")
{
	unsigned sink = 13;
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == static_cast<size_t>(-1));
	CHECK(!ptr->is_set());
	REQUIRE(sink == 13);
	REQUIRE(ptr->to_str() == "13");
}


TEST_CASE("Knob rejects negative values", "[argparse::knob]")
{
	unsigned sink = 13;
	string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	arguments.push_back("-47");
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == static_cast<size_t>(-1));
	CHECK(!ptr->is_set());
	REQUIRE(sink == 13);
	REQUIRE(ptr->to_str() == "13");
}


TEST_CASE("Knob accepts zero", "[argparse::knob]")
{
	unsigned sink = 13;
	string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	arguments.push_back("0");
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == 0);
	REQUIRE(ptr->to_str() == "0");
}


TEST_CASE("Knob accepts unsigned upper bound", "[argparse::knob]")
{
	unsigned sink = 13;
	string_vec arguments;
	consumer_autoptr ptr(new argparse::knob(&sink));
	arguments.push_back("4294967295");
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == 4294967295);
	REQUIRE(ptr->to_str() == "4294967295");
}


TEST_CASE("Knob rejects past unsigned upper bound", "[argparse::knob]")
{
    unsigned sink = 13;
    string_vec arguments;
    consumer_autoptr ptr(new argparse::knob(&sink));
    arguments.push_back("4294967296");
    CHECK(!ptr->is_set());
    REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == static_cast<size_t>(-1));
    CHECK(!ptr->is_set());
}


TEST_CASE("Knob rejects trailing garbage", "[argparse::knob]")
{
    unsigned sink = 13;
    string_vec arguments;
    consumer_autoptr ptr(new argparse::knob(&sink));
    arguments.push_back("7913w");
    CHECK(!ptr->is_set());
    REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == static_cast<size_t>(-1));
    CHECK(!ptr->is_set());
}


///////////////////////////////////////////////////////////////////////////////
// floaty_knob -- double

TEST_CASE("Floaty knob defaults", "[argparse::floaty_knob]")
{
	double sink = 0;
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "");
	REQUIRE(ptr->help() == "");
	REQUIRE(ptr->to_str() == "0");
}


TEST_CASE("Floaty knob treats NaN as unset", "[argparse::floaty_knob]")
{
	double sink = std::numeric_limits<double>::quiet_NaN();
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "");
	REQUIRE(ptr->help() == "");
	REQUIRE(ptr->to_str() == "<not set>");
}


TEST_CASE("Floaty knob requires sink", "[argparse::floaty_knob]")
{
	REQUIRE_THROWS_AS(argparse::floaty_knob(nullptr), assert_failed);
}


TEST_CASE("Floaty knob with values", "[argparse::floaty_knob]")
{
	double sink = 3.142;
	consumer_autoptr ptr(new argparse::floaty_knob(&sink, "a metavar", "help!"));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->metavar() == "a metavar");
	REQUIRE(ptr->help() == "help!");
	REQUIRE(ptr->to_str() == "3.142");
}


TEST_CASE("Floaty knob consumes one argument", "[argparse::floaty_knob]")
{
	double sink = 47.0;
	string_vec arguments;
	arguments.push_back("-19.84");
	arguments.push_back("bar");
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == 1);
	CHECK(ptr->is_set());
	REQUIRE(sink == -19.84);
	REQUIRE(ptr->to_str() == "-19.84");
}


TEST_CASE("Floaty knob does not consume empty list of arguments", "[argparse::floaty_knob]")
{
	double sink = 13;
	const string_vec arguments;
	consumer_autoptr ptr(new argparse::floaty_knob(&sink));
	CHECK(!ptr->is_set());
	REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == static_cast<size_t>(-1));
	CHECK(!ptr->is_set());
	REQUIRE(sink == 13);
	REQUIRE(ptr->to_str() == "13");
}


TEST_CASE("Floaty knob rejects trailing garbage", "[argparse::floaty_knob]")
{
    double sink = 47.0;
    string_vec arguments;
    arguments.push_back("-19.84wat");
    consumer_autoptr ptr(new argparse::floaty_knob(&sink));
    CHECK(!ptr->is_set());
    REQUIRE(ptr->consume(arguments.begin(), arguments.end()) == static_cast<size_t>(-1));
    CHECK(!ptr->is_set());
}


///////////////////////////////////////////////////////////////////////////////
// parser

} // namespace ar
