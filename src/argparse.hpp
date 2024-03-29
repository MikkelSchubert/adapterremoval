/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#pragma once

#include <algorithm> // for max
#include <map>       // for map
#include <memory>    // for shared_ptr, unique_ptr
#include <stddef.h>  // for size_t
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector

#include "commontypes.hpp" // for string_vec_citer, string_vec

namespace adapterremoval {

namespace argparse {

class argument;
class argument_parser;

class sink;
class bool_sink;
class double_sink;
class uint_sink;
class str_sink;
class vec_sink;

using argument_ptr = std::shared_ptr<argument>;

//! Parse results for command-line arguments
enum class parse_result
{
  //! Terminate now (e.g. --version or --help used)
  exit,
  //! Error occurred parsing arguments / invalid combination of args
  error,
  //! No errors parsing commandline arguments
  ok
};

/**
 * Simple type-safe parsing of command-line options.
 *
 * To handle a an argument, an object of type argument is assigned
 * to the argparse::parser using the [] operator. For example, the parse and
 * save an integer value to a variable when the user supplies the argument
 * '--example', the following is done:
 *
 * int target = 0; // Must be initialized!
 * argparse::parser argparser(...);
 * argparser
 *  .add_knob("--example", &target)
 *  .metavar("N")
 *  .help("the number of examples");
 * argparser.parse_args(argc, argv);
 *
 * Aliases can be created for command-line arguments simply by assigning the
 * same parser to multiple keys; for example, to alias '--example' as '-e':
 *
 * argparser["-e"] = argparser["--example"];
 *
 * The class automatically handles the following arguments:
 *  --help, --h, -help, -h: Displays the program help
 *  --version, -v: Displays the program name and version string
 *
 * In both cases, the parse_args function returns false.
 *
 * The assigned consumer_ptrs are owned by and freed by the argparse::parser
 * upon destruction of the object. Pointers assigned multiple times (i.e. when
 * used with aliases) are only freed once.
 */
class parser
{
public:
  parser();

  /** Sets the name used in --help and --version messages */
  void set_name(const std::string& name);
  /** Sets the version string used in --help and --version messages */
  void set_version(const std::string& version);
  /** Sets the preamble text used in --help */
  void set_preamble(const std::string& text);
  /** Sets the license text used in --licenses */
  void set_licenses(const std::string& text);

  /** Parses a set of command-line options as passed to main(argc, argv). */
  parse_result parse_args(int argc, char const* const* argv);

  /** Returns true if the option with the given key has been set. */
  bool is_set(const std::string& key) const;
  /** Returns the value associated with the argument as a string. */
  std::string to_str(const std::string& key) const;

  /** Add argument with metavar. By default this takes no values. */
  argument& add(const std::string& name, const std::string& metavar = "");

  /** Add a blank line between the previous and the next command. */
  void add_separator();

  /** Add a blank line and a header between the previous and next command. */
  void add_header(const std::string& header);

  /** Helper function; prints the program name and version string. */
  void print_version() const;
  /** Helper functions; prints the full set of help-text. */
  void print_help() const;
  /** Helper functions; formats and prints the licenses. */
  void print_licenses() const;

  /** Set the maximum terminal width. */
  void set_terminal_width(unsigned w);

  //! Copy construction not supported
  parser(const parser&) = delete;
  //! Assignment not supported
  parser& operator=(const parser&) = delete;

private:
  void update_argument_map();

  argument_ptr find_argument(const std::string& key);

  struct argument_entry
  {
    std::string header;
    argument_ptr argument;
  };

  std::vector<argument_entry> m_args;
  std::map<std::string, argument_ptr, std::less<>> m_keys;

  //! Name of the program
  std::string m_name;
  //! Version string for the program (excluding the name)
  std::string m_version;
  //! Preamble text for the program.
  std::string m_preamble;
  //! Licenses for the program.
  std::string m_licenses;
  //! Maximum terminal width used for printing help messages
  unsigned m_terminal_width;
};

/**
 * Base class for arguments;
 *
 * Each consumer must implement the consume function, which takes iterators to
 * the arguments following the key for this parser (i.e. not including the
 * --option). These then consume zero or more values, returning the number
 * thus consumed, or (size_t)-1 if the values were missing or invalid.
 */
class argument
{
public:
  struct argument_key
  {
    std::string name;
    bool deprecated;
  };

  argument(const std::string& key, const std::string& metavar = "");

  /** Returns true if the consumer has consumed a value. */
  bool is_set() const;
  /** Returns true if the argument is deprecated and hidden. */
  bool is_deprecated() const;

  /** Returns the canonical argument key. */
  const std::string& key() const;
  /** Returns the short argument key; may be an empty string. */
  const std::string& short_key() const;
  /** Returns long, short, and deprecated argument keys. */
  string_vec keys() const;
  /** Returns true if this key is a deprecated alias for this argument. */
  bool is_deprecated_alias(const std::string& key) const;
  /** Returns the metavariable. May be an empty string. */
  const std::string& metavar() const;
  /** Returns help string with %default replaced with the current value. */
  std::string help() const;

  /** Indicates the minimum number of values taken by this argument */
  size_t min_values() const;
  /** Indicates the maximum number of values taken by this argument */
  size_t max_values() const;

  /** Options that MUST be specified along with this argument. */
  const string_vec& depends_on() const;
  /** Options that must NOT be specified along with this argument. */
  const string_vec& conflicts_with() const;

  /** Returns the value associated with the argument as a string. */
  std::string to_str() const;

  /** Set the metavar for this argument. */
  argument& metavar(const std::string& metavar);
  /** Set help string for this argument. */
  argument& help(const std::string& alias);
  /** Create a short form of the argument. */
  argument& abbreviation(char key);
  /** Create deprecated alias for the argument. */
  argument& deprecated_alias(const std::string& alias);
  /** The argument is deprecated / not to be printed by -h/--help. */
  argument& deprecated();

  /** Option `key` MUST be specified along with this argument. */
  argument& depends_on(const std::string& key);
  /** Option `key` must NOT be specified along with this argument. */
  argument& conflicts_with(const std::string& key);

  bool_sink& bind_bool(bool* sink);
  uint_sink& bind_uint(unsigned* sink);
  double_sink& bind_double(double* sink);
  str_sink& bind_str(std::string* sink);
  vec_sink& bind_vec(string_vec* sink);

  /** Parse the next arguments, returning the number of items parsed or -1. */
  size_t parse(string_vec_citer start, const string_vec_citer& end);

  //! Copy construction not supported
  argument(const argument&) = delete;
  //! Assignment not supported
  argument& operator=(const argument&) = delete;

private:
  //! Number of times the argument has been specified
  unsigned m_times_set;
  //! Default sink value
  bool m_default_sink;
  //! Indicates if the argument is deprecated / hidden.
  bool m_deprecated;
  //! Deprecated keys (long and short) for this argument
  string_vec m_deprecated_keys;

  //! The long, canonical argument key
  std::string m_key_long;
  //! An optional, short argument key
  std::string m_key_short;

  //! Optional metavar (defaults to uppercase `m_name` without dashes)
  std::string m_metavar;
  //! Help string; the string '%default' will be replaced with the current value
  std::string m_help;

  //! This argument must be specified along with these arguments.
  string_vec m_depends_on;
  //! This argument cannot be specified along with these arguments.
  string_vec m_conflicts_with;

  std::unique_ptr<sink> m_sink;
};

class sink
{
public:
  /** Creates a sink that takes exactly `n_values` values */
  explicit sink(size_t n_values);
  /** Creates a sink that takes between min and max values (incl.) */
  sink(size_t min_values, size_t max_values);

  virtual ~sink() = default;

  /** Returns the value associated with the consumer as a string. **/
  virtual std::string to_str() const = 0;

  /** See argument::consume */
  virtual size_t consume(string_vec_citer start,
                         const string_vec_citer& end) = 0;

  /** Indicates if the sink has been supplied with a default value. */
  virtual bool has_default() const;

  /** Returns the list of valid choices, if any, formatted as strings */
  virtual string_vec choices() const;

  /** Indicates the minimum number of values taken by this sink */
  size_t min_values() const;
  /** Indicates the maximum number of values taken by this sink */
  size_t max_values() const;

protected:
  //! Indicates if the sink has been supplied with a default value
  bool m_has_default;
  //! The minimum number of values taken by this sink
  size_t m_min_values;
  //! The maximum number of values taken by this sink
  size_t m_max_values;

private:
  //! Copy construction not supported
  sink(const sink&) = delete;
  //! Assignment not supported
  sink& operator=(const sink&) = delete;
};

class bool_sink : public sink
{
public:
  explicit bool_sink(bool* sink);

  std::string to_str() const override;
  size_t consume(string_vec_citer, const string_vec_citer& end) override;

private:
  //! Copy construction not supported
  bool_sink(const bool_sink&) = delete;
  //! Assignment not supported
  bool_sink& operator=(const bool_sink&) = delete;

  bool* m_sink;
};

class uint_sink : public sink
{
public:
  explicit uint_sink(unsigned* sink);

  uint_sink& with_default(unsigned value);

  std::string to_str() const override;
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

private:
  //! Copy construction not supported
  uint_sink(const uint_sink&) = delete;
  //! Assignment not supported
  uint_sink& operator=(const uint_sink&) = delete;

  unsigned* m_sink;
};

class double_sink : public sink
{
public:
  explicit double_sink(double* sink);

  double_sink& with_default(double value);

  std::string to_str() const override;
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

private:
  //! Copy construction not supported
  double_sink(const double_sink&) = delete;
  //! Assignment not supported
  double_sink& operator=(const double_sink&) = delete;

  double* m_sink;
};

class str_sink : public sink
{
public:
  explicit str_sink(std::string* sink);

  str_sink& with_default(const char* value);
  str_sink& with_default(const std::string& value);
  str_sink& with_choices(const string_vec& choices);

  std::string to_str() const override;
  string_vec choices() const override;
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

private:
  //! Copy construction not supported
  str_sink(const str_sink&) = delete;
  //! Assignment not supported
  str_sink& operator=(const str_sink&) = delete;

  std::string* m_sink;
  string_vec m_choices;
};

class vec_sink : public sink
{
public:
  explicit vec_sink(string_vec* sink);

  /** The minimum number of values expected on the command-line (default 1) */
  vec_sink& with_min_values(size_t n);
  /** The maximum number of values expected on the command-line (default inf) */
  vec_sink& with_max_values(size_t n);

  std::string to_str() const override;
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

private:
  //! Copy construction not supported
  vec_sink(const vec_sink&) = delete;
  //! Assignment not supported
  vec_sink& operator=(const vec_sink&) = delete;

  string_vec* m_sink;
};

} // namespace argparse

} // namespace adapterremoval
