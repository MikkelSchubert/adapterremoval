// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp" // for string_vec_citer, string_vec
#include <cstddef>         // for size_t
#include <cstdint>         // for uint32_t
#include <functional>      // for less, function
#include <iosfwd>          // for ostream
#include <map>             // for map
#include <memory>          // for shared_ptr, unique_ptr
#include <string>          // for string
#include <string_view>     // for string_view
#include <vector>          // for vector

namespace adapterremoval::argparse {

class argument;

class sink;
class bool_sink;
class double_sink;
class u32_sink;
class str_sink;
class vec_sink;

using argument_ptr = std::shared_ptr<argument>;
using preprocess_func = std::function<void(std::string&)>;

//! Parse results for command-line arguments
enum class parse_result
{
  //! Terminate now (e.g. --version or --help used)
  exit,
  //! Error occurred parsing arguments / invalid combination of args
  error,
  //! No errors parsing command-line arguments
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
  ~parser() = default;

  /** Sets the name used in --help and --version messages */
  void set_name(const std::string_view& name);
  /** Sets the version string used in --help and --version messages */
  void set_version(const std::string_view& version);
  /** Sets the preamble text used in --help */
  void set_preamble(const std::string_view& text);
  /** Sets the license text used in --licenses */
  void set_licenses(const std::string_view& text);

  /** Parses a set of command-line options as passed to main(argc, argv). */
  parse_result parse_args(const string_vec& args);

  /** Returns true if the option with the given key has been set. */
  [[nodiscard]] bool is_set(std::string_view key) const;
  /** Returns the value associated with the argument as a string. */
  [[nodiscard]] std::string value(std::string_view key) const;

  /** Add argument with metavar. By default this takes no values. */
  argument& add(std::string_view name, std::string_view metavar = {});

  /** Add a blank line between the previous and the next command. */
  void add_separator();

  /** Add a blank line and a header between the previous and next command. */
  void add_header(std::string_view header);

  /** Helper function; prints the program name and version string. */
  void print_version() const;
  /** Helper functions; prints the full set of help-text. */
  void print_help() const;
  /** Helper functions; formats and prints the licenses. */
  void print_licenses() const;

  /** Set the maximum terminal width. */
  void set_terminal_width(unsigned w);

  parser(const parser&) = delete;
  parser(parser&&) = delete;
  parser& operator=(const parser&) = delete;
  parser& operator=(parser&&) = delete;

private:
  /** Collect and validate all arguments */
  void update_argument_map();

  /** Return matching argument or log lexiographically similar arguments */
  argument_ptr find_argument(std::string_view key);

  struct argument_entry
  {
    std::string header{};
    argument_ptr argument{};
  };

  std::vector<argument_entry> m_args{};
  std::map<std::string, argument_ptr, std::less<>> m_keys{};

  //! Name of the program
  std::string m_name{};
  //! Version string for the program (excluding the name)
  std::string m_version{};
  //! Preamble text for the program.
  std::string m_preamble{};
  //! Licenses for the program.
  std::string m_licenses{};
  //! Maximum terminal width used for printing help messages
  unsigned m_terminal_width = 100;
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
  static constexpr size_t parsing_failed = static_cast<size_t>(-1);
  static constexpr size_t invalid_choice = static_cast<size_t>(-2);

  struct argument_key
  {
    std::string name;
    bool deprecated;
  };

  explicit argument(std::string_view key, std::string_view metavar = {});
  ~argument() = default;

  /** Returns true if the consumer has consumed a value. */
  [[nodiscard]] bool is_set() const { return m_times_set; }

  /** Returns true if the argument is deprecated. */
  [[nodiscard]] bool is_deprecated() const { return m_deprecated; }

  /** Returns true if the argument has been completely removed. */
  [[nodiscard]] bool is_removed() const { return m_removed; }

  /** Returns true if the argument is hidden. */
  [[nodiscard]] bool is_hidden() const { return m_hidden; }

  /** Returns the canonical argument key. */
  [[nodiscard]] const std::string& key() const { return m_key_long; }

  /** Returns the short argument key; may be an empty string. */
  [[nodiscard]] const std::string& short_key() const { return m_key_short; }

  /** Returns long, short, and deprecated argument keys. */
  [[nodiscard]] string_vec keys() const;
  /** Returns true if this key is a deprecated alias for this argument. */
  [[nodiscard]] bool is_deprecated_alias(std::string_view key) const;

  /** Returns the meta-variable. May be an empty string. */
  [[nodiscard]] const std::string& metavar() const { return m_metavar; }

  /** Returns help string with %default replaced with the default (if any). */
  [[nodiscard]] std::string help() const;

  /** Indicates the minimum number of values taken by this argument */
  [[nodiscard]] size_t min_values() const;
  /** Indicates the maximum number of values taken by this argument */
  [[nodiscard]] size_t max_values() const;

  /** Options that MUST be specified along with this argument. */
  [[nodiscard]] const string_vec& depends_on() const { return m_depends_on; }

  /** Options that must NOT be specified along with this argument. */
  [[nodiscard]] const string_vec& conflicts_with() const
  {
    return m_conflicts_with;
  }

  /** Returns the value associated with the argument as a string. */
  [[nodiscard]] std::string value() const;
  /** Returns the default value associated as a string. */
  [[nodiscard]] std::string default_value() const;

  /** Set the metavar for this argument. */
  argument& metavar(std::string_view metavar);
  /** Set help string for this argument. */
  argument& help(std::string_view text);
  /** Create a short form of the argument. */
  argument& abbreviation(char key);
  /** Create deprecated alias for the argument. */
  argument& deprecated_alias(std::string_view key);
  /** The argument is deprecated. Implies `hidden()` */
  argument& deprecated();
  /** The argument is deprecated. Implies `hidden()` */
  argument& removed();
  /** The argument will not be printed by -h/--help */
  argument& hidden();

  /** Option `key` MUST be specified along with this argument. */
  argument& depends_on(std::string key);
  /** Option `key` must NOT be specified along with this argument. */
  argument& conflicts_with(std::string key);

  bool_sink& bind_bool(bool* ptr);
  u32_sink& bind_u32(uint32_t* ptr);
  double_sink& bind_double(double* ptr);
  str_sink& bind_str(std::string* ptr);
  vec_sink& bind_vec(string_vec* ptr);

  /** Parse the next arguments, returning the number of items parsed or -1. */
  size_t parse(string_vec_citer start, const string_vec_citer& end);

  argument(const argument&) = delete;
  argument(argument&&) = delete;
  argument& operator=(const argument&) = delete;
  argument& operator=(argument&&) = delete;

private:
  //! Number of times the argument has been specified
  unsigned m_times_set{};
  //! Default sink value
  bool m_default_sink{};
  //! Indicates if the argument is deprecated
  bool m_deprecated{};
  //! Indicates if the argument has been removed
  bool m_removed{};
  //! Deprecated keys (long and short) for this argument
  string_vec m_deprecated_keys{};
  //! Indicates if the argument is hidden
  bool m_hidden{};

  //! The long, canonical argument key
  std::string m_key_long{};
  //! An optional, short argument key
  std::string m_key_short{};

  //! Optional metavar (defaults to uppercase `m_name` without dashes)
  std::string m_metavar{};
  //! Help string; the string '%default' will be replaced with the current value
  std::string m_help{};

  //! This argument must be specified along with these arguments.
  string_vec m_depends_on{};
  //! This argument cannot be specified along with these arguments.
  string_vec m_conflicts_with{};

  std::unique_ptr<sink> m_sink{};
};

class sink
{
public:
  /** Creates a sink that takes exactly `n_values` values */
  explicit sink(size_t n_values);
  /** Creates a sink that takes between min and max values (incl.) */
  sink(size_t min_values, size_t max_values);

  virtual ~sink() = default;

  /** Returns the current argument value as a string. **/
  [[nodiscard]] virtual std::string value() const = 0;
  /** Returns string-representation of the default value */
  [[nodiscard]] virtual std::string default_value() const = 0;

  /** Indicates if the sink has been supplied with a default value. */
  [[nodiscard]] virtual bool has_default() const { return m_has_default; }

  /** See argument::consume */
  virtual size_t consume(string_vec_citer start,
                         const string_vec_citer& end) = 0;

  /** Returns the list of valid choices, if any, formatted as strings */
  [[nodiscard]] virtual string_vec choices() const { return {}; };

  /** Indicates the minimum number of values taken by this sink */
  [[nodiscard]] size_t min_values() const { return m_min_values; };

  /** Indicates the maximum number of values taken by this sink */
  [[nodiscard]] size_t max_values() const { return m_max_values; };

  /** Sets pre-processor function used before validating input  */
  sink& with_preprocessor(preprocess_func func)
  {
    m_preprocess = func;
    return *this;
  }

  sink(const sink&) = delete;
  sink(sink&&) = delete;
  sink& operator=(const sink&) = delete;
  sink& operator=(sink&&) = delete;

protected:
  void set_has_default() { m_has_default = true; }

  void set_min_values(size_t n) { m_min_values = n; }

  void set_max_values(size_t n) { m_max_values = n; }

  /** Preprocess the value if a preprocessor was set */
  std::string preprocess(std::string value) const;

private:
  //! Indicates if the sink has been supplied with a default value
  bool m_has_default{};
  //! The minimum number of values taken by this sink
  size_t m_min_values{};
  //! The maximum number of values taken by this sink
  size_t m_max_values{};
  //! Function used to pre-process the user supplied arguments
  preprocess_func m_preprocess{};
};

class bool_sink : public sink
{
public:
  explicit bool_sink(bool* ptr);
  ~bool_sink() override = default;

  std::string value() const override;
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

  /** Not implemented */
  [[nodiscard]] std::string default_value() const override;

  bool_sink(const bool_sink&) = delete;
  bool_sink(bool_sink&&) = delete;
  bool_sink& operator=(const bool_sink&) = delete;
  bool_sink& operator=(bool_sink&&) = delete;

private:
  bool* m_sink;
  //! Sink variable used if no sink was supplied
  bool m_fallback_sink = false;
};

class u32_sink : public sink
{
public:
  explicit u32_sink(uint32_t* ptr);
  ~u32_sink() override = default;

  u32_sink& with_default(uint32_t value);
  [[nodiscard]] std::string default_value() const override;

  /** Set minimum allowed value (inclusive) */
  u32_sink& with_minimum(uint32_t value);
  /** Set maximum allowed value (inclusive) */
  u32_sink& with_maximum(uint32_t value);

  [[nodiscard]] std::string value() const override;
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

  u32_sink(const u32_sink&) = delete;
  u32_sink(u32_sink&&) = delete;
  u32_sink& operator=(const u32_sink&) = delete;
  u32_sink& operator=(u32_sink&&) = delete;

private:
  uint32_t* m_sink = nullptr;
  //! Default value used for -h/--help output
  uint32_t m_default{};
  //! Minimum allowed value (inclusive)
  uint32_t m_minimum;
  //! Maximum allowed value (inclusive)
  uint32_t m_maximum;
};

class double_sink : public sink
{
public:
  explicit double_sink(double* ptr);
  ~double_sink() override = default;

  double_sink& with_default(double value);
  [[nodiscard]] std::string default_value() const override;

  /** Set minimum allowed value (inclusive) */
  double_sink& with_minimum(double value);
  /** Set maximum allowed value (inclusive) */
  double_sink& with_maximum(double value);

  [[nodiscard]] std::string value() const override;
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

  double_sink(const double_sink&) = delete;
  double_sink(double_sink&&) = delete;
  double_sink& operator=(const double_sink&) = delete;
  double_sink& operator=(double_sink&&) = delete;

private:
  double* m_sink = nullptr;
  //! Default value used for -h/--help output
  double m_default{};
  //! Minimum allowed value (inclusive)
  double m_minimum;
  //! Maximum allowed value (inclusive)
  double m_maximum;
};

class str_sink : public sink
{
public:
  explicit str_sink(std::string* ptr);
  ~str_sink() override = default;

  str_sink& with_default(std::string_view value);
  str_sink& with_choices(const string_vec& choices);
  [[nodiscard]] std::string default_value() const override;

  /** Implicit argument if the user does not supply one */
  str_sink& with_implicit_argument(std::string_view value);

  [[nodiscard]] std::string value() const override;

  [[nodiscard]] string_vec choices() const override { return m_choices; }

  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

  str_sink(const str_sink&) = delete;
  str_sink(str_sink&&) = delete;
  str_sink& operator=(const str_sink&) = delete;
  str_sink& operator=(str_sink&&) = delete;

private:
  std::string* m_sink = nullptr;
  string_vec m_choices{};
  //! Default value used for -h/--help output
  std::string m_default{};
  //! Specifies if a default argument was provided
  bool m_has_implicit_argument{};
  //! Default argument if no value is supplied by the user
  std::string m_implicit_argument{};
  //! Sink variable used if no sink was supplied
  std::string m_fallback_sink{};
};

class vec_sink : public sink
{
public:
  explicit vec_sink(string_vec* ptr);
  ~vec_sink() override = default;

  /** The minimum number of values expected on the command-line (default 1) */
  vec_sink& with_min_values(size_t n);
  /** The maximum number of values expected on the command-line (default inf) */
  vec_sink& with_max_values(size_t n);

  /** Returns string representation of user-provided values */
  [[nodiscard]] std::string value() const override;
  /** Not implemented */
  [[nodiscard]] std::string default_value() const override;

  /** Consumes the provided values, using the specified pre-processor, if any */
  size_t consume(string_vec_citer start, const string_vec_citer& end) override;

  vec_sink(const vec_sink&) = delete;
  vec_sink(vec_sink&&) = delete;
  vec_sink& operator=(const vec_sink&) = delete;
  vec_sink& operator=(vec_sink&&) = delete;

private:
  string_vec* m_sink = nullptr;
};

/** Stream operator for debugging output */
std::ostream&
operator<<(std::ostream& os, const parse_result& value);

} // namespace adapterremoval::argparse