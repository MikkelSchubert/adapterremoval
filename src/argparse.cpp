/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <algorithm>   // for min, copy, max, replace
#include <limits>      // for numeric_limits
#include <memory>      // for make_unique, make_shared
#include <set>         // for set
#include <sstream>     // for istringstream, ostringstream
#include <stdexcept>   // for invalid_argument
#include <sys/ioctl.h> // for ioctl, winsize, TIOCGWINSZ
#include <unistd.h>    // for STDERR_FILENO

#include "argparse.hpp" // header
#include "debug.hpp"    // for AR_REQUIRE
#include "logging.hpp"  // for log
#include "strutils.hpp" // for cli_formatter, str_to_unsigned, toupper

namespace adapterremoval {

namespace argparse {

const size_t parsing_failed = static_cast<size_t>(-1);

/** Returns the number of columns available in the terminal. */
size_t
get_terminal_columns()
{
  struct winsize params;
  if (ioctl(STDERR_FILENO, TIOCGWINSZ, &params)) {
    // Default to 80 columns if the parameters could not be retrieved.
    return 80;
  }

  return std::min<size_t>(100, std::max<size_t>(80, params.ws_col));
}

/** Detect similar arguments based on prefixes or max edit distance. */
bool
is_similar_argument(const std::string& user,
                    const std::string& ref,
                    size_t max_distance)
{
  const auto overlap = std::min(user.size(), ref.size());
  if (overlap == user.size() && user == ref.substr(0, overlap)) {
    return true;
  }

  const auto diff = std::max(user.size(), ref.size()) - overlap;
  if (diff <= max_distance && levenshtein(user, ref) <= max_distance) {
    return true;
  }

  return false;
}

bool
to_double(const std::string& value, double& out)
{
  std::istringstream stream(value);
  if (!(stream >> out)) {
    return false;
  }

  // Failing on trailing, non-numerical values
  char trailing;
  if (stream >> trailing) {
    return false;
  }

  return true;
}

std::string
escape(const std::string& s)
{
  std::string out;
  out.append("\"");
  for (auto c : s) {
    switch (c) {
      case '\\':
        out.push_back('\\');
        break;
      case '\"':
        out.append("\\\"");
        break;
      default:
        out.push_back(c);
    }
  }

  out.append("\"");

  return out;
}

///////////////////////////////////////////////////////////////////////////////

parser::parser(const std::string& name,
               const std::string& version,
               const std::string& help)
  : m_args()
  , m_keys()
  , m_name(name)
  , m_version(version)
  , m_help(help)
  , m_terminal_width(100)
{
  add_header("OPTIONS:");

  // Built-in arguments
  add("--help").abbreviation('h').help("Display this message.");
  add("--version").abbreviation('v').help("Print the version string.");
  add_separator();
}

parse_result
parser::parse_args(int argc, char const* const* argv)
{
  update_argument_map();
  const string_vec argvec(argv + 1, argv + argc);

  string_vec_citer it = argvec.begin();
  while (it != argvec.end()) {
    const auto argument = find_argument(*it);
    if (argument) {
      const size_t consumed = argument->parse(it, argvec.end());

      if (consumed == parsing_failed) {
        return parse_result::error;
      }

      it += static_cast<string_vec::iterator::difference_type>(consumed);
      AR_REQUIRE(it <= argvec.end());
    } else {
      return parse_result::error;
    }
  }

  parse_result result = parse_result::ok;
  for (const auto& arg : m_args) {
    if (arg.argument && arg.argument->is_set()) {
      const auto& key = arg.argument->key();

      for (const auto& requirement : arg.argument->depends_on()) {
        if (!is_set(requirement)) {
          result = parse_result::error;
          log::error() << "Option " << requirement << " is required when "
                       << "using option " << key;
        }
      }

      for (const auto& prohibited : arg.argument->conflicts_with()) {
        if (is_set(prohibited)) {
          result = parse_result::error;
          log::error() << "Option " << prohibited
                       << " cannot be used together with option " << key;
        }
      }
    }
  }

  if (is_set("--help")) {
    print_help();
    return parse_result::exit;
  } else if (is_set("--version")) {
    print_version();
    return parse_result::exit;
  }

  return result;
}

bool
parser::is_set(const std::string& key) const
{
  const auto it = m_keys.find(key);
  AR_REQUIRE(it != m_keys.end());

  return it->second->is_set();
}

argument&
parser::add(const std::string& name, const std::string& metavar)
{
  auto ptr = std::make_shared<argument>(name, metavar);
  m_args.push_back({ std::string(), ptr });

  return *ptr;
}

void
parser::add_separator()
{
  m_args.push_back({ std::string(), argument_ptr() });
}

void
parser::add_header(const std::string& header)
{
  add_separator();
  m_args.push_back({ header, argument_ptr() });
}

void
parser::print_version() const
{
  log::cerr() << m_name << " " << m_version << "\n";
}

void
parser::print_help() const
{
  print_version();

  auto cerr = log::cerr();
  cerr << "\n" << m_help;

  string_vec signatures;

  size_t indentation = 0;
  for (const auto& entry : m_args) {
    if (entry.argument && !entry.argument->is_deprecated()) {
      const auto& arg = *entry.argument;

      std::ostringstream ss;
      if (arg.short_key().empty()) {
        ss << "   " << arg.key();
      } else {
        ss << "   " << arg.short_key() << ", " << arg.key();
      }

      for (size_t i = 0; i < arg.min_values(); ++i) {
        ss << " <" << arg.metavar() << ">";
      }

      if (arg.max_values() == std::numeric_limits<size_t>::max()) {
        ss << " [" << arg.metavar() << ", ...]";
      } else {
        for (size_t i = arg.min_values(); i < arg.max_values(); ++i) {
          ss << " [" << arg.metavar() << "]";
        }
      }

      indentation = std::max<size_t>(indentation, ss.str().size() + 3);
      signatures.push_back(ss.str());
    } else {
      signatures.push_back(std::string());
    }
  }

  cli_formatter fmt;
  fmt.set_indent(indentation);
  fmt.set_indent_first_line(false);
  fmt.set_column_width(m_terminal_width - indentation);

  for (size_t i = 0; i < m_args.size(); i++) {
    const auto& entry = m_args.at(i);
    if (entry.argument) {
      if (!entry.argument->is_deprecated()) {
        const auto& arg = *entry.argument;
        const auto signature = signatures.at(i);

        cerr << signature;

        const std::string help = arg.help();
        if (!help.empty()) {
          // Format into columns and indent lines (except the first line)
          cerr << std::string(indentation - signature.length(), ' ')
               << fmt.format(help);
        }

        cerr << "\n";
      }
    } else {
      cerr << entry.header << "\n";
    }
  }
}

void
parser::set_terminal_width(unsigned w)
{
  m_terminal_width = w;
}

void
parser::update_argument_map()
{
  m_keys.clear();

  for (auto& it : m_args) {
    if (it.argument) {
      for (const auto& key : it.argument->keys()) {
        const auto result = m_keys.emplace(key, it.argument);
        AR_REQUIRE(result.second);
      }
    }
  }

  bool any_errors = false;
  for (auto& it : m_args) {
    if (it.argument) {
      for (const auto& key : it.argument->conflicts_with()) {
        if (!m_keys.count(key)) {
          any_errors = true;
          log::error() << it.argument->key() << " conflicts with "
                       << "unknown command-line option " << key;
        }
      }

      for (const auto& key : it.argument->depends_on()) {
        if (!m_keys.count(key)) {
          any_errors = true;
          log::error() << it.argument->key() << " requires "
                       << "unknown command-line option " << key;
        }
      }
    }
  }

  AR_REQUIRE(!any_errors, "bugs in argument parsing");
}

argument_ptr
parser::find_argument(const std::string& key)
{
  auto it = m_keys.find(key);
  if (it != m_keys.end()) {
    return it->second;
  }

  if (key.size() > 1 && key.front() == '-' && key.back() != '-') {
    string_vec candidates;
    const size_t max_distance = 1 + key.size() / 4;

    for (const auto& arg : m_args) {
      if (arg.argument) {
        for (const auto& name : arg.argument->keys()) {

          if (is_similar_argument(key, name, max_distance)) {
            candidates.push_back(name);
          }
        }
      }
    }

    auto error = log::error();
    error << "Unknown argument '" << key << "'\n";

    std::sort(candidates.begin(), candidates.end());
    if (!candidates.empty()) {
      error << "\n    Did you mean\n";

      for (const auto& candidate : candidates) {
        error << "      " << candidate << "\n";
      }
    }
  } else {
    log::error() << "Unexpected positional argument '" << key << "'";
  }

  return argument_ptr();
}

///////////////////////////////////////////////////////////////////////////////

argument::argument(const std::string& key, const std::string& metavar)
  : m_times_set()
  , m_default_sink()
  , m_deprecated()
  , m_deprecated_keys()
  , m_key_long(key)
  , m_key_short()
  , m_metavar(metavar)
  , m_help()
  , m_depends_on()
  , m_conflicts_with()
  , m_sink(std::make_unique<bool_sink>(&m_default_sink))
{
  AR_REQUIRE(key.size() && key.at(0) == '-');
}

argument&
argument::help(const std::string& text)
{
  m_help = text;

  return *this;
}

bool
argument::is_set() const
{
  return m_times_set;
}

bool
argument::is_deprecated() const
{
  return m_deprecated;
}

const std::string&
argument::key() const
{
  return m_key_long;
}

const std::string&
argument::short_key() const
{
  return m_key_short;
}

string_vec
argument::keys() const
{
  string_vec keys = m_deprecated_keys;
  keys.push_back(m_key_long);

  if (!m_key_short.empty()) {
    keys.push_back(m_key_short);
  }

  return keys;
}

const std::string&
argument::metavar() const
{
  return m_metavar;
}

std::string
argument::help() const
{
  // Append string representation of current (default) value
  std::ostringstream ss(m_help);
  ss << m_help;

  const auto choices = m_sink->choices();
  if (!choices.empty()) {
    ss << ". Possible values are ";

    for (size_t i = 0; i < choices.size(); ++i) {
      if (i && i == choices.size() - 1) {
        ss << ", and ";
      } else if (i) {
        ss << ", ";
      }

      ss << choices.at(i);
    }
  }

  if (m_sink->has_default() && m_help.find("[default:") == std::string::npos) {
    ss << " [default: " << to_str() << "]";
  }

  return ss.str();
}

size_t
argument::min_values() const
{
  return m_sink->min_values();
}

size_t
argument::max_values() const
{
  return m_sink->max_values();
}

const string_vec&
argument::depends_on() const
{
  return m_depends_on;
}

const string_vec&
argument::conflicts_with() const
{
  return m_conflicts_with;
}

std::string
argument::to_str() const
{
  return m_sink->to_str();
}

template<typename A, typename B>
A&
bind(std::unique_ptr<sink>& ptr, B* sink)
{
  ptr = std::make_unique<A>(sink);

  return static_cast<A&>(*ptr);
}

bool_sink&
argument::bind_bool(bool* sink)
{
  return bind<bool_sink>(m_sink, sink);
}

uint_sink&
argument::bind_uint(unsigned* sink)
{
  return bind<uint_sink>(m_sink, sink);
}

double_sink&
argument::bind_double(double* sink)
{
  return bind<double_sink>(m_sink, sink);
}

str_sink&
argument::bind_str(std::string* sink)
{
  return bind<str_sink>(m_sink, sink);
}

vec_sink&
argument::bind_vec(string_vec* sink)
{
  return bind<vec_sink>(m_sink, sink);
}

argument&
argument::abbreviation(char key)
{
  m_key_short.clear();
  m_key_short.push_back('-');
  m_key_short.push_back(key);

  return *this;
}

argument&
argument::deprecated_alias(const std::string& key)
{
  AR_REQUIRE(key.size() && key.at(0) == '-');
  m_deprecated_keys.emplace_back(key);

  return *this;
}

argument&
argument::deprecated()
{
  m_deprecated = true;

  return *this;
}

argument&
argument::depends_on(const std::string& key)
{
  AR_REQUIRE(key.size() && key.at(0) == '-');
  m_depends_on.push_back(key);

  return *this;
}

argument&
argument::conflicts_with(const std::string& key)
{
  AR_REQUIRE(key.size() && key.at(0) == '-');
  m_conflicts_with.push_back(key);

  return *this;
}

void
n_args_error(const std::string& key,
             size_t limit,
             const char* relation,
             size_t n)
{
  auto out = log::error();

  out << "Command-line argument " << key << " takes" << relation << " " << limit
      << " value";

  if (limit != 1) {
    out << "s";
  }

  if (n == 1) {
    out << ", but 1 value was provided!";
  } else {
    out << ", but " << n << " values were provided!";
  }
}

size_t
argument::parse(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(start != end);

  bool is_deprecated = is_deprecated_alias(*start);
  AR_REQUIRE(is_deprecated || *start == m_key_long ||
             (m_key_short.size() && *start == m_key_short));

  if (m_deprecated) {
    log::warn() << "Option " << *start << " is deprecated and will "
                << "be removed in the future.";
  } else if (is_deprecated) {
    log::warn() << "Option " << *start << " is deprecated and will "
                << "be removed in the future. Please use " << key()
                << " instead.";
  }

  if (m_times_set == 1) {
    log::warn() << "Command-line option " << key()
                << " has been specified more than once.";
  }

  double numeric_sink = 0;
  auto end_of_values = start + 1;
  for (; end_of_values != end; ++end_of_values) {
    if (end_of_values->size() > 1 && end_of_values->front() == '-' &&
        // Avoid confusing numeric values for command-line arguments
        !to_double(*end_of_values, numeric_sink)) {
      break;
    }
  }

  const auto min_values = m_sink->min_values();
  const auto max_values = m_sink->max_values();

  AR_REQUIRE(start < end_of_values);
  auto n_values = static_cast<size_t>(end_of_values - start - 1);

  if (n_values != min_values && min_values == max_values) {
    n_args_error(*start, min_values, "", n_values);

    return parsing_failed;
  } else if (n_values < min_values) {
    n_args_error(*start, min_values, " at least", n_values);

    return parsing_failed;
  } else if (n_values > max_values) {
    n_args_error(*start, max_values, " at most", n_values);

    return parsing_failed;
  }

  const auto result = m_sink->consume(start + 1, end_of_values);
  if (result == parsing_failed) {
    log::error() << "Invalid value for " << *start << ": "
                 << escape(*(start + 1));

    return result;
  }

  m_times_set++;
  return result + 1;
}

bool
argument::is_deprecated_alias(const std::string& key) const
{
  return std::find(m_deprecated_keys.begin(), m_deprecated_keys.end(), key) !=
         m_deprecated_keys.end();
}

///////////////////////////////////////////////////////////////////////////////
// sink

sink::sink(size_t n_values)
  : sink(n_values, n_values)
{
}

sink::sink(size_t min_values, size_t max_values)
  : m_has_default()
  , m_min_values(min_values)
  , m_max_values(max_values)
{
}

sink::~sink() {}

bool
sink::has_default() const
{
  return m_has_default;
}

string_vec
sink::choices() const
{
  return string_vec();
}

size_t
sink::min_values() const
{
  return m_min_values;
}

size_t
sink::max_values() const
{
  return m_max_values;
}

///////////////////////////////////////////////////////////////////////////////
// bool_sink

bool_sink::bool_sink(bool* ptr)
  : sink(0)
  , m_sink(ptr)
{
  AR_REQUIRE(ptr);

  *m_sink = false;
}

std::string
bool_sink::to_str() const
{
  return *m_sink ? "on" : "off";
}

size_t
bool_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(start == end);
  *m_sink = true;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// uint_sink

uint_sink::uint_sink(unsigned* ptr)
  : sink(1)
  , m_sink(ptr)
{
  AR_REQUIRE(ptr);

  *m_sink = 0;
}

uint_sink&
uint_sink::with_default(unsigned value)
{
  m_has_default = true;
  *m_sink = value;

  return *this;
}

std::string
uint_sink::to_str() const
{
  return std::to_string(*m_sink);
}

size_t
uint_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(end - start == 1);

  try {
    *m_sink = str_to_unsigned(*start);

    return 1;
  } catch (const std::invalid_argument&) {
    return parsing_failed;
  }
}

///////////////////////////////////////////////////////////////////////////////
// double_sink

double_sink::double_sink(double* ptr)
  : sink(1)
  , m_sink(ptr)
{
  AR_REQUIRE(ptr);

  *m_sink = 0.0;
}

double_sink&
double_sink::with_default(double value)
{
  m_has_default = true;
  *m_sink = value;

  return *this;
}

std::string
double_sink::to_str() const
{
  auto s = std::to_string(*m_sink);
  s.erase(s.find_last_not_of('0') + 1, std::string::npos);
  s.erase(s.find_last_not_of('.') + 1, std::string::npos);

  return s;
}

size_t
double_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(end - start == 1);

  double value = 0;
  if (!to_double(*start, value)) {
    return parsing_failed;
  }

  *m_sink = value;
  return 1;
}

///////////////////////////////////////////////////////////////////////////////
// str_sink

str_sink::str_sink(std::string* ptr)
  : sink(1)
  , m_sink(ptr)
  , m_choices()
{
  AR_REQUIRE(ptr);

  m_sink->clear();
}

str_sink&
str_sink::with_default(const char* value)
{
  m_has_default = true;
  *m_sink = value;

  return *this;
}

str_sink&
str_sink::with_choices(const string_vec& value)
{
  m_choices = value;

  return *this;
}

std::string
str_sink::to_str() const
{
  return escape(*m_sink);
}

string_vec
str_sink::choices() const
{
  return m_choices;
}

size_t
str_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(end - start == 1);

  if (m_choices.empty()) {
    *m_sink = *start;
    return 1;
  } else {
    const auto value = tolower(*start);
    for (const auto& it : m_choices) {
      if (value == tolower(it)) {
        *m_sink = it;
        return 1;
      }
    }

    return parsing_failed;
  }
}

vec_sink::vec_sink(string_vec* ptr)
  : sink(1, std::numeric_limits<size_t>::max())
  , m_sink(ptr)
{
  AR_REQUIRE(ptr);

  m_sink->clear();
}

vec_sink&
vec_sink::max_values(size_t n)
{
  AR_REQUIRE(n >= m_min_values);
  m_max_values = n;

  return *this;
}

size_t
vec_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(end - start >= 1);

  m_sink->assign(start, end);

  return static_cast<size_t>(end - start);
}

std::string
vec_sink::to_str() const
{
  std::string output;

  for (const auto& s : *m_sink) {
    if (!output.empty()) {
      output.push_back(';');
    }

    output.append(escape(s));
  }

  return output;
}

} // namespace argparse

} // namespace adapterremoval
