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
#include "argparse.hpp"
#include "debug.hpp"    // for AR_REQUIRE
#include "logging.hpp"  // for log_stream, error, cerr, warn
#include "strutils.hpp" // for string_vec, shell_escape, to_lower
#include <algorithm>    // for max, copy, find, min, sort
#include <limits>       // for numeric_limits
#include <memory>       // for __shared_ptr_access, unique_ptr, share...
#include <sstream>      // for operator<<, basic_ostream, ostringstream
#include <stdexcept>    // for invalid_argument
#include <utility>      // for pair

namespace adapterremoval {

namespace argparse {

namespace {

const size_t parsing_failed = static_cast<size_t>(-1);
const size_t invalid_choice = static_cast<size_t>(-2);

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

  return (diff <= max_distance && levenshtein(user, ref) <= max_distance);
}

} // namespace

///////////////////////////////////////////////////////////////////////////////

parser::parser()
{
  add_header("OPTIONS:");

  // Built-in arguments
  add("--help").abbreviation('h').help("Display this message.");
  add("--version").abbreviation('v').help("Print the version string.");
  add("--licenses").help("Print licenses for this software.");
  add_separator();
}

void
parser::set_name(const std::string& name)
{
  m_name = name;
}

void
parser::set_version(const std::string& version)
{
  m_version = version;
}

void
parser::set_preamble(const std::string& text)
{
  m_preamble = text;
}

void
parser::set_licenses(const std::string& text)
{
  m_licenses = text;
}

parse_result
parser::parse_args(const string_vec& args)
{
  AR_REQUIRE(!args.empty());
  update_argument_map();

  for (auto it = args.begin() + 1; it != args.end();) {
    const auto argument = find_argument(*it);
    if (argument) {
      const size_t consumed = argument->parse(it, args.end());

      if (consumed == parsing_failed) {
        return parse_result::error;
      }

      it += static_cast<string_vec::iterator::difference_type>(consumed);
      AR_REQUIRE(it <= args.end());
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
  } else if (is_set("--licenses")) {
    print_licenses();
    return parse_result::exit;
  }

  return result;
}

bool
parser::is_set(const std::string& key) const
{
  const auto it = m_keys.find(key);
  AR_REQUIRE(it != m_keys.end(), shell_escape(key));

  return it->second->is_set();
}

std::string
parser::value(const std::string& key) const
{
  const auto it = m_keys.find(key);
  AR_REQUIRE(it != m_keys.end(), shell_escape(key));

  return it->second->value();
}

argument&
parser::add(const std::string& name, const std::string& metavar)
{
  AR_REQUIRE(starts_with(name, "--"), name);
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
  if (!m_preamble.empty()) {
    cli_formatter fmt;
    fmt.set_indent(0);
    fmt.set_indent_first_line(false);
    fmt.set_column_width(m_terminal_width);

    cerr << "\n";
    for (auto line : split_lines(m_preamble)) {
      if (starts_with(line, " ")) {
        cerr << line << "\n";
      } else {
        cerr << fmt.format(line) << "\n";
      }
    }
  }

  string_vec signatures;

  size_t indentation = 0;
  for (const auto& entry : m_args) {
    if (entry.argument && !entry.argument->is_hidden()) {
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
      signatures.emplace_back(ss.str());
    } else {
      signatures.emplace_back();
    }
  }

  cli_formatter fmt;
  fmt.set_indent(indentation);
  fmt.set_indent_first_line(false);
  fmt.set_column_width(m_terminal_width - indentation);

  for (size_t i = 0; i < m_args.size(); i++) {
    const auto& entry = m_args.at(i);
    if (entry.argument) {
      if (!entry.argument->is_hidden()) {
        const auto& arg = *entry.argument;
        const auto& signature = signatures.at(i);

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
parser::print_licenses() const
{
  auto cerr = log::cerr();

  for (const auto& block : split_lines(m_licenses)) {
    if (block.empty()) {
      cerr << "\n";
    } else {
      size_t indentation = block.find_first_not_of(' ');
      if (indentation == std::string::npos) {
        indentation = 0;
      }

      size_t ljust = 0;
      if (!block.empty() && block.at(indentation) == '*') {
        ljust = 2;
      } else if (block.size() > indentation + 1 &&
                 block.at(indentation + 1) == '.') {
        ljust = 3;
      }

      const auto width = 80 - indentation;
      for (const auto& line : wrap_text(block, width, ljust)) {
        cerr << std::string(indentation, ' ') << line << "\n";
      }
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
        AR_REQUIRE(result.second, shell_escape(key));
      }
    }
  }

  bool any_errors = false;
  for (const auto& it : m_args) {
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
    error << "Unknown argument '" << key << "'";

    if (!candidates.empty()) {
      std::sort(candidates.begin(), candidates.end());

      error << ". Did you mean " << candidates.front();
      for (size_t i = 1; i < candidates.size() - 1; ++i) {
        error << ", " << candidates.at(i);
      }

      if (candidates.size() > 1) {
        error << " or " << candidates.back();
      }

      error << "?";
    }
  } else {
    log::error() << "Unexpected positional argument '" << key << "'";
  }

  return {};
}

///////////////////////////////////////////////////////////////////////////////

argument::argument(const std::string& key, std::string metavar)
  : m_key_long(key)
  , m_metavar(std::move(metavar))
  , m_sink(std::make_unique<bool_sink>(&m_default_sink))
{
  AR_REQUIRE(!key.empty() && key.at(0) == '-', shell_escape(key));
}

argument&
argument::help(const std::string& text)
{
  m_help = text;

  return *this;
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

std::string
argument::help() const
{
  // Append string representation of current (default) value
  std::ostringstream ss(m_help);
  ss << m_help;

  const auto choices = m_sink->choices();
  if (!choices.empty()) {
    ss << ". Possible values are " << join_text(choices, ", ", ", and ");
  }

  if (m_sink->has_default() && m_help.find("[default:") == std::string::npos) {
    ss << " [default: " << default_value() << "]";
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

std::string
argument::value() const
{
  return m_sink->value();
}

std::string
argument::default_value() const
{
  return m_sink->default_value();
}

template<typename A, typename B>
A&
bind(std::unique_ptr<sink>& ptr, B* sink)
{
  ptr = std::make_unique<A>(sink);

  return static_cast<A&>(*ptr);
}

bool_sink&
argument::bind_bool(bool* ptr)
{
  return bind<bool_sink>(m_sink, ptr);
}

u32_sink&
argument::bind_u32(unsigned* ptr)
{
  return bind<u32_sink>(m_sink, ptr);
}

double_sink&
argument::bind_double(double* ptr)
{
  return bind<double_sink>(m_sink, ptr);
}

str_sink&
argument::bind_str(std::string* ptr)
{
  return bind<str_sink>(m_sink, ptr);
}

vec_sink&
argument::bind_vec(string_vec* ptr)
{
  return bind<vec_sink>(m_sink, ptr);
}

argument&
argument::abbreviation(char key)
{
  // To avoid ambiguity only lower-case alphabetical arguments are allowed
  AR_REQUIRE(key >= 'a' && key <= 'z', std::string{ key });
  m_key_short.clear();
  m_key_short.push_back('-');
  m_key_short.push_back(key);

  return *this;
}

argument&
argument::deprecated_alias(const std::string& key)
{
  AR_REQUIRE(!key.empty() && key.at(0) == '-', shell_escape(key));
  m_deprecated_keys.emplace_back(key);

  return *this;
}

argument&
argument::deprecated()
{
  m_deprecated = true;

  return hidden();
}

argument&
argument::hidden()
{
  m_hidden = true;

  return *this;
}

argument&
argument::depends_on(const std::string& key)
{
  AR_REQUIRE(!key.empty() && key.at(0) == '-', shell_escape(key));
  m_depends_on.push_back(key);

  return *this;
}

argument&
argument::conflicts_with(const std::string& key)
{
  AR_REQUIRE(!key.empty() && key.at(0) == '-', shell_escape(key));
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

  const bool deprecated_alias = is_deprecated_alias(*start);
  AR_REQUIRE(deprecated_alias || *start == m_key_long ||
             (!m_key_short.empty() && *start == m_key_short));

  if (m_deprecated) {
    log::warn() << "Option " << *start << " is deprecated and will "
                << "be removed in the future.";
  } else if (deprecated_alias) {
    log::warn() << "Option " << *start << " has been renamed to " << key()
                << ". Support for the old name will be removed in the future.";
  }

  if (m_times_set == 1) {
    log::warn() << "Command-line option " << key()
                << " has been specified more than once.";
  }

  auto end_of_values = start + 1;
  for (; end_of_values != end; ++end_of_values) {
    if (end_of_values->size() > 1 && end_of_values->front() == '-') {
      // Avoid confusing numeric values for command-line arguments
      try {
        str_to_double(*end_of_values);
      } catch (const std::invalid_argument&) {
        break;
      }
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

  size_t result = parsing_failed;
  std::string error_message;

  try {
    result = m_sink->consume(start + 1, end_of_values);
  } catch (const std::invalid_argument& error) {
    error_message = error.what();
  }

  if (result == parsing_failed) {
    auto error = log::error();

    error << "Invalid command-line argument " << *start;
    for (auto it = start + 1; it != end_of_values; ++it) {
      error << " " << shell_escape(*it);
    }

    if (!error_message.empty()) {
      error << ": " << error_message;
    }

    return result;
  } else if (result == invalid_choice) {
    auto error = log::error();

    error << "Invalid command-line argument " << *start;
    for (auto it = start + 1; it != end_of_values; ++it) {
      error << " " << shell_escape(*it);
    }

    error << ". Valid values for " << *start << " are "
          << join_text(m_sink->choices(), ", ", ", and ");
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
  : m_min_values(min_values)
  , m_max_values(max_values)
{
}

std::string
sink::default_value() const
{
  AR_FAIL("sink::default_value not implemented");
}

std::string
sink::preprocess(std::string value) const
{
  if (m_preprocess) {
    m_preprocess(value);
  }

  return value;
}

///////////////////////////////////////////////////////////////////////////////
// bool_sink

bool_sink::bool_sink(bool* ptr)
  : sink(0)
  , m_sink(ptr ? ptr : &m_fallback_sink)
{
  *m_sink = false;
}

std::string
bool_sink::value() const
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
// u32_sink

u32_sink::u32_sink(uint32_t* ptr)
  : sink(1)
  , m_sink(ptr)
  , m_minimum(std::numeric_limits<decltype(m_minimum)>::lowest())
  , m_maximum(std::numeric_limits<decltype(m_maximum)>::max())
{
  AR_REQUIRE(ptr);

  *m_sink = m_default;
}

u32_sink&
u32_sink::with_default(uint32_t value)
{
  AR_REQUIRE(value >= m_minimum && value <= m_maximum);
  set_has_default();
  *m_sink = m_default = value;

  return *this;
}

std::string
u32_sink::value() const
{
  return std::to_string(*m_sink);
}

std::string
u32_sink::default_value() const
{
  AR_REQUIRE(has_default());
  return std::to_string(m_default);
}

u32_sink&
u32_sink::with_minimum(uint32_t value)
{
  AR_REQUIRE(m_default >= value && value <= m_maximum);
  m_minimum = value;
  return *this;
}

u32_sink&
u32_sink::with_maximum(uint32_t value)
{
  AR_REQUIRE(m_default <= value && value >= m_minimum);
  m_maximum = value;
  return *this;
}

size_t
u32_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(end - start == 1);

  const auto value = str_to_u32(preprocess(*start));
  if (value < m_minimum) {
    throw std::invalid_argument("value must be at least " +
                                std::to_string(m_minimum));
  } else if (value > m_maximum) {
    throw std::invalid_argument("value must be at most " +
                                std::to_string(m_maximum));
  }

  *m_sink = value;
  return 1;
}

///////////////////////////////////////////////////////////////////////////////
// double_sink

double_sink::double_sink(double* ptr)
  : sink(1)
  , m_sink(ptr)
  , m_minimum(std::numeric_limits<decltype(m_minimum)>::lowest())
  , m_maximum(std::numeric_limits<decltype(m_maximum)>::max())
{
  AR_REQUIRE(ptr);
  *m_sink = m_default;
}

double_sink&
double_sink::with_default(double value)
{
  AR_REQUIRE(value >= m_minimum && value <= m_maximum);
  set_has_default();
  *m_sink = m_default = value;

  return *this;
}

namespace {

std::string
double_to_string(double value)
{
  auto s = std::to_string(value);
  s.erase(s.find_last_not_of('0') + 1, std::string::npos);
  s.erase(s.find_last_not_of('.') + 1, std::string::npos);

  return s;
}

} // namespace

std::string
double_sink::default_value() const
{
  AR_REQUIRE(has_default());
  return double_to_string(m_default);
}

double_sink&
double_sink::with_minimum(double value)
{
  AR_REQUIRE(m_default >= value && value >= m_minimum);
  m_minimum = value;

  return *this;
}

double_sink&
double_sink::with_maximum(double value)
{
  AR_REQUIRE(m_default <= value && value >= m_minimum);
  m_maximum = value;
  return *this;
}

std::string
double_sink::value() const
{
  return double_to_string(*m_sink);
}

size_t
double_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(end - start == 1);

  const auto value = str_to_double(*start);
  if (value < m_minimum) {
    throw std::invalid_argument("value must be at least " +
                                double_to_string(m_minimum));
  } else if (value > m_maximum) {
    throw std::invalid_argument("value must be at most " +
                                double_to_string(m_maximum));
  }

  *m_sink = value;
  return 1;
}

///////////////////////////////////////////////////////////////////////////////
// str_sink

str_sink::str_sink(std::string* ptr)
  : sink(1)
  , m_sink(ptr ? ptr : &m_fallback_sink)
{
  AR_REQUIRE(m_sink);

  *m_sink = m_default;
}

str_sink&
str_sink::with_default(const char* value)
{
  set_has_default();
  *m_sink = m_default = value;

  return *this;
}

str_sink&
str_sink::with_default(const std::string& value)
{
  set_has_default();
  *m_sink = m_default = value;

  return *this;
}

str_sink&
str_sink::with_choices(const string_vec& choices)
{
  m_choices = choices;

  return *this;
}

std::string
str_sink::value() const
{
  return *m_sink;
}

std::string
str_sink::default_value() const
{
  AR_REQUIRE(has_default());
  return m_default;
}

size_t
str_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(end - start == 1);

  if (m_choices.empty()) {
    *m_sink = preprocess(*start);
    return 1;
  } else {
    const auto value = preprocess(*start);
    const auto choice = to_lower(value);
    for (const auto& it : m_choices) {
      if (choice == to_lower(it)) {
        *m_sink = it;
        return 1;
      }
    }

    return invalid_choice;
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
vec_sink::with_min_values(size_t n)
{
  AR_REQUIRE(n <= max_values());
  set_min_values(n);

  return *this;
}

vec_sink&
vec_sink::with_max_values(size_t n)
{
  AR_REQUIRE(n >= min_values());
  set_max_values(n);

  return *this;
}

size_t
vec_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_REQUIRE(static_cast<size_t>(end - start) >= min_values());
  AR_REQUIRE(static_cast<size_t>(end - start) <= max_values());
  for (auto it = start; it != end; ++it) {
    m_sink->emplace_back(preprocess(*it));
  }

  return static_cast<size_t>(end - start);
}

std::string
vec_sink::value() const
{
  std::string output;

  for (const auto& s : *m_sink) {
    if (!output.empty()) {
      output.push_back(';');
    }

    output.append(shell_escape(s));
  }

  return output;
}

} // namespace argparse

} // namespace adapterremoval
