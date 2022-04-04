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
#include <algorithm>   // for min, copy, max, replace
#include <iomanip>     // for operator<<, setw
#include <iostream>    // for operator<<, basic_ostream, endl, cerr, ostream
#include <limits>      // for numeric_limits
#include <set>         // for set
#include <sstream>     // for stringstream
#include <stdexcept>   // for invalid_argument
#include <sys/ioctl.h> // for ioctl, winsize, TIOCGWINSZ
#include <unistd.h>    // for STDERR_FILENO

#include "argparse.hpp"
#include "debug.hpp"    // for AR_DEBUG_ASSERT
#include "strutils.hpp" // for cli_formatter, str_to_unsigned, toupper

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

  return std::min<size_t>(88, std::max<size_t>(50, params.ws_col));
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
  , m_stream(&std::cerr)
  , m_terminal_width(get_terminal_columns())
{
  add_header("OPTIONS:");

  // Built-in arguments
  add("-h").alias("--help").help("Display this message.");
  add("-v").alias("--version").help("Print the version string.");
  add_separator();
}

parse_result
parser::parse_args(int argc, char const* const* argv)
{
  update_argument_map();
  const string_vec argvec(argv + 1, argv + argc);

  string_vec_citer it = argvec.begin();
  while (it != argvec.end()) {
    auto argument = find_argument(*it);
    if (argument) {
      const size_t consumed = argument->parse(it, argvec.end());

      if (consumed == parsing_failed) {
        return parse_result::error;
      }

      it += static_cast<string_vec::iterator::difference_type>(consumed);
    } else {
      return parse_result::error;
    }
  }

  parse_result result = parse_result::ok;
  for (const auto& it : m_args) {
    if (it.argument && it.argument->is_set()) {
      const auto& key = it.argument->key();

      for (const auto& requirement : it.argument->requires()) {
        if (!is_set(requirement)) {
          result = parse_result::error;
          *m_stream << "ERROR: Option " << requirement << " is required when "
                    << "using option " << key << std::endl;
        }
      }

      for (const auto& prohibited : it.argument->conflicts()) {
        if (is_set(prohibited)) {
          result = parse_result::error;
          *m_stream << "ERROR: Option " << prohibited
                    << " cannot be used together with option " << key
                    << std::endl;
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
  AR_DEBUG_ASSERT(it != m_keys.end());

  return it->second->is_set();
}

argument&
parser::add(const std::string& name, const std::string& metavar)
{
  argument_ptr ptr = argument_ptr(new argument(name, metavar));
  m_args.push_back({ std::string(), ptr });
  ptr->set_ostream(m_stream);

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
  *m_stream << m_name << " " << m_version << std::endl;
}

void
parser::print_help() const
{
  print_version();
  *m_stream << "\n" << m_help;

  const size_t indentation = 4;

  cli_formatter fmt;
  fmt.set_indent(indentation);
  fmt.set_column_width(m_terminal_width - indentation);

  for (const auto& entry : m_args) {
    const auto& ptr = entry.argument;

    if (!ptr) {
      *m_stream << entry.header << "\n";
      continue;
    } else if (ptr->is_deprecated()) {
      continue;
    }

    bool first = true;
    *m_stream << "  ";
    for (const auto& key : ptr->keys()) {
      if (!key.deprecated) {
        *m_stream << (first ? "" : ", ") << key.name;
        first = false;
      }
    };

    if (!ptr->metavar().empty()) {
      *m_stream << " " << ptr->metavar();
    }

    *m_stream << "\n";

    const std::string help = ptr->help();
    if (!help.empty()) {
      // Format into columns and indent lines (except the first line)
      *m_stream << fmt.format(help) << "\n";
    }
  }

  m_stream->flush();
}

void
parser::set_ostream(std::ostream* stream)
{
  AR_DEBUG_ASSERT(stream);

  m_stream = stream;

  for (auto& it : m_args) {
    if (it.argument) {
      it.argument->set_ostream(stream);
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
        const auto result = m_keys.emplace(key.name, it.argument);
        AR_DEBUG_ASSERT(result.second);
      }
    }
  }

  bool any_errors = false;
  for (auto& it : m_args) {
    if (it.argument) {
      for (const auto& key : it.argument->conflicts()) {
        if (!m_keys.count(key)) {
          any_errors = true;
          *m_stream << "ERROR: " << it.argument->key() << " conflicts with "
                    << "unknown command-line option " << key << std::endl;
        }
      }

      for (const auto& key : it.argument->requires()) {
        if (!m_keys.count(key)) {
          any_errors = true;
          *m_stream << "ERROR: " << it.argument->key() << " requires "
                    << "unknown command-line option " << key << std::endl;
        }
      }
    }
  }

  if (any_errors) {
    AR_DEBUG_FAIL("bugs in argument parsing");
  }
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
        for (const auto& it : arg.argument->keys()) {
          if (!it.deprecated &&
              is_similar_argument(key, it.name, max_distance)) {
            candidates.push_back(it.name);
          }
        }
      }
    }

    *m_stream << "ERROR: Unknown argument '" << key << "'" << std::endl;

    std::sort(candidates.begin(), candidates.end());
    if (!candidates.empty()) {
      *m_stream << "\n    Did you mean" << std::endl;

      for (const auto& key : candidates) {
        *m_stream << "      " << key << std::endl;
      }
    }
  } else {
    *m_stream << "ERROR: Unexpected positional argument '" << key << "'"
              << std::endl;
  }

  return argument_ptr();
}

///////////////////////////////////////////////////////////////////////////////

argument::argument(const std::string& key, const std::string& metavar)
  : m_times_set()
  , m_default_sink()
  , m_deprecated()
  , m_keys()
  , m_metavar(metavar)
  , m_help()
  , m_requires()
  , m_conflicts()
  , m_sink(new bool_sink(&m_default_sink))
  , m_stream(&std::cerr)
{
  AR_DEBUG_ASSERT(key.size() && key.at(0) == '-');

  m_keys.push_back({ key, false });
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
  AR_DEBUG_ASSERT(m_keys.size());
  return m_keys.front().name;
}

const std::vector<argument::argument_key>&
argument::keys() const
{
  return m_keys;
}

const std::string&
argument::metavar() const
{
  return m_metavar;
}

std::string
argument::help() const
{
  if (m_sink->has_default() && m_help.find("[default:") == std::string::npos) {
    // Append string representation of current (default) value
    std::stringstream ss;
    ss << m_help << " [default: " << to_str() << "]";

    return ss.str();
  } else {
    return m_help;
  }
}

const string_vec& argument::requires() const
{
  return m_requires;
}

const string_vec&
argument::conflicts() const
{
  return m_conflicts;
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
  auto* argsink = new A(sink);
  ptr.reset(argsink);

  return *argsink;
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
argument::alias(const std::string& key)
{
  m_keys.push_back({ key, false });

  return *this;
}

argument&
argument::deprecated_alias(const std::string& key)
{
  m_keys.push_back({ key, true });

  return *this;
}

argument&
argument::deprecated()
{
  m_deprecated = true;

  return *this;
}

argument& argument::requires(const std::string& key)
{
  m_requires.push_back(key);

  return *this;
}

argument&
argument::conflicts(const std::string& key)
{
  m_conflicts.push_back(key);

  return *this;
}

void
n_args_error(std::ostream& out,
             const std::string& key,
             size_t limit,
             const char* relation,
             size_t n)
{
  out << "ERROR: Command-line argument " << key << " takes" << relation << " "
      << limit << " value";

  if (limit != 1) {
    out << "s";
  }

  if (n == 1) {
    out << ", but 1 value was provided!";
  } else {
    out << ", but " << n << " values were provided!";
  }

  out << std::endl;
}

size_t
argument::parse(string_vec_citer start, const string_vec_citer& end)
{
  AR_DEBUG_ASSERT(start != end);
  for (const auto& key_ : m_keys) {
    if (*start == key_.name) {
      if (m_deprecated) {
        *m_stream << "WARNING: Option " << *start << " is deprecated and will "
                  << "be removed in the future." << std::endl;
      } else if (key_.deprecated) {
        *m_stream << "WARNING: Option " << *start << " is deprecated and will "
                  << "be removed in the future. Please use " << key()
                  << " instead." << std::endl;
      }

      if (m_times_set == 1) {
        *m_stream << "WARNING: Command-line option " << key()
                  << " has been specified more than once." << std::endl;
      }

      auto end_of_values = start + 1;
      for (; end_of_values != end; ++end_of_values) {
        if (!end_of_values->empty() && end_of_values->front() == '-') {
          break;
        }
      }

      AR_DEBUG_ASSERT(start < end_of_values);
      const auto n_values = static_cast<size_t>(end_of_values - start - 1);
      const auto min_values = m_sink->min_values();
      const auto max_values = m_sink->max_values();

      if (n_values != min_values && min_values == max_values) {
        n_args_error(*m_stream, *start, min_values, "", n_values);

        return parsing_failed;
      } else if (n_values < min_values) {
        n_args_error(*m_stream, *start, min_values, " at least", n_values);

        return parsing_failed;
      } else if (n_values > max_values) {
        n_args_error(*m_stream, *start, max_values, " at most", n_values);

        return parsing_failed;
      }

      const auto result = m_sink->consume(start + 1, end_of_values);
      if (result == parsing_failed) {
        *m_stream << "ERROR: Invalid value for " << *start << ": "
                  << escape(*(start + 1)) << std::endl;

        return result;
      }

      m_times_set++;
      return result + 1;
    }
  }

  return parsing_failed;
}

void
argument::set_ostream(std::ostream* stream)
{
  AR_DEBUG_ASSERT(stream);

  m_stream = stream;
}

///////////////////////////////////////////////////////////////////////////////
// sink

sink::sink(size_t n_values)
  : sink(n_values, n_values)
{}

sink::sink(size_t min_values, size_t max_values)
  : m_has_default()
  , m_min_values(min_values)
  , m_max_values(max_values)
{}

sink::~sink() {}

bool
sink::has_default() const
{
  return m_has_default;
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
  AR_DEBUG_ASSERT(ptr);

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
  AR_DEBUG_ASSERT(start == end);
  *m_sink = true;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// uint_sink

uint_sink::uint_sink(unsigned* ptr)
  : sink(1)
  , m_sink(ptr)
{
  AR_DEBUG_ASSERT(ptr);

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
  AR_DEBUG_ASSERT(end - start == 1);

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
  AR_DEBUG_ASSERT(ptr);

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
  AR_DEBUG_ASSERT(end - start == 1);

  double value = 0;
  std::stringstream stream(*start);
  if (!(stream >> value)) {
    return parsing_failed;
  }

  // Failing on trailing, non-numerical values
  std::string trailing;
  if (stream >> trailing) {
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
{
  AR_DEBUG_ASSERT(ptr);

  m_sink->clear();
}

str_sink&
str_sink::with_default(const char* value)
{
  m_has_default = true;
  *m_sink = value;

  return *this;
}

std::string
str_sink::to_str() const
{
  return escape(*m_sink);
}

size_t
str_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_DEBUG_ASSERT(end - start == 1);

  *m_sink = *start;

  return 1;
}

vec_sink::vec_sink(string_vec* ptr)
  : sink(1, std::numeric_limits<size_t>::max())
  , m_sink(ptr)
{
  AR_DEBUG_ASSERT(ptr);

  m_sink->clear();
}

size_t
vec_sink::consume(string_vec_citer start, const string_vec_citer& end)
{
  AR_DEBUG_ASSERT(end - start >= 1);

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
