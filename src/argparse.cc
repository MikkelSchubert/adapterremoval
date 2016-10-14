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
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <limits>
#include <set>

#include <sys/types.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "argparse.h"
#include "debug.h"
#include "strutils.h"


namespace ar
{
namespace argparse
{

typedef std::set<consumer_ptr> consumer_set;


/** Returns the number of columns available in the terminal. */
size_t get_terminal_columns()
{
    struct winsize params;
    if (ioctl(STDERR_FILENO, TIOCGWINSZ, &params)) {
        // Default to 80 columns if the parameters could not be retrieved.
        return 80;
    }

    return std::min<size_t>(120, std::max<size_t>(80, params.ws_col));
}


parser::parser(const std::string& name,
               const std::string& version,
               const std::string& help)
    : m_keys()
    , m_parsers()
    , m_name(name)
    , m_version(version)
    , m_help(help)
{
    // Built-in arguments (aliases are not shown!)
    (*this)["--help"] = new flag(NULL, "Display this message.");
    create_alias("--help", "-help");
    create_alias("--help", "-h");

    (*this)["--version"] = new flag(NULL, "Print the version string.");
    create_alias("--version", "-version");
    create_alias("--version", "-v");

    add_seperator();
}


parser::~parser()
{
    // Collect set of unique pointers, to handle parsers assigned to multiple keys.
    consumer_set pointers;
    for (consumer_map::iterator it = m_parsers.begin(); it != m_parsers.end(); ++it) {
        pointers.insert(it->second);
    }

    for (consumer_set::iterator it = pointers.begin(); it != pointers.end(); ++it) {
        delete *it;
    }
}


consumer_ptr& parser::operator[](const std::string& key)
{
    AR_DEBUG_ASSERT(!key.empty());
    if (m_parsers.find(key) == m_parsers.end()) {
        m_keys.push_back(key_pair(true, key));
    }

    return m_parsers[key];
}


const consumer_ptr& parser::at(const std::string& key) const
{
    return m_parsers.at(key);
}


parse_result parser::parse_args(int argc, char* argv[])
{
    const string_vec argvec(argv + 1, argv + argc);
    string_vec_citer it = argvec.begin();
    while (it != argvec.end()) {
        consumer_map::iterator parser = m_parsers.end();
        if (find_argument(parser, *it)) {
            if (parser->second->is_set()) {
                std::cerr << "WARNING: Command-line option " << parser->first
                          << " has been specified more than once."
                          << std::endl;
            }

            const size_t consumed = parser->second->consume(++it, argvec.end());

            if (consumed == static_cast<size_t>(-1)) {
                if (it != argvec.end()) {
                    std::cerr << "ERROR: Invalid value for " << *(it - 1) << ": '"
                              << *it << "'; aborting ..." << std::endl;
                } else {
                    std::cerr << "ERROR: No value supplied for " << *(it - 1)
                              << "; aborting ..." << std::endl;
                }

                return pr_error;
            }

            it += static_cast<consumer_map::iterator::difference_type>(consumed);
        } else {
            return pr_error;
        }
    }

    if (is_set("--help")) {
        print_help();
        return pr_exit;
    } else if (is_set("--version")) {
        print_version();
        return pr_exit;
    }

    return pr_ok;
}


bool parser::is_set(const std::string& key) const
{
    const consumer_map::const_iterator it = m_parsers.find(key);
    AR_DEBUG_ASSERT(it != m_parsers.end());

    return it->second->is_set();
}


void parser::add_seperator()
{
    m_keys.push_back(key_pair(false, std::string()));
}


void parser::add_header(const std::string& header)
{
    add_seperator();
    m_keys.push_back(key_pair(false, header));
}


void parser::create_alias(const std::string& key, const std::string& alias)
{
    AR_DEBUG_ASSERT(m_parsers.find(alias) == m_parsers.end());

    m_parsers[alias] = m_parsers.at(key);
}


void parser::print_version() const
{
    std::cerr << m_name << " " << m_version << std::endl;
}


void parser::print_help() const
{
    print_version();
    std::cerr <<"\n" << m_help << "\n\n";

    print_arguments(m_keys);
}


void parser::print_arguments(const key_pair_vec& keys) const
{
    typedef key_pair_vec::const_iterator key_pair_citer;

    size_t indentation = 0;
    for (key_pair_citer it = keys.begin(); it != keys.end(); ++it) {
        if (it->first) {
            const consumer_ptr ptr = m_parsers.at(it->second);
            const std::string metavar = get_metavar_str(ptr, it->second);
            size_t current_len = it->second.length();

            if (!metavar.empty()) {
                current_len += metavar.length();
            }

            // For simplicity, we always include the space for the metavar
            indentation = std::max<size_t>(indentation, current_len + 1);
        }
    }

    // indentation + 4 space before description
    indentation = 2 + indentation + 4;

    std::cerr << std::left << std::setw(indentation)
                  << "Arguments:" << "Description:\n";

    cli_formatter fmt;
    fmt.set_ljust(2);  // Indent subsequent lines two spaces
    fmt.set_indent(indentation);
    fmt.set_indent_first_line(false);
    fmt.set_column_width(get_terminal_columns() - indentation - 3);

    for (key_pair_citer it = keys.begin(); it != keys.end(); ++it) {
        if (!it->first) {
            std::cerr << it->second << "\n";
            continue;
        }

        const consumer_ptr ptr = m_parsers.at(it->second);
        if (ptr->help() == "HIDDEN") {
            continue;
        }

        const std::string metavar = get_metavar_str(ptr, it->second);
        std::cerr << std::left << std::setw(indentation)
                  << ("  " + it->second + " " + metavar);

        std::string value = ptr->help();
        if (!value.empty()) {
            // Replace "%default" with string representation of current value.
            size_t index = value.find("%default");
            if (index != std::string::npos) {
                const std::string default_value = ptr->to_str();
                value.replace(index, 8, default_value);
            }

            // Format into columns and indent lines (except the first line)
            std::cerr << fmt.format(value) << "\n";
        }
    }

    std::cerr << std::endl;
}


bool parser::find_argument(consumer_map::iterator& it, const std::string& str)
{
    it = m_parsers.find(str);
    if (it != m_parsers.end()) {
        return true;
    }

    // Locate partial arguments by finding arguments with 'str' as prefix
    if (str != "-" && str != "--") {
        key_pair_vec matches;
        consumer_map::iterator cit = m_parsers.begin();
        for (; cit != m_parsers.end(); ++cit) {
            if (cit->first.substr(0, str.size()) == str) {
                matches.push_back(key_pair(true, cit->first));
            }
        }

        if (matches.size() == 1) {
            it = m_parsers.find(matches.front().second);
            return true;
        } else if (matches.size() > 1) {
            std::cerr << "ERROR: Ambiguous argument '" << str << "'; "
                      << "Candidate arguments are\n\n";

            print_arguments(matches);

            return false;
        }
    }

    std::cerr << "ERROR: Unknown argument: '" << str << "'" << std::endl;

    return false;
}


std::string parser::get_metavar_str(const consumer_ptr ptr,
                                    const std::string& key) const
{
    if (!ptr->metavar().empty()) {
        return ptr->metavar();
    }

    std::string metavar = key;
    metavar.erase(0, metavar.find_first_not_of("-"));
    std::replace(metavar.begin(), metavar.end(), '-', '_');

    return toupper(metavar);
}


///////////////////////////////////////////////////////////////////////////////

consumer_base::consumer_base(const std::string& metavar,
                             const std::string& help)
    : m_value_set(false)
    , m_metavar(metavar)
    , m_help(help)
{
}


consumer_base::~consumer_base()
{
}


bool consumer_base::is_set() const
{
    return m_value_set;
}


const std::string& consumer_base::metavar() const
{
    return m_metavar;
}


const std::string& consumer_base::help() const
{
    return m_help;
}



///////////////////////////////////////////////////////////////////////////////

flag::flag(bool* value, const std::string& help)
    // Non-empty metavar to avoid auto-generated metavar
    : consumer_base(" ", help)
    , m_ptr(value)
{
    m_value_set = m_ptr ? *value : false;
}


size_t flag::consume(string_vec_citer, const string_vec_citer&)
{
    if (m_ptr) {
        *m_ptr = true;
    }
    m_value_set = true;

    return 0;
}


std::string flag::to_str() const
{
    return std::string(m_value_set ? "on" : "off");
}


///////////////////////////////////////////////////////////////////////////////

any::any(std::string* value, const std::string& metavar, const std::string& help)
    : consumer_base(metavar, help)
    , m_ptr(value)
    , m_sink()
{
}


size_t any::consume(string_vec_citer start, const string_vec_citer& end)
{
    if (start != end) {
        m_value_set = true;
        if (m_ptr) {
            *m_ptr = *start;
        } else {
            m_sink = *start;
        }

        return 1;
    }

    return static_cast<size_t>(-1);
}


std::string any::to_str() const
{
    const std::string& result = m_ptr ? *m_ptr : m_sink;
    if (result.empty()) {
        return "<not set>";
    } else {
        return result;
    }
}


///////////////////////////////////////////////////////////////////////////////

knob::knob(unsigned* value, const std::string& metavar, const std::string& help)
    : consumer_base(metavar, help)
    , m_ptr(value)
{
    AR_DEBUG_ASSERT(m_ptr);
}


size_t knob::consume(string_vec_citer start, const string_vec_citer& end)
{
    if (start != end) {
        std::stringstream stream(*start);
        int64_t temp = 0;

        if (!(stream >> temp)) {
            return static_cast<size_t>(-1);
        }

        // Failing on trailing, non-numerical values
        std::string trailing;
        if (stream >> trailing) {
            return static_cast<size_t>(-1);
        }

        if (temp >= 0 && temp <= std::numeric_limits<unsigned>::max()) {
            *m_ptr = static_cast<unsigned>(temp);
            m_value_set = true;
            return 1;
        }
    }

    return static_cast<size_t>(-1);
}


std::string knob::to_str() const
{
    std::stringstream stream;
    stream << *m_ptr;
    return stream.str();
}


///////////////////////////////////////////////////////////////////////////////

floaty_knob::floaty_knob(double* value, const std::string& metavar, const std::string& help)
    : consumer_base(metavar, help)
    , m_ptr(value)
{
    AR_DEBUG_ASSERT(m_ptr);
}


size_t floaty_knob::consume(string_vec_citer start, const string_vec_citer& end)
{
    if (start != end) {
        std::stringstream stream(*start);
        if (!(stream >> *m_ptr)) {
            return static_cast<size_t>(-1);
        }

        // Failing on trailing, non-numerical values
        std::string trailing;
        if (stream >> trailing) {
            return static_cast<size_t>(-1);
        }

        m_value_set = true;
        return 1;
    }

    return static_cast<size_t>(-1);

}


std::string floaty_knob::to_str() const
{
    std::stringstream stream;
    stream << *m_ptr;
    return stream.str();
}

} // namespace argparse
} // namespace ar
