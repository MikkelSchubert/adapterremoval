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
#include <iomanip>
#include <limits>
#include <set>

#include <sys/ioctl.h>
#include <unistd.h>

#include "argparse.h"

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


std::string toupper(const std::string& str)
{
    std::string news = str;
    for(size_t i = 0; i < news.length(); ++i) {
        const char current = news.at(i);
        if (current >= 'a' && current <= 'z') {
            news.at(i) -= 32;
        }

    }

    return news;
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
    if (key.empty()) {
        throw std::invalid_argument("empty key");
    } else if (m_parsers.find(key) == m_parsers.end()) {
        m_keys.push_back(key);
    }

	return m_parsers[key];
}


const consumer_ptr& parser::at(const std::string& key) const
{
	return m_parsers.at(key);
}


parse_result parser::parse_args(int argc, char* argv[])
{
    const StringVec argvec(argv + 1, argv + argc);
    StringVecConstIter it = argvec.begin();
    while (it != argvec.end()) {
        consumer_map::iterator parser = m_parsers.find(*it);
        if (parser != m_parsers.end()) {
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
            std::cerr << "ERROR: Unknown argument: '" << *it << "'; aborting ..."
                      << std::endl;

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
	return m_parsers.at(key)->is_set();
}


void parser::add_seperator()
{
    m_keys.push_back(std::string());
}


void parser::create_alias(const std::string& key, const std::string& alias)
{
    if (m_parsers.find(alias) != m_parsers.end()) {
        throw std::invalid_argument("Attemping to overwrite existing command-"
                                    "line argument with new alias.");
    }

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

    size_t ljust = 0;
    for (StringVecConstIter it = m_keys.begin(); it != m_keys.end(); ++it) {
        if (it->empty()) {
            continue;
        }

        const consumer_ptr ptr = m_parsers.at(*it);
        const std::string metavar = get_metavar_str(ptr, *it);
        size_t current_len = it->length();

        if (!metavar.empty()) {
            current_len += metavar.length();
        }

        // For simplicity, we always include the space for the metavar
        ljust = std::max<size_t>(ljust, current_len + 1);
    }

    // indentation + ljust + 4 space before description
    ljust = 2 + ljust + 4;

    std::cerr << std::left << std::setw(ljust)
              << "Arguments:" << "Description:\n";

    const size_t max_columns = get_terminal_columns() - 2;
    for (StringVecConstIter it = m_keys.begin(); it != m_keys.end(); ++it) {
        if (it->empty()) {
            std::cerr << "\n";
            continue;
        }

        const consumer_ptr ptr = m_parsers.at(*it);
        if (ptr->help() == "HIDDEN") {
            continue;
        }

        const std::string metavar = get_metavar_str(ptr, *it);
        std::cerr << std::left << std::setw(ljust)
                  << ("  " + *it + " " + metavar);

        std::string value = ptr->help();
        if (!value.empty()) {
            // Replace "%default" with string representation of current value.
            size_t index = value.find("%default");
            if (index != std::string::npos) {
                const std::string default_value = ptr->to_str();
                value.replace(index, 8, default_value);
            }

            index = 0;
            size_t current_width = 0;
            size_t current_ljust = ljust;
            while (index < value.length()) {
                const size_t endpos = value.find(' ', index);
                const std::string substr = value.substr(index, endpos - index);

                if (current_width) {
                    if (current_ljust + current_width + 1 + substr.length() > max_columns) {
                        // Add another 2 characters of indentation
                        current_ljust = ljust + 2;
                        std::cerr << "\n" << std::left
                                  << std::setw(current_ljust)
                                  << "" << substr;
                        current_width = substr.length();
                    } else {
                        std::cerr << " " << substr;
                        current_width += substr.length() + 1;
                    }
                } else {
                    std::cerr << substr;
                    current_width += substr.length();
                }

                index += substr.length() + 1;
            } while (index < value.length());

            std::cerr << "\n";
        }

    }

    std::cerr.flush();
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


size_t flag::consume(StringVecConstIter, const StringVecConstIter&)
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


size_t any::consume(StringVecConstIter start, const StringVecConstIter& end)
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
    if (m_ptr) {
        return *m_ptr;
    }

    return m_sink;
}


///////////////////////////////////////////////////////////////////////////////

knob::knob(unsigned* value, const std::string& metavar, const std::string& help)
    : consumer_base(metavar, help)
    , m_ptr(value)
{
    if (!m_ptr) {
        throw std::invalid_argument("sink is null");
    }
}


size_t knob::consume(StringVecConstIter start, const StringVecConstIter& end)
{
    if (start != end) {
        std::stringstream stream(*start);
        signed long temp = 0;

        stream >> temp;
        if (!stream.fail()) {
            if (temp >= 0 && temp <= std::numeric_limits<unsigned>::max()) {
                *m_ptr = static_cast<unsigned>(temp);
                m_value_set = true;
                return 1;
            }
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
    if (!m_ptr) {
        throw std::invalid_argument("sink is null");
    }
}


size_t floaty_knob::consume(StringVecConstIter start, const StringVecConstIter& end)
{
    if (start != end) {
        std::stringstream stream(*start);
        stream >> *m_ptr;
        if (!stream.fail()) {
            m_value_set = true;
            return 1;
        }
    }

    return static_cast<size_t>(-1);

}


std::string floaty_knob::to_str() const
{
    std::stringstream stream;
    stream << *m_ptr;
    return stream.str();
}

}
