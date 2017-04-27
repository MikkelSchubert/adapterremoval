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
#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <sstream>
#include <vector>
#include <map>

#include "commontypes.h"

namespace ar
{
namespace argparse
{

class consumer_base;
typedef consumer_base* consumer_ptr;
typedef std::map<std::string, consumer_ptr> consumer_map;


//! Parse results for command-line arguments
enum parse_result {
    //! Terminate now (e.g. --version or --help used)
    pr_exit,
    //! Error occurred parsing arguments / invalid combination of args
    pr_error,
    //! No errors parsing commandline arguments
    pr_ok
};


/**
 * Simple type-safe parsing of command-line options.
 *
 * To handle a an argument, an object of type consumer_base is assigned
 * to the argparse::parser using the [] operator. For example, the parse and
 * save an integer value to a variable when the user supplies the argument
 * '--example', the following is done:
 *
 * int target = 0; // Must be initialized!
 * argparse::parser argparser(...);
 * argparser["--example"] = new argparse::knob(&target);
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
    /**
     * Arguments:
     *   name: Name of the program; used in --version and --help.
     *   version: Version string (excluding the name); is used by the
     *              arguments --help and --version.
     *   help: General help string used by --help; should not include the
     *           parameters themselves, as this is done automatically.
     */
    parser(const std::string& name,
           const std::string& version,
           const std::string& help);

    /** Deletes all (unique) parsers assigned to the set. */
    ~parser();

    /** Parses a set of command-line options as passed to main(argc, argv). */
    parse_result parse_args(int argc, char* argv[]);

    /** Returns true if the option with the given key has been set. */
    bool is_set(const std::string& key) const;


    /**
     * Returns a reference to the pointer for the parser with the given key.
     *
     * If an entry does not exist for the given key, it is created upon access.
     */
    consumer_ptr& operator[](const std::string& key);

    /**
     * Returns a reference to the pointer for the parser with the given key.
     *
     * If an entry does not exist for the given key, out_of_range is thrown.
     */
    const consumer_ptr& at(const std::string& key) const;


    /** Add a blank line between the previous and the next command. */
    void add_seperator();

    /** Add a blank line and a header between the previous and next command. */
    void add_header(const std::string& header);


    /** Create alias for key with name alias. */
    void create_alias(const std::string& key, const std::string& alias);


    /** Helper function; prints the program name and version string. */
    void print_version() const;
    /** Helper functions; prints the full set of help-text. */
    void print_help() const;

private:
    typedef std::pair<bool, std::string> key_pair;
    typedef std::vector<key_pair> key_pair_vec;

    /** Pretty-print the listed arguments. */
    void print_arguments(const key_pair_vec& keys) const;

    /**
     * Attempt to find argument by similarity.
     *
     * @param it If a match is found, is set to iterator of matching parser,
                 otherwise it is set to m_parsers.end().
     * @param str String containing expected argument.
     * @return True if a match is found, false otherwise.
     *
     * If an exact match or a single partial match is found, 'it' is set and
     * true is returned; otherwise it is set to m_parsers.end() and false is
     * returned. In the case of multiple partial matches, the help string for
     * each candidate is printed.
     */
    bool find_argument(consumer_map::iterator& it, const std::string& str);

    /** Generate metavar from argument, namely uppercase without dashes. */
    std::string get_metavar_str(const consumer_ptr, const std::string&) const;

    //! Vector of keys (command-line options), tracking the order of addition.
    key_pair_vec m_keys;
    //! Map of keys (command-line args) to parser pointers; multiple
    //! keys may be associated with the same pointer.
    consumer_map m_parsers;

    //! Name of the program
    std::string m_name;
    //! Version string for the program (excluding the name)
    std::string m_version;
    //! Help text for the program.
    std::string m_help;
};


/**
 * Base class for argument parsers;
 *
 * Each consumer must implement the consume function, which takes iterators to
 * the arguments following the key for this parser (i.e. not including the
 * --option). These then consume zero or more values, returning the number
 * thus consumed, or (size_t)-1 if the values were missing or invalid.
 */
class consumer_base
{
public:
    /**
     * Base constructor; sets various values used when printing --help.
     *
     *  metavar - Used to represent the input value; if empty,
     *            argparse::parser will use the current key associated with
     *            the parser to generate a metavar.
     *  help    - Help string; the value %default is replaced with the default
     *            value.
     */
    consumer_base(const std::string& metavar, const std::string& help);

    /* Destructor; does nothing in the base class. */
    virtual ~consumer_base();

    /**
     * Attempts to consume a value specified on the command-line; returns the
     * number of values consumed (if any), or -1 if parsing failed (e.g. due to
     * the specified value being of the wrong type).
     *
     * Parameters:
     *   start - Iterator pointing to the value following the command-line
     *           argument, if any remain to be consumed.
     *   end - Iterator past-the-end of the list of command-line arguments.
     */
    virtual size_t consume(string_vec_citer start, const string_vec_citer& end) = 0;

    /** Returns true if the consumer has consumed a value. **/
    virtual bool is_set() const;

    /** Returns the metavariable associated with the consumer. **/
    virtual const std::string& metavar() const;

    /** Returns the help string associated with the consumer. **/
    virtual const std::string& help() const;

    /** Returns the value associated with the consumer as a string. **/
    virtual std::string to_str() const = 0;

protected:
    //! Should be set to true if a value has been consumed in a derived class.
    bool m_value_set;

private:
    //! Not implemented
    consumer_base(const consumer_base&);
    //! Not implemented
    consumer_base& operator=(const consumer_base&);

    //! Stores the metavar associated with the consumer
    std::string m_metavar;
    //! Stores the optional description of default behavior.
    std::string m_help;
};


/**
 * Consumer for boolean values (i.e. flags).
 *
 * Unlike typical consumers, this consumer does not expected a value associated
 * with the command-line argument, but instead sets the associated value to
 * true if the command-line argument is specified one or more times.
 */
class flag : public consumer_base
{
public:
    /**
     * See consumer_base::consumer_base
     *
     * Unlike the base constructor, this class does not take a 'metavar', as
     * no values are consumed during parsing.
     */
    flag(bool* sink = NULL, const std::string& help = "");

    /** See consumer_base::consume */
    virtual size_t consume(string_vec_citer start, const string_vec_citer& end);

    /** See consumer_base::to_str */
    virtual std::string to_str() const;

private:
    //! Not implemented
    flag(const flag&);
    //! Not implemented
    flag& operator=(const flag&);

    //! Optional pointer to storage for boolean value; if NULL, m_value is used.
    bool* m_ptr;
};



/**
 * Consumer for string values (filenames, etc.).
 */
class any : public consumer_base
{
public:
    /**
     * See consumer_base::consumer_base
     */
    any(std::string* sink = NULL, const std::string& metavar = "", const std::string& help = "");

    /** See consumer_base::consume */
    virtual size_t consume(string_vec_citer start, const string_vec_citer& end);

    /** See consumer_base::to_str */
    virtual std::string to_str() const;

private:
    //! Not implemented
    any(const any&);
    //! Not implemented
    any& operator=(const any&);

    //! Optional pointer to storage for string value; if NULL, m_value is used.
    std::string* m_ptr;
    //! Value sink used if a pointer to a sink is not provided.
    std::string m_sink;
};


/**
 * Consumer for multiple string values; consumes values until another option
 * (value starting with '-') is encountered.
 */
class many : public consumer_base
{
public:
    /**
     * See consumer_base::consumer_base
     */
    many(string_vec* sink = NULL, const std::string& metavar = "", const std::string& help = "");

    /** See consumer_base::consume */
    virtual size_t consume(string_vec_citer start, const string_vec_citer& end);

    /** See consumer_base::to_str */
    virtual std::string to_str() const;

private:
    //! Not implemented
    many(const many&);
    //! Not implemented
    many& operator=(const many&);

    //! Optional pointer to storage for string value; if NULL, m_value is used.
    string_vec* m_ptr;
    //! Value sink used if a pointer to a sink is not provided.
    string_vec m_sink;
};


/**
 * Consumer for unsigned integer values.
 *
 * Signed values are rejected. On 32 bit systems, the range of values is
 * limited to 0 .. 2^31 - 1, on 64 bit systems the range is 0 .. 2^64 - 1.
 */
class knob : public consumer_base
{
public:
    /**
     * See consumer_base::consumer_base; a sink must be set.
     */
    knob(unsigned* sink, const std::string& metavar = "", const std::string& help = "");

    /** See consumer_base::consume */
    virtual size_t consume(string_vec_citer start, const string_vec_citer& end);

    /** See consumer_base::to_str */
    virtual std::string to_str() const;

private:
    //! Not implemented
    knob(const knob&);
    //! Not implemented
    knob& operator=(const knob&);

    //! Pointer to storage for unsigned value (required).
    unsigned* m_ptr;
};



/**
 * Consumer for floating point values (doubles).
 */
class floaty_knob : public consumer_base
{
public:
    /**
     * See consumer_base::consumer_base; a sink must be set.
     */
    floaty_knob(double* sink, const std::string& metavar = "", const std::string& help = "");

    /** See consumer_base::consume */
    virtual size_t consume(string_vec_citer start, const string_vec_citer& end);

    /** See consumer_base::to_str */
    virtual std::string to_str() const;

private:
    //! Not implemented
    floaty_knob(const floaty_knob&);
    //! Not implemented
    floaty_knob& operator=(const floaty_knob&);

    //! Pointer to storage for unsigned value (required).
    double* m_ptr;
};

} // namespace argparse
} // namespace ar

#endif
