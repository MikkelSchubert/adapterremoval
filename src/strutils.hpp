/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#ifndef STRUTILS_H
#define STRUTILS_H

#include <string>
#include <sstream>
#include <iomanip>

namespace ar
{

const size_t DEFAULT_MAX_COLUMNS = 78;
const size_t DEFAULT_INDENTATION = 4;


/** Uppercases letters in the range a-z */
std::string toupper(const std::string& str);


/** Split text by newlines and add fixed indentation following newlines. */
std::string indent_lines(const std::string& lines, size_t identation = DEFAULT_INDENTATION);


/**
 * Formats text into fixed-width columns.
 *
 * @param value Text representing a single paragraph to be formatted.
 * @param max_width Maximum width of output lines in characters.
 * @param ljust Indent lines after the first line by this amount of characters.
 *
 * Note that all whitespace in the input string is consumed; output words are
 * seperated by a single space, and the terminal line does not end with a
 * newline.
 */
std::string columnize_text(const std::string& value, size_t max_width = DEFAULT_MAX_COLUMNS, size_t ljust = 0);


/**
 * Wrapper around 'indent_lines' and 'columnize_text'.
 *
 * Defaults to 78 character columns, with 4 character indentation, and
 * no ljust, corresponding to function defaults. Unlike simply combining
 * 'indent_lines' and 'columnize_text', the 'format' function preserves
 * newlines, allowing paragraphs to be printed easily.
 */
class cli_formatter
{
public:
    /** Creates formatter using default parameters (see above). */
    cli_formatter();

    /** Include (or not) indentation on the first line of output. */
    cli_formatter& set_indent_first_line(bool value);
    /** Maximum number of columns in each line in output. */
    cli_formatter& set_column_width(size_t value);
    /** Number of spaces to indent line 2+ in each paragraph. */
    cli_formatter& set_ljust(size_t value);
    /** Fixed indentation for each line (see also set_indent_first_line). */
    cli_formatter& set_indent(size_t value);

    /** Formats string using the current settings. */
    std::string format(const std::string& value) const;

    /** Format string using default parameters. */
    static std::string fmt(const std::string& value);

    /**
     * Format string using default parameters, but include prefix on first line
     * and indent subsequent lines using the width of the prefix.
     */
    static std::string fmt(const std::string& prefix, const std::string& value);

private:
    //! Specifices whether or not to indent the first line of output.
    bool m_indent_first;
    //! Number of spaces to indent the 2+ lines in each paragraph.
    size_t m_ljust;
    //! The maximum number of columns (in characters).
    size_t m_columns;
    //! The number of spaces to indent each line (see also m_indent_first_line)
    size_t m_indentation;
};

} // namespace ar

#endif
