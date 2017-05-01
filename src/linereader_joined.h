/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2017 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#ifndef LINEREADER_JOINED_H
#define LINEREADER_JOINED_H

#include <memory>

#include "commontypes.h"
#include "linereader.h"


namespace ar
{


/**
 * Multi-file line-reader
 *
 * Wrapper around line_reader that automatically reads through one or more
 * files, returning the content as a contiguous stream of lines. No assumptions
 * are made about the format of the individual files.
 */
class joined_line_readers : public line_reader_base
{
public:
    /** Creates line-reader over multiple files in the specified order. */
    joined_line_readers(const string_vec& filenames);

    /** Closes any still open files. */
    ~joined_line_readers();

    /**
     * Reads a line from the currently open file; if EOF is encountered, the
     * currently open file is closed and the next file is opened. Returns true
     * if a line was successfully read, or false if no files remain.
     */
    bool getline(std::string& dst);

private:
    /**
     * Open the next file, removes it from the queue, and returns true; returns
     * false if no files remain to be processed.
     */
    bool open_next_file();

    //! Not implemented
    joined_line_readers(const joined_line_readers&);
    //! Not implemented
    joined_line_readers& operator=(const joined_line_readers&);

    //! Files left to read; stored in reverse order.
    string_vec m_filenames;
    //! Currently open file, if any.
    std::unique_ptr<line_reader> m_reader;
};


} // namespace ar

#endif
