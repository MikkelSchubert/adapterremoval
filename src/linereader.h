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
#ifndef GZFILE_H
#define GZFILE_H

#include <string>
#include <zlib.h>


/** Represents errors during basic IO. */
class io_error : public std::ios_base::failure
{
public:
    io_error(const std::string& message);
};


/** Represents errors during GZip (de)compression. */
class gzip_error : public io_error
{
public:
    gzip_error(const std::string& message);
};


/**
 * Simple line reader.
 *
 * Currently reads
 *  - uncompressed files
 *  - gzip compressed files
 *
 * Errors are reported using either 'io_error' or 'gzip_error'.
 */
class line_reader
{
public:
    /** Basic constructor; call 'open' afterwards. */
    line_reader();

    /** Directly opens file; see 'open'. */
    line_reader(const std::string& fpath);

    /** Closes the file, if still open. */
    ~line_reader();

    /** Opens file; throws io_error or gzip_error on failure. */
    void open(const std::string& fpath);

    /** Reads a lien into dst, returning false on EOF. */
    bool getline(std::string& dst);

    /** Closes the file, if still open. */
    void close();

    /** Returns true if the file has been closed. */
    bool is_open() const;

    /** Returns true if a read after EOF has been attempted. */
    bool eof() const;

private:
    //! Not implemented
    line_reader(const line_reader&);
    //! Not implemented
    line_reader& operator=(const line_reader&);

    //! Refills 'm_buffer' and sets 'm_buffer_ptr' and 'm_buffer_end'.
    void refill_buffer();

    //! GZip file object; currently this takes care of reading gz files
    gzFile m_file;

    //! Pointer to buffer of decompressed data.
    char* m_buffer;
    //! Pointer to current location in input buffer.
    char* m_buffer_ptr;
    //! Pointer to end of current buffer.
    char* m_buffer_end;

    //! Indicates if a read across the EOF has been attempted.
    bool m_eof;
};


#endif
