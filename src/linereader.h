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
#include <cstdio>

#ifdef AR_GZIP_SUPPORT
#include <zlib.h>
#endif

#ifdef AR_BZIP2_SUPPORT
#include <bzlib.h>
#endif

namespace ar
{

/** Represents errors during basic IO. */
class io_error : public std::ios_base::failure
{
public:
    io_error(const std::string& message, int error_number = 0);
};


/** Represents errors during GZip (de)compression. */
class gzip_error : public io_error
{
public:
    gzip_error(const std::string& message, const char* gzip_msg = NULL);
};


/** Represents errors during BZip2 (de)compression. */
class bzip2_error : public io_error
{
public:
    bzip2_error(const std::string& message);
};


/** Base-class for line reading; used by recievers. */
class line_reader_base
{
public:
    /** Does nothing. */
    line_reader_base();

    /** Closes the file, if still open. */
    virtual ~line_reader_base();

    /** Reads a lien into dst, returning false on EOF. */
    virtual bool getline(std::string& dst) = 0;
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
class line_reader : public line_reader_base
{
public:
    /** Constructor; opens file and throws on errors. */
    line_reader(const std::string& fpath);

    /** Closes the file, if still open. */
    ~line_reader();

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
    void refill_buffers();

    //! Raw file used to read input.
    FILE* m_file;
    /** Refills 'm_raw_buffer'; sets 'm_raw_buffer_ptr' and 'm_raw_buffer_end'. */
    void refill_raw_buffer();
    /** Points 'm_buffer' and other points to corresponding 'm_raw_buffer's. */
    void refill_buffers_uncompressed();


#ifdef AR_GZIP_SUPPORT
    //! GZip stream pointer; used if input it detected to be gzip compressed.
    z_stream* m_gzip_stream;
#endif

    /** Returns true if the raw buffer contains gzip'd data. */
    bool identify_gzip() const;
    /** Initializes gzip stream and output buffers. */
    void initialize_buffers_gzip();
    /** Refills 'm_buffer' from compressed data; may refill raw buffers. */
    void refill_buffers_gzip();
    /** Closes gzip buffers and frees assosiated memory. */
    void close_buffers_gzip();


#ifdef AR_BZIP2_SUPPORT
    //! GZip stream pointer; used if input it detected to be gzip compressed.
    bz_stream* m_bzip2_stream;
#endif

    /** Returns true if the raw buffer contains bzip2'd data. */
    bool identify_bzip2() const;
    /** Initializes bzip2 stream and output buffers. */
    void initialize_buffers_bzip2();
    /** Refills 'm_buffer' from compressed data; may refill raw buffers. */
    void refill_buffers_bzip2();
    /** Closes gzip2 buffers and frees assosiated memory. */
    void close_buffers_bzip2();


    //! Pointer to buffer of decompressed data.
    char* m_buffer;
    //! Pointer to current location in input buffer.
    char* m_buffer_ptr;
    //! Pointer to end of current buffer.
    char* m_buffer_end;

    //! Pointer to buffer of raw data.
    char* m_raw_buffer;
    //! Pointer to end of current raw buffer.
    char* m_raw_buffer_end;

    //! Indicates if a read across the EOF has been attempted.
    bool m_eof;
};


///////////////////////////////////////////////////////////////////////////////

inline line_reader_base::line_reader_base()
{
}


inline line_reader_base::~line_reader_base()
{
}

} // namespace ar

#endif
