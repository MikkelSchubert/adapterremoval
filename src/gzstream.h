// ============================================================================
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : gzstream.h
// Revision      : $Revision: 1.5 $
// Revision_date : $Date: 2002/04/26 23:30:15 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
//
// Standard streambuf implementation following Nicolai Josuttis, "The
// Standard C++ Library".
// ============================================================================

#ifndef GZSTREAM_H
#define GZSTREAM_H

#include <iostream>
#include <fstream>

#include <zlib.h>


namespace gzip
{

// ----------------------------------------------------------------------------
// Internal classes to implement gzstream. See below for user classes.
// ----------------------------------------------------------------------------

class gzstreambuf : public std::streambuf
{
private:
    static const int bufferSize = 47+256;    // size of data buff
    // totals 512 bytes under g++ for igzstream at the end.

    gzFile           file;               // file handle for compressed file
    char             buffer[bufferSize]; // data buffer
    bool             opened;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();

public:
    gzstreambuf();
    bool is_open() const;
    gzstreambuf* open(const char* name, int open_mode, int level);
    gzstreambuf* close();
    ~gzstreambuf();

    virtual int     overflow(int c = EOF);
    virtual int     underflow();
    virtual int     sync();
};


class gzstreambase : virtual public std::ios
{
protected:
    gzstreambuf buf;
public:
    gzstreambase();
    gzstreambase(const char* name, int open_mode, int level);
    ~gzstreambase();
    void open(const char* name, int open_mode, int level);
    void close();
    gzstreambuf* rdbuf();
};

// ----------------------------------------------------------------------------
// User classes. Use igzstream and ogzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the gz*
// function interface of the zlib. Files are compatible with gzip compression.
// ----------------------------------------------------------------------------

class igzstream : public gzstreambase, public std::istream
{
public:
    igzstream();
    igzstream(const char* name, int open_mode = std::ios::in);
    gzstreambuf* rdbuf();
    void open(const char* name, int open_mode = std::ios::in);
    bool is_open() const;
};


class ogzstream : public gzstreambase, public std::ostream
{
public:
    ogzstream();
    ogzstream(const char* name, int mode = std::ios::out, int level = Z_DEFAULT_COMPRESSION);

    gzstreambuf* rdbuf();
    void open(const char* name, int open_mode = std::ios::out, int level = Z_DEFAULT_COMPRESSION);
    bool is_open() const;
};


///////////////////////////////////////////////////////////////////////////////
// Implementations

inline gzstreambuf::gzstreambuf()
  : opened(false)
{
    setp( buffer, buffer + (bufferSize-1));
    setg( buffer + 4,     // beginning of putback area
          buffer + 4,     // read position
          buffer + 4);    // end position
    // ASSERT: both input & output capabilities will not be used together
}


inline bool gzstreambuf::is_open() const
{
    return opened;
}


inline gzstreambuf::~gzstreambuf()
{
    close();
}


inline gzstreambase::gzstreambase()
{
    init(&buf);
}


inline gzstreambuf* gzstreambase::rdbuf()
{
    return &buf;
}


inline igzstream::igzstream()
  : std::istream( &buf)
{
}


inline igzstream::igzstream(const char* name, int open_mode)
  : gzstreambase(name, open_mode, Z_DEFAULT_COMPRESSION)
  , std::istream(&buf)
{
}


inline gzstreambuf* igzstream::rdbuf()
{
    return gzstreambase::rdbuf();
}


inline void igzstream::open(const char* name, int open_mode)
{
    gzstreambase::open(name, open_mode, Z_DEFAULT_COMPRESSION);
}


inline bool igzstream::is_open() const
{
    return buf.is_open();
}


inline ogzstream::ogzstream()
  : std::ostream(&buf)
{

}


inline ogzstream::ogzstream(const char* name, int mode, int level)
  : gzstreambase(name, mode, level)
  , std::ostream(&buf)
{
}


inline gzstreambuf* ogzstream::rdbuf()
{
    return gzstreambase::rdbuf();
}


inline void ogzstream::open(const char* name, int open_mode, int level)
{
    gzstreambase::open(name, open_mode, level);
}

inline bool ogzstream::is_open() const
{
    return buf.is_open();
}

}

#endif
