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
#ifndef FASTQ_IO_H
#define FASTQ_IO_H

#include <vector>
#include <memory>

#include <zlib.h>

#ifdef AR_BZIP2_SUPPORT
#include <bzlib.h>
#endif


#include "commontypes.h"
#include "fastq.h"
#include "scheduler.h"
#include "timer.h"
#include "linereader.h"


class userconfig;

typedef std::pair<size_t, unsigned char*> buffer_pair;
typedef std::vector<buffer_pair> buffer_vec;


/**
 * Container object for raw reads.
 */
class fastq_file_chunk : public analytical_chunk
{
public:
    /** Create chunk representing lines starting at line offset (1-based). */
    fastq_file_chunk(size_t offset_);

    //! Indicates that EOF has been reached.
    bool eof;

    //! The line-number offset from which the lines start
    size_t offset;

    //! Lines read from the mate 1 and mate 2 files
    std::vector<string_vec> mates;
};


/**
 * Container object for processed reads.
 */
class fastq_output_chunk : public analytical_chunk
{
public:
    /** Constructor; does nothing. */
    fastq_output_chunk(bool eof);

    /** Destructor; frees buffers. */
    ~fastq_output_chunk();

    /** Add FASTQ read, accounting for one or more input reads. */
    void add(const fastq_encoding& encoding, const fastq& read, size_t count = 1);

    //! Indicates that EOF has been reached.
    const bool eof;

private:
    friend class gzip_paired_fastq;
    friend class bzip2_paired_fastq;
    friend class write_paired_fastq;

    //! The number of reads used to generate this chunk; may differ from the
    //! the number of reads, in the case of collapsed reads.
    size_t count;

    //! Lines read from the mate 1 and mate 2 files
    string_vec reads;

    //! Buffers of compressed lines
    buffer_vec buffers;
};


/**
 * Simple file reading step.
 *
 * Reads from either the mate 1 or the mate 2 file, storing the reads in the
 * mates variable of a fastq_file_chunk, using the index corresponding to
 * either rt_mate_1 or rt_mate_2. Once the EOF has been reached, a single
 * empty of lines will be returned.
 *
 * The class will re-use existing fastq_file_chunk objects passed to the
 * 'process' function, resizing the list of lines as nessesary to match the
 * number of lines read.
 */
class read_paired_fastq : public analytical_step
{
public:
    /**
     * Constructor.
     *
     * @param config User settings; needed for 'open_ifstream'.
     * @param mate Either rt_mate_1 or rt_mate_2; other values throw.
     *
     * Opens the input file corresponding to the specified mate.
     */
    read_paired_fastq(const userconfig& config, read_type mate, size_t next_step);

    /** Reads N lines from the input file and saves them in an fastq_file_chunk. */
    virtual chunk_list process(analytical_chunk* chunk);

private:
    static std::string get_filename(const userconfig& config, read_type mate);

    //! Current line in the input file (1-based)
    size_t m_line_offset;
    //! Pointer to iostream opened using userconfig::open_ifstream
    line_reader m_io_input;
    //! Read type; either rt_mate_1 or rt_mate_2.
    const read_type m_type;
    //! The analytical step following this step
    const size_t m_next_step;
};


#ifdef AR_BZIP2_SUPPORT
/**
 * BZip2 compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class bzip2_paired_fastq : public analytical_step
{
public:
    /** Constructor; 'next_step' sets the destination of compressed chunks. */
    bzip2_paired_fastq(const userconfig& config, size_t next_step);

    /** Destructor; frees BZip2 stream. */
    virtual ~bzip2_paired_fastq();

    /** Compresses input lines, saving compressed chunks to chunk->buffers. */
    virtual chunk_list process(analytical_chunk* chunk);

private:
    //! The analytical step following this step
    const size_t m_next_step;

    //! BZip2 stream object
    bz_stream m_stream;
};

#endif


#ifdef AR_GZIP_SUPPORT
/**
 * GZip compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class gzip_paired_fastq : public analytical_step
{
public:
    /** Constructor; 'next_step' sets the destination of compressed chunks. */
    gzip_paired_fastq(const userconfig& config, size_t next_step);

    /** Destructor; frees GZip stream. */
    virtual ~gzip_paired_fastq();

    /** Compresses input lines, saving compressed chunks to chunk->buffers. */
    virtual chunk_list process(analytical_chunk* chunk);

private:
    //! The analytical step following this step
    const size_t m_next_step;

    //! GZip stream object
    z_stream m_stream;
};
#endif


/**
 * Simple file reading step.
 *
 * The 'process' function takes a fastq_file_chunk object and writes the lines
 * at the offset corresponding to the 'type' argument to the corresponding
 * output file. The list of lines is cleared upon writing.
 */
class write_paired_fastq : public analytical_step
{
public:
    /**
     * Constructor.
     *
     * @param config User settings.
     * @param read_type The type of reads to write.
     *
     * Based on the read-type specified, and SE / PE mode, the corresponding
     * output file is opened
     */
    write_paired_fastq(const userconfig& config, read_type type);

    /** Destructor; closes output file. */
    ~write_paired_fastq();

    /** Writes the reads of the type specified in the constructor. */
    virtual chunk_list process(analytical_chunk* chunk);

    /** Flushes the output file and prints progress report (if enabled). */
    virtual void finalize();

private:
    //! The read type written by this instance.
    const read_type m_type;
    //! Pointer to output file opened using userconfig::open_with_default_filename.
    std::auto_ptr<std::ostream> m_output;

    //! Specifies if progress reports are to be printed
    bool m_progress;
};


#endif
