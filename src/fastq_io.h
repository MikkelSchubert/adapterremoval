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
#include <fstream>

#include <zlib.h>

#ifdef AR_BZIP2_SUPPORT
#include <bzlib.h>
#endif


#include "commontypes.h"
#include "fastq.h"
#include "scheduler.h"
#include "timer.h"
#include "linereader.h"
#include "strutils.h"

namespace ar
{

class userconfig;
class fastq_read_chunk;
class fastq_output_chunk;

typedef std::unique_ptr<fastq_output_chunk> output_chunk_ptr;
typedef std::unique_ptr<fastq_read_chunk> read_chunk_ptr;
typedef std::pair<size_t, unsigned char*> buffer_pair;
typedef std::vector<buffer_pair> buffer_vec;


//! Number of FASTQ records to read for each data-chunk
const size_t FASTQ_CHUNK_SIZE = 2 * 1024;

#if defined(AR_GZIP_SUPPORT) || defined(AR_BZIP2_SUPPORT)
//! Size of compressed chunks used to transport compressed data
const size_t FASTQ_COMPRESSED_CHUNK = 40 * 1024;
#endif




/**
 * Container object for (demultiplexed) reads.
 */
class fastq_read_chunk : public analytical_chunk
{
public:
    /** Create chunk representing lines starting at line offset (1-based). */
    fastq_read_chunk(bool eof_ = false);

    //! Indicates that EOF has been reached.
    bool eof;

    //! Lines read from the mate 1 files
    fastq_vec reads_1;
    //! Lines read from the mate 2 files
    fastq_vec reads_2;
};


/**
 * Container object for processed reads.
 */
class fastq_output_chunk : public analytical_chunk
{
public:
    /** Constructor; does nothing. */
    fastq_output_chunk(bool eof_ = false);

    /** Destructor; frees buffers. */
    ~fastq_output_chunk();

    /** Add FASTQ read, accounting for one or more input reads. */
    void add(const fastq_encoding& encoding, const fastq& read, size_t count = 1);

    //! Indicates that EOF has been reached.
    bool eof;

    //! The number of reads used to generate this chunk; may differ from the
    //! the number of reads, in the case of collapsed reads.
    size_t count;

private:
    friend class gzip_fastq;
    friend class bzip2_fastq;
    friend class write_fastq;

    //! Lines read from the mate 1 and mate 2 files
    string_vec reads;

    //! Buffers of compressed lines
    buffer_vec buffers;
};


/**
 * Simple file reading step.
 *
 * Reads from a single FASTQ file, storing the reads in a fastq_file_chunk.
 * Once the EOF has been reached, a single empty chunk will be returned,
 * marked using the 'eof' property.
 */
class read_single_fastq : public analytical_step
{
public:
    /**
     * Constructor.
     *
     * @param encoding FASTQ encoding for reading quality scores.
     * @param filename Path to FASTQ file containing mate 1 / 2 reads.
     * @param next_step ID of analytical step to which data is forwarded.
     *
     * Opens the input file corresponding to the specified mate.
     */
    read_single_fastq(const fastq_encoding* encoding,
                      const std::string& filename,
                      size_t next_step);

    /** Reads N lines from the input file and saves them in an fastq_read_chunk. */
    virtual chunk_vec process(analytical_chunk* chunk);

    /** Finalizer; checks that all input has been processed. */
    virtual void finalize();

private:
    //! Not implemented
    read_single_fastq(const read_single_fastq&);
    //! Not implemented
    read_single_fastq& operator=(const read_single_fastq&);

    //! Encoding used to parse FASTQ reads.
    const fastq_encoding* m_encoding;
    //! Current line in the input file (1-based)
    size_t m_line_offset;
    //! Line reader used to read raw / gzip'd / bzip2'd FASTQ files.
    line_reader m_io_input;
    //! The analytical step following this step
    const size_t m_next_step;
    //! Used to track whether an EOF block has been received.
    bool m_eof;
};


/**
 * Simple file reading step.
 *
 * Reads from the mate 1 and the mate 2 files, storing the reads in a
 * fastq_file_chunk. Once the EOF has been reached, a single empty chunk will
 * be returned, marked using the 'eof' property.
 */
class read_paired_fastq : public analytical_step
{
public:
    /**
     * Constructor.
     */
    read_paired_fastq(const fastq_encoding* encoding,
                      const std::string& filename_1,
                      const std::string& filename_2,
                      size_t next_step);

    /** Reads N lines from the input file and saves them in an fastq_file_chunk. */
    virtual chunk_vec process(analytical_chunk* chunk);

    /** Finalizer; checks that all input has been processed. */
    virtual void finalize();

private:
    //! Not implemented
    read_paired_fastq(const read_paired_fastq&);
    //! Not implemented
    read_paired_fastq& operator=(const read_paired_fastq&);

    //! Encoding used to parse FASTQ reads.
    const fastq_encoding* m_encoding;
    //! Current line in the input file (1-based)
    size_t m_line_offset;
    //! Line reader used to read raw / gzip'd / bzip2'd FASTQ files.
    line_reader m_io_input_1;
    //! Line reader used to read raw / gzip'd / bzip2'd FASTQ files.
    line_reader m_io_input_2;
    //! The analytical step following this step
    const size_t m_next_step;
    //! Used to track whether an EOF block has been received.
    bool m_eof;
};


/**
 * Simple file reading step.
 *
 * Reads from an FASTQ file containing interleaved FASTQ reads, storing the
 * reads in a fastq_file_chunk. Once the EOF has been reached, a single empty
 * chunk will be returned, marked using the 'eof' property.
 */
class read_interleaved_fastq : public analytical_step
{
public:
    /**
     * Constructor.
     */
    read_interleaved_fastq(const fastq_encoding* encoding,
                           const std::string& filename,
                           size_t next_step);

    /** Reads N lines from the input file and saves them in an fastq_file_chunk. */
    virtual chunk_vec process(analytical_chunk* chunk);

    /** Finalizer; checks that all input has been processed. */
    virtual void finalize();

private:
    //! Not implemented
    read_interleaved_fastq(const read_interleaved_fastq&);
    //! Not implemented
    read_interleaved_fastq& operator=(const read_interleaved_fastq&);

    //! Encoding used to parse FASTQ reads.
    const fastq_encoding* m_encoding;
    //! Current line in the input file (1-based)
    size_t m_line_offset;
    //! Line reader used to read raw / gzip'd / bzip2'd FASTQ files.
    line_reader m_io_input;
    //! The analytical step following this step
    const size_t m_next_step;
    //! Used to track whether an EOF block has been received.
    bool m_eof;
};



#ifdef AR_BZIP2_SUPPORT
/**
 * BZip2 compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class bzip2_fastq : public analytical_step
{
public:
    /** Constructor; 'next_step' sets the destination of compressed chunks. */
    bzip2_fastq(const userconfig& config, size_t next_step);

    /** Compresses input lines, saving compressed chunks to chunk->buffers. */
    virtual chunk_vec process(analytical_chunk* chunk);

    /** Checks that all input has been processed and frees stream. */
    virtual void finalize();

private:
    //! Not implemented
    bzip2_fastq(const bzip2_fastq&);
    //! Not implemented
    bzip2_fastq& operator=(const bzip2_fastq&);

    //! N reads which did not result in an output chunk
    size_t m_buffered_reads;
    //! The analytical step following this step
    const size_t m_next_step;
    //! BZip2 stream object
    bz_stream m_stream;
    //! Used to track whether an EOF block has been received.
    bool m_eof;
};

#endif


#ifdef AR_GZIP_SUPPORT
/**
 * GZip compression step; takes any lines in the input chunk, compresses them,
 * and adds them to the buffer list of the chunk, before forwarding it. */
class gzip_fastq : public analytical_step
{
public:
    /** Constructor; 'next_step' sets the destination of compressed chunks. */
    gzip_fastq(const userconfig& config, size_t next_step);

    /** Compresses input lines, saving compressed chunks to chunk->buffers. */
    virtual chunk_vec process(analytical_chunk* chunk);

    /** Checks that all input has been processed and frees stream. */
    virtual void finalize();

private:
    //! Not implemented
    gzip_fastq(const gzip_fastq&);
    //! Not implemented
    gzip_fastq& operator=(const gzip_fastq&);

    //! N reads which did not result in an output chunk
    size_t m_buffered_reads;
    //! The analytical step following this step
    const size_t m_next_step;
    //! GZip stream object
    z_stream m_stream;
    //! Used to track whether an EOF block has been received.
    bool m_eof;
};
#endif


/**
 * Simple file reading step.
 *
 * The 'process' function takes a fastq_file_chunk object and writes the lines
 * at the offset corresponding to the 'type' argument to the corresponding
 * output file. The list of lines is cleared upon writing.
 */
class write_fastq : public analytical_step
{
public:
    /**
     * Constructor.
     *
     * @param filename Filename to which FASTQ reads are written.
     *
     * Based on the read-type specified, and SE / PE mode, the corresponding
     * output file is opened
     */
    write_fastq(const std::string& filename);

    /** Writes the reads of the type specified in the constructor. */
    virtual chunk_vec process(analytical_chunk* chunk);

    /** Flushes the output file and prints progress report (if enabled). */
    virtual void finalize();

private:
    //! Pointer to output file opened using userconfig::open_with_default_filename.
    std::ofstream m_output;

    //! Used to track whether an EOF block has been received.
    bool m_eof;
};

} // namespace ar

#endif
