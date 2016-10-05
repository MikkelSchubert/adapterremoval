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
#ifndef TRIMMED_READS_H
#define TRIMMED_READS_H

#include <string>

#include "fastq_io.h"
#include "scheduler.h"
#include "statistics.h"


class fastq;
class userconfig;


namespace ar
{


//! Enum representing the possible states of read processing
enum read_status
{
    //! Read passed all checks, and should be written to the main output file
    PASSED,
    //! Read failed one or more checks, and should be discarded; this may
    //! include stripping the sequence and qualities, and flagging the read
    FAILED,
};


/**
 * Helper-class for directing trimmed reads to the next step in the pipeline.
 */
class trimmed_reads
{
public:
    /**
     * Constructor.
     *
     * @param config Global User-config instance; must outlive instance.
     * @param config offset The file-offset for the reads being processed.
     * @param config eof If true, this chunk of reads are at the EOF.
     */
    trimmed_reads(const userconfig& config, size_t offset, bool eof);

    /**
     * Encodes and caches the specified mate 1 read.
     *
     * @param read Processed FASTQ read.
     * @param state FAILED or PASSED; determines destination file.
     * @param read_count Number of actual reads represented by 'read'; for
     *                   merged sequences this may be > 1.
     *
     * Note that 'read' may be modified (truncated) depending on user options.
     */
    void add_mate_1_read(fastq& read, read_status state, size_t read_count = 1);
    /**
     * Encodes and caches the specified mate 2 read.
     *
     * @param read Processed FASTQ read.
     * @param state FAILED or PASSED; determines destination file.
     * @param read_count Number of actual reads represented by 'read'; for
     *                   merged sequences this may be > 1.
     *
     * Note that 'read' may be modified (truncated) depending on user options.
     */
    void add_mate_2_read(fastq& read, read_status state, size_t read_count = 1);

    /**
     * Encodes and caches the specified pair of reads.
     *
     * @param read_1 Processed mate 1 FASTQ read.
     * @param state_1 FAILED or PASSED; determines destination file for mate 1.
     * @param read_2 Processed mate 2 FASTQ read.
     * @param state_2 FAILED or PASSED; determines destination file for mate 2.
     *
     * Note that 'read' may be modified (truncated) depending on user options.
     */
    void add_pe_reads(fastq& read_1, read_status state_1,
                      fastq& read_2, read_status state_2);

    /**
     * Encodes and caches the specified collapsed read.
     *
     * @param read Processed FASTQ read.
     * @param state FAILED or PASSED; determines destination file.
     * @param read_count Number of actual reads represented by 'read'; for
     *                   merged sequences this may be > 1.
     *
     * Note that 'read' may be modified (truncated) depending on user options.
     */
    void add_collapsed_read(fastq& read, read_status, size_t read_count = 1);
    /**
     * Encodes and caches the specified collapsed, truncated read.
     *
     * @param read Processed FASTQ read.
     * @param state FAILED or PASSED; determines destination file.
     * @param read_count Number of actual reads represented by 'read'; for
     *                   merged sequences this may be > 1.
     *
     * Note that 'read' may be modified (truncated) depending on user options.
     */
    void add_collapsed_truncated_read(fastq& read, read_status state, size_t read_count = 1);

    /** Returns vector of chunks from all cached reads. */
    chunk_vec finalize();

private:
    /*
     * Helper function; assigns a given read to a cache depending on state and
     * user settings, it's state (state_1) and the state of its mate (state_2).
     *
     * If interleaved output is enabled, reads are typically (depending on
     * state, etc.) written to 'regular', and are otherwise written to
     * 'interleaved'.
     */
    void distribute_read(output_chunk_ptr& regular,
                         output_chunk_ptr& interleaved,
                         fastq& read,
                         read_status state_1,
                         read_status state_2,
                         size_t read_count = 1);

    //! User configuration; must outlive instance.
    const userconfig& m_config;
    //! Output-encoding used to write reads.
    const fastq_encoding& m_encoding;

    //! The offset of this chunk of reads.
    size_t m_offset;

    //! Pointer to cached mate 1 reads.
    output_chunk_ptr m_mate_1;
    //! Pointer to cached mate 2 reads; may be NULL.
    output_chunk_ptr m_mate_2;
    //! Pointer to cached singleton reads; may be NULL.
    output_chunk_ptr m_singleton;
    //! Pointer to cached collapsed reads; may be NULL.
    output_chunk_ptr m_collapsed;
    //! Pointer to cached collapsed, truncated reads; may be NULL.
    output_chunk_ptr m_collapsed_truncated;
    //! Pointer to cached discarded reads; may be NULL.
    output_chunk_ptr m_discarded;
};


} // namespace ar

#endif
