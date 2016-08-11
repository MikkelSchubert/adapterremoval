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
    trimmed_reads(const userconfig& config, size_t offset, bool eof);

    void add_mate_1_read(fastq& read, read_status state, size_t read_count = 1);
    void add_mate_2_read(fastq& read, read_status state, size_t read_count = 1);

    /**
     *
     *
     *
     *
     *
     * Note that each read is assumed to represent one read (read_count = 1).
     */
    void add_pe_reads(fastq& read_1, read_status state_1,
                      fastq& read_2, read_status state_2);

    void add_collapsed_read(fastq& read, read_status, size_t read_count = 1);
    void add_collapsed_truncated_read(fastq& read, read_status state, size_t read_count = 1);

    chunk_vec finalize();

private:
    void distribute_read(output_chunk_ptr& regular,
                         output_chunk_ptr& interleaved,
                         fastq& read,
                         read_status state_1,
                         read_status state_2,
                         size_t read_count = 1);

    const userconfig& m_config;
    const fastq_encoding& m_encoding;

    //! The offset 
    size_t m_offset;

    output_chunk_ptr m_mate_1;
    output_chunk_ptr m_mate_2;
    output_chunk_ptr m_singleton;
    output_chunk_ptr m_collapsed;
    output_chunk_ptr m_collapsed_truncated;
    output_chunk_ptr m_discarded;
};


} // namespace ar

#endif
