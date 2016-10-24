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
#include "trimmed_reads.h"
#include "userconfig.h"



namespace ar
{

inline void add_chunk(chunk_vec& chunks, size_t target, output_chunk_ptr chunk)
{
    if (chunk.get()) {
        chunks.push_back(chunk_pair(target, std::move(chunk)));
    }
}


trimmed_reads::trimmed_reads(const userconfig& config, size_t offset, bool eof)
    : m_config(config)
    , m_encoding(*config.quality_output_fmt)
    , m_offset(offset)
    , m_mate_1()
    , m_mate_2()
    , m_singleton()
    , m_collapsed()
    , m_collapsed_truncated()
    , m_discarded()
{
    m_mate_1.reset(new fastq_output_chunk(eof));
    if (config.paired_ended_mode && !config.interleaved_output) {
        m_mate_2.reset(new fastq_output_chunk(eof));
    }

    if (!config.combined_output) {
        m_discarded.reset(new fastq_output_chunk(eof));

        if (config.paired_ended_mode) {
            m_singleton.reset(new fastq_output_chunk(eof));
        }

        if (config.collapse) {
            m_collapsed.reset(new fastq_output_chunk(eof));
            m_collapsed_truncated.reset(new fastq_output_chunk(eof));
        }
    }
}


void trimmed_reads::add_mate_1_read(fastq& read, read_status state,
                                    size_t read_count)
{
    // Single end reads always go into the mate 1 file or the discarded file
    distribute_read(m_mate_1, m_mate_1, read, state, PASSED, read_count);
}


void trimmed_reads::add_mate_2_read(fastq& read, read_status state,
                                    size_t read_count)
{
    // Single end reads always go into the mate 2 file or the discarded file
    distribute_read(m_mate_2, m_mate_1, read, state, PASSED, read_count);
}


void trimmed_reads::add_pe_reads(fastq& read_1, read_status state_1,
                                 fastq& read_2, read_status state_2)
{
    distribute_read(m_mate_1, m_mate_1, read_1, state_1, state_2);
    distribute_read(m_mate_2, m_mate_1, read_2, state_2, state_1);
}


void trimmed_reads::add_collapsed_read(fastq& read,
                                       read_status state,
                                       size_t read_count)
{
    output_chunk_ptr& destination = m_config.combined_output ? m_mate_1 : m_collapsed;

    // Collapsed reads may go into the mate 1, mate 2, or discard file
    distribute_read(destination, destination, read, state, PASSED, read_count);
}


void trimmed_reads::add_collapsed_truncated_read(fastq& read,
                                                 read_status state,
                                                 size_t read_count)
{
    output_chunk_ptr& destination = m_config.combined_output ? m_mate_1 : m_collapsed_truncated;

    // Collapsed tr. reads may go into the mate 1, mate 2, or discard file
    distribute_read(destination, destination, read, state, PASSED, read_count);
}


chunk_vec trimmed_reads::finalize()
{
    chunk_vec chunks;

    add_chunk(chunks, m_offset + ai_write_mate_1, std::move(m_mate_1));
    add_chunk(chunks, m_offset + ai_write_mate_2, std::move(m_mate_2));
    add_chunk(chunks, m_offset + ai_write_singleton, std::move(m_singleton));
    add_chunk(chunks, m_offset + ai_write_collapsed, std::move(m_collapsed));
    add_chunk(chunks, m_offset + ai_write_collapsed_truncated, std::move(m_collapsed_truncated));
    add_chunk(chunks, m_offset + ai_write_discarded, std::move(m_discarded));

    return chunks;
}


void trimmed_reads::distribute_read(output_chunk_ptr& regular,
                                    output_chunk_ptr& interleaved,
                                    fastq& read,
                                    read_status state_1,
                                    read_status state_2,
                                    size_t read_count)
{
    if (state_1 == PASSED) {
        if (state_2 == PASSED || m_config.combined_output) {
            if (m_config.interleaved_output) {
                interleaved->add(m_encoding, read, read_count);
            } else {
                regular->add(m_encoding, read, read_count);
            }
        } else {
            m_singleton->add(m_encoding, read, read_count);
        }
    } else if (m_config.combined_output) {
        read.discard();

        if (m_config.interleaved_output) {
            interleaved->add(m_encoding, read, read_count);
        } else {
            regular->add(m_encoding, read, read_count);
        }
    } else {
        m_discarded->add(m_encoding, read, read_count);
    }
}


} // namespace ar
