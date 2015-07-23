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
#include <stdexcept>
#include <iostream>
#include <cerrno>
#include <cstring>

#include "fastq_io.h"
#include "userconfig.h"


//! Number of lines to read for each data-chunk
const size_t CHUNK_SIZE = 4 * 1000;


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'fastq_file_chunk'

fastq_file_chunk::fastq_file_chunk(size_t offset_)
  : offset(offset_)
  , mates(2)
  , output(rt_max)
{
    typedef std::vector<string_vec>::iterator iter;

    for (iter it = mates.begin(); it != mates.end(); ++it) {
        it->reserve(CHUNK_SIZE);
    }

    for (iter it = output.begin(); it != output.end(); ++it) {
        it->reserve(CHUNK_SIZE);
    }
}


///////////////////////////////////////////////////////////////////////////////
// Implementations for 'read_paired_fastq'

read_paired_fastq::read_paired_fastq(const userconfig& config, read_type mate)
  : analytical_step(analytical_step::ordered, true)
  , m_line_offset(1)
  , m_io_input()
  , m_type(mate)
{
    if (m_type == rt_mate_1) {
        m_io_input = config.open_ifstream(config.input_file_1);
    } else if (m_type == rt_mate_2) {
        m_io_input = config.open_ifstream(config.input_file_2);
    } else {
        throw std::invalid_argument("mate must be rt_mate_1 or rt_mate_2");
    }
}


analytical_chunk* read_paired_fastq::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_file_chunk> file_chunk;
    if (chunk) {
        file_chunk.reset(dynamic_cast<fastq_file_chunk*>(chunk));
    } else {
        file_chunk.reset(new fastq_file_chunk(m_line_offset));
    }

    if (!m_io_input.get()) {
        return NULL;
    }

    string_vec& lines = file_chunk->mates.at(m_type);
    if (lines.size() != CHUNK_SIZE) {
        lines.resize(CHUNK_SIZE);
    }

    string_vec::iterator it = lines.begin();
    try {
        for (; it != lines.end(); ++it) {
            if (!std::getline(*m_io_input, *it)) {
                lines.resize(it - lines.begin());
                break;
            }
        }
    } catch (const std::ios_base::failure&) {
        const size_t offset = m_line_offset + (it - lines.begin());

        std::cerr << "Error reading FASTQ file at line "
                  << offset << "); aborting:\n"
                  << "    " << std::strerror(errno)
                  << std::endl;
        throw;
    }

    if (lines.empty()) {
        // Close the file, but still return a chunk to catch unbalanced files,
        // and to allow downstream nodes to flush (if needed)
        m_io_input.reset();
    }

    m_line_offset += lines.size();

    return file_chunk.release();
}


write_paired_fastq::write_paired_fastq(const userconfig& config, read_type type, bool progress)
  : analytical_step(analytical_step::ordered, true)
  , m_type(type)
  , m_progress(progress)
  , m_output()
  , m_timer("pairs", config.quiet)
{
    switch (m_type) {
        case rt_mate_1:
            if (config.paired_ended_mode) {
                m_output = config.open_with_default_filename("--output1", ".pair1.truncated");
            } else {
                m_output = config.open_with_default_filename("--output1", ".truncated");
            }

            break;

        case rt_mate_2:
            m_output = config.open_with_default_filename("--output2", ".pair2.truncated");
            break;

        case rt_singleton:
            m_output = config.open_with_default_filename("--singleton", ".singleton.truncated");
            break;

        case rt_collapsed:
            m_output = config.open_with_default_filename("--outputcollapsed", ".collapsed");
            break;

        case rt_collapsed_truncated:
            m_output = config.open_with_default_filename("--outputcollapsedtruncated", ".collapsed.truncated");
            break;

        case rt_discarded:
            m_output = config.open_with_default_filename("--discarded", ".discarded");
            break;

        default:
            throw std::invalid_argument("invalid read-type in write_paired_fastq constructor");
    }
}


write_paired_fastq::~write_paired_fastq()
{
}


analytical_chunk* write_paired_fastq::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_file_chunk> file_chunk(dynamic_cast<fastq_file_chunk*>(chunk));
    std::vector<string_vec>& lines = file_chunk->output;

    write_lines(m_output, lines.at(m_type));

    if (m_progress) {
        print_locker lock;
        m_timer.increment(file_chunk->mates.at(0).size() / 4);
    }

    return file_chunk.release();
}


void write_paired_fastq::write_lines(std::auto_ptr<std::ostream>& file, string_vec& lines)
{
    if (lines.empty()) {
        return;
    } else if (!file.get()) {
        throw std::logic_error("write_paired_fastq::write_lines: attempted to write to closed FASTQ file");
    }

    std::ostream& stream = *file;

    for (string_vec::iterator it = lines.begin(); it != lines.end(); ++it) {
        stream << *it;
    }

    lines.clear();
}


void write_paired_fastq::finalize()
{
    if (m_output.get()) {
        m_output->flush();
    }

    if (m_progress) {
        print_locker lock;
        m_timer.finalize();
    }
}
