/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2016 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "debug.h"
#include "demultiplex.h"
#include "fastq.h"
#include "fastq_io.h"
#include "main.h"
//#include "strutils.h"
#include "userconfig.h"


namespace ar
{

//! Implemented in main_adapter_rm.cc
void add_write_step(const userconfig& config, scheduler& sch, size_t offset,
                    const std::string& name, analytical_step* step);


void write_demultiplex_statistics(std::ofstream& output,
                                  const userconfig& config,
                                  const demultiplex_reads* step)
{
    const demux_statistics stats = step->statistics();
    const size_t total = stats.total();

    output.precision(3);
    output << std::fixed << std::setw(3)
           << "\n\n[Demultiplexing statistics]"
           << "\nName\tBarcode_1\tBarcode_2\tHits\tFraction\n"
           << "unidentified\tNA\tNA\t" << stats.unidentified << "\t"
           << stats.unidentified / static_cast<double>(total) << "\n"
           << "ambiguous\tNA\tNA\t" << stats.ambiguous << "\t"
           << stats.ambiguous / static_cast<double>(total) << "\n";

    const fastq_pair_vec barcodes = config.adapters.get_barcodes();
    for (size_t nth = 0; nth < barcodes.size(); ++nth) {
        const fastq_pair& current = barcodes.at(nth);

        output << config.adapters.get_sample_name(nth) << "\t"
               << current.first.sequence() << "\t";
        if (current.second.length()) {
            output << current.second.sequence() << "\t";
        } else {
            output << "*\t";
        }

        output << stats.barcodes.at(nth) << "\t"
               << stats.barcodes.at(nth) / static_cast<double>(total)
               << "\n";
    }

    output << "*\t*\t*\t" << total << "\t" << 1.0 << std::endl;
}


bool write_demultiplex_settings(const userconfig& config,
                                const demultiplex_reads* step,
                                int nth = -1)
{
    const std::string filename \
        = config.get_output_filename(nth == -1 ? "demux_stats" : "--settings", nth);

    try {
        std::ofstream output(filename.c_str(), std::ofstream::out);

        if (!output.is_open()) {
            std::string message = std::string("Failed to open file '") + filename + "': ";
            throw std::ofstream::failure(message + std::strerror(errno));
        }

        output.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        output << NAME << " " << VERSION
                 << "\nDemultiplexing of ";

        if (config.adapters.barcode_count()) {
            if (config.adapters.get_barcodes().front().second.length()) {
                output << "double-indexed ";
            } else {
                output << "single-indexed ";
            }
        }

        if (config.paired_ended_mode) {
            if (config.interleaved_input) {
                output << "interleaved ";
            }

            output << "paired-end reads";
        } else {
            output << "single-end reads";
        }

        output << "\n\n\n[Demultiplexing]"
               << "\nMaximum mismatches (total): " << config.barcode_mm;

        if (config.paired_ended_mode) {
            output << "\nMaximum mate 1 mismatches: " << config.barcode_mm_r1;
            output << "\nMaximum mate 2 mismatches: " << config.barcode_mm_r2;
        }

        output << "\n\n\n[Demultiplexing samples]"
               << "\nName\tBarcode_1\tBarcode_2\n";

        const fastq_pair_vec barcodes = config.adapters.get_barcodes();
        for (size_t idx = 0; idx < barcodes.size(); ++idx) {
            output << config.adapters.get_sample_name(idx);
            if (static_cast<size_t>(nth) == idx) {
                output << "*";
            }

            const fastq_pair& current = barcodes.at(idx);
            output << "\t" << current.first.sequence();

            if (current.second.length()) {
                output << "\t" << current.second.sequence() << "\n";
            } else {
                output << "\t*\n";
            }
        }

        output << "\n\n[Adapter sequences]";
        if (nth == -1) {
            const fastq_pair_vec adapters = config.adapters.get_raw_adapters();
            size_t adapter_id = 0;
            for (fastq_pair_vec::const_iterator it = adapters.begin(); it != adapters.end(); ++it, ++adapter_id) {
                output << "\nAdapter1[" << adapter_id + 1 << "]: " << it->first.sequence();

                fastq adapter_2 = it->second;
                adapter_2.reverse_complement();
                output << "\nAdapter2[" << adapter_id + 1 << "]: " << adapter_2.sequence() << "\n";
            }

            write_demultiplex_statistics(output, config, step);
        } else {
            const string_pair_vec adapters = config.adapters.get_pretty_adapter_set(nth);
            size_t adapter_id = 0;
            for (string_pair_vec::const_iterator it = adapters.begin(); it != adapters.end(); ++it, ++adapter_id) {
                output << "\nAdapter1[" << adapter_id + 1 << "]: " << it->first;
                output << "\nAdapter2[" << adapter_id + 1 << "]: " << it->second << "\n";
            }
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error writing settings file; aborting:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;
        return false;
    }

    return true;
}


class se_demultiplexed_reads_processor : public analytical_step
{
public:
    se_demultiplexed_reads_processor(const userconfig& config, size_t nth)
      : analytical_step(analytical_step::unordered)
      , m_config(config)
      , m_nth(nth)
    {
    }

    chunk_vec process(analytical_chunk* chunk)
    {
        const size_t offset = m_nth * ai_analyses_offset;
        read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));

        output_chunk_ptr encoded_reads(new fastq_output_chunk(read_chunk->eof));

        for (auto const& read: read_chunk->reads_1) {
            encoded_reads->add(*m_config.quality_output_fmt, read);
        }

        chunk_vec chunks;
        chunks.push_back(chunk_pair(offset + ai_write_mate_1, std::move(encoded_reads)));

        return chunks;
    }

private:
    const userconfig& m_config;
    const size_t m_nth;
};


class pe_demultiplexed_reads_processor : public analytical_step
{
public:
    pe_demultiplexed_reads_processor(const userconfig& config, size_t nth)
      : analytical_step(analytical_step::unordered)
      , m_config(config)
      , m_nth(nth)
    {
    }

    chunk_vec process(analytical_chunk* chunk)
    {
        const size_t offset = m_nth * ai_analyses_offset;
        read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
        AR_DEBUG_ASSERT(read_chunk->reads_1.size() == read_chunk->reads_2.size());

        output_chunk_ptr encoded_reads_1(new fastq_output_chunk(read_chunk->eof));
        output_chunk_ptr encoded_reads_2;
        if (!m_config.interleaved_output) {
            encoded_reads_2.reset(new fastq_output_chunk(read_chunk->eof));
        }

        fastq_vec::iterator it_1 = read_chunk->reads_1.begin();
        fastq_vec::iterator it_2 = read_chunk->reads_2.begin();
        while (it_1 != read_chunk->reads_1.end()) {
            fastq read_1 = *it_1++;
            fastq read_2 = *it_2++;

            encoded_reads_1->add(*m_config.quality_output_fmt, read_1);

            if (m_config.interleaved_output) {
                encoded_reads_1->add(*m_config.quality_output_fmt, read_2);
            } else {
                encoded_reads_2->add(*m_config.quality_output_fmt, read_2);
            }
        }

        chunk_vec chunks;
        chunks.push_back(chunk_pair(offset + ai_write_mate_1, std::move(encoded_reads_1)));
        if (!m_config.interleaved_output) {
            chunks.push_back(chunk_pair(offset + ai_write_mate_2, std::move(encoded_reads_2)));
        }

        return chunks;
    }

private:
    const userconfig& m_config;
    const size_t m_nth;
};


int demultiplex_sequences_se(const userconfig& config)
{
    std::cerr << "Demultiplexing single ended reads ..." << std::endl;

    scheduler sch;
    demultiplex_reads* demultiplexer = NULL;

    try {
        // Step 1: Read input file
        sch.add_step(ai_read_fastq, "read_fastq",
                     new read_single_fastq(config.quality_input_fmt.get(),
                                           config.input_files_1,
                                           ai_demultiplex));

        // Step 2: Parse and demultiplex reads based on single or double indices
        sch.add_step(ai_demultiplex, "demultiplex_se",
                     demultiplexer = new demultiplex_se_reads(&config));

        add_write_step(config, sch, ai_write_unidentified_1, "unidentified",
                       new write_fastq(config.get_output_filename("demux_unknown")));

        // Step 3 - N: Trim and write demultiplexed reads
        for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
            const size_t offset = nth * ai_analyses_offset;
            const std::string& sample = config.adapters.get_sample_name(nth);

            sch.add_step(offset + ai_trim_se, "process_se_" + sample,
                         new se_demultiplexed_reads_processor(config, nth));

            add_write_step(config, sch, offset + ai_write_mate_1, sample + "_fastq",
                           new write_fastq(config.get_output_filename("--output1", nth)));
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;
        return 1;
    }

    if (!sch.run(config.max_threads)) {
        return 1;
    } else if (!write_demultiplex_settings(config, demultiplexer)) {
        return 1;
    }

    for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
        if (!write_demultiplex_settings(config, demultiplexer, nth)) {
            return 1;
        }
    }

    return 0;
}


int demultiplex_sequences_pe(const userconfig& config)
{
    std::cerr << "Demultiplexing paired end reads ..." << std::endl;

    scheduler sch;
    demultiplex_reads* demultiplexer = NULL;

    try {
        // Step 1: Read input file
        if (config.interleaved_input) {
            sch.add_step(ai_read_fastq, "read_interleaved_fastq",
                         new read_interleaved_fastq(config.quality_input_fmt.get(),
                                                    config.input_files_1,
                                                    ai_demultiplex));
        } else {
            sch.add_step(ai_read_fastq, "read_paired_fastq",
                         new read_paired_fastq(config.quality_input_fmt.get(),
                                               config.input_files_1,
                                               config.input_files_2,
                                               ai_demultiplex));
        }

        // Step 2: Parse and demultiplex reads based on single or double indices
        sch.add_step(ai_demultiplex, "demultiplex_pe",
                     demultiplexer = new demultiplex_pe_reads(&config));

        add_write_step(config, sch, ai_write_unidentified_1, "unidentified_mate_1",
                       new write_fastq(config.get_output_filename("demux_unknown", 1)));

        if (!config.interleaved_output) {
            add_write_step(config, sch, ai_write_unidentified_2, "unidentified_mate_2",
                           new write_fastq(config.get_output_filename("demux_unknown", 2)));
        }

        // Step 3 - N: Write demultiplexed reads
        for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
            const size_t offset = nth * ai_analyses_offset;
            const std::string& sample = config.adapters.get_sample_name(nth);

            sch.add_step(offset + ai_trim_pe, "process_pe_" + sample,
                         new pe_demultiplexed_reads_processor(config, nth));

            add_write_step(config, sch, offset + ai_write_mate_1, sample + "_mate_1",
                           new write_fastq(config.get_output_filename("--output1", nth)));

            if (!config.interleaved_output) {
                add_write_step(config, sch, offset + ai_write_mate_2, sample + "_mate_2",
                               new write_fastq(config.get_output_filename("--output2", nth)));
            }
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;
        return 1;
    }

    if (!sch.run(config.max_threads)) {
        return 1;
    } else if (!write_demultiplex_settings(config, demultiplexer)) {
        return 1;
    }

    for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
        if (!write_demultiplex_settings(config, demultiplexer, nth)) {
            return 1;
        }
    }

    return 0;
}


int demultiplex_sequences(const userconfig& config)
{
    if (config.paired_ended_mode) {
        return demultiplex_sequences_pe(config);
    } else {
        return demultiplex_sequences_se(config);
    }
}

} // namespace ar
