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

#include "debug.hpp"
#include "demultiplexing.hpp"
#include "fastq.hpp"
#include "fastq_io.hpp"
#include "main.hpp"
#include "reports.hpp"
#include "trimming.hpp"
#include "userconfig.hpp"

//! Implemented in main_adapter_rm.cpp
void
add_write_step(const userconfig& config,
               scheduler& sch,
               size_t offset,
               const std::string& name,
               analytical_step* step);

bool
demux_write_report(const userconfig& config,
                   const std::vector<reads_processor*>& processors,
                   ar_statistics& stats)
{
  for (auto ptr : processors) {
    stats.trimming.push_back(*ptr->get_final_statistics());
  }

  return write_report(config, stats);
}

class se_demultiplexed_reads_processor : public reads_processor
{
public:
  se_demultiplexed_reads_processor(const userconfig& config, size_t nth)
    : reads_processor(config, nth)
  {}

  chunk_vec process(analytical_chunk* chunk)
  {
    const size_t offset = m_nth * ai_analyses_offset;
    read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));

    auto stats = m_stats.acquire();
    output_chunk_ptr encoded_reads(new fastq_output_chunk(read_chunk->eof));

    for (const auto& read : read_chunk->reads_1) {
      stats->read_1.process(read);
      encoded_reads->add(*m_config.quality_output_fmt, read);
    }

    chunk_vec chunks;
    chunks.push_back(
      chunk_pair(offset + ai_write_mate_1, std::move(encoded_reads)));

    m_stats.release(stats);

    return chunks;
  }
};

class pe_demultiplexed_reads_processor : public reads_processor
{
public:
  pe_demultiplexed_reads_processor(const userconfig& config, size_t nth)
    : reads_processor(config, nth)
  {}

  chunk_vec process(analytical_chunk* chunk)
  {
    const size_t offset = m_nth * ai_analyses_offset;
    read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
    AR_DEBUG_ASSERT(read_chunk->reads_1.size() == read_chunk->reads_2.size());

    statistics_ptr stats = m_stats.acquire();
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

      stats->read_1.process(read_1);
      stats->read_2.process(read_2);

      encoded_reads_1->add(*m_config.quality_output_fmt, read_1);

      if (m_config.interleaved_output) {
        encoded_reads_1->add(*m_config.quality_output_fmt, read_2);
      } else {
        encoded_reads_2->add(*m_config.quality_output_fmt, read_2);
      }
    }

    chunk_vec chunks;
    chunks.push_back(
      chunk_pair(offset + ai_write_mate_1, std::move(encoded_reads_1)));
    if (!m_config.interleaved_output) {
      chunks.push_back(
        chunk_pair(offset + ai_write_mate_2, std::move(encoded_reads_2)));
    }

    m_stats.release(stats);

    return chunks;
  }
};

int
demultiplex_sequences_se(const userconfig& config)
{
  std::cerr << "Demultiplexing single ended reads ..." << std::endl;

  scheduler sch;
  std::vector<reads_processor*> processors;
  ar_statistics stats;

  try {
    // Step 1: Read input file
    sch.add_step(ai_read_fastq,
                 "read_fastq",
                 new read_single_fastq(config.quality_input_fmt.get(),
                                       config.input_files_1,
                                       ai_demultiplex,
                                       &stats.input_1));

    // Step 2: Parse and demultiplex reads based on single or double indices
    sch.add_step(ai_demultiplex,
                 "demultiplex_se",
                 new demultiplex_se_reads(&config, &stats.demultiplexing));

    add_write_step(
      config,
      sch,
      ai_write_unidentified_1,
      "unidentified",
      new write_fastq(config.get_output_filename("demux_unknown", 1)));

    // Step 3 - N: Trim and write demultiplexed reads
    for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
      const size_t offset = nth * ai_analyses_offset;
      const std::string& sample = config.adapters.get_sample_name(nth);

      processors.push_back(new se_demultiplexed_reads_processor(config, nth));
      sch.add_step(
        offset + ai_trim_se, "process_se_" + sample, processors.back());

      add_write_step(
        config,
        sch,
        offset + ai_write_mate_1,
        sample + "_fastq",
        new write_fastq(config.get_output_filename("--output1", nth)));
    }
  } catch (const std::ios_base::failure& error) {
    std::cerr << "IO error opening file; aborting:\n"
              << cli_formatter::fmt(error.what()) << std::endl;
    return 1;
  }

  if (!sch.run(config.max_threads)) {
    return 1;
  } else if (!demux_write_report(config, processors, stats)) {
    return 1;
  }

  return 0;
}

int
demultiplex_sequences_pe(const userconfig& config)
{
  std::cerr << "Demultiplexing paired end reads ..." << std::endl;

  scheduler sch;
  std::vector<reads_processor*> processors;
  ar_statistics stats;

  try {
    // Step 1: Read input file
    if (config.interleaved_input) {
      sch.add_step(ai_read_fastq,
                   "read_interleaved_fastq",
                   new read_interleaved_fastq(config.quality_input_fmt.get(),
                                              config.input_files_1,
                                              ai_demultiplex,
                                              &stats.input_1,
                                              &stats.input_2));
    } else {
      sch.add_step(ai_read_fastq,
                   "read_paired_fastq",
                   new read_paired_fastq(config.quality_input_fmt.get(),
                                         config.input_files_1,
                                         config.input_files_2,
                                         ai_demultiplex,
                                         &stats.input_1,
                                         &stats.input_2));
    }

    // Step 2: Parse and demultiplex reads based on single or double indices
    sch.add_step(ai_demultiplex,
                 "demultiplex_pe",
                 new demultiplex_pe_reads(&config, &stats.demultiplexing));

    add_write_step(
      config,
      sch,
      ai_write_unidentified_1,
      "unidentified_mate_1",
      new write_fastq(config.get_output_filename("demux_unknown", 1)));

    if (!config.interleaved_output) {
      add_write_step(
        config,
        sch,
        ai_write_unidentified_2,
        "unidentified_mate_2",
        new write_fastq(config.get_output_filename("demux_unknown", 2)));
    }

    // Step 3 - N: Write demultiplexed reads
    for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
      const size_t offset = nth * ai_analyses_offset;
      const std::string& sample = config.adapters.get_sample_name(nth);

      processors.push_back(new pe_demultiplexed_reads_processor(config, nth));
      sch.add_step(
        offset + ai_trim_pe, "process_pe_" + sample, processors.back());

      add_write_step(
        config,
        sch,
        offset + ai_write_mate_1,
        sample + "_mate_1",
        new write_fastq(config.get_output_filename("--output1", nth)));

      if (!config.interleaved_output) {
        add_write_step(
          config,
          sch,
          offset + ai_write_mate_2,
          sample + "_mate_2",
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
  } else if (!demux_write_report(config, processors, stats)) {
    return 1;
  }

  return 0;
}

int
demultiplex_sequences(const userconfig& config)
{
  if (config.paired_ended_mode) {
    return demultiplex_sequences_pe(config);
  } else {
    return demultiplex_sequences_se(config);
  }
}
