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
#include <algorithm> // for copy, max
#include <cstring>   // for size_t
#include <iostream>  // for operator<<, endl, basic_ostream, cerr
#include <limits>    // for numeric_limits
#include <memory>    // for unique_ptr
#include <string>    // for string, operator+
#include <vector>    // for vector, vector<>::iterator

#include "adapterset.hpp"     // for adapter_set
#include "commontypes.hpp"    // for fastq_vec, read_type, read_type::mate_1
#include "debug.hpp"          // for AR_DEBUG_ASSERT
#include "demultiplexing.hpp" // for post_demux_steps, demultiplex_pe_reads
#include "fastq_io.hpp"       // for fastq_read_chunk, read_chunk_ptr, read...
#include "reports.hpp"        // for write_report
#include "scheduler.hpp"      // for scheduler, threadstate, analytical_chunk
#include "statistics.hpp"     // for trimming_statistics, ar_statistics
#include "trimming.hpp"       // for trimmed_reads, reads_processor
#include "userconfig.hpp"     // for userconfig, output_files, output_sampl...

class fastq;

//! Implemented in main_adapter_rm.cpp
size_t
add_write_step(const userconfig& config,
               scheduler& sch,
               const std::string& name,
               const std::string& filename);

class se_demuxed_processor : public reads_processor
{
public:
  se_demuxed_processor(const userconfig& config,
                       const output_sample_files& output,
                       size_t nth)
    : reads_processor(config, output, nth)
  {}

  chunk_vec process(analytical_chunk* chunk)
  {
    read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));

    statistics_ptr stats = m_stats.acquire();
    trimmed_reads chunks(m_config, m_output, read_chunk->eof);

    for (auto& read : read_chunk->reads_1) {
      stats->read_1.process(read);
      chunks.add(read, read_type::mate_1);
    }

    m_stats.release(stats);

    return chunks.finalize();
  }
};

class pe_demuxed_processor : public reads_processor
{
public:
  pe_demuxed_processor(const userconfig& config,
                       const output_sample_files& output,
                       size_t nth)
    : reads_processor(config, output, nth)
  {}

  chunk_vec process(analytical_chunk* chunk)
  {
    read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
    AR_DEBUG_ASSERT(read_chunk->reads_1.size() == read_chunk->reads_2.size());

    statistics_ptr stats = m_stats.acquire();
    trimmed_reads chunks(m_config, m_output, read_chunk->eof);

    fastq_vec::iterator it_1 = read_chunk->reads_1.begin();
    fastq_vec::iterator it_2 = read_chunk->reads_2.begin();
    while (it_1 != read_chunk->reads_1.end()) {
      fastq& read_1 = *it_1++;
      fastq& read_2 = *it_2++;

      stats->read_1.process(read_1);
      stats->read_2.process(read_2);

      chunks.add(read_1, read_type::mate_1);
      chunks.add(read_2, read_type::mate_2);
    }

    m_stats.release(stats);

    return chunks.finalize();
  }
};

int
demultiplex_sequences(const userconfig& config)
{
  std::cerr << "Demultiplexing reads ..." << std::endl;

  scheduler sch;
  std::vector<reads_processor*> processors;
  ar_statistics stats(config.report_sample_rate);

  auto out_files = config.get_output_filenames();

  post_demux_steps steps;

  // Step 3 - N: Trim and write (demultiplexed) reads
  for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
    const std::string& sample = config.adapters.get_sample_name(nth);

    auto& mapping = out_files.samples.at(nth);
    for (const auto& filename : mapping.filenames) {
      mapping.steps.push_back(add_write_step(config, sch, sample, filename));
    }

    if (config.paired_ended_mode) {
      processors.push_back(new pe_demuxed_processor(config, mapping, nth));
    } else {
      processors.push_back(new se_demuxed_processor(config, mapping, nth));
    }

    steps.samples.push_back(sch.add_step("post_" + sample, processors.back()));
  }

  size_t processing_step = std::numeric_limits<size_t>::max();

  // Step 2: Parse and demultiplex reads based on single or double indices
  if (config.adapters.barcode_count()) {
    steps.unidentified_1 = add_write_step(
      config, sch, "unidentified_mate_1", out_files.unidentified_1);

    if (config.paired_ended_mode && !config.interleaved_output) {
      steps.unidentified_2 = add_write_step(
        config, sch, "unidentified_mate_2", out_files.unidentified_2);
    }

    if (config.paired_ended_mode) {
      processing_step = sch.add_step(
        "demultiplex",
        new demultiplex_pe_reads(config, steps, &stats.demultiplexing));

    } else {
      processing_step = sch.add_step(
        "demultiplex",
        new demultiplex_se_reads(config, steps, &stats.demultiplexing));
    }
  } else {
    processing_step = steps.samples.back();
  }

  // Step 1: Read input file
  sch.add_step("read_fastq",
               new read_fastq(config.quality_input_fmt.get(),
                              config.input_files_1,
                              config.input_files_2,
                              processing_step,
                              config.interleaved_input,
                              &stats.input_1,
                              &stats.input_2));

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  for (auto ptr : processors) {
    stats.trimming.push_back(*ptr->get_final_statistics());
  }

  return !write_report(config, stats, out_files.settings);
}
