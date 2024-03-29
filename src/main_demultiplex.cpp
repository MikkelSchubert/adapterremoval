/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2016 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <limits>    // for numeric_limits
#include <memory>    // for unique_ptr, make_unique
#include <string>    // for string, operator+
#include <vector>    // for vector, vector<>::iterator

#include "adapterset.hpp"     // for adapter_set
#include "commontypes.hpp"    // for fastq_vec, read_type, read_type::mate_1
#include "debug.hpp"          // for AR_REQUIRE
#include "demultiplexing.hpp" // for post_demux_steps, demultiplex_pe_reads
#include "fastq_io.hpp"       // for fastq_read_chunk, read_chunk_ptr, read...
#include "reports.hpp"        // for write_report
#include "scheduler.hpp"      // for scheduler, threadstate, analytical_chunk
#include "statistics.hpp"     // for trimming_statistics, ar_statistics
#include "trimming.hpp"       // for trimmed_reads, reads_processor
#include "userconfig.hpp"     // for userconfig, output_files, output_sampl...

namespace adapterremoval {

class fastq;

//! Implemented in main_adapter_rm.cpp
size_t
add_write_step(scheduler& sch,
               const userconfig& config,
               const std::string& filename);

class se_demuxed_processor : public reads_processor
{
public:
  se_demuxed_processor(const userconfig& config,
                       const output_sample_files& output,
                       const size_t nth,
                       trim_stats_ptr sink)
    : reads_processor(config, output, nth, sink)
  {
  }

  ~se_demuxed_processor() override = default;

  chunk_vec process(chunk_ptr chunk) override
  {
    auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);

    auto stats = m_stats.acquire();
    trimmed_reads chunks(m_output, read_chunk.eof);

    for (auto& read : read_chunk.reads_1) {
      stats->read_1->process(read);
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
                       const size_t nth,
                       trim_stats_ptr sink)
    : reads_processor(config, output, nth, sink)
  {
  }

  ~pe_demuxed_processor() override = default;

  chunk_vec process(chunk_ptr chunk) override
  {
    auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
    AR_REQUIRE(read_chunk.reads_1.size() == read_chunk.reads_2.size());

    auto stats = m_stats.acquire();
    trimmed_reads chunks(m_output, read_chunk.eof);

    auto it_1 = read_chunk.reads_1.begin();
    auto it_2 = read_chunk.reads_2.begin();
    while (it_1 != read_chunk.reads_1.end()) {
      fastq& read_1 = *it_1++;
      fastq& read_2 = *it_2++;

      stats->read_1->process(read_1);
      stats->read_2->process(read_2);

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
  scheduler sch;

  statistics stats = statistics_builder()
                       .sample_rate(config.report_sample_rate)
                       .estimate_duplication(config.report_duplication)
                       .demultiplexing(config.adapters.barcode_count())
                       .initialize();

  auto out_files = config.get_output_filenames();

  post_demux_steps steps;

  // Step 4 - N: Trim and write (demultiplexed) reads
  for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
    auto& mapping = out_files.samples.at(nth);
    for (const auto& filename : mapping.filenames()) {
      mapping.push_pipeline_step(add_write_step(sch, config, filename));
    }

    stats.trimming.push_back(std::make_shared<trimming_statistics>());
    if (config.paired_ended_mode) {
      steps.samples.push_back(sch.add<pe_demuxed_processor>(
        config, mapping, nth, stats.trimming.back()));
    } else {
      steps.samples.push_back(sch.add<se_demuxed_processor>(
        config, mapping, nth, stats.trimming.back()));
    }
  }

  size_t processing_step = std::numeric_limits<size_t>::max();

  // Step 3: Parse and demultiplex reads based on single or double indices
  if (config.adapters.barcode_count()) {
    if (out_files.unidentified_1 != DEV_NULL) {
      steps.unidentified_1 =
        add_write_step(sch, config, out_files.unidentified_1);
    }

    if (config.paired_ended_mode && !config.interleaved_output &&
        out_files.unidentified_2 != DEV_NULL) {
      steps.unidentified_2 =
        add_write_step(sch, config, out_files.unidentified_2);
    }

    if (config.paired_ended_mode) {
      processing_step =
        sch.add<demultiplex_pe_reads>(config, steps, stats.demultiplexing);

    } else {
      processing_step =
        sch.add<demultiplex_se_reads>(config, steps, stats.demultiplexing);
    }
  } else {
    processing_step = steps.samples.back();
  }

  // Step 2: Post-process, validate, and collect statistics on FASTQ reads
  const size_t postproc_step =
    sch.add<post_process_fastq>(processing_step, stats, config.io_encoding);

  // Step 1: Read input file(s)
  sch.add<read_fastq>(config, postproc_step);

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  if (!write_json_report(config, stats, out_files.settings_json)) {
    return 1;
  }

  return !write_html_report(config, stats, out_files.settings_html);
}

} // namespace adapterremoval
