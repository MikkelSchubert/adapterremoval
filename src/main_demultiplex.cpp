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
#include "adapterset.hpp"     // for adapter_set
#include "commontypes.hpp"    // for fastq_vec, read_type, read_type::mate_1
#include "debug.hpp"          // for AR_REQUIRE
#include "demultiplexing.hpp" // for demultiplex_pe_reads, demultiplex_se_r...
#include "fastq.hpp"          // for fastq
#include "fastq_io.hpp"       // for read_chunk, post_process_fastq
#include "output.hpp"         // for output_files, processed reads
#include "reports.hpp"        // for write_html_report, write_json_report
#include "scheduler.hpp"      // for scheduler, threadstate, analytical_chunk
#include "serializer.hpp"     // for fastq_flags
#include "simd.hpp"           // for size_t
#include "statistics.hpp"     // for trim_stats_ptr, trimming_statistics
#include "trimming.hpp"       // for  reads_processor
#include "userconfig.hpp"     // for userconfig, output_files, DEV_NULL
#include <algorithm>          // for max
#include <cstring>            // for size_t
#include <limits>             // for numeric_limits
#include <memory>             // for unique_ptr, __shared_ptr_access, make_...
#include <string>             // for operator!=, basic_string, string
#include <vector>             // for vector

namespace adapterremoval {

class se_demuxed_processor : public reads_processor
{
public:
  se_demuxed_processor(const userconfig& config,
                       const sample_output_files& output,
                       const size_t nth,
                       trim_stats_ptr sink)
    : reads_processor(config, output, nth, sink)
  {
  }

  ~se_demuxed_processor() override = default;

  chunk_vec process(chunk_ptr chunk) override
  {
    AR_REQUIRE(chunk);
    auto stats = m_stats.acquire();
    processed_reads chunks{ m_output, chunk->first };

    for (auto& read : chunk->reads_1) {
      stats->read_1->process(read);
      chunks.add(read, read_type::mate_1, fastq_flags::se);
    }

    m_stats.release(stats);

    return chunks.finalize(chunk->eof);
  }

  se_demuxed_processor(const se_demuxed_processor&) = delete;
  se_demuxed_processor(se_demuxed_processor&&) = delete;
  se_demuxed_processor& operator=(const se_demuxed_processor&) = delete;
  se_demuxed_processor& operator=(se_demuxed_processor&&) = delete;
};

class pe_demuxed_processor : public reads_processor
{
public:
  pe_demuxed_processor(const userconfig& config,
                       const sample_output_files& output,
                       const size_t nth,
                       trim_stats_ptr sink)
    : reads_processor(config, output, nth, sink)
  {
  }

  ~pe_demuxed_processor() override = default;

  chunk_vec process(chunk_ptr chunk) override
  {
    AR_REQUIRE(chunk);
    AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size());

    auto stats = m_stats.acquire();
    processed_reads chunks{ m_output, chunk->first };

    auto it_1 = chunk->reads_1.begin();
    auto it_2 = chunk->reads_2.begin();
    while (it_1 != chunk->reads_1.end()) {
      fastq& read_1 = *it_1++;
      fastq& read_2 = *it_2++;

      stats->read_1->process(read_1);
      stats->read_2->process(read_2);

      chunks.add(read_1, read_type::mate_1, fastq_flags::pe_1);
      chunks.add(read_2, read_type::mate_2, fastq_flags::pe_2);
    }

    m_stats.release(stats);

    return chunks.finalize(chunk->eof);
  }

  pe_demuxed_processor(const pe_demuxed_processor&) = delete;
  pe_demuxed_processor(pe_demuxed_processor&&) = delete;
  pe_demuxed_processor& operator=(const pe_demuxed_processor&) = delete;
  pe_demuxed_processor& operator=(pe_demuxed_processor&&) = delete;
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

  auto output = config.get_output_filenames();
  output.add_write_steps(sch, config);

  post_demux_steps steps{ output };

  // Step 4 - N: Trim and write (demultiplexed) reads
  for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
    stats.trimming.push_back(std::make_shared<trimming_statistics>());

    if (config.paired_ended_mode) {
      steps.samples.push_back(sch.add<pe_demuxed_processor>(
        config, output.get_sample(nth), nth, stats.trimming.back()));
    } else {
      steps.samples.push_back(sch.add<se_demuxed_processor>(
        config, output.get_sample(nth), nth, stats.trimming.back()));
    }
  }

  // Step 3: Parse and demultiplex reads based on single or double indices
  size_t processing_step = std::numeric_limits<size_t>::max();
  if (config.is_demultiplexing_enabled()) {
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
    sch.add<post_process_fastq>(config, processing_step, stats);

  // Step 1: Read input file(s)
  read_fastq::add_steps(sch, config, postproc_step);

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  if (!write_json_report(config, stats, output.settings_json)) {
    return 1;
  }

  return !write_html_report(config, stats, output.settings_html);
}

} // namespace adapterremoval
