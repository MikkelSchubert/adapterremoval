/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "debug.hpp"          // for AR_REQUIRE
#include "demultiplexing.hpp" // for demultiplex_pe_reads, demultiplex_se_r...
#include "fastq_io.hpp"       // for gzip_split_fastq, post_process_fastq
#include "output.hpp"         // for outpuT_file, DEV_NULL
#include "reports.hpp"        // for write_html_report, write_json_report
#include "scheduler.hpp"      // for scheduler
#include "sequence_sets.hpp"  // for adapter_set
#include "simd.hpp"           // for size_t
#include "statistics.hpp"     // for trim_stats_ptr, trimming_statistics
#include "strutils.hpp"       // for ends_with, to_lower, string_vec
#include "trimming.hpp"       // for pe_reads_processor, se_reads_processor
#include "userconfig.hpp"     // for userconfig, output_files, DEV_NULL
#include <algorithm>          // for max
#include <cstring>            // for size_t
#include <limits>             // for numeric_limits
#include <memory>             // for make_shared
#include <string>             // for operator!=, basic_string, string
#include <vector>             // for vector

namespace adapterremoval {

int
remove_adapter_sequences(const userconfig& config)
{
  scheduler sch;

  statistics stats = statistics_builder()
                       .sample_rate(config.report_sample_rate)
                       .estimate_duplication(config.report_duplication)
                       .demultiplexing(config.adapters.barcode_count())
                       .initialize();

  auto output = config.get_output_filenames();
  // Add write steps for demultiplexed and per-samples output
  output.add_write_steps(sch, config);

  post_demux_steps steps;

  // Step 4 - N: Trim and write (demultiplexed) reads
  for (size_t nth = 0; nth < output.samples().size(); ++nth) {
    stats.trimming.push_back(std::make_shared<trimming_statistics>());

    if (!config.is_adapter_trimming_enabled()) {
      steps.samples.push_back(sch.add<process_demultiplexed>(
        config, output.get_sample(nth), stats.trimming.back()));
    } else if (config.paired_ended_mode) {
      steps.samples.push_back(sch.add<pe_reads_processor>(
        config, output.get_sample(nth), nth, stats.trimming.back()));
    } else {
      steps.samples.push_back(sch.add<se_reads_processor>(
        config, output.get_sample(nth), nth, stats.trimming.back()));
    }
  }

  // Step 3: Parse and demultiplex reads based on single or double indices
  size_t processing_step = std::numeric_limits<size_t>::max();
  if (config.is_demultiplexing_enabled()) {
    // Statistics and serialization of unidentified reads
    steps.unidentified =
      sch.add<processes_unidentified>(config, output, stats.demultiplexing);

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
