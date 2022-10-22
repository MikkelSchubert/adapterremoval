/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <algorithm> // for max
#include <cstring>   // for size_t
#include <limits>    // for numeric_limits
#include <memory>    // for unique_ptr, make_unique
#include <string>    // for operator+, string
#include <vector>    // for vector

#include "adapterset.hpp"     // for adapter_set
#include "commontypes.hpp"    // for string_vec
#include "demultiplexing.hpp" // for post_demux_steps, demultiplex_pe_reads
#include "fastq_io.hpp"       // for gzip_fastq, gzip_split_fastq, read_fastq
#include "reports.hpp"        // for write_report
#include "scheduler.hpp"      // for scheduler
#include "statistics.hpp"     // for trimming_statistics, ar_statistics
#include "strutils.hpp"       // for ends_with
#include "trimming.hpp"       // for pe_reads_processor, reads_processor
#include "userconfig.hpp"     // for userconfig, output_files, output_sampl...

namespace adapterremoval {

size_t
add_write_step(scheduler& sch,
               const userconfig& config,
               const std::string& filename)
{
  AR_REQUIRE(filename != DEV_NULL);
  size_t step_id = sch.add<write_fastq>(config, filename);

  if (config.gzip || ends_with(tolower(filename), ".gz")) {
    step_id = sch.add<gzip_split_fastq>(config, filename, step_id);
    step_id = sch.add<split_fastq>(config, filename, step_id);
  }

  return step_id;
}

int
remove_adapter_sequences(const userconfig& config)
{
  scheduler sch;

  statistics stats = statistics_builder()
                       .sample_rate(config.report_sample_rate)
                       .estimate_duplication(config.report_duplication)
                       .demultiplexing(config.adapters.barcode_count())
                       .initialize();

  auto out_files = config.get_output_filenames();

  post_demux_steps steps;
  size_t processing_step = std::numeric_limits<size_t>::max();

  // Step 4 - N: Trim and write (demultiplexed) reads
  for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
    auto& mapping = out_files.samples.at(nth);

    for (const auto& filename : mapping.filenames()) {
      mapping.push_pipeline_step(add_write_step(sch, config, filename));
    }

    stats.trimming.push_back(std::make_shared<trimming_statistics>());

    if (config.paired_ended_mode) {
      steps.samples.push_back(sch.add<pe_reads_processor>(
        config, mapping, nth, stats.trimming.back()));
    } else {
      steps.samples.push_back(sch.add<se_reads_processor>(
        config, mapping, nth, stats.trimming.back()));
    }
  }

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
    sch.add<post_process_fastq>(processing_step, stats);

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
