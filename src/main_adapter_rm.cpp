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
#include <algorithm> // for max
#include <cstring>   // for size_t
#include <iostream>  // for operator<<, endl, basic_ostream, cerr
#include <limits>    // for numeric_limits
#include <memory>    // for unique_ptr
#include <string>    // for operator+, string
#include <vector>    // for vector

#include "adapterset.hpp"     // for adapter_set
#include "commontypes.hpp"    // for string_vec
#include "demultiplexing.hpp" // for post_demux_steps, demultiplex_pe_reads
#include "fastq_io.hpp"       // for gzip_fastq, gzip_split_fastq, read_fastq
#include "reports.hpp"        // for write_report
#include "scheduler.hpp"      // for scheduler
#include "statistics.hpp"     // for trimming_statistics, ar_statistics
#include "trimming.hpp"       // for pe_reads_processor, reads_processor
#include "userconfig.hpp"     // for userconfig, output_files, output_sampl...

size_t
add_write_step(const userconfig& config,
               scheduler& sch,
               const std::string& name,
               const std::string& filename)
{
  size_t step_id = sch.add_step("write_" + name, new write_fastq(filename));

  if (config.gzip_stream) {
    step_id = sch.add_step("gzip_" + name, new gzip_fastq(config, step_id));
  } else if (config.gzip) {
    step_id =
      sch.add_step("gzip_" + name, new gzip_split_fastq(config, step_id));

    step_id = sch.add_step("split_" + name, new split_fastq(step_id));
  }

  return step_id;
}

int
remove_adapter_sequences(const userconfig& config)
{
  std::cerr << "Trimming reads" << std::endl;

  scheduler sch;
  std::vector<reads_processor*> processors;

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
    const std::string& sample = config.adapters.get_sample_name(nth);
    auto& mapping = out_files.samples.at(nth);

    for (const auto& filename : mapping.filenames) {
      mapping.steps.push_back(add_write_step(config, sch, sample, filename));
    }

    if (config.paired_ended_mode) {
      processors.push_back(new pe_reads_processor(config, mapping, nth));
    } else {
      processors.push_back(new se_reads_processor(config, mapping, nth));
    }

    steps.samples.push_back(sch.add_step("trim_" + sample, processors.back()));
  }

  // Step 3: Parse and demultiplex reads based on single or double indices
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
        new demultiplex_pe_reads(config, steps, stats.demultiplexing));

    } else {
      processing_step = sch.add_step(
        "demultiplex",
        new demultiplex_se_reads(config, steps, stats.demultiplexing));
    }
  } else {
    processing_step = steps.samples.back();
  }

  // Step 2: Post-process, validate, and collect statistics on FASTQ reads
  const size_t postproc_step =
    sch.add_step("post_process_fastq",
                 new post_process_fastq(config, processing_step, &stats));

  // Step 1: Read input file(s)
  sch.add_step("read_fastq", new read_fastq(config, postproc_step));

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  for (auto ptr : processors) {
    stats.trimming.push_back(ptr->get_final_statistics());
  }

  if (!write_json_report(config, stats, out_files.settings_json)) {
    return 1;
  }

  return !write_html_report(config, stats, out_files.settings_html);
}
