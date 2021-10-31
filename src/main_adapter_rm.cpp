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
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "debug.hpp"
#include "demultiplexing.hpp"
#include "fastq.hpp"
#include "reports.hpp"
#include "trimmed_reads.hpp"
#include "trimming.hpp"
#include "userconfig.hpp"

void
add_write_step(const userconfig& config,
               scheduler& sch,
               size_t offset,
               const std::string& name,
               analytical_step* step)
{
  if (config.gzip) {
    if (!config.gzip_stream) {
      sch.add_step(
        offset, "split_" + name, new split_fastq(offset + ai_split_offset));

      sch.add_step(offset + ai_split_offset,
                   "gzip_" + name,
                   new gzip_split_fastq(config, offset + ai_zip_offset));
    } else {
      sch.add_step(
        offset, "gzip_" + name, new gzip_fastq(config, offset + ai_zip_offset));
    }

    sch.add_step(offset + ai_zip_offset, "write_gzip_" + name, step);
  } else {
    sch.add_step(offset, "write_" + name, step);
  }
}

int
remove_adapter_sequences(const userconfig& config)
{
  std::cerr << "Trimming reads ..." << std::endl;

  scheduler sch;
  std::vector<reads_processor*> processors;
  ar_statistics stats(config.report_sample_rate);

  const auto out_files = config.get_output_filenames();

  // Step 1: Read input file
  const size_t next_step =
    config.adapters.barcode_count() ? ai_demultiplex : ai_analyses_offset;
  sch.add_step(ai_read_fastq,
               "read_fastq",
               new read_fastq(config.quality_input_fmt.get(),
                              config.input_files_1,
                              config.input_files_2,
                              next_step,
                              config.interleaved_input,
                              &stats.input_1,
                              &stats.input_2));

  if (config.adapters.barcode_count()) {
    // Step 2: Parse and demultiplex reads based on single or double indices
    if (config.paired_ended_mode) {
      sch.add_step(ai_demultiplex,
                   "demultiplex",
                   new demultiplex_pe_reads(&config, &stats.demultiplexing));

    } else {
      sch.add_step(ai_demultiplex,
                   "demultiplex",
                   new demultiplex_se_reads(&config, &stats.demultiplexing));
    }

    add_write_step(config,
                   sch,
                   ai_write_unidentified_1,
                   "unidentified_mate_1",
                   new write_fastq(out_files.unidentified_1));

    if (config.paired_ended_mode && !config.interleaved_output) {
      add_write_step(config,
                     sch,
                     ai_write_unidentified_2,
                     "unidentified_mate_2",
                     new write_fastq(out_files.unidentified_2));
    }
  }

  // Step 3 - N: Trim and write demultiplexed reads
  for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
    const size_t offset = (nth + 1) * ai_analyses_offset;
    const std::string& sample = config.adapters.get_sample_name(nth);
    const auto& filemap = out_files.samples.at(nth);
    const auto& filenames = filemap.filenames;

    if (config.paired_ended_mode) {
      processors.push_back(new pe_reads_processor(config, nth));
    } else {
      processors.push_back(new se_reads_processor(config, nth));
    }

    sch.add_step(offset + ai_trim_se, "trim_se_" + sample, processors.back());

    add_write_step(config,
                   sch,
                   offset + ai_write_mate_1,
                   sample + "_mate_1",
                   new write_fastq(filenames.at(filemap.output_1)));

    if (config.paired_ended_mode && !config.interleaved_output) {
      add_write_step(config,
                     sch,
                     offset + ai_write_mate_2,
                     sample + "_mate_2",
                     new write_fastq(filenames.at(filemap.output_2)));
    }

    if (!config.combined_output) {
      add_write_step(config,
                     sch,
                     offset + ai_write_discarded,
                     sample + "_discarded",
                     new write_fastq(filenames.at(filemap.discarded)));

      if (config.paired_ended_mode) {
        add_write_step(config,
                       sch,
                       offset + ai_write_singleton,
                       sample + "_singleton",
                       new write_fastq(filenames.at(filemap.singleton)));
      }

      if (config.collapse) {
        add_write_step(config,
                       sch,
                       offset + ai_write_collapsed,
                       sample + "_collapsed",
                       new write_fastq(filenames.at(filemap.merged)));
      }
    }
  }

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  for (auto ptr : processors) {
    stats.trimming.push_back(*ptr->get_final_statistics());
  }

  return !write_report(config, stats, out_files.settings);
}
