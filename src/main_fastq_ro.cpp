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
#include <cstring> // for size_t

#include "fastq_io.hpp"   // for read_fastq
#include "reports.hpp"    // for write_report
#include "scheduler.hpp"  // for scheduler
#include "statistics.hpp" // for ar_statistics
#include "userconfig.hpp" // for userconfig, output_files

namespace adapterremoval {

class reads_sink : public analytical_step
{
public:
  reads_sink()
    : analytical_step(processing_order::unordered, "reads_sink")
  {
  }

  virtual ~reads_sink() override;

  chunk_vec process(chunk_ptr) override { return chunk_vec(); }
};

int
fastq_report_only(const userconfig& config)
{
  scheduler sch;

  statistics stats = statistics_builder()
                       .sample_rate(config.report_sample_rate)
                       .estimate_duplication(config.report_duplication)
                       .initialize();

  // Discard all written reads
  size_t sink_step = sch.add<reads_sink>();

  // Step 2: Post-processing, validate, and collect statistics on FASTQ reads
  const size_t postproc_step =
    sch.add<post_process_fastq>(sink_step, stats, config.io_encoding);

  // Step 1: Read input file(s)
  sch.add<read_fastq>(config, postproc_step);

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  const auto out_files = config.get_output_filenames();
  if (!write_json_report(config, stats, out_files.settings_json)) {
    return 1;
  }

  return !write_html_report(config, stats, out_files.settings_html);
}

// Out-of-line definition to make -Wweak-vtables happy
reads_sink::~reads_sink() {}

} // namespace adapterremoval
