// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_selector.hpp" // for select_adapter
#include "commontypes.hpp"      // for adapter_selection
#include "demultiplexing.hpp"   // for demultiplex_pe_reads, demultiplex_se_r...
#include "fastq_io.hpp"         // for gzip_split_fastq, post_process_fastq
#include "main.hpp"             // declarations
#include "output.hpp"           // for outpuT_file, DEV_NULL
#include "reports.hpp"          // for write_html_report, write_json_report
#include "scheduler.hpp"        // for scheduler
#include "statistics.hpp"       // for trim_stats_ptr, trimming_statistics
#include "trimming.hpp"         // for pe_reads_processor, se_reads_processor
#include "userconfig.hpp"       // for userconfig, output_files, DEV_NULL
#include <cstring>              // for size_t
#include <limits>               // for numeric_limits
#include <memory>               // for make_shared

namespace adapterremoval {

int
remove_adapter_sequences(const userconfig& config)
{
  scheduler sch;
  statistics stats = statistics_builder()
                       .sample_rate(config.report_sample_rate)
                       .estimate_duplication(config.report_duplication)
                       .demultiplexing(config.samples.get_reader()->size())
                       .initialize();

  auto output = config.get_output_filenames();
  // Add write steps for demultiplexed and per-samples output
  output.add_write_steps(sch, config);

  post_demux_steps steps;

  // Step 4 - N: Trim and write (demultiplexed) reads
  for (size_t nth = 0; nth < output.samples().size(); ++nth) {
    stats.trimming.push_back(std::make_shared<trimming_statistics>());

    if (!config.is_adapter_trimming_enabled()) {
      steps.samples.push_back(
        sch.add<process_demultiplexed>(config,
                                       output.get_sample(nth),
                                       nth,
                                       stats.trimming.back()));
    } else if (config.paired_ended_mode) {
      steps.samples.push_back(
        sch.add<pe_reads_processor>(config,
                                    output.get_sample(nth),
                                    nth,
                                    stats.trimming.back()));
    } else {
      steps.samples.push_back(
        sch.add<se_reads_processor>(config,
                                    output.get_sample(nth),
                                    nth,
                                    stats.trimming.back()));
    }
  }

  // Step 3: Parse and demultiplex reads based on single or double indices
  size_t step = std::numeric_limits<size_t>::max();
  if (config.is_demultiplexing_enabled()) {
    // Statistics and serialization of unidentified reads
    steps.unidentified =
      sch.add<processes_unidentified>(config, output, stats.demultiplexing);

    if (config.paired_ended_mode) {
      step = sch.add<demultiplex_pe_reads>(config, steps, stats.demultiplexing);
    } else {
      step = sch.add<demultiplex_se_reads>(config, steps, stats.demultiplexing);
    }
  } else {
    step = steps.samples.back();
  }

  if (config.adapter_selection_strategy == adapter_selection::automatic) {
    const auto database = config.known_adapters();

    step = sch.add<adapter_finalizer>(database,
                                      config.samples,
                                      config.adapter_fallback_strategy,
                                      step);
    step = sch.add<adapter_selector>(database,
                                     config.simd,
                                     config.mismatch_threshold,
                                     step,
                                     config.max_threads);
    step = sch.add<adapter_preselector>(step);
  }

  // Step 2: Post-process, validate, and collect statistics on FASTQ reads
  step = sch.add<post_process_fastq>(config, step, stats);

  // Step 1: Read input file(s)
  read_fastq::add_steps(sch, config, step, stats);

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  if (!write_json_report(config, stats, output.settings_json)) {
    return 1;
  }

  return !write_html_report(config, stats, output.settings_html);
}

} // namespace adapterremoval
