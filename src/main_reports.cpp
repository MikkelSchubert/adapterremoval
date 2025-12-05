// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_id.hpp"       // for adapter_id_statistics
#include "adapter_selector.hpp" // for adapter_selector, ...
#include "alignment.hpp"        // for extract_adapter_sequences, ...
#include "debug.hpp"            // for AR_REQUIRE
#include "fastq.hpp"            // for ACGTN, fastq, ACGT, ACGT:...
#include "fastq_io.hpp"         // for read_fastq, read_chunk
#include "main.hpp"             // declarations
#include "output.hpp"           // for output_files
#include "reports.hpp"          // for write_html_report, write_json_report
#include "scheduler.hpp"        // for threadstate, scheduler, analytical_step
#include "sequence_sets.hpp"    // for adapter_set
#include "simd_selector.hpp"    // for simd_selector
#include "statistics.hpp"       // for trimming_statistics
#include "userconfig.hpp"       // for userconfig
#include <cstddef>              // for size_t
#include <memory>               // for unique_ptr, __shared_ptr_access, ...
#include <string>               // for string, operator<<, char_traits
#include <vector>               // for vector

namespace adapterremoval {

class reads_sink : public analytical_step
{
public:
  reads_sink()
    : analytical_step(processing_order::unordered, "reads_sink")
  {
  }

  chunk_vec process(chunk_ptr /* chunk */) override { return {}; }
};

class adapter_identification : public analytical_step
{
public:
  explicit adapter_identification(const userconfig& config, statistics& stats)
    : analytical_step(processing_order::unordered, "adapter_identification")
    , m_config(config)
    , m_sink_id(stats.adapter_id)
  {
    AR_REQUIRE(!stats.trimming.empty());
    AR_REQUIRE(stats.adapter_id);
    m_sink_trim = stats.trimming.back();

    // FIXME: Shouldn't be adapter specific
    m_stats_id.emplace_back_n(m_config.max_threads,
                              m_sink_id->adapter1.max_length());
    m_stats_ins.emplace_back_n(m_config.max_threads);
  }

  ~adapter_identification() override = default;

  chunk_vec process(chunk_ptr data) override
  {
    auto chunk = dynamic_cast_unique<fastq_chunk>(data);
    AR_REQUIRE(chunk);

    const adapter_set adapters = { { "", "" } };

    auto aligner = sequence_aligner(adapters,
                                    *m_config.simd.get_reader(),
                                    m_config.mismatch_threshold);
    aligner.set_merge_threshold(m_config.merge_threshold);

    auto stats_id = m_stats_id.acquire();
    auto stats_ins = m_stats_ins.acquire();

    AR_REQUIRE(chunk->reads_1.size() == chunk->reads_2.size());
    auto it_1 = chunk->reads_1.begin();
    auto it_2 = chunk->reads_2.begin();

    for (; it_1 != chunk->reads_1.end(); ++it_1, ++it_2) {
      auto& read_1 = *it_1;
      auto& read_2 = *it_2;

      // Reverse complement to match the orientation of read1
      read_2.reverse_complement();

      const auto alignment =
        aligner.align_paired_end(read_1, read_2, m_config.shift);

      if (alignment.type() >= alignment_type::good) {
        stats_id->aligned_pairs++;
        if (alignment.type() == alignment_type::mergeable) {
          const size_t insert_size = alignment.insert_size(read_1, read_2);
          stats_ins->insert_sizes.resize_up_to(insert_size + 1);
          stats_ins->insert_sizes.inc(insert_size);

          if (alignment.extract_adapter_sequences(read_1, read_2)) {
            stats_id->pairs_with_adapters++;

            stats_id->adapter1.process(read_1.sequence());
            read_2.reverse_complement();
            stats_id->adapter2.process(read_2.sequence());
          }
        }
      }
    }

    m_stats_id.release(std::move(stats_id));
    m_stats_ins.release(std::move(stats_ins));

    return {};
  }

  void finalize() override
  {
    m_stats_id.merge_into(*m_sink_id);
    m_stats_ins.merge_into(*m_sink_trim);
  }

  adapter_identification(const adapter_identification&) = delete;
  adapter_identification(adapter_identification&&) = delete;
  adapter_identification& operator=(const adapter_identification&) = delete;
  adapter_identification& operator=(adapter_identification&&) = delete;

private:
  const userconfig& m_config;

  threadstate<adapter_id_statistics> m_stats_id{};
  threadstate<trimming_statistics> m_stats_ins{};

  trim_stats_ptr m_sink_trim{};
  adapter_id_stats_ptr m_sink_id{};
};

int
generate_reports(const userconfig& config)
{
  scheduler sch;

  // FIXME: Length picked based on known sequences
  const auto max_adapter_length = config.paired_ended_mode ? 42 : 0;
  statistics stats = statistics_builder()
                       .sample_rate(config.report_sample_rate)
                       .estimate_duplication(config.report_duplication)
                       .adapter_identification(max_adapter_length)
                       .initialize();
  // FIXME: Required for insert size statistics
  stats.trimming.push_back(std::make_shared<trimming_statistics>());

  // Step 5: Identify adapters from pair-wise alignments and infer insert sizes
  size_t step = std::numeric_limits<size_t>::max();
  if (config.paired_ended_mode) {
    // Attempt to identify adapters through pair-wise alignments
    step = sch.add<adapter_identification>(config, stats);
  } else {
    // Discard all written reads
    step = sch.add<reads_sink>();
  }

  // Step 4: Determine optional SIMD instruction set for alignments
  if (config.simd_auto_select) {
    step = sch.add<simd_selector>(config.samples,
                                  config.simd,
                                  config.mismatch_threshold,
                                  config.shift,
                                  step);
  }

  // Step 3: Attempt to identify adapters based on known sequences
  if (config.adapter_selection_strategy == adapter_selection::automatic) {
    const auto database = config.known_adapters();

    step = sch.add<adapter_finalizer>(database,
                                      config.samples,
                                      config.adapter_fallback_strategy,
                                      step);
    step = sch.add<adapter_selector>(database,
                                     *config.simd.get_reader(),
                                     config.mismatch_threshold,
                                     step,
                                     config.max_threads);
    step = sch.add<adapter_preselector>(step);
  }

  // Step 2: Post-processing, validate, and collect statistics on FASTQ reads
  step = sch.add<post_process_fastq>(config, step, stats);

  // Step 1: Read input file(s)
  read_fastq::add_steps(sch, config, step, stats);

  if (!sch.run(config.max_threads)) {
    return 1;
  }

  const auto out_files = config.get_output_filenames();
  if (!write_json_report(config, stats, out_files.settings_json)) {
    return 1;
  }

  return !write_html_report(config, stats, out_files.settings_html);
}

} // namespace adapterremoval
