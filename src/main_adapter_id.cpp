/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "adapter_id.hpp" // for adapter_id_statistics
#include "alignment.hpp"  // for extract_adapter_sequences, sequence_aligner
#include "debug.hpp"      // for AR_REQUIRE
#include "fastq.hpp"      // for ACGTN, fastq, fastq_pair_vec, ACGT, ACGT:...
#include "fastq_io.hpp"   // for read_fastq, read_chunk
#include "reports.hpp"    // for write_html_report, write_json_report
#include "scheduler.hpp"  // for threadstate, scheduler, analytical_step
#include "userconfig.hpp" // for userconfig
#include <cstddef>        // for size_t
#include <string>         // for string, operator<<, char_traits

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

  chunk_vec process(chunk_ptr chunk) override
  {
    AR_REQUIRE(chunk);

    const fastq empty_adapter("dummy", "", "");
    fastq_pair_vec adapters;
    adapters.emplace_back(empty_adapter, empty_adapter);

    auto aligner = sequence_aligner(adapters, m_config.simd);

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

      if (m_config.is_good_alignment(alignment)) {
        stats_id->aligned_pairs++;
        if (m_config.can_merge_alignment(alignment)) {
          const size_t insert_size = alignment.insert_size(read_1, read_2);
          stats_ins->insert_sizes.resize_up_to(insert_size + 1);
          stats_ins->insert_sizes.inc(insert_size);
          stats_ins->overlapping_reads += 2;

          if (extract_adapter_sequences(alignment, read_1, read_2)) {
            stats_id->pairs_with_adapters++;

            stats_id->adapter1.process(read_1.sequence());
            read_2.reverse_complement();
            stats_id->adapter2.process(read_2.sequence());
          }
        }
      }
    }

    m_stats_id.release(stats_id);
    m_stats_ins.release(stats_ins);

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
identify_adapter_sequences(const userconfig& config)
{
  scheduler sch;

  statistics stats = statistics_builder()
                       .sample_rate(config.report_sample_rate)
                       .estimate_duplication(config.report_duplication)
                       // FIXME: Length picked based on known sequences
                       .adapter_identification(42)
                       .initialize();
  // FIXME: Required for insert size statistics
  stats.trimming.push_back(std::make_shared<trimming_statistics>());

  // Step 3:
  size_t final_step;
  if (config.paired_ended_mode) {
    // Attempt to identify adapters through pair-wise alignments
    final_step = sch.add<adapter_identification>(config, stats);
  } else {
    // Discard all written reads
    final_step = sch.add<reads_sink>();
  }

  // Step 2: Post-processing, validate, and collect statistics on FASTQ reads
  const size_t postproc_step =
    sch.add<post_process_fastq>(config, final_step, stats);

  // Step 1: Read input file(s)
  read_fastq::add_steps(sch, config, postproc_step);

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
