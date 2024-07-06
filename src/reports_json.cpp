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
#include "adapterset.hpp"  // for adapter_set
#include "commontypes.hpp" // for read_type, read_type::mate_1, read_typ...
#include "counts.hpp"      // for counts, counts_tmpl, indexed_count
#include "debug.hpp"       // for AR_FAIL
#include "fastq.hpp"       // for ACGT, fastq_pair_vec, fastq, ACGT::values
#include "json.hpp"        // for json_dict, json_dict_ptr, json_list
#include "logging.hpp"     // for log_stream, error
#include "main.hpp"        // for NAME, VERSION
#include "reports.hpp"     // for write_json_report
#include "simd.hpp"        // for size_t
#include "statistics.hpp"  // for fastq_stats_ptr, trimming_statistics
#include "strutils.hpp"    // for string_vec, to_lower, indent_lines
#include "userconfig.hpp"  // for userconfig, output_files, output_sampl...
#include <algorithm>       // for max
#include <array>           // for array
#include <cstring>         // for size_t, strerror
#include <errno.h>         // for errno
#include <fstream>         // for ofstream, ios_base::failure, operator|
#include <memory>          // for __shared_ptr_access, shared_ptr, make_...
#include <string>          // for basic_string, string, operator+, char_...
#include <utility>         // for pair
#include <vector>          // for vector

namespace adapterremoval {

//! Trimming statistics
struct feature_stats
{
  //! Processing stage name
  std::string key;
  //! Whether or not this step is enabled by command-line options
  bool enabled;
  //! Number of reads/bases trimmed/filtered
  reads_and_bases count;
};

void
write_report_meta(const userconfig& config, json_dict& report)
{
  const auto meta = report.dict("meta");

  meta->str("version", NAME + " " + VERSION);
  meta->str_vec("command", config.args);
  meta->f64("runtime", config.runtime());
}

void
write_report_summary_stats(const json_dict_ptr& json,
                           const std::vector<fastq_stats_ptr>& stats)
{
  size_t n_reads = 0;
  size_t n_reads_s = 0;
  size_t n_bases = 0;
  size_t n_a = 0;
  size_t n_c = 0;
  size_t n_g = 0;
  size_t n_t = 0;
  size_t n_n = 0;
  size_t n_q20 = 0;
  size_t n_q30 = 0;

  for (const auto& it : stats) {
    n_reads += it->number_of_output_reads();
    n_reads_s += it->number_of_sampled_reads();
    n_bases += it->length_dist().product();
    // The following stats are all (potentially) based on a subset of reads
    n_a += it->nucleotides_pos('A').sum();
    n_c += it->nucleotides_pos('C').sum();
    n_g += it->nucleotides_pos('G').sum();
    n_t += it->nucleotides_pos('T').sum();
    n_n += it->nucleotides_pos('N').sum();
    n_q20 += it->quality_dist().sum(20);
    n_q30 += it->quality_dist().sum(30);
  }

  const auto n_bases_s = n_a + n_c + n_g + n_t + n_n;

  json->i64("reads", n_reads);
  json->i64("bases", n_bases);
  json->f64("mean_length", static_cast<double>(n_bases) / n_reads);
  json->i64("reads_sampled", n_reads_s);
  json->f64("q20_rate", static_cast<double>(n_q20) / n_bases_s);
  json->f64("q30_rate", static_cast<double>(n_q30) / n_bases_s);
  json->f64("uncalled_rate", static_cast<double>(n_n) / n_bases_s);
  json->f64("gc_content", static_cast<double>(n_g + n_c) / n_bases_s);
}

void
write_report_counts(const json_dict_ptr& json,
                    const std::vector<feature_stats>& stats)
{
  for (const auto& it : stats) {
    if (it.enabled) {
      const auto dict = json->inline_dict(it.key);

      dict->i64("reads", it.count.reads());
      dict->i64("bases", it.count.bases());
    }
  }
}

//! Poly-X trimming statistics
struct poly_x_stats
{
  //! Processing stage name
  std::string key;
  //! X trimmed
  std::string nucleotides;
  //! Number of reads trimmed
  indexed_count<ACGT> reads;
  //! Number of bases trimmed
  indexed_count<ACGT> bases;
};

void
write_report_poly_x(const json_dict_ptr& json,
                    const std::vector<poly_x_stats>& stats)
{
  for (const auto& it : stats) {
    if (!it.nucleotides.empty()) {
      const auto dict = json->dict(it.key);
      for (const auto nuc : it.nucleotides) {
        const auto nuc_stats = dict->inline_dict(to_lower(std::string(1, nuc)));
        nuc_stats->i64("reads", it.reads.get(nuc));
        nuc_stats->i64("bases", it.bases.get(nuc));
      }
    }
  }
}

void
write_report_trimming(const userconfig& config,
                      const json_dict_ptr& json,
                      const trimming_statistics& totals,
                      const fastq_pair_vec& adapters)
{
  if (!config.is_adapter_trimming_enabled()) {
    json->null("adapter_trimming");
    json->null("quality_trimming");
    json->null("poly_x_trimming");
    json->null("filtering");
    return;
  }

  {
    const auto trimming = json->dict("adapter_trimming");
    const auto adapter_list = trimming->list("adapters");

    for (size_t i = 0; i < adapters.size(); ++i) {
      const auto adapter = adapter_list->dict();

      adapter->str("sequence_1", adapters.at(i).first.sequence());
      adapter->str("sequence_2", adapters.at(i).second.sequence());
      adapter->i64("reads", totals.adapter_trimmed_reads.get(i));
      adapter->i64("bases", totals.adapter_trimmed_bases.get(i));
    }

    trimming->i64("overlapping_reads", totals.overlapping_reads);
    if (config.paired_ended_mode) {
      trimming->i64_vec("insert_sizes", totals.insert_sizes);
    } else {
      trimming->null("insert_sizes");
    }
  }

  write_report_counts(json->dict("quality_trimming"),
                      { { "terminal_pre",
                          config.is_terminal_base_pre_trimming_enabled(),
                          totals.terminal_pre_trimmed },
                        { "terminal_post",
                          config.is_terminal_base_post_trimming_enabled(),
                          totals.terminal_post_trimmed },
                        { "low_quality",
                          config.is_low_quality_trimming_enabled(),
                          totals.low_quality_trimmed } });

  write_report_poly_x(json->dict("poly_x_trimming"),
                      { { "pre",
                          config.pre_trim_poly_x,
                          totals.poly_x_pre_trimmed_reads,
                          totals.poly_x_pre_trimmed_bases },
                        { "post",
                          config.post_trim_poly_x,
                          totals.poly_x_post_trimmed_reads,
                          totals.poly_x_post_trimmed_bases } });

  write_report_counts(json->dict("filtering"),
                      { { "min_length",
                          config.is_short_read_filtering_enabled(),
                          totals.filtered_min_length },
                        { "max_length",
                          config.is_long_read_filtering_enabled(),
                          totals.filtered_max_length },
                        { "ambiguous_bases",
                          config.is_ambiguous_base_filtering_enabled(),
                          totals.filtered_ambiguous },
                        { "low_complexity",
                          config.is_low_complexity_filtering_enabled(),
                          totals.filtered_low_complexity } });
}

void
write_report_summary(const userconfig& config,
                     json_dict& report,
                     const statistics& stats)
{
  const bool demux_only = config.run_type == ar_command::demultiplex_only;

  const auto summary = report.dict("summary");

  write_report_summary_stats(summary->dict("input"),
                             { stats.input_1, stats.input_2 });

  if (config.adapters.barcode_count()) {
    const auto summary_demux = summary->dict("demultiplexing");
    const auto& demux = *stats.demultiplexing;
    const auto total = demux.total();

    summary_demux->i64("total_reads", total);
    summary_demux->i64("assigned_reads",
                       total - demux.unidentified - demux.ambiguous);
    summary_demux->i64("ambiguous_reads", demux.ambiguous);
    summary_demux->i64("unassigned_reads", demux.unidentified);

    const auto samples = summary_demux->dict("samples");
    for (size_t i = 0; i < demux.barcodes.size(); ++i) {
      samples->i64(config.adapters.get_sample_name(i), demux.barcodes.at(i));
    }
  } else {
    summary->null("demultiplexing");
  }

  {
    trimming_statistics totals;
    for (const auto& it : stats.trimming) {
      totals += *it;
    }

    write_report_trimming(
      config, summary, totals, config.adapters.get_raw_adapters());
  }

  if (config.run_type == ar_command::report_only) {
    summary->null("output");
  } else {
    const auto output = summary->dict("output");

    std::vector<fastq_stats_ptr> passed;
    for (const auto& it : stats.trimming) {
      passed.push_back(it->read_1);
      passed.push_back(it->read_2);
      passed.push_back(it->merged);
    }

    write_report_summary_stats(output->dict("passed"), passed);

    if (config.adapters.barcode_count()) {
      write_report_summary_stats(
        output->dict("unidentified"),
        { stats.demultiplexing->unidentified_stats_1,
          stats.demultiplexing->unidentified_stats_2 });
    } else {
      output->null("unidentified");
    }

    if (demux_only) {
      output->null("discarded");
    } else {
      std::vector<fastq_stats_ptr> discarded;
      for (const auto& it : stats.trimming) {
        discarded.push_back(it->discarded);
      }

      write_report_summary_stats(output->dict("discarded"), discarded);
    }
  }
}

/** Helper struct used to simplify writing of multiple io sections. */
struct io_section
{
  io_section(read_type rtype,
             const fastq_stats_ptr& stats,
             const output_sample_files& sample_files)
    : io_section(rtype, stats, string_vec())
  {
    const auto offset = sample_files.offset(rtype);
    if (offset != output_sample_files::disabled) {
      m_filenames.push_back(sample_files.filenames().at(offset));
    }
  }

  io_section(read_type rtype,
             const fastq_stats_ptr& stats,
             const string_vec& filenames)
    : m_read_type(rtype)
    , m_stats(stats)
    , m_filenames(filenames)
  {
  }

  const char* name() const
  {
    switch (m_read_type) {
      case read_type::mate_1:
        return "read1";
      case read_type::mate_2:
        return "read2";
      case read_type::merged:
        return "merged";
      case read_type::singleton:
        return "singleton";
      case read_type::discarded:
        return "discarded";
      case read_type::unidentified_1:
        return "unidentified_1";
      case read_type::unidentified_2:
        return "unidentified_2";
      case read_type::max:
        AR_FAIL("unsupported read type");
      default:
        AR_FAIL("invalid read type");
    }
  }

  void write_to_if(const json_dict_ptr& json, bool enabled = true) const
  {
    if (!enabled) {
      json->null(name());
      return;
    }

    const auto section = json->dict(name());
    if (m_filenames.empty()) {
      section->null("filenames");
    } else {
      section->str_vec("filenames", m_filenames);
    }

    section->i64("input_reads", m_stats->number_of_input_reads());
    section->i64("output_reads", m_stats->number_of_output_reads());
    section->i64("reads_sampled", m_stats->number_of_sampled_reads());
    section->i64_vec("lengths", m_stats->length_dist());

    if (m_stats->length_dist().product()) {
      const auto total_bases = m_stats->nucleotides_pos();
      const auto total_quality = m_stats->qualities_pos();

      {
        const auto quality_curves = section->dict("quality_curves");

        for (const auto nuc : ACGT::values) {
          const auto nucleotides = m_stats->nucleotides_pos(nuc);
          const auto quality = m_stats->qualities_pos(nuc);

          quality_curves->f64_vec(std::string(1, to_lower(nuc)),
                                  quality / nucleotides);
        }

        quality_curves->f64_vec("mean", total_quality / total_bases);
      }

      {
        const auto content_curves = section->dict("content_curves");

        for (const auto nuc : ACGTN::values) {
          const auto bases = m_stats->nucleotides_pos(nuc);

          content_curves->f64_vec(std::string(1, to_lower(nuc)),
                                  bases / total_bases);
        }
      }

      const auto quality_dist = m_stats->quality_dist().trim();
      section->i64_vec("quality_scores", quality_dist);
      section->f64_vec("gc_content", m_stats->gc_content());

      // Currently only for input 1/2
      const auto dup_stats = m_stats->duplication();

      if (dup_stats) {
        // Must be enabled, but key is always written for applicable files
        if (dup_stats->max_unique()) {
          const auto duplication = section->dict("duplication");

          const auto summary = dup_stats->summarize();
          duplication->str_vec("labels", summary.labels);
          duplication->f64_vec("unique_sequences", summary.unique_sequences);
          duplication->f64_vec("total_sequences", summary.total_sequences);
          duplication->f64("unique_frac", summary.unique_frac);
        } else {
          section->null("duplication");
        }
      }
    } else {
      section->null("quality_curves");
      section->null("content_curves");
      section->null("quality_scores");
      section->null("gc_content");

      // Currently only for input 1/2

      if (m_stats->duplication()) {
        section->null("duplication");
      }
    }
  }

private:
  const read_type m_read_type;

  const fastq_stats_ptr m_stats;
  string_vec m_filenames;
};

void
write_report_input(const userconfig& config,
                   json_dict& report,
                   const statistics& stats)
{
  const auto input = report.dict("input");
  const auto mate_2_filenames =
    config.interleaved_input ? config.input_files_1 : config.input_files_2;

  io_section(read_type::mate_1, stats.input_1, config.input_files_1)
    .write_to_if(input);
  io_section(read_type::mate_2, stats.input_2, mate_2_filenames)
    .write_to_if(input, config.paired_ended_mode);
}

void
write_report_demultiplexing(const userconfig& config,
                            json_dict& report,
                            const statistics& sample_stats)
{
  const bool demux_only = config.run_type == ar_command::demultiplex_only;
  const auto out_files = config.get_output_filenames();

  if (config.adapters.barcode_count()) {
    const auto demultiplexing = report.dict("demultiplexing");

    const auto& demux = *sample_stats.demultiplexing;
    size_t assigned_reads = 0;
    for (size_t it : demux.barcodes) {
      assigned_reads += it;
    }

    demultiplexing->i64("assigned_reads", assigned_reads);
    demultiplexing->i64("ambiguous_reads", demux.ambiguous);
    demultiplexing->i64("unassigned_reads", demux.unidentified);

    const auto samples = demultiplexing->dict("samples");
    for (size_t i = 0; i < demux.barcodes.size(); ++i) {
      const auto sample = samples->dict(config.adapters.get_sample_name(i));
      const auto& stats = *sample_stats.trimming.at(i);
      const auto& files = out_files.samples.at(i);

      sample->i64("reads", demux.barcodes.at(i));

      write_report_trimming(
        config, sample, stats, config.adapters.get_adapter_set(i));

      const auto output = sample->dict("output");
      io_section(read_type::mate_1, stats.read_1, files)
        .write_to_if(output, true);

      io_section(read_type::mate_2, stats.read_2, files)
        .write_to_if(output, config.paired_ended_mode);

      io_section(read_type::singleton, stats.singleton, files)
        .write_to_if(output, config.paired_ended_mode && !demux_only);

      io_section(read_type::merged, stats.merged, files)
        .write_to_if(output, config.is_read_merging_enabled());

      io_section(read_type::discarded, stats.discarded, files)
        .write_to_if(output, !demux_only);
    }
  } else {
    report.null("demultiplexing");
  }
}

string_vec
collect_files(const output_files& files, read_type rtype)
{
  string_vec filenames;

  for (const auto& sample_files : files.samples) {
    const auto offset = sample_files.offset(rtype);
    if (offset != output_sample_files::disabled) {
      filenames.push_back(sample_files.filenames().at(offset));
    }
  }

  return filenames;
}

void
write_report_output(const userconfig& config,
                    json_dict& report,
                    const statistics& stats)
{
  if (config.run_type == ar_command::report_only) {
    report.null("output");
    return;
  }

  auto output_1 = std::make_shared<fastq_statistics>();
  auto output_2 = std::make_shared<fastq_statistics>();
  auto merged = std::make_shared<fastq_statistics>();
  auto singleton = std::make_shared<fastq_statistics>();
  auto discarded = std::make_shared<fastq_statistics>();

  for (const auto& it : stats.trimming) {
    *output_1 += *it->read_1;
    *output_2 += *it->read_2;
    *merged += *it->merged;
    *singleton += *it->singleton;
    *discarded += *it->discarded;
  }

  const auto out_files = config.get_output_filenames();
  const auto mate_1_files = collect_files(out_files, read_type::mate_1);
  const auto mate_2_files = collect_files(out_files, read_type::mate_2);
  const auto merged_files = collect_files(out_files, read_type::merged);
  const auto singleton_files = collect_files(out_files, read_type::singleton);
  const auto discarded_files = collect_files(out_files, read_type::discarded);

  const bool demux_only = config.run_type == ar_command::demultiplex_only;

  const auto output = report.dict("output");
  io_section(read_type::mate_1, output_1, mate_1_files)
    .write_to_if(output, true);
  io_section(read_type::mate_2, output_2, mate_2_files)
    .write_to_if(output, config.paired_ended_mode);
  io_section(read_type::merged, merged, merged_files)
    .write_to_if(output, config.is_read_merging_enabled());

  io_section(read_type::unidentified_1,
             stats.demultiplexing->unidentified_stats_1,
             { out_files.unidentified_1 })
    .write_to_if(output, config.adapters.barcode_count());
  io_section(read_type::unidentified_2,
             stats.demultiplexing->unidentified_stats_2,
             { out_files.unidentified_2 })
    .write_to_if(output,
                 config.adapters.barcode_count() && config.paired_ended_mode);

  io_section(read_type::singleton, singleton, singleton_files)
    .write_to_if(output, config.paired_ended_mode && !demux_only);

  io_section(read_type::discarded, discarded, discarded_files)
    .write_to_if(output, !demux_only);
}

bool
write_json_report(const userconfig& config,
                  const statistics& stats,
                  const std::string& filename)
{
  if (filename == DEV_NULL) {
    // User disabled the report
    return true;
  }

  try {
    std::ofstream output(filename, std::ofstream::out);
    if (!output.is_open()) {
      throw std::ofstream::failure(std::strerror(errno));
    }

    output.exceptions(std::ofstream::failbit | std::ofstream::badbit);

    {
      json_dict report;
      write_report_meta(config, report);
      write_report_summary(config, report, stats);
      write_report_input(config, report, stats);
      write_report_demultiplexing(config, report, stats);
      write_report_output(config, report, stats);
      report.write(output);
    }

    output << std::endl;
  } catch (const std::ios_base::failure& error) {
    log::error() << "Error writing JSON report to '" << filename << "':\n"
                 << indent_lines(error.what());
    return false;
  }

  return true;
}

} // namespace adapterremoval
