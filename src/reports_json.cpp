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
#include <cstring>  // for size_t, strerror
#include <errno.h>  // for errno
#include <fstream>  // for ofstream
#include <iostream> // for ofstream, operator<<, basic_ostream, endl
#include <string>   // for operator+, string, operator<<
#include <vector>   // for vector

#include "adapterset.hpp" // for adapter_set
#include "counts.hpp"     // for counts, counts_tmpl
#include "debug.hpp"      // for AR_FAIL
#include "fastq.hpp"      // for fastq_pair_vec, IDX_TO_ACGT, fastq
#include "json.hpp"       // for json_writer, json_section
#include "main.hpp"       // for NAME, VERSION
#include "statistics.hpp" // for fastq_statistics, trimming_statistics, ar_...
#include "strutils.hpp"   // for cli_formatter
#include "userconfig.hpp" // for userconfig, ar_command, ar_command::demult...
#include "utilities.hpp"  // for make_shared

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
    n_reads += it->number_of_input_reads();
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
write_report_trimming(const userconfig& config,
                      const json_dict_ptr& json,
                      const trimming_statistics& totals,
                      const fastq_pair_vec& adapters)
{
  if (config.run_type == ar_command::demultiplex_sequences ||
      config.run_type == ar_command::report_only) {
    json->null("trimming_and_filtering");
    return;
  }

  const auto trimming = json->dict("trimming_and_filtering");
  const auto adapter_list = trimming->list("adapter_sequences");

  for (size_t i = 0; i < adapters.size(); ++i) {
    const auto adapter = adapter_list->dict();

    adapter->str("adapter_sequence_1", adapters.at(i).first.sequence());
    adapter->str("adapter_sequence_2", adapters.at(i).second.sequence());
    adapter->i64("adapter_trimmed_reads", totals.adapter_trimmed_reads.get(i));
    adapter->i64("adapter_trimmed_bases", totals.adapter_trimmed_bases.get(i));
  }

  trimming->i64("overlapping_reads", totals.overlapping_reads);
  trimming->i64("terminal_bases_trimmed", totals.terminal_bases_trimmed);
  trimming->i64("low_quality_trimmed_reads", totals.low_quality_trimmed_reads);
  trimming->i64("low_quality_trimmed_bases", totals.low_quality_trimmed_bases);
  trimming->i64("filtered_min_length_reads", totals.filtered_min_length_reads);
  trimming->i64("filtered_min_length_bases", totals.filtered_min_length_bases);
  trimming->i64("filtered_max_length_reads", totals.filtered_max_length_reads);
  trimming->i64("filtered_max_length_bases", totals.filtered_max_length_bases);
  trimming->i64("filtered_ambiguous_reads", totals.filtered_ambiguous_reads);
  trimming->i64("filtered_ambiguous_bases", totals.filtered_ambiguous_bases);
  trimming->i64("filtered_low_complexity_reads",
                totals.filtered_low_complexity_reads);
  trimming->i64("filtered_low_complexity_bases",
                totals.filtered_low_complexity_bases);

  if (config.paired_ended_mode) {
    trimming->i64_vec("insert_sizes", totals.insert_sizes);
  } else {
    trimming->null("insert_sizes");
  }
}

void
write_report_summary(const userconfig& config,
                     json_dict& report,
                     const statistics& stats)
{
  const bool demux_only = config.run_type == ar_command::demultiplex_sequences;

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
      summary_demux->i64(config.adapters.get_sample_name(i),
                         demux.barcodes.at(i));
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
      m_filenames.push_back(sample_files.filenames.at(offset));
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
      default:
        AR_FAIL("unknown read type");
    }
  }

  void write_to_if(const json_dict_ptr& json, bool enabled = true)
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

    const auto total_bases = m_stats->nucleotides_pos();
    const auto total_quality = m_stats->qualities_pos();

    {
      const auto quality_curves = section->dict("quality_curves");

      for (const auto nuc : ACGT) {
        const auto& nucleotides = m_stats->nucleotides_pos(nuc);
        const auto& quality = m_stats->qualities_pos(nuc);

        quality_curves->f64_vec(std::string(1, tolower(nuc)),
                                quality / nucleotides);
      }

      quality_curves->f64_vec("mean", total_quality / total_bases);
    }

    {
      const auto content_curves = section->dict("content_curves");

      for (const auto nuc : ACGTN) {
        const auto& bases = m_stats->nucleotides_pos(nuc);

        content_curves->f64_vec(std::string(1, tolower(nuc)),
                                bases / total_bases);
      }

      content_curves->f64_vec("gc",
                              m_stats->nucleotides_gc_pos() / total_bases);
    }

    const auto quality_dist = m_stats->quality_dist().trim();
    section->i64_vec("quality_scores", quality_dist);
    section->i64_vec("gc_content", m_stats->gc_content());

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
  const bool demux_only = config.run_type == ar_command::demultiplex_sequences;
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
        .write_to_if(output, config.paired_ended_mode);

      io_section(read_type::merged, stats.merged, files)
        .write_to_if(output, config.merge && !demux_only);

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
      filenames.push_back(sample_files.filenames.at(offset));
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

  fastq_stats_ptr output_1 = make_shared<fastq_statistics>();
  fastq_stats_ptr output_2 = make_shared<fastq_statistics>();
  fastq_stats_ptr merged = make_shared<fastq_statistics>();
  fastq_stats_ptr singleton = make_shared<fastq_statistics>();
  fastq_stats_ptr discarded = make_shared<fastq_statistics>();

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

  const bool demux_only = config.run_type == ar_command::demultiplex_sequences;

  const auto output = report.dict("output");
  io_section(read_type::mate_1, output_1, mate_1_files)
    .write_to_if(output, true);
  io_section(read_type::mate_2, output_2, mate_2_files)
    .write_to_if(output, config.paired_ended_mode);
  io_section(read_type::merged, merged, merged_files)
    .write_to_if(output, config.merge && !demux_only);

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
    .write_to_if(output, !demux_only);

  io_section(read_type::discarded, discarded, discarded_files)
    .write_to_if(output, !demux_only);
}

bool
write_json_report(const userconfig& config,
                  const statistics& stats,
                  const std::string& filename)
{
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
    std::cerr << "Error writing JSON report to '" << filename << "':\n"
              << cli_formatter::fmt(error.what()) << std::endl;
    return false;
  }

  return true;
}
