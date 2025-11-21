// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_id.hpp"    // for consensus_adapter_stats
#include "commontypes.hpp"   // for read_file, read_file::mate_1, read_typ...
#include "counts.hpp"        // for counts, counts_tmpl, indexed_count
#include "debug.hpp"         // for AR_REQUIRE
#include "errors.hpp"        // for io_error
#include "fastq.hpp"         // for ACGT, fastq, ACGT::values
#include "fastq_enc.hpp"     // for ACGT, ACGTN
#include "json.hpp"          // for json_dict, json_dict_ptr, json_list
#include "logging.hpp"       // for log_stream, error
#include "managed_io.hpp"    // for managed_writer
#include "output.hpp"        // for sample_output_files
#include "reports.hpp"       // for write_json_report
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for adapter_set
#include "statistics.hpp"    // for fastq_stats_ptr, trimming_statistics
#include "strutils.hpp"      // for string_vec, to_lower, indent_lines
#include "timeutils.hpp"     // for format_time, start_time
#include "userconfig.hpp"    // for userconfig, output_files, output_sampl...
#include "version.hpp"       // for version, short_version
#include <cstdint>           // for uint64_t
#include <cstring>           // for size_t, strerror
#include <memory>            // for __shared_ptr_access, shared_ptr, make_...
#include <sstream>           // for basic_ostringstream, basic_ostream, bas...
#include <string>            // for basic_string, string, operator+, char_...
#include <string_view>       // for string_view, operator!=, operator==
#include <utility>           // for pair
#include <vector>            // for vector

namespace adapterremoval {

namespace {

////////////////////////////////////////////////////////////////////////////////
// Meta data

void
write_report_meta(const userconfig& config, json_dict& report)
{
  const auto meta = report.dict("meta");

  meta->str("version", program::long_version());
  meta->str_vec("command", config.args);
  meta->f64("runtime", config.runtime());
  meta->str("timestamp", format_time(start_time(), "%FT%T%z"));
}

////////////////////////////////////////////////////////////////////////////////
// Summary statistics

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

  json->u64("reads", n_reads);
  json->u64("bases", n_bases);
  json->f64("mean_length",
            static_cast<double>(n_bases) / static_cast<double>(n_reads));
  json->u64("reads_sampled", n_reads_s);
  json->f64("q20_rate",
            static_cast<double>(n_q20) / static_cast<double>(n_bases_s));
  json->f64("q30_rate",
            static_cast<double>(n_q30) / static_cast<double>(n_bases_s));
  json->f64("uncalled_rate",
            static_cast<double>(n_n) / static_cast<double>(n_bases_s));
  json->f64("gc_content",
            static_cast<double>(n_g + n_c) / static_cast<double>(n_bases_s));
}

void
write_report_summary(const userconfig& config,
                     json_dict& report,
                     const statistics& stats)
{
  const auto summary = report.dict("summary");

  write_report_summary_stats(summary->dict("input"),
                             { stats.input_1, stats.input_2 });

  if (config.run_type == ar_command::report_only) {
    summary->null("output");
  } else {
    const auto output = summary->dict("output");

    std::vector<fastq_stats_ptr> passed;
    for (const auto& it : stats.trimming) {
      passed.push_back(it->read_1);
      passed.push_back(it->read_2);
      passed.push_back(it->merged);
      passed.push_back(it->singleton);

      // Discarded reads are excluded, even if saved
    }

    write_report_summary_stats(summary->dict("output"), passed);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Input

/** Helper struct used to simplify writing of multiple io sections. */
struct io_section
{
  io_section(read_file rtype,
             std::string name,
             fastq_stats_ptr stats,
             const sample_output_files& sample_files)
    : io_section(std::move(name), std::move(stats), {})
  {
    const auto offset = sample_files.offset(rtype);
    if (offset != sample_output_files::disabled) {
      m_filenames.push_back(sample_files.filename(offset));
    }
  }

  io_section(std::string name, fastq_stats_ptr stats, string_vec filenames)
    : m_stats(std::move(stats))
    , m_name(std::move(name))
    , m_filenames(std::move(filenames))
  {
  }

  void write_to_if(const json_dict_ptr& json, bool enabled = true) const
  {
    if (!enabled) {
      json->null(m_name);
      return;
    }

    const auto section = json->dict(m_name);
    if (m_filenames.empty()) {
      section->null("filenames");
    } else {
      section->str_vec("filenames", m_filenames);
    }

    section->u64("input_reads", m_stats->number_of_input_reads());
    section->u64("output_reads", m_stats->number_of_output_reads());
    section->u64("reads_sampled", m_stats->number_of_sampled_reads());
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

          // FIXME: Should be raw counts instead of fractions
          content_curves->f64_vec(std::string(1, to_lower(nuc)),
                                  bases / total_bases);
        }
      }

      const auto quality_dist = m_stats->quality_dist().trim();
      section->i64_vec("quality_scores", quality_dist);
      section->f64_vec("gc_content", m_stats->gc_content());
    } else {
      section->null("quality_curves");
      section->null("content_curves");
      section->null("quality_scores");
      section->null("gc_content");
    }
  }

private:
  const fastq_stats_ptr m_stats;
  //! Name of the section
  const std::string m_name;
  //! Filenames (if any) generated for this file type
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

  io_section("read1", stats.input_1, config.input_files_1).write_to_if(input);
  io_section("read2", stats.input_2, mate_2_filenames)
    .write_to_if(input, config.paired_ended_mode);
}

////////////////////////////////////////////////////////////////////////////////
// Demultiplexing

void
write_report_demultiplexing(const userconfig& config,
                            const sample_set& samples,
                            json_dict& report,
                            const statistics& sample_stats)
{
  const bool demux_only = config.run_type == ar_command::demultiplex_only;
  const auto out_files = config.get_output_filenames();

  if (config.is_demultiplexing_enabled()) {
    const auto demultiplexing = report.dict("demultiplexing");

    const auto& demux = *sample_stats.demultiplexing;
    size_t assigned_reads = 0;
    for (const auto& it : demux.samples) {
      assigned_reads += it.sum();
    }

    demultiplexing->u64("assigned_reads", assigned_reads);
    demultiplexing->u64("ambiguous_reads", demux.ambiguous);
    demultiplexing->u64("unassigned_reads", demux.unidentified);

    const auto sample_dict = demultiplexing->dict("samples");
    for (size_t i = 0; i < demux.samples.size(); ++i) {
      const auto sample = sample_dict->dict(samples.at(i).name());
      const auto& stats = *sample_stats.trimming.at(i);
      const auto& files = out_files.get_sample(i);

      const auto& barcodes = samples.at(i);
      const auto barcode_list = sample->list("barcodes");
      for (size_t j = 0; j < barcodes.size(); ++j) {
        const auto it = barcodes.at(j);

        const auto dict = barcode_list->inline_dict();
        dict->str("barcode1", it.barcode_1.as_string());
        dict->str("barcode2", it.barcode_2.as_string());

        switch (it.orientation) {
          case barcode_orientation::unspecified:
            dict->null("orientation");
            break;
          case barcode_orientation::forward:
            dict->str("orientation", "forward");
            break;
          case barcode_orientation::reverse:
            dict->str("orientation", "reverse");
            break;
          default:
            AR_FAIL("invalid barcode orientation");
        }

        dict->u64("reads", demux.samples.at(i).get(j));
      }

      const auto output = sample->dict("output");
      io_section(read_file::mate_1, "read1", stats.read_1, files)
        .write_to_if(output, true);

      io_section(read_file::mate_2, "read2", stats.read_2, files)
        .write_to_if(output, config.paired_ended_mode);

      io_section(read_file::singleton, "singleton", stats.singleton, files)
        .write_to_if(output,
                     config.paired_ended_mode && !demux_only &&
                       config.is_any_filtering_enabled());

      io_section(read_file::merged, "merged", stats.merged, files)
        .write_to_if(output, config.is_read_merging_enabled());

      io_section(read_file::discarded, "discarded", stats.discarded, files)
        .write_to_if(output, !demux_only && config.is_any_filtering_enabled());
    }
  } else {
    report.null("demultiplexing");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Processing

//! The kind of action performed by a processing step
enum class feature_type
{
  //! Quality or complexity filtering
  filter,
  //! Read merging and error correction
  merge,
  //! Adapter, quality, or poly-X trimming
  trim,
};

std::string_view
feature_name(const feature_type action)
{
  switch (action) {
    case feature_type::filter:
      return "filter";
    case feature_type::merge:
      return "merge";
    case feature_type::trim:
      return "trim";
    default:
      AR_FAIL("invalid processing step type");
  }
}

//! Basic processing/filtering statistics
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
write_report_count(const json_list_ptr& json,
                   const feature_type action,
                   const feature_stats& it)
{
  if (it.enabled) {
    const auto dict = json->inline_dict();

    dict->str("step", it.key);
    dict->str("action", feature_name(action));
    dict->u64("reads", it.count.reads());
    dict->u64("bases", it.count.bases());
  }
}

void
write_report_counts(const json_list_ptr& json,
                    const feature_type action,
                    const std::vector<feature_stats>& stats)
{
  for (const auto& it : stats) {
    write_report_count(json, action, it);
  }
}

void
write_report_poly_x(json_dict& json,
                    const std::string_view step,
                    const std::string_view nucleotides,
                    const indexed_count<ACGT>& reads,
                    const indexed_count<ACGT>& bases)
{
  json.str("step", step);
  json.str("action", "trim");
  json.u64("reads", reads.sum());
  json.u64("bases", bases.sum());

  const auto dict = json.dict("x");
  for (const auto nuc : nucleotides) {
    const auto nuc_stats = dict->inline_dict(std::string(1, to_lower(nuc)));

    nuc_stats->u64("reads", reads.get(nuc));
    nuc_stats->u64("bases", bases.get(nuc));
  }
}

void
write_report_processing(const userconfig& config,
                        const sample_set& samples,
                        json_dict_ptr report,
                        const statistics& stats)
{
  if (!config.is_adapter_trimming_enabled()) {
    report->null("processing");
    return;
  }

  trimming_statistics totals;
  for (const auto& it : stats.trimming) {
    totals += *it;
  }

  auto json = report->list("processing");
  write_report_count(json,
                     feature_type::trim,
                     { "terminal_pre",
                       config.is_terminal_base_pre_trimming_enabled(),
                       totals.terminal_pre_trimmed });

  if (config.is_poly_x_tail_pre_trimming_enabled()) {
    write_report_poly_x(*json->dict(),
                        "poly_x_pre",
                        config.pre_trim_poly_x,
                        totals.poly_x_pre_trimmed_reads,
                        totals.poly_x_pre_trimmed_bases);
  }

  {
    AR_REQUIRE(totals.adapter_trimmed_reads.size() ==
               totals.adapter_trimmed_bases.size());

    uint64_t reads = 0;
    uint64_t bases = 0;
    for (size_t i = 0; i < totals.adapter_trimmed_reads.size(); ++i) {
      reads += totals.adapter_trimmed_reads.get(i);
      bases += totals.adapter_trimmed_bases.get(i);
    }

    const auto dict = json->dict();

    dict->str("step", "adapters");
    dict->str("action", "trim");
    dict->u64("reads", reads);
    dict->u64("bases", bases);

    const auto adapters = samples.adapters().to_read_orientation();
    const auto adapter_list = dict->list("adapter_list");
    for (size_t i = 0; i < adapters.size(); ++i) {
      const auto adapter = adapter_list->inline_dict();

      adapter->str("adapter1", adapters.at(i).first.as_string());
      adapter->str("adapter2", adapters.at(i).second.as_string());
      adapter->u64("reads", totals.adapter_trimmed_reads.get(i));
      adapter->u64("bases", totals.adapter_trimmed_bases.get(i));
    }
  }

  if (config.is_read_merging_enabled()) {
    const auto dict = json->inline_dict();

    dict->str("step", "merging");
    dict->str("action", feature_name(feature_type::merge));
    dict->u64("reads", totals.reads_merged.reads());
    dict->u64("bases", totals.reads_merged.bases());
    dict->u64("corrected", totals.mismatches_resolved);
    dict->u64("masked", totals.mismatches_unresolved);
    dict->u64("unmasked", totals.ns_resolved);
  }

  write_report_count(json,
                     feature_type::trim,
                     { "terminal_post",
                       config.is_terminal_base_post_trimming_enabled(),
                       totals.terminal_post_trimmed });

  if (config.is_poly_x_tail_post_trimming_enabled()) {
    write_report_poly_x(*json->dict(),
                        "poly_x_post",
                        config.post_trim_poly_x,
                        totals.poly_x_post_trimmed_reads,
                        totals.poly_x_post_trimmed_bases);
  }

  write_report_count(json,
                     feature_type::trim,
                     { "low_quality",
                       config.is_low_quality_trimming_enabled(),
                       totals.low_quality_trimmed });

  // Filtering is (currently) performed after trimming
  write_report_counts(json,
                      feature_type::filter,
                      { { "min_length",
                          config.is_short_read_filtering_enabled(),
                          totals.filtered_min_length },
                        { "max_length",
                          config.is_long_read_filtering_enabled(),
                          totals.filtered_max_length },
                        { "ambiguous_bases",
                          config.is_ambiguous_base_filtering_enabled(),
                          totals.filtered_ambiguous },
                        { "mean_quality",
                          config.is_mean_quality_filtering_enabled(),
                          totals.filtered_mean_quality },
                        { "low_complexity",
                          config.is_low_complexity_filtering_enabled(),
                          totals.filtered_low_complexity } });
}

////////////////////////////////////////////////////////////////////////////////
// Analyses

void
write_report_duplication(json_dict_ptr json,
                         const std::string_view key,
                         const duplication_stats_ptr& stats)
{
  AR_REQUIRE(stats);
  auto duplication = json->dict(key);

  const auto summary = stats->summarize();
  duplication->str_vec("labels", summary.labels);
  duplication->f64_vec("unique_sequences", summary.unique_sequences);
  duplication->f64_vec("total_sequences", summary.total_sequences);
  duplication->f64("unique_frac", summary.unique_frac);
}

void
write_report_consensus_adapter(json_dict_ptr json,
                               const std::string_view key,
                               const consensus_adapter_stats& stats)
{
  auto dict = json->dict(key);

  const auto adapter = stats.summarize();
  dict->str("consensus", adapter.adapter().sequence());
  dict->str("qualities", adapter.adapter().qualities());

  auto kmer_dict = dict->dict("kmers");
  for (const auto& it : adapter.top_kmers()) {
    kmer_dict->u64(it.first, it.second);
  }
}

void
write_report_analyses(const userconfig& config,
                      json_dict_ptr json,
                      const statistics& stats)
{
  json = json->dict("analyses");

  if (config.report_duplication) {
    auto dict = json->dict("duplication");

    write_report_duplication(dict, "read1", stats.duplication_1);
    if (config.paired_ended_mode) {
      write_report_duplication(dict, "read2", stats.duplication_2);
    } else {
      dict->null("read2");
    }
  } else {
    json->null("duplication");
  }

  if (config.paired_ended_mode) {
    counts total_insert_sizes;
    for (const auto& it : stats.trimming) {
      total_insert_sizes += it->insert_sizes;
    }

    json->i64_vec("insert_sizes", total_insert_sizes);
  } else {
    json->null("insert_sizes");
  }

  if (stats.adapter_id) {
    auto consensus = json->dict("consensus_adapters");

    consensus->u64("aligned_pairs", stats.adapter_id->aligned_pairs);
    consensus->u64("pairs_with_adapters",
                   stats.adapter_id->pairs_with_adapters);

    write_report_consensus_adapter(consensus,
                                   "read1",
                                   stats.adapter_id->adapter1);
    write_report_consensus_adapter(consensus,
                                   "read2",
                                   stats.adapter_id->adapter2);
  } else {
    json->null("consensus_adapters");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Output

string_vec
collect_files(const output_files& files, read_file rtype)
{
  string_vec filenames;

  for (const auto& sample_files : files.samples()) {
    const auto offset = sample_files.offset(rtype);
    if (offset != sample_output_files::disabled) {
      filenames.push_back(sample_files.filename(offset));
    }
  }

  return filenames;
}

string_vec
filter_output_file(output_file file)
{
  string_vec out;
  if (file.name != DEV_NULL) {
    out.emplace_back(file.name);
  }

  return out;
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
  const auto mate_1_files = collect_files(out_files, read_file::mate_1);
  const auto mate_2_files = collect_files(out_files, read_file::mate_2);
  const auto merged_files = collect_files(out_files, read_file::merged);
  const auto singleton_files = collect_files(out_files, read_file::singleton);
  const auto discarded_files = collect_files(out_files, read_file::discarded);

  const bool demux_only = config.run_type == ar_command::demultiplex_only;

  const auto output = report.dict("output");
  io_section("read1", output_1, mate_1_files).write_to_if(output, true);
  io_section("read2", output_2, mate_2_files)
    .write_to_if(output, config.paired_ended_mode);
  io_section("merged", merged, merged_files)
    .write_to_if(output, config.is_read_merging_enabled());

  io_section("unidentified1",
             stats.demultiplexing->unidentified_stats_1,
             filter_output_file(out_files.unidentified_1))
    .write_to_if(output, config.is_demultiplexing_enabled());
  io_section("unidentified2",
             stats.demultiplexing->unidentified_stats_2,
             filter_output_file(out_files.unidentified_2))
    .write_to_if(output,
                 config.is_demultiplexing_enabled() &&
                   config.paired_ended_mode);

  io_section("singleton", singleton, singleton_files)
    .write_to_if(output,
                 config.paired_ended_mode && !demux_only &&
                   config.is_any_filtering_enabled());

  io_section("discarded", discarded, discarded_files)
    .write_to_if(output, !demux_only && config.is_any_filtering_enabled());
}

} // namespace

bool
write_json_report(const userconfig& config,
                  const statistics& stats,
                  std::string_view filename)
{
  if (filename == DEV_NULL) {
    // User disabled the report
    return true;
  }

  std::ostringstream output;

  {
    std::ostringstream url;
    url << "https://MikkelSchubert.github.io/adapterremoval/schemas/"
        << program::short_version() << ".json";

    auto samples = config.samples.get_reader();
    json_dict_ptr report = std::make_shared<json_dict>();
    report->str("$schema", url.str());

    write_report_meta(config, *report);
    write_report_summary(config, *report, stats);
    write_report_input(config, *report, stats);
    write_report_demultiplexing(config, *samples, *report, stats);
    write_report_processing(config, *samples, report, stats);
    write_report_analyses(config, report, stats);
    write_report_output(config, *report, stats);
    report->write(output);
  }

  output << std::endl;

  try {
    managed_writer writer{ std::string{ filename } };
    writer.write(output.str());
    writer.close();
  } catch (const io_error& error) {
    log::error() << "Error writing JSON report to '" << filename << "':\n"
                 << indent_lines(error.what());
    return false;
  }

  return true;
}

} // namespace adapterremoval
