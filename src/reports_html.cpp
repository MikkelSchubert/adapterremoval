/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <cstring> // for size_t, strerror
#include <errno.h> // for errno
#include <fstream> // for ofstream
#include <iomanip> // for setprecision
#include <memory>  // for make_shared, make_unique
#include <sstream> // for ostringstream
#include <string>  // for operator+, string, operator<<
#include <vector>  // for vector

#include "adapterset.hpp"            // for adapter_set
#include "counts.hpp"                // for counts, counts_tmpl
#include "debug.hpp"                 // for AR_FAIL
#include "fastq.hpp"                 // for fastq_pair_vec, IDX_TO_ACGT, fastq
#include "json.hpp"                  // for json_writer, json_section
#include "logging.hpp"               // for log
#include "main.hpp"                  // for NAME, VERSION
#include "reports_template_html.hpp" // for template strings
#include "statistics.hpp"            // for fastq_statistics, ...
#include "strutils.hpp"              // for cli_formatter, ...
#include "userconfig.hpp"            // for userconfig, ...

namespace adapterremoval {

using fastq_stats_vec = std::vector<fastq_stats_ptr>;
using template_ptr = std::unique_ptr<html_template>;

//! Size chosen to allow fitting two pages side-by-side on a 1920 width display
const char* FIGURE_WIDTH = "736";
//! Per figure width for two-column facet figures; approximate
const char* FACET_WIDTH_2 = "351";
//! Per figure width for one-column facet figures; approximate
const char* FACET_WIDTH_1 = FIGURE_WIDTH;

////////////////////////////////////////////////////////////////////////////////

namespace {

/** Escapes a string that needs to be embedded in a JS */
std::string
json_encode(const std::string& s)
{
  return json_token::from_str(s)->to_string();
}

/** JSON escaped string */
std::string operator""_json(const char* s, size_t length)
{
  return json_encode(std::string(s, length));
}

std::string
runtime_to_str(double seconds)
{
  std::ostringstream ss;

  if (seconds >= 3600.0) {
    ss << static_cast<size_t>(seconds / 3600.0) << " "
       << (seconds >= 7200.0 ? "hours, " : "hour, ") << std::setw(2);
  }

  if (seconds >= 60.0) {
    size_t minutes = static_cast<size_t>(std::fmod(seconds, 3600.0) / 60.0);
    ss << minutes << " "
       << ((!minutes || minutes >= 120.0) ? "minutes" : "minute") << ", and "
       << std::setw(4);
  }

  ss << std::fixed << std::setprecision(1) << std::fmod(seconds, 60.0)
     << " seconds";

  return ss.str();
}

std::string
mean_of_counts(const counts& count)
{
  auto reads = count.sum();
  auto bases = count.product();

  if (!reads) {
    return "NA";
  }

  if (bases % reads == 0) {
    return std::to_string(bases / reads);
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(1)
     << (bases / static_cast<double>(reads));

  return ss.str();
}

/**
 * VEGA-lite will omit plots if there are no values; this function therefore
 * ensures that at least one value is written for a given measurement.
 */
template<typename T>
counts_tmpl<T>
require_values(counts_tmpl<T> r, T fallback = T())
{
  if (r.size()) {
    return r;
  }

  return counts_tmpl<T>({ fallback });
}

std::string
format_average_bases(const reads_and_bases& counts)
{
  return format_fraction(counts.bases(), counts.reads());
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

class io_summary_writer
{
public:
  enum class io
  {
    input,
    output
  };

  io_summary_writer(const std::string& title, io type)
    : m_writer()
    , m_type(type)
  {
    m_writer.set_title(title);
  }

  void write(std::ofstream& output) { m_writer.write(output); }

  void add_column(const std::string& title, const fastq_statistics& stats)
  {
    m_writer.add_columns(title);

    const auto n_reads = (m_type == io::input) ? stats.number_of_input_reads()
                                               : stats.number_of_output_reads();
    m_writer.add_n_reads(format_rough_number(n_reads));

    m_writer.add_n_bases(format_rough_number(stats.length_dist().product()));

    m_writer.add_lengths(mean_of_counts(stats.length_dist()));

    auto total = stats.quality_dist().sum();
    m_writer.add_q30(format_percentage(stats.quality_dist().sum(30), total));
    m_writer.add_q20(format_percentage(stats.quality_dist().sum(20), total));
    m_writer.add_ns(format_percentage(stats.nucleotides_pos('N').sum(), total));
    m_writer.add_gc(format_percentage(stats.nucleotides_gc_pos().sum(), total));
  }

private:
  html_summary_io m_writer;
  io m_type;
};

std::string
build_base_qualities(const fastq_stats_vec& reads, const string_vec& names)
{
  json_list qualities;

  for (size_t i = 0; i < reads.size(); ++i) {
    const auto& stats = *reads.at(i);

    auto total_quality = stats.qualities_pos();
    auto total_bases = stats.nucleotides_pos();

    for (const auto nuc : ACGT::values) {
      const auto nucleotides = stats.nucleotides_pos(nuc);
      const auto quality = stats.qualities_pos(nuc);

      auto dict = qualities.dict();
      dict->str("read", names.at(i));
      dict->i64("offset", 1);
      dict->str("group", std::string(1, ::toupper(nuc)));
      dict->f64_vec("y", quality / nucleotides);
    }

    auto dict = qualities.dict();
    dict->str("read", names.at(i));
    dict->i64("offset", 1);
    dict->str("group", "Mean");

    // Ensure that values get written, to prevent the plot being omitted
    dict->f64_vec("y", require_values(total_quality / total_bases));
  }

  return qualities.to_string();
}

std::string
build_quality_distribution(const fastq_stats_vec& reads,
                           const string_vec& names)
{
  json_list data;

  for (size_t i = 0; i < reads.size(); ++i) {
    const auto& stats = reads.at(i);
    auto count = stats->quality_dist().trim();
    // A max that should give a uniform look to most data
    count.resize_up_to(44);

    const auto m = data.dict();
    m->str("group", names.at(i));
    m->i64("offset", 0);
    m->f64_vec("y", count.normalize());
  }

  return data.to_string();
}

std::string
build_base_content(const fastq_stats_vec& reads, const string_vec& names)
{
  json_list content;

  for (size_t i = 0; i < reads.size(); ++i) {
    const auto& stats = *reads.at(i);

    auto total_bases = stats.nucleotides_pos();

    for (const auto nuc : ACGTN::values) {
      const auto bases = stats.nucleotides_pos(nuc);

      const auto dict = content.dict();
      dict->str("read", names.at(i));
      dict->i64("offset", 1);
      dict->str("group", std::string(1, nuc));

      // Ensure that values get written, to prevent the plot being omitted
      dict->f64_vec("y", require_values(bases / total_bases));
    }

    {
      const auto gc_content = stats.nucleotides_gc_pos();
      auto dict = content.dict();
      dict->str("read", names.at(i));
      dict->i64("offset", 1);
      dict->str("group", "GC");

      // Ensure that values get written, to prevent the plot being omitted
      dict->f64_vec("y", require_values(gc_content / total_bases));
    }
  }

  return content.to_string();
}

////////////////////////////////////////////////////////////////////////////////
// Main sections

void
write_html_summary_section(const userconfig& config,
                           const statistics& stats,
                           std::ofstream& output)
{
  html_head().set_name(NAME).set_version(VERSION).write(output);

  html_body_start().write(output);

  // Basic information about the executable / call
  {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::ostringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");

    html_summary()
      .set_date_and_time(ss.str())
      .set_version(VERSION)
      .set_command(shell_escape_command(config.args))
      .set_runtime(runtime_to_str(config.runtime()))
      .write(output);
  }

  fastq_statistics output_1;
  fastq_statistics output_2;
  fastq_statistics merged;
  fastq_statistics singleton;
  fastq_statistics discarded;

  for (const auto& it : stats.trimming) {
    output_1 += *it->read_1;
    output_2 += *it->read_2;
    merged += *it->merged;
    singleton += *it->singleton;
    discarded += *it->discarded;
  }

  if (config.paired_ended_mode) {
    // Summary statistics for input files
    {
      fastq_statistics totals;
      totals += *stats.input_1;
      totals += *stats.input_2;

      io_summary_writer summary("Input", io_summary_writer::io::input);

      if (config.paired_ended_mode) {
        summary.add_column("Summary", totals);
        summary.add_column("File 1", *stats.input_1);
        summary.add_column("File 2", *stats.input_2);
      }

      summary.write(output);
    }

    // Summary statistics for output files
    {
      fastq_statistics totals;
      totals += output_1;
      totals += output_2;
      totals += merged;
      totals += singleton;
      // discarded reads not counted in the output
      // totals += discarded;

      io_summary_writer summary("Output", io_summary_writer::io::output);
      summary.add_column("Passed*", totals);

      if (config.paired_ended_mode) {
        summary.add_column("File 1", output_1);
        summary.add_column("File 2", output_2);

        if (config.merge) {
          summary.add_column("Merged", merged);
        }

        summary.add_column("Singleton", singleton);
      }

      summary.add_column("Discarded*", discarded);
      summary.write(output);

      // Note regarding passed / discarded reads
      html_output_note_pe().write(output);
    }
  } else {
    io_summary_writer summary("Input/Output summary",
                              io_summary_writer::io::input);

    summary.add_column("Input", *stats.input_1);
    summary.add_column("Output", output_1);
    summary.add_column("Discarded*", discarded);

    summary.write(output);

    // Note regarding discarded reads in output
    html_output_note_se().write(output);
  }
}

//! Trimming statistics
struct trimming_stats
{
  //! Processing stage relative to adapter trimming (pre, X, post)
  std::string stage;
  //! Row label 1 (step)
  std::string label_1;
  //! Row label 1 (sub-step)
  std::string label_2;
  //! Whether or not this step is enabled by command-line options
  bool enabled;
  //! Number of reads/bases trimmed/filtered
  reads_and_bases count;
};

void
write_html_trimming_stats(std::ofstream& output,
                          const std::vector<trimming_stats>& stats)
{
  size_t n_enabled = 0;
  for (const auto& it : stats) {
    n_enabled += it.enabled;
  }

  html_summary_trimming_head().write(output);

  if (!n_enabled) {
    return;
  }

  std::string previous_stage;
  std::string previous_label_1;

  reads_and_bases total;
  for (const auto& it : stats) {
    if (it.enabled) {
      const auto label_1 = it.label_1 == previous_label_1 ? "" : it.label_1;
      const auto stage = it.stage == previous_stage ? "" : it.stage;

      previous_stage = it.stage;
      previous_label_1 = it.label_1;

      html_summary_trimming_row()
        .set_stage(stage)
        .set_label_1(label_1)
        .set_label_2(it.label_2)
        .set_reads(format_rough_number(it.count.reads()))
        .set_bases(format_rough_number(it.count.bases()))
        .set_avg_bases(format_average_bases(it.count))
        .write(output);

      total += it.count;
    }
  }

  if (n_enabled > 1) {
    // FIXME: Total is wrong for trimming, since reads get counted N times
    html_summary_trimming_row()
      .set_stage("")
      .set_label_1("")
      .set_label_2("")
      .set_reads("FIXME")
      .set_bases(format_rough_number(total.bases()))
      .set_avg_bases("FIXME")
      .write(output);
  }

  html_summary_trimming_tail()
    .set_n_enabled(std::to_string(n_enabled))
    .set_n_total(std::to_string(stats.size()))
    .write(output);
}

//! Filtering statistics
struct filtering_stats
{
  //! Filtering step
  std::string label;
  //! Whether or not this step is enabled by command-line options
  bool enabled;
  //! Number of reads/bases trimmed/filtered
  reads_and_bases count;
};

void
write_html_filtering_stats(std::ofstream& output,
                           const std::vector<filtering_stats>& stats)
{
  size_t n_enabled = 0;
  for (const auto& it : stats) {
    n_enabled += it.enabled;
  }

  html_summary_filtering_head().write(output);

  if (!n_enabled) {
    return;
  }

  reads_and_bases total;
  for (const auto& it : stats) {
    if (it.enabled) {
      html_summary_filtering_row()
        .set_label(it.label)
        .set_reads(format_rough_number(it.count.reads()))
        .set_bases(format_rough_number(it.count.bases()))
        .set_avg_bases(format_average_bases(it.count))
        .write(output);

      total += it.count;
    }
  }

  if (n_enabled > 1) {
    html_summary_filtering_row()
      .set_label("")
      .set_reads(format_rough_number(total.reads()))
      .set_bases(format_rough_number(total.bases()))
      .set_avg_bases(format_average_bases(total))
      .write(output);
  }

  html_summary_filtering_tail()
    .set_n_enabled(std::to_string(n_enabled))
    .set_n_total(std::to_string(stats.size()))
    .write(output);
}

void
write_html_processing_section(const userconfig& config,
                              const statistics& stats,
                              std::ofstream& output)
{
  if (config.run_type == ar_command::demultiplex_sequences ||
      config.run_type == ar_command::report_only) {
    return;
  }

  trimming_statistics totals;
  for (const auto& it : stats.trimming) {
    totals += *it;
  }

  uint64_t adapter_reads = 0;
  uint64_t adapter_bases = 0;

  for (size_t i = 0; i < config.adapters.adapter_count(); ++i) {
    adapter_reads += totals.adapter_trimmed_reads.get(i);
    adapter_bases += totals.adapter_trimmed_bases.get(i);
  }

  // Trimming steps prior to adapter trimming
  std::vector<trimming_stats> trimming = {
    { "Pre",
      "Terminal bases",
      std::string(),
      config.is_terminal_base_pre_trimming_enabled(),
      totals.terminal_pre_trimmed },
  };

  for (const auto nucleotide : config.pre_trim_poly_x) {
    trimming.push_back(
      { "Pre",
        "Poly-X tails",
        std::string(1, nucleotide),
        true,
        reads_and_bases(totals.poly_x_pre_trimmed_reads.get(nucleotide),
                        totals.poly_x_pre_trimmed_bases.get(nucleotide)) });
  }

  trimming.push_back({ "Main",
                       "Adapters",
                       "-",
                       config.is_adapter_trimming_enabled(),
                       reads_and_bases(adapter_reads, adapter_bases) });

  trimming.push_back({ "Post",
                       "Terminal bases",
                       "-",
                       config.is_terminal_base_post_trimming_enabled(),
                       totals.terminal_post_trimmed });

  for (const auto nucleotide : config.post_trim_poly_x) {
    trimming.push_back(
      { "Post",
        "Poly-X tails",
        std::string(1, nucleotide),
        true,
        reads_and_bases(totals.poly_x_post_trimmed_reads.get(nucleotide),
                        totals.poly_x_post_trimmed_bases.get(nucleotide)) });
  }

  trimming.push_back({ "Post",
                       "Low quality bases",
                       "-",
                       config.is_low_quality_trimming_enabled(),
                       totals.low_quality_trimmed });

  write_html_trimming_stats(output, trimming);

  write_html_filtering_stats(output,
                             { { "Short reads",
                                 config.is_short_read_filtering_enabled(),
                                 totals.filtered_min_length },
                               { "Long reads",
                                 config.is_long_read_filtering_enabled(),
                                 totals.filtered_max_length },
                               { "Ambiguous reads",
                                 config.is_ambiguous_base_filtering_enabled(),
                                 totals.filtered_ambiguous },
                               { "Low complexity reads",
                                 config.is_low_complexity_filtering_enabled(),
                                 totals.filtered_low_complexity } });
}

void
write_html_section_title(const std::string& title, std::ofstream& output)
{
  html_h2_tag().set_title(title).write(output);
}

void
write_html_io_section(const userconfig& config,
                      const fastq_stats_vec& statistics,
                      const string_vec& names,
                      std::ofstream& output,
                      const fastq_stats_ptr& merged = fastq_stats_ptr())
{
  AR_REQUIRE(statistics.size() == names.size());

  const char* dynamic_width =
    config.paired_ended_mode || merged ? FACET_WIDTH_2 : FACET_WIDTH_1;

  html_facet_line_plot()
    .set_title("Position quality distribution"_json)
    .set_x_axis(config.merge && merged ? "null" : "Position"_json)
    .set_y_axis("Phred score"_json)
    .set_width(dynamic_width)
    .set_values(build_base_qualities(statistics, names))
    .write(output);

  if (config.merge && merged) {
    html_line_plot()
      .set_title(" "_json)
      .set_sub_title("Merged"_json)
      .set_x_axis("Position"_json)
      .set_y_axis("Phred score"_json)
      .set_width(FIGURE_WIDTH)
      .set_values(build_base_qualities({ merged }, { "merged" }))
      .write(output);
  }

  html_facet_line_plot()
    .set_title("Nucleotide content"_json)
    .set_x_axis(config.merge && merged ? "null" : "Position"_json)
    .set_y_axis("Frequency"_json)
    .set_width(dynamic_width)
    .set_values(build_base_content(statistics, names))
    .write(output);

  if (config.merge && merged) {
    html_line_plot()
      .set_title(" "_json)
      .set_sub_title("Merged"_json)
      .set_title_anchor("start"_json)
      .set_x_axis("Position"_json)
      .set_y_axis("Frequency"_json)
      .set_width(FIGURE_WIDTH)
      .set_values(build_base_content({ merged }, { "merged" }))
      .write(output);
  }

  html_line_plot()
    .set_title("Quality score distribution"_json)
    .set_title_anchor("start"_json)
    .set_x_axis("Phred score"_json)
    .set_y_axis("Frequency"_json)
    .set_width(FIGURE_WIDTH)
    .set_values(build_quality_distribution(statistics, names))
    .write(output);

  {
    json_list data;

    for (size_t i = 0; i < statistics.size(); ++i) {
      const auto m = data.dict();
      m->str("group", names.at(i));
      m->i64("offset", 0);
      m->f64_vec("y", statistics.at(i)->gc_content().normalize());
    }

    html_line_plot()
      .set_title("GC Content"_json)
      .set_title_anchor("start"_json)
      .set_x_axis("%GC"_json)
      .set_y_axis("Frequency"_json)
      .set_width(FIGURE_WIDTH)
      .set_values(data.to_string())
      .write(output);
  }
}

void
write_html_input_section(const userconfig& config,
                         const statistics& stats,
                         std::ofstream& output)
{
  fastq_stats_vec stats_vec = { stats.input_1 };
  string_vec names = { "File 1" };

  if (config.paired_ended_mode) {
    stats_vec.push_back(stats.input_2);
    names.push_back("File 2");
  }

  write_html_section_title("Input", output);

  // Insert size distribution
  if (config.paired_ended_mode) {
    counts insert_sizes;
    for (const auto& it : stats.trimming) {
      insert_sizes += it->insert_sizes;
    }

    json_list samples;
    const auto sample = samples.dict();
    sample->str("group", "insert_sizes");
    sample->i64("offset", 0);
    sample->f64_vec("y", insert_sizes.normalize());

    // FIXME: Specify "identified reads" when in demultiplexing mode and
    // correct format_percentage to merged / n_identified.
    std::ostringstream ss;
    ss << "Insert sizes inferred for "
       << format_percentage(insert_sizes.sum(),
                            stats.input_1->number_of_input_reads())
       << "% of reads";

    html_line_plot()
      .set_title("Insert-size distribution"_json)
      .set_title_anchor("start"_json)
      .set_x_axis("Insert size"_json)
      .set_y_axis("Frequency"_json)
      .set_legend("null")
      .set_width(FIGURE_WIDTH)
      .set_sub_title(json_encode(ss.str()))
      .set_values(samples.to_string())
      .write(output);
  }

  write_html_io_section(config, stats_vec, names, output);
}

void
write_html_demultiplexing_section(const userconfig& config,
                                  const statistics& stats,
                                  std::ofstream& output)

{
  if (!config.adapters.barcode_count()) {
    return;
  }

  write_html_section_title("Demultiplexing", output);

  json_list data;

  const size_t input_reads = stats.input_1->number_of_input_reads() +
                             stats.input_2->number_of_input_reads();

  for (size_t i = 0; i < config.adapters.barcode_count(); ++i) {
    auto m = data.dict();
    m->str("x", config.adapters.get_sample_name(i));

    if (input_reads) {
      m->f64("y", (100.0 * stats.demultiplexing->barcodes.at(i)) / input_reads);
    } else {
      m->null("y");
    }
  }

  html_bar_plot()
    .set_title("Samples identified"_json)
    .set_x_axis("Samples"_json)
    .set_y_axis("Percent"_json)
    .set_width(FIGURE_WIDTH)
    .set_values(data.to_string())
    .write(output);

  html_demultiplexing_head().write(output);

  {
    const size_t unidentified = stats.demultiplexing->unidentified;

    fastq_statistics total;
    total += *stats.demultiplexing->unidentified_stats_1;
    total += *stats.demultiplexing->unidentified_stats_2;

    const auto output_reads = total.length_dist().sum();
    const auto output_bp = total.nucleotides_pos().sum();

    html_demultiplexing_row()
      .set_n("")
      .set_barcode_1("")
      .set_barcode_2("")
      .set_name("<b>Unidentified</b>")
      .set_pct(format_percentage(unidentified, input_reads, 2))
      .set_reads(format_rough_number(output_reads))
      .set_bp(format_rough_number(output_bp))
      .set_length(mean_of_counts(total.length_dist()))
      .set_gc(format_percentage(total.nucleotides_gc_pos().sum(), output_bp))
      .write(output);
  }

  const auto barcodes = config.adapters.get_barcodes();

  for (size_t i = 0; i < config.adapters.barcode_count(); ++i) {
    const auto& sample = *stats.trimming.at(i);

    fastq_statistics total;

    total += *sample.read_1;
    total += *sample.read_2;
    total += *sample.merged;
    total += *sample.singleton;
    // Not included in overview:
    // total += *sample.discarded;

    const auto output_reads = total.length_dist().sum();
    const auto output_bp = total.nucleotides_pos().sum();

    html_demultiplexing_row()
      .set_n(std::to_string(i + 1))
      .set_barcode_1(barcodes.at(i).first.sequence())
      .set_barcode_2(barcodes.at(i).second.sequence())
      .set_name(config.adapters.get_sample_name(i))
      .set_pct(
        format_percentage(stats.demultiplexing->barcodes.at(i), input_reads, 2))
      .set_reads(format_rough_number(output_reads))
      .set_bp(format_rough_number(output_bp))
      .set_length(mean_of_counts(total.length_dist()))
      .set_gc(format_percentage(total.nucleotides_gc_pos().sum(), output_bp))
      .write(output);
  }

  html_demultiplexing_tail().write(output);
}

void
write_html_output_section(const userconfig& config,
                          const statistics& stats,
                          std::ofstream& output)

{
  fastq_stats_vec stats_vec;
  string_vec names;

  fastq_stats_ptr merged = std::make_shared<fastq_statistics>();

  {
    fastq_stats_ptr output_1 = std::make_shared<fastq_statistics>();
    fastq_stats_ptr output_2 = std::make_shared<fastq_statistics>();
    fastq_stats_ptr singleton = std::make_shared<fastq_statistics>();
    fastq_stats_ptr discarded = std::make_shared<fastq_statistics>();

    for (const auto& it : stats.trimming) {
      *output_1 += *it->read_1;
      *output_2 += *it->read_2;
      *merged += *it->merged;
      *singleton += *it->singleton;
      *discarded += *it->discarded;
    }

    stats_vec.push_back(output_1);
    names.push_back("Output 1");

    if (config.paired_ended_mode) {
      stats_vec.push_back(output_2);
      stats_vec.push_back(singleton);

      names.push_back("Output 2");
      names.push_back("Singleton");
    }

    stats_vec.push_back(discarded);
    names.push_back("Discarded");
  }

  write_html_section_title("Output", output);
  write_html_io_section(config, stats_vec, names, output, merged);
}

////////////////////////////////////////////////////////////////////////////////

bool
write_html_report(const userconfig& config,
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

    write_html_summary_section(config, stats, output);
    write_html_processing_section(config, stats, output);
    write_html_input_section(config, stats, output);
    write_html_demultiplexing_section(config, stats, output);
    write_html_output_section(config, stats, output);

    html_body_end().write(output);
  } catch (const std::ios_base::failure& error) {
    log::error() << "Error writing JSON report to '" << filename << "':\n"
                 << indent_lines(error.what());
    return false;
  }

  return true;
}

} // namespace adapterremoval
