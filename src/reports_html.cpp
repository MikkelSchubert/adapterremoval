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
#include "adapter_id.hpp" // for adapter_id_statistics
#include "adapterset.hpp" // for adapter_set
#include "counts.hpp"     // for counts, indexed_count, counts_tmpl
#include "debug.hpp"      // for AR_REQUIRE
#include "fastq.hpp"      // for ACGT, ACGT::values, fastq, ACGTN
#include "json.hpp"       // for json_dict, json_list, json_ptr
#include "logging.hpp"    // for log_stream, error
#include "main.hpp"       // for VERSION, NAME
#include "managed_io.hpp"
#include "reports.hpp"               // for write_html_report
#include "reports_template_html.hpp" // for html_line_plot, html_demultiple...
#include "simd.hpp"                  // for size_t
#include "statistics.hpp"            // for fastq_stats_ptr, fastq_statistics
#include "strutils.hpp"              // for format_percentage, format_rough...
#include "userconfig.hpp"            // for userconfig, ar_command, DEV_NULL
#include <algorithm>                 // for max
#include <array>                     // for array
#include <cctype>                    // for toupper
#include <cerrno>                    // for errno
#include <chrono>                    // for system_clock
#include <cmath>                     // for fmod
#include <cstdint>                   // for uint64_t
#include <cstring>                   // for size_t, strerror
#include <ctime>                     // for localtime
#include <fstream>                   // for operator<<, ofstream, basic_ost...
#include <iomanip>                   // for operator<<, setprecision, setw
#include <memory>                    // for __shared_ptr_access, shared_ptr
#include <sstream>                   // for ostringstream
#include <string>                    // for string, operator==, to_string
#include <utility>                   // for pair
#include <vector>                    // for vector

namespace adapterremoval {

using fastq_stats_vec = std::vector<fastq_stats_ptr>;
using template_ptr = std::unique_ptr<html_template>;

//! Size chosen to allow fitting two pages side-by-side on a 1920 width display
const char* const FIGURE_WIDTH = "736";
//! Per figure width for two-column facet figures; approximate
const char* const FACET_WIDTH_2 = "351";
//! Per figure width for one-column facet figures; approximate
const char* const FACET_WIDTH_1 = FIGURE_WIDTH;

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
    auto minutes = static_cast<size_t>(std::fmod(seconds, 3600.0) / 60.0);
    ss << minutes << " "
       << ((!minutes || minutes >= 120) ? "minutes" : "minute") << ", and "
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
  return format_fraction(counts.bases(), counts.reads(), 1);
}

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

  void write(std::ostream& output) { m_writer.write(output); }

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

    for (const auto nucleotide : ACGT::values) {
      const auto nucleotides = stats.nucleotides_pos(nucleotide);
      const auto quality = stats.qualities_pos(nucleotide);

      auto dict = qualities.dict();
      dict->str("read", names.at(i));
      dict->i64("offset", 1);
      dict->str("group", std::string(1, ::toupper(nucleotide)));
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

    for (const auto nucleotide : ACGTN::values) {
      const auto bases = stats.nucleotides_pos(nucleotide);

      const auto dict = content.dict();
      dict->str("read", names.at(i));
      dict->i64("offset", 1);
      dict->str("group", std::string(1, nucleotide));

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
write_html_sampling_note(const userconfig& config,
                         const std::string& label,
                         const fastq_statistics& stats,
                         std::ostream& output)
{
  if (config.report_sample_rate < 1.0) {
    html_sampling_note()
      .set_label(label)
      .set_reads(format_rough_number((stats.number_of_sampled_reads())))
      .set_pct(format_percentage(stats.number_of_sampled_reads(),
                                 stats.number_of_input_reads()))
      .write(output);
  }
}

void
write_html_summary_section(const userconfig& config,
                           const statistics& stats,
                           std::ostream& output)
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

      write_html_sampling_note(config, "input", totals, output);
    }

    // Summary statistics for output files
    if (config.run_type != ar_command::report_only) {
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

        if (config.is_read_merging_enabled()) {
          summary.add_column("Merged", merged);
        }

        summary.add_column("Singleton", singleton);
      }

      summary.add_column("Discarded*", discarded);
      summary.write(output);

      write_html_sampling_note(config, "output", totals, output);

      // Note regarding passed / discarded reads
      html_output_footnote()
        .set_symbol("*")
        .set_text("The <b>Passed</b> column includes all read types except for "
                  "<b>Discarded</b> reads.")
        .write(output);
    }
  } else if (config.run_type == ar_command::report_only) {
    io_summary_writer summary("Input summary", io_summary_writer::io::input);
    summary.add_column("Input", *stats.input_1);
    summary.write(output);

    write_html_sampling_note(config, "input", *stats.input_1, output);
  } else {
    io_summary_writer summary("Input/Output summary",
                              io_summary_writer::io::input);

    summary.add_column("Input", *stats.input_1);
    summary.add_column("Output", output_1);
    summary.add_column("Discarded*", discarded);
    summary.write(output);

    fastq_statistics totals;
    totals += *stats.input_1;
    totals += output_1;

    write_html_sampling_note(config, "input/output", totals, output);

    // Note regarding discarded reads in output
    html_output_footnote()
      .set_symbol("*")
      .set_text(
        "<b>Discarded</b> reads are not included in the <b>Output</b> column.")
      .write(output);
  }
}

//! Trimming statistics
struct trimming_stats
{
  size_t id;
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
write_html_trimming_stats(std::ostream& output,
                          const std::vector<trimming_stats>& stats,
                          const reads_and_bases& totals)
{
  size_t n_processing_steps = 0;
  size_t n_processing_steps_on = 0;
  size_t n_filtering_steps = 0;
  size_t n_filtering_steps_on = 0;

  size_t last_id = -1;
  size_t last_enabled = -1;
  for (const auto& it : stats) {
    if (it.id != last_id) {
      if (it.stage == "Processing") {
        n_processing_steps++;
      } else if (it.stage == "Filtering") {
        n_filtering_steps++;
      }

      last_id = it.id;
    }

    if (it.enabled && it.id != last_enabled) {
      if (it.stage == "Processing") {
        n_processing_steps_on++;
      } else if (it.stage == "Filtering") {
        n_filtering_steps_on++;
      }

      last_enabled = it.id;
    }
  }

  html_summary_trimming_head().write(output);

  std::string previous_stage;
  std::string previous_label_1;

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
        .set_pct_reads(format_percentage(it.count.reads(), totals.reads()))
        .set_bases(format_rough_number(it.count.bases()))
        .set_pct_bases(format_percentage(it.count.bases(), totals.bases()))
        .set_avg_bases(format_average_bases(it.count))
        .write(output);
    }
  }

  html_summary_trimming_tail()
    .set_n_enabled_filt(std::to_string(n_filtering_steps_on))
    .set_n_total_filt(std::to_string(n_filtering_steps))
    .set_n_enabled_proc(std::to_string(n_processing_steps_on))
    .set_n_total_proc(std::to_string(n_processing_steps))
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

reads_and_bases
summarize_input(const fastq_stats_ptr& ptr)
{
  const auto n_bases = ptr->length_dist().product();
  AR_REQUIRE(n_bases >= 0);

  return { ptr->number_of_input_reads(), static_cast<uint64_t>(n_bases) };
}

void
build_polyx_trimming_rows(std::vector<trimming_stats>& out,
                          const std::string& polyx_nucleotides,
                          const indexed_count<ACGT>& reads,
                          const indexed_count<ACGT>& bases,
                          const size_t id)
{
  for (const auto nucleotide : ACGT::values) {
    out.push_back(
      { id,
        "Processing",
        "Poly-X tails",
        std::string(1, nucleotide),
        polyx_nucleotides.find(nucleotide) != std::string::npos,
        reads_and_bases(reads.get(nucleotide), bases.get(nucleotide)) });
  }

  out.push_back({ id,
                  "Processing",
                  "Poly-X tails",
                  "*",
                  polyx_nucleotides.size() > 1,
                  reads_and_bases(reads.sum(), bases.sum()) });
}

void
write_html_processing_section(const userconfig& config,
                              const statistics& stats,
                              std::ostream& output)
{
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

  const auto total_input =
    summarize_input(stats.input_1) + summarize_input(stats.input_2);

  reads_and_bases total_output;
  for (const auto& it : stats.trimming) {
    total_output += summarize_input(it->read_1);
    total_output += summarize_input(it->read_2);
    total_output += summarize_input(it->singleton);
    total_output += summarize_input(it->merged);
  }

  // Trimming steps prior to adapter trimming
  size_t step_id = 0;
  std::vector<trimming_stats> trimming = {
    { step_id++, "Input", "Raw reads", "-", true, total_input },
    { step_id++,
      "Processing",
      "Terminal bases",
      "-",
      config.is_terminal_base_pre_trimming_enabled(),
      totals.terminal_pre_trimmed },
  };

  build_polyx_trimming_rows(trimming,
                            config.pre_trim_poly_x,
                            totals.poly_x_pre_trimmed_reads,
                            totals.poly_x_pre_trimmed_bases,
                            step_id++);

  trimming.push_back({ step_id++,
                       "Processing",
                       "Adapters",
                       "-",
                       config.is_adapter_trimming_enabled(),
                       reads_and_bases(adapter_reads, adapter_bases) });

  trimming.push_back({ step_id++,
                       "Processing",
                       "Merging",
                       "-",
                       config.is_read_merging_enabled(),
                       totals.reads_merged });

  trimming.push_back({ step_id++,
                       "Processing",
                       "Terminal bases",
                       "-",
                       config.is_terminal_base_post_trimming_enabled(),
                       totals.terminal_post_trimmed });

  build_polyx_trimming_rows(trimming,
                            config.post_trim_poly_x,
                            totals.poly_x_post_trimmed_reads,
                            totals.poly_x_post_trimmed_bases,
                            step_id++);

  trimming.push_back({ step_id++,
                       "Processing",
                       "Low quality bases",
                       "-",
                       config.is_low_quality_trimming_enabled(),
                       totals.low_quality_trimmed });

  trimming.push_back({ step_id++,
                       "Filtering",
                       "Short reads",
                       "-",
                       config.is_short_read_filtering_enabled(),
                       totals.filtered_min_length });

  trimming.push_back({ step_id++,
                       "Filtering",
                       "Long reads",
                       "-",
                       config.is_long_read_filtering_enabled(),
                       totals.filtered_max_length });
  trimming.push_back({ step_id++,
                       "Filtering",
                       "Ambiguous bases",
                       "-",
                       config.is_ambiguous_base_filtering_enabled(),
                       totals.filtered_ambiguous });
  trimming.push_back({ step_id++,
                       "Filtering",
                       "Low complexity reads",
                       "-",
                       config.is_low_complexity_filtering_enabled(),
                       totals.filtered_low_complexity });

  trimming.push_back(
    { step_id++, "Output", "Filtered reads", "-", true, total_output });

  write_html_trimming_stats(output, trimming, total_input);
}

void
write_html_section_title(const std::string& title, std::ostream& output)
{
  html_h2_tag().set_title(title).write(output);
}

void
write_html_io_section(const userconfig& config,
                      const fastq_stats_vec& statistics,
                      const string_vec& names,
                      std::ostream& output,
                      const fastq_stats_ptr& merged = fastq_stats_ptr())
{
  AR_REQUIRE(statistics.size() == names.size());

  const char* dynamic_width =
    config.paired_ended_mode || merged ? FACET_WIDTH_2 : FACET_WIDTH_1;

  html_facet_line_plot()
    .set_title("Position quality distribution"_json)
    .set_x_axis(config.is_read_merging_enabled() && merged ? "null"
                                                           : "Position"_json)
    .set_y_axis("Phred score"_json)
    .set_width(dynamic_width)
    .set_values(build_base_qualities(statistics, names))
    .write(output);

  if (config.is_read_merging_enabled() && merged) {
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
    .set_x_axis(config.is_read_merging_enabled() && merged ? "null"
                                                           : "Position"_json)
    .set_y_axis("Frequency"_json)
    .set_width(dynamic_width)
    .set_values(build_base_content(statistics, names))
    .write(output);

  if (config.is_read_merging_enabled() && merged) {
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
                         std::ostream& output)
{
  fastq_stats_vec stats_vec = { stats.input_1 };
  string_vec names = { "File 1" };

  if (config.paired_ended_mode) {
    stats_vec.push_back(stats.input_2);
    names.emplace_back("File 2");
  }

  write_html_section_title("Input", output);

  write_html_io_section(config, stats_vec, names, output);
}

void
write_html_analyses_section(const userconfig& config,
                            const statistics& stats,
                            std::ostream& output)

{
  write_html_section_title("Analyses", output);

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

    if (config.run_type == ar_command::report_only) {
      html_output_note()
        .set_text(
          "Insert size distribution inferred using adapter-free alignments.")
        .write(output);
    }
  }

  // Consensus adapter sequence inference
  if (config.paired_ended_mode && config.run_type == ar_command::report_only) {
    AR_REQUIRE(stats.adapter_id);

    const auto adapter_1 = stats.adapter_id->adapter1.summarize();
    const auto adapter_2 = stats.adapter_id->adapter2.summarize();

    // Consensus adapter sequences
    {
      const auto reference_adapters =
        config.adapters.get_raw_adapters().front();
      const auto& reference_adapter_1 = reference_adapters.first.sequence();
      const auto& reference_adapter_2 = reference_adapters.second.sequence();

      html_consensus_adapter_head()
        .set_overlapping_pairs(
          format_rough_number(stats.adapter_id->aligned_pairs))
        .set_pairs_with_adapters(
          format_rough_number(stats.adapter_id->pairs_with_adapters))
        .write(output);

      html_consensus_adapter_table()
        .set_name_1("--adapter1")
        .set_reference_1(reference_adapter_1)
        .set_alignment_1(adapter_1.compare_with(reference_adapter_1))
        .set_consensus_1(adapter_1.adapter().sequence())
        .set_qualities_1(adapter_1.adapter().qualities())
        .set_name_2("--adapter2")
        .set_reference_2(reference_adapter_2)
        .set_alignment_2(adapter_2.compare_with(reference_adapter_2))
        .set_consensus_2(adapter_2.adapter().sequence())
        .set_qualities_2(adapter_2.adapter().qualities())
        .write(output);
    }

    // Top N most common 5' kmers in adapter fragments
    {
      const auto& top_kmers_1 = adapter_1.top_kmers();
      const auto& top_kmers_2 = adapter_2.top_kmers();

      html_consensus_adapter_kmer_head()
        .set_n_kmers(std::to_string(consensus_adapter_stats::top_n_kmers))
        .set_kmer_length(std::to_string(consensus_adapter_stats::kmer_length))
        .write(output);

      const auto kmers = std::max(top_kmers_1.size(), top_kmers_2.size());
      for (size_t i = 0; i < kmers; ++i) {
        html_consensus_adapter_kmer_row row;
        row.set_index(std::to_string(i + 1));

        if (top_kmers_1.size() > i) {
          const auto& kmer = top_kmers_1.at(i);

          row.set_kmer_1(kmer.first)
            .set_count_1(format_rough_number(kmer.second))
            .set_pct_1(format_percentage(kmer.second, adapter_1.total_kmers()));
        } else {
          row.set_kmer_1("").set_count_1("").set_pct_1("");
        }

        if (top_kmers_2.size() > i) {
          const auto& kmer = top_kmers_2.at(i);

          row.set_kmer_2(kmer.first)
            .set_count_2(format_rough_number(kmer.second))
            .set_pct_2(format_percentage(kmer.second, adapter_2.total_kmers()));
        } else {
          row.set_kmer_2("").set_count_2("").set_pct_2("");
        }

        row.write(output);
      }

      html_consensus_adapter_kmer_tail().write(output);
    }
  }
}

void
write_html_demultiplexing_section(const userconfig& config,
                                  const statistics& stats,
                                  std::ostream& output)

{
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
                          std::ostream& output)

{
  fastq_stats_vec stats_vec;
  string_vec names;

  auto merged = std::make_shared<fastq_statistics>();

  {
    auto output_1 = std::make_shared<fastq_statistics>();
    auto output_2 = std::make_shared<fastq_statistics>();
    auto singleton = std::make_shared<fastq_statistics>();
    auto discarded = std::make_shared<fastq_statistics>();

    for (const auto& it : stats.trimming) {
      *output_1 += *it->read_1;
      *output_2 += *it->read_2;
      *merged += *it->merged;
      *singleton += *it->singleton;
      *discarded += *it->discarded;
    }

    stats_vec.push_back(output_1);
    names.emplace_back("Output 1");

    if (config.paired_ended_mode) {
      stats_vec.push_back(output_2);
      stats_vec.push_back(singleton);

      names.emplace_back("Output 2");
      names.emplace_back("Singleton");
    }

    stats_vec.push_back(discarded);
    names.emplace_back("Discarded");
  }

  write_html_section_title("Output", output);
  write_html_io_section(config, stats_vec, names, output, merged);
}

} // namespace

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

  std::ostringstream output;

  write_html_summary_section(config, stats, output);

  if (config.run_type != ar_command::demultiplex_only &&
      config.run_type != ar_command::report_only) {
    write_html_processing_section(config, stats, output);
  }

  write_html_input_section(config, stats, output);

  if (config.run_type == ar_command::report_only) {
    write_html_analyses_section(config, stats, output);
  }

  if (config.adapters.barcode_count()) {
    write_html_demultiplexing_section(config, stats, output);
  }

  if (config.run_type != ar_command::report_only) {
    write_html_output_section(config, stats, output);
  }

  html_body_end().write(output);

  try {
    managed_writer writer{ filename };
    writer.write(output.str());
    writer.close();
  } catch (const std::ios_base::failure& error) {
    log::error() << "Error writing JSON report to '" << filename << "':\n"
                 << indent_lines(error.what());
    return false;
  }

  return true;
}

} // namespace adapterremoval
