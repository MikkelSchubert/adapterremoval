// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_id.hpp"            // for adapter_id_statistics
#include "counts.hpp"                // for counts, indexed_count, counts_tmpl
#include "debug.hpp"                 // for AR_REQUIRE
#include "errors.hpp"                // for io_error
#include "fastq.hpp"                 // for ACGT, ACGT::values, fastq, ACGTN
#include "json.hpp"                  // for json_dict, json_list, json_ptr
#include "logging.hpp"               // for log_stream, error
#include "main.hpp"                  // for VERSION, NAME
#include "managed_io.hpp"            // for managed_io
#include "output.hpp"                // for DEV_NULL, output_files
#include "reports.hpp"               // for write_html_report
#include "reports_template_html.hpp" // for html_frequency_plot, html_demultiple...
#include "sequence_sets.hpp"         // for adapter_set
#include "simd.hpp"                  // for size_t
#include "statistics.hpp"            // for fastq_stats_ptr, fastq_statistics
#include "strutils.hpp"              // for format_percentage, format_rough...
#include "userconfig.hpp"            // for userconfig, ar_command, DEV_NULL
#include <algorithm>                 // for max
#include <cctype>                    // for toupper
#include <cerrno>                    // for errno
#include <cmath>                     // for fmod
#include <cstdint>                   // for uint64_t
#include <cstring>                   // for size_t, strerror
#include <iomanip>                   // for operator<<, setprecision, setw
#include <memory>                    // for __shared_ptr_access, shared_ptr
#include <sstream>                   // for ostringstream
#include <string>                    // for string, operator==, to_string
#include <string_view>               // for string_view
#include <utility>                   // for pair
#include <vector>                    // for vector

namespace adapterremoval {

namespace {

using fastq_stats_vec = std::vector<fastq_stats_ptr>;
using template_ptr = std::unique_ptr<html_template>;

//! Size chosen to allow fitting two pages side-by-side on a 1920 width display
const char* const FIGURE_WIDTH = "736";
//! Per figure width for two-column facet figures; approximate
const char* const FACET_WIDTH_2 = "351";
//! Per figure width for one-column facet figures; approximate
const char* const FACET_WIDTH_1 = FIGURE_WIDTH;

////////////////////////////////////////////////////////////////////////////////

/** Escapes a string that needs to be embedded in a JS */
std::string
json_encode(const std::string& s)
{
  return json_token::from_str(s)->to_string();
}

/** JSON escaped string */
std::string
operator""_json(const char* s, size_t length)
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
mean_of_bp_counts(const counts& count)
{
  auto reads = count.sum();
  auto bases = count.product();

  if (!reads) {
    return "NA";
  }

  if (bases % reads == 0) {
    return std::to_string(bases / reads) + " bp";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(1)
     << (bases / static_cast<double>(reads)) << " bp";

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
  const auto reads = counts.reads();

  if (reads) {
    return format_fraction(counts.bases(), reads, 1) + " bp";
  } else {
    return "NA";
  }
}

std::string
orientation_to_label(const sample_sequences& it)
{
  switch (it.orientation) {
    case barcode_orientation::unspecified:
      return {};
    case barcode_orientation::forward:
      return "+";
    case barcode_orientation::reverse:
      return "-";
    default:
      AR_FAIL("invalid barcode orientation");
  }
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

  io_summary_writer(std::ostream& output, const io type)
    : m_output(output)
    , m_type(type)

  {
  }

  void write_head(const std::string& title, const std::string& href)
  {
    html_summary_io_head().set_title(title).set_href(href).write(m_output);
  }

  void write_row(const std::string& title, const fastq_statistics& stats)
  {
    const auto n_reads = (m_type == io::input) ? stats.number_of_input_reads()
                                               : stats.number_of_output_reads();
    const auto total = stats.quality_dist().sum();

    html_summary_io_row()
      .set_name(title)
      .set_n_reads(format_rough_number(n_reads))
      .set_n_bases(format_rough_number(stats.length_dist().product()))
      .set_lengths(mean_of_bp_counts(stats.length_dist()))
      .set_q30(format_percentage(stats.quality_dist().sum(30), total))
      .set_q20(format_percentage(stats.quality_dist().sum(20), total))
      .set_ns(format_percentage(stats.nucleotides_pos('N').sum(), total))
      .set_gc(format_percentage(stats.nucleotides_gc_pos().sum(), total))
      .write(m_output);
  }

  void write_tail() { html_summary_io_tail().write(m_output); }

private:
  std::ostream& m_output;
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
    m->i64_vec("y", count);
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
  html_head().set_title(config.report_title).write(output);

  html_body_start().set_title(config.report_title).write(output);

  // Basic information about the executable / call
  {
    html_summary()
      .set_date_and_time(userconfig::start_time)
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

      io_summary_writer summary(output, io_summary_writer::io::input);
      summary.write_head("Input", "summary-input");
      if (config.paired_ended_mode) {
        summary.write_row("Summary", totals);
        summary.write_row("File 1", *stats.input_1);
        summary.write_row("File 2", *stats.input_2);
      }
      summary.write_tail();

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

      io_summary_writer summary{ output, io_summary_writer::io::output };
      summary.write_head("Output", "summary-output");
      summary.write_row("Passed*", totals);
      if (config.paired_ended_mode) {
        summary.write_row("File 1", output_1);
        summary.write_row("File 2", output_2);

        if (config.is_read_merging_enabled()) {
          summary.write_row("Merged", merged);
        }

        if (config.is_any_filtering_enabled()) {
          summary.write_row("Singleton", singleton);
        }
      }

      if (config.is_any_filtering_enabled()) {
        summary.write_row("Discarded*", discarded);
      }
      summary.write_tail();

      write_html_sampling_note(config, "output", totals, output);

      // Note regarding passed / discarded reads
      html_output_footnote()
        .set_symbol("*")
        .set_html("The <b>Passed</b> column includes all read types except "
                  "for <b>Discarded</b> reads.")
        .write(output);
    }
  } else if (config.run_type == ar_command::report_only) {
    io_summary_writer summary{ output, io_summary_writer::io::input };
    summary.write_head("Input summary", "summary-input");
    summary.write_row("Input", *stats.input_1);
    summary.write_tail();

    write_html_sampling_note(config, "input", *stats.input_1, output);
  }

  else {
    io_summary_writer summary{ output, io_summary_writer::io::input };
    summary.write_head("Input/Output summary", "summary-input-output");
    summary.write_row("Input", *stats.input_1);
    summary.write_row("Output", output_1);
    if (config.is_any_filtering_enabled()) {
      summary.write_row("Discarded*", discarded);
    }
    summary.write_tail();

    fastq_statistics totals;
    totals += *stats.input_1;
    totals += output_1;

    write_html_sampling_note(config, "input/output", totals, output);

    if (config.is_any_filtering_enabled()) {
      // Note regarding discarded reads in output
      html_output_footnote()
        .set_symbol("*")
        .set_html("<b>Discarded</b> reads are not included in the "
                  "<b>Output</b> column.")
        .write(output);
    }
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

  return reads_and_bases{ ptr->number_of_input_reads(),
                          static_cast<uint64_t>(n_bases) };
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

  for (size_t i = 0; i < config.samples.adapters().size(); ++i) {
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
                       "Mean quality",
                       "-",
                       config.is_mean_quality_filtering_enabled(),
                       totals.filtered_mean_quality });
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
  html_h2_tag().set_title(title).set_href(to_lower(title)).write(output);
}

void
write_html_io_section(const userconfig& config,
                      std::ostream& output,
                      const std::string& title,
                      fastq_stats_vec statistics,
                      string_vec names,
                      const fastq_stats_ptr& merged = fastq_stats_ptr())
{
  AR_REQUIRE(statistics.size() == names.size());

  write_html_section_title(title, output);

  const char* dynamic_width =
    config.paired_ended_mode || merged ? FACET_WIDTH_2 : FACET_WIDTH_1;

  html_plot_title()
    .set_href(to_lower(title) + "-position-qualities")
    .set_title("Position quality distribution")
    .write(output);
  html_facet_line_plot()
    .set_x_axis(config.is_read_merging_enabled() && merged ? "null"
                                                           : "Position"_json)
    .set_y_axis("Phred score"_json)
    .set_width(dynamic_width)
    .set_values(build_base_qualities(statistics, names))
    .write(output);

  if (config.is_read_merging_enabled() && merged) {
    html_facet_line_plot()
      .set_x_axis("Position"_json)
      .set_y_axis("Phred score"_json)
      .set_width(FIGURE_WIDTH)
      .set_values(build_base_qualities({ merged }, { "Merged" }))
      .write(output);
  }

  html_plot_title()
    .set_href(to_lower(title) + "-nucleotide-content")
    .set_title("Nucleotide content")
    .write(output);
  html_facet_line_plot()
    .set_x_axis(config.is_read_merging_enabled() && merged ? "null"
                                                           : "Position"_json)
    .set_y_axis("Frequency"_json)
    .set_width(dynamic_width)
    .set_values(build_base_content(statistics, names))
    .write(output);

  if (config.is_read_merging_enabled() && merged) {
    html_facet_line_plot()
      .set_x_axis("Position"_json)
      .set_y_axis("Frequency"_json)
      .set_width(FIGURE_WIDTH)
      .set_values(build_base_content({ merged }, { "Merged" }))
      .write(output);

    // Subsequent plots should include merged reads
    names.push_back("Merged");
    statistics.push_back(merged);
  }

  html_plot_title()
    .set_href(to_lower(title) + "-quality-scores")
    .set_title("Quality score distribution")
    .write(output);
  html_frequency_plot()
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
      m->f64_vec("y", statistics.at(i)->gc_content());
    }

    html_plot_title()
      .set_href(to_lower(title) + "-gc-content")
      .set_title("GC Content")
      .write(output);
    html_frequency_plot()
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

  write_html_io_section(config,
                        output,
                        "Input",
                        std::move(stats_vec),
                        std::move(names));
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
    sample->i64_vec("y", insert_sizes);

    // FIXME: Specify "identified reads" when in demultiplexing mode and
    // correct format_percentage to merged / n_identified.
    std::ostringstream ss;
    ss << "Insert sizes inferred for "
       << format_percentage(insert_sizes.sum(),
                            stats.input_1->number_of_input_reads())
       << " of reads";

    html_plot_title()
      .set_href("analyses-insert-sizes")
      .set_title("Insert-size distribution")
      .write(output);
    html_plot_sub_title().set_sub_title(ss.str()).write(output);
    html_frequency_plot()
      .set_x_axis("Insert size"_json)
      .set_y_axis("Frequency"_json)
      .set_legend("null")
      .set_width(FIGURE_WIDTH)
      .set_values(samples.to_string())
      .write(output);

    if (config.run_type == ar_command::report_only) {
      html_output_note()
        .set_text(
          "Insert size distribution inferred using adapter-free alignments.")
        .write(output);
    }
  }

  if (config.report_duplication) {
    AR_REQUIRE(stats.duplication_1 && stats.duplication_2);
    const auto dupes_1 = stats.duplication_1->summarize();
    const auto dupes_2 = stats.duplication_2->summarize();
    const auto mean_uniq_frac = (dupes_1.unique_frac + dupes_2.unique_frac) / 2;

    const auto to_percent = [](double value) {
      std::ostringstream os;
      os << std::fixed << std::setprecision(1) << (value * 100.0) << " %";
      return os.str();
    };

    html_duplication_head().write(output);
    if (config.paired_ended_mode) {
      html_duplication_body_pe()
        .set_pct_unique(to_percent(mean_uniq_frac))
        .set_pct_unique_1(to_percent(dupes_1.unique_frac))
        .set_pct_unique_2(to_percent(dupes_2.unique_frac))
        .write(output);
    } else {
      html_duplication_body_se()
        .set_pct_unique(to_percent(dupes_1.unique_frac))
        .write(output);
    }

    const auto add_line = [](json_list& list,
                             std::string_view read,
                             std::string_view group,
                             const std::vector<std::string>& labels,
                             const rates& values) {
      AR_REQUIRE(labels.size() == values.size());
      for (size_t i = 0; i < labels.size(); ++i) {
        auto dict = list.dict();
        dict->str("read", read);
        dict->str("group", group);
        dict->str("x", labels.at(i));
        dict->f64("y", values.get(i));
      }
    };

    json_list data;
    const auto add_lines = [add_line, &data](const decltype(dupes_1)& s,
                                             std::string_view label) {
      add_line(data, label, "All", s.labels, s.total_sequences);
      add_line(data, label, "Unique", s.labels, s.unique_sequences);
    };

    add_lines(dupes_1, "File 1");
    if (config.paired_ended_mode) {
      add_lines(dupes_2, "File 2");
    }

    html_duplication_plot()
      .set_width(config.paired_ended_mode ? FACET_WIDTH_2 : FACET_WIDTH_1)
      .set_values(data.to_string())
      .write(output);
  }

  // Consensus adapter sequence inference
  if (config.paired_ended_mode && config.run_type == ar_command::report_only) {
    AR_REQUIRE(stats.adapter_id);

    const auto adapter_1 = stats.adapter_id->adapter1.summarize();
    const auto adapter_2 = stats.adapter_id->adapter2.summarize();

    // Consensus adapter sequences
    {
      const auto reference_adapters =
        config.samples.adapters().to_read_orientation().front();
      std::string reference_adapter_1{ reference_adapters.first };
      std::string reference_adapter_2{ reference_adapters.second };

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
        }

        if (top_kmers_2.size() > i) {
          const auto& kmer = top_kmers_2.at(i);

          row.set_kmer_2(kmer.first)
            .set_count_2(format_rough_number(kmer.second))
            .set_pct_2(format_percentage(kmer.second, adapter_2.total_kmers()));
        }

        row.write(output);
      }

      html_consensus_adapter_kmer_tail().write(output);
    }
  }
}

void
write_html_demultiplexing_barplot(const userconfig& config,
                                  const statistics& stats,
                                  std::ostream& output)
{
  json_list data;

  const size_t input_reads = stats.input_1->number_of_input_reads() +
                             stats.input_2->number_of_input_reads();

  for (size_t i = 0; i < config.samples.size(); ++i) {
    const auto& sample = config.samples.at(i);

    for (size_t j = 0; j < sample.size(); ++j) {
      auto count = stats.demultiplexing->samples.at(i).get(j);

      const auto& sequences = sample.at(j);
      std::string key{ sequences.barcode_1 };
      if (!sequences.barcode_2.empty()) {
        key.push_back('-');
        key.append(sequences.barcode_2);
      }

      auto m = data.dict();
      m->i64("n", j + 1);
      m->str("barcodes", key);

      if (sequences.orientation != barcode_orientation::unspecified) {
        m->str("orientation", orientation_to_label(sequences));
      }

      m->str("sample", sample.name());

      if (input_reads) {
        m->f64("pct", (100.0 * count) / input_reads);
      } else {
        m->null("pct");
      }
    }
  }

  html_plot_title()
    .set_href("demux-samples")
    .set_title("Samples identified")
    .write(output);
  html_bar_plot()
    .set_x_axis("Samples"_json)
    .set_y_axis("Percent"_json)
    .set_width(FIGURE_WIDTH)
    .set_values(data.to_string())
    .write(output);
}

void
write_html_demultiplexing_table(const userconfig& config,
                                const statistics& stats,
                                std::ostream& output,
                                const bool multiple_barcodes,
                                const bool mixed_orientation)
{
  const size_t input_reads = stats.input_1->number_of_input_reads() +
                             stats.input_2->number_of_input_reads();

  html_demultiplexing_table_head()
    .set_orientation(mixed_orientation ? "<th></th>" : "")
    .write(output);

  {
    const size_t unidentified = stats.demultiplexing->unidentified;

    fastq_statistics total;
    total += *stats.demultiplexing->unidentified_stats_1;
    total += *stats.demultiplexing->unidentified_stats_2;

    const auto output_reads = total.length_dist().sum();
    const auto output_bp = total.nucleotides_pos().sum();

    html_demultiplexing_row()
      .set_name("<b>Unidentified</b>")
      .set_sample_pct(format_percentage(unidentified, input_reads, 2))
      .set_reads(format_rough_number(output_reads))
      .set_bp(format_rough_number(output_bp))
      .set_length(mean_of_bp_counts(total.length_dist()))
      .set_gc(format_percentage(total.nucleotides_gc_pos().sum(), output_bp))
      .set_orientation(mixed_orientation ? "<td></td>" : "")
      .write(output);
  }

  size_t sample_idx = 0;
  for (const auto& sample : config.samples) {
    const auto& output_stats = *stats.trimming.at(sample_idx);
    const auto& barcode_counts = stats.demultiplexing->samples.at(sample_idx);
    const auto sample_reads = barcode_counts.sum();

    fastq_statistics total;

    total += *output_stats.read_1;
    total += *output_stats.read_2;
    total += *output_stats.merged;
    total += *output_stats.singleton;
    // Not included in overview:
    // total += *sample.discarded;

    const auto output_reads = total.length_dist().sum();
    const auto output_bp = total.nucleotides_pos().sum();

    html_demultiplexing_row row;
    if (sample.size() < 2) {
      const auto& it = sample.at(0);
      row.set_barcode_1(std::string{ it.barcode_1 })
        .set_barcode_2(std::string{ it.barcode_2 });

      if (mixed_orientation) {
        row.set_orientation("<td>" + orientation_to_label(it) + "</td>");
      }
    } else {
      const auto cell = "<i>" + std::to_string(sample.size()) + " barcodes</i>";
      row.set_barcode_1(cell).set_barcode_2(cell);

      if (mixed_orientation) {
        row.set_orientation("<td></td>");
      }
    }

    row.set_n(std::to_string(sample_idx + 1))
      .set_name(sample.name())
      .set_sample_pct(format_percentage(sample_reads, input_reads, 2))
      .set_reads(format_rough_number(output_reads))
      .set_bp(format_rough_number(output_bp))
      .set_length(mean_of_bp_counts(total.length_dist()))
      .set_gc(format_percentage(total.nucleotides_gc_pos().sum(), output_bp))
      .write(output);

    if (sample.size() > 1) {
      const auto total = barcode_counts.sum();

      for (size_t j = 0; j < sample.size(); j++) {
        const auto& it = sample.at(j);
        const auto count = barcode_counts.get(j);

        html_demultiplexing_barcode_row row;
        row.set_barcode_1(std::string{ it.barcode_1 })
          .set_barcode_2(std::string{ it.barcode_2 })
          .set_barcode_pct_row(format_percentage(count, total, 2));

        if (mixed_orientation) {
          row.set_orientation("<td>" + orientation_to_label(it) + "</td>");
        }

        row.write(output);
      }
    }

    ++sample_idx;
  }

  html_demultiplexing_table_tail().write(output);

  if (multiple_barcodes || mixed_orientation) {
    html_demultiplexing_toggle().write(output);
  }
}

void
write_html_demultiplexing_section(const userconfig& config,
                                  const statistics& stats,
                                  std::ostream& output)

{
  bool multiple_barcodes = false;
  bool mixed_orientation = false;
  for (const auto& sample : config.samples) {
    multiple_barcodes |= sample.size() > 1;
    for (const auto& it : sample) {
      mixed_orientation |= it.orientation != barcode_orientation::unspecified;
    }
  }

  write_html_section_title("Demultiplexing", output);
  html_demultiplexing_head().write(output);
  write_html_demultiplexing_barplot(config, stats, output);
  write_html_demultiplexing_table(config,
                                  stats,
                                  output,
                                  multiple_barcodes,
                                  mixed_orientation);
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
      names.emplace_back("Output 2");

      if (config.is_any_filtering_enabled()) {
        stats_vec.push_back(singleton);
        names.emplace_back("Singleton");
      }
    }

    if (config.is_any_filtering_enabled()) {
      stats_vec.push_back(discarded);
      names.emplace_back("Discarded");
    }
  }

  write_html_io_section(config,
                        output,
                        "Output",
                        std::move(stats_vec),
                        std::move(names),
                        merged);
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

  if (config.paired_ended_mode || config.report_duplication ||
      config.run_type == ar_command::report_only) {
    write_html_analyses_section(config, stats, output);
  }

  if (config.is_demultiplexing_enabled()) {
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
  } catch (const io_error& error) {
    log::error() << "Error writing JSON report to '" << filename << "':\n"
                 << indent_lines(error.what());
    return false;
  }

  return true;
}

} // namespace adapterremoval
