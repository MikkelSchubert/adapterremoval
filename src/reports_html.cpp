/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <iomanip>  // for setprecision
#include <iostream> // for ofstream, operator<<, basic_ostream, endl
#include <memory>   // for make_shared, make_unique
#include <sstream>  // for stringstream
#include <string>   // for operator+, string, operator<<
#include <vector>   // for vector

#include "adapterset.hpp"            // for adapter_set
#include "counts.hpp"                // for counts, counts_tmpl
#include "debug.hpp"                 // for AR_FAIL
#include "fastq.hpp"                 // for fastq_pair_vec, IDX_TO_ACGT, fastq
#include "json.hpp"                  // for json_writer, json_section
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
  std::stringstream ss;

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

  std::stringstream ss;
  ss << std::fixed << std::setprecision(1)
     << (bases / static_cast<double>(reads));

  return ss.str();
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

    const auto total_quality = stats.qualities_pos();
    auto total_bases = stats.nucleotides_pos();
    // Ensure that values get written, to prevent the plot being omitted
    total_bases.resize_up_to(1);

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
    dict->f64_vec("y", total_quality / total_bases);
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
    // Ensure that values get written, to prevent the plot being omitted
    total_bases.resize_up_to(1);

    for (const auto nuc : ACGTN::values) {
      const auto bases = stats.nucleotides_pos(nuc);

      const auto dict = content.dict();
      dict->str("read", names.at(i));
      dict->i64("offset", 1);
      dict->str("group", std::string(1, nuc));
      dict->f64_vec("y", bases / total_bases);
    }

    {
      auto dict = content.dict();
      dict->str("read", names.at(i));
      dict->i64("offset", 1);
      dict->str("group", "GC");
      dict->f64_vec("y", stats.nucleotides_gc_pos() / total_bases);
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

    std::stringstream ss;
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

  html_summary_processing section;

  section.set_adapter_trimming_reads(format_rough_number(adapter_reads));
  section.set_adapter_trimming_bases(format_rough_number(adapter_bases));
  section.set_adapter_trimming_avg(
    format_fraction(adapter_bases, adapter_reads));

  // FIXME:
  // section.set_overlapping_reads(format_rough_number(totals.overlapping_reads));
  // section.set_terminal_bases_trimmed(
  //   format_rough_number(totals.terminal_bases_trimmed));
  section.set_low_quality_reads(
    format_rough_number(totals.low_quality_trimmed_reads));
  section.set_low_quality_bases(
    format_rough_number(totals.low_quality_trimmed_bases));
  section.set_low_quality_avg(format_fraction(
    totals.low_quality_trimmed_bases, totals.low_quality_trimmed_reads));

  section.set_min_length_reads(
    format_rough_number(totals.filtered_min_length_reads));
  section.set_min_length_bases(
    format_rough_number(totals.filtered_min_length_bases));
  section.set_min_length_avg(format_fraction(totals.filtered_min_length_bases,
                                             totals.filtered_min_length_reads));

  section.set_max_length_reads(
    format_rough_number(totals.filtered_max_length_reads));
  section.set_max_length_bases(
    format_rough_number(totals.filtered_max_length_bases));
  section.set_max_length_avg(format_fraction(totals.filtered_max_length_bases,
                                             totals.filtered_max_length_reads));

  section.set_ambiguous_reads(
    format_rough_number(totals.filtered_ambiguous_reads));
  section.set_ambiguous_bases(
    format_rough_number(totals.filtered_ambiguous_bases));
  section.set_ambiguous_avg(format_fraction(totals.filtered_ambiguous_bases,
                                            totals.filtered_ambiguous_reads));

  section.set_complexity_reads(
    format_rough_number(totals.filtered_low_complexity_reads));
  section.set_complexity_bases(
    format_rough_number(totals.filtered_low_complexity_bases));
  section.set_complexity_avg(
    format_fraction(totals.filtered_low_complexity_bases,
                    totals.filtered_low_complexity_reads));

  section.write(output);
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

    // FIXME: Specify "identified reads" when in demultiplexing mode and correct
    //        format_percentage to merged / n_identified.
    std::stringstream ss;
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
    std::cerr << "Error writing JSON report to '" << filename << "':\n"
              << cli_formatter::fmt(error.what()) << std::endl;
    return false;
  }

  return true;
}

} // namespace adapterremoval
