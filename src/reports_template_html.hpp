/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Stinus Lindgreen - stinus@binf.ku.dk            *
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
#pragma once

#include <fstream>
#include <string>
#include <vector>

namespace adapterremoval {

class html_template
{
public:
  html_template();
  virtual ~html_template();
  virtual void write(std::ofstream& out) = 0;
};

class html_head : public html_template
{
public:
  html_head();
  virtual ~html_head() override;

  html_head(const html_head&) = delete;
  html_head& operator=(const html_head&) = delete;

  html_head& set_name(const std::string& value);
  html_head& set_version(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_name;
  bool m_name_is_set;
  std::string m_version;
  bool m_version_is_set;
};

class html_body_start : public html_template
{
public:
  html_body_start();
  virtual ~html_body_start() override;

  html_body_start(const html_body_start&) = delete;
  html_body_start& operator=(const html_body_start&) = delete;

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
};

class html_summary : public html_template
{
public:
  html_summary();
  virtual ~html_summary() override;

  html_summary(const html_summary&) = delete;
  html_summary& operator=(const html_summary&) = delete;

  html_summary& set_command(const std::string& value);
  html_summary& set_date_and_time(const std::string& value);
  html_summary& set_runtime(const std::string& value);
  html_summary& set_version(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_command;
  bool m_command_is_set;
  std::string m_date_and_time;
  bool m_date_and_time_is_set;
  std::string m_runtime;
  bool m_runtime_is_set;
  std::string m_version;
  bool m_version_is_set;
};

class html_summary_io : public html_template
{
public:
  html_summary_io();
  virtual ~html_summary_io() override;

  html_summary_io(const html_summary_io&) = delete;
  html_summary_io& operator=(const html_summary_io&) = delete;

  html_summary_io& add_columns(const std::string& value);
  html_summary_io& add_gc(const std::string& value);
  html_summary_io& add_lengths(const std::string& value);
  html_summary_io& add_n_bases(const std::string& value);
  html_summary_io& add_n_reads(const std::string& value);
  html_summary_io& add_ns(const std::string& value);
  html_summary_io& add_q20(const std::string& value);
  html_summary_io& add_q30(const std::string& value);
  html_summary_io& set_title(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::vector<std::string> m_columns;
  bool m_columns_is_set;
  std::vector<std::string> m_gc;
  bool m_gc_is_set;
  std::vector<std::string> m_lengths;
  bool m_lengths_is_set;
  std::vector<std::string> m_n_bases;
  bool m_n_bases_is_set;
  std::vector<std::string> m_n_reads;
  bool m_n_reads_is_set;
  std::vector<std::string> m_ns;
  bool m_ns_is_set;
  std::vector<std::string> m_q20;
  bool m_q20_is_set;
  std::vector<std::string> m_q30;
  bool m_q30_is_set;
  std::string m_title;
  bool m_title_is_set;
};

class html_output_note_pe : public html_template
{
public:
  html_output_note_pe();
  virtual ~html_output_note_pe() override;

  html_output_note_pe(const html_output_note_pe&) = delete;
  html_output_note_pe& operator=(const html_output_note_pe&) = delete;

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
};

class html_output_note_se : public html_template
{
public:
  html_output_note_se();
  virtual ~html_output_note_se() override;

  html_output_note_se(const html_output_note_se&) = delete;
  html_output_note_se& operator=(const html_output_note_se&) = delete;

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
};

class html_summary_processing : public html_template
{
public:
  html_summary_processing();
  virtual ~html_summary_processing() override;

  html_summary_processing(const html_summary_processing&) = delete;
  html_summary_processing& operator=(const html_summary_processing&) = delete;

  html_summary_processing& set_adapter_trimming_avg(const std::string& value);
  html_summary_processing& set_adapter_trimming_bases(const std::string& value);
  html_summary_processing& set_adapter_trimming_reads(const std::string& value);
  html_summary_processing& set_ambiguous_avg(const std::string& value);
  html_summary_processing& set_ambiguous_bases(const std::string& value);
  html_summary_processing& set_ambiguous_reads(const std::string& value);
  html_summary_processing& set_complexity_avg(const std::string& value);
  html_summary_processing& set_complexity_bases(const std::string& value);
  html_summary_processing& set_complexity_reads(const std::string& value);
  html_summary_processing& set_low_quality_avg(const std::string& value);
  html_summary_processing& set_low_quality_bases(const std::string& value);
  html_summary_processing& set_low_quality_reads(const std::string& value);
  html_summary_processing& set_max_length_avg(const std::string& value);
  html_summary_processing& set_max_length_bases(const std::string& value);
  html_summary_processing& set_max_length_reads(const std::string& value);
  html_summary_processing& set_min_length_avg(const std::string& value);
  html_summary_processing& set_min_length_bases(const std::string& value);
  html_summary_processing& set_min_length_reads(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_adapter_trimming_avg;
  bool m_adapter_trimming_avg_is_set;
  std::string m_adapter_trimming_bases;
  bool m_adapter_trimming_bases_is_set;
  std::string m_adapter_trimming_reads;
  bool m_adapter_trimming_reads_is_set;
  std::string m_ambiguous_avg;
  bool m_ambiguous_avg_is_set;
  std::string m_ambiguous_bases;
  bool m_ambiguous_bases_is_set;
  std::string m_ambiguous_reads;
  bool m_ambiguous_reads_is_set;
  std::string m_complexity_avg;
  bool m_complexity_avg_is_set;
  std::string m_complexity_bases;
  bool m_complexity_bases_is_set;
  std::string m_complexity_reads;
  bool m_complexity_reads_is_set;
  std::string m_low_quality_avg;
  bool m_low_quality_avg_is_set;
  std::string m_low_quality_bases;
  bool m_low_quality_bases_is_set;
  std::string m_low_quality_reads;
  bool m_low_quality_reads_is_set;
  std::string m_max_length_avg;
  bool m_max_length_avg_is_set;
  std::string m_max_length_bases;
  bool m_max_length_bases_is_set;
  std::string m_max_length_reads;
  bool m_max_length_reads_is_set;
  std::string m_min_length_avg;
  bool m_min_length_avg_is_set;
  std::string m_min_length_bases;
  bool m_min_length_bases_is_set;
  std::string m_min_length_reads;
  bool m_min_length_reads_is_set;
};

class html_h2_tag : public html_template
{
public:
  html_h2_tag();
  virtual ~html_h2_tag() override;

  html_h2_tag(const html_h2_tag&) = delete;
  html_h2_tag& operator=(const html_h2_tag&) = delete;

  html_h2_tag& set_title(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_title;
  bool m_title_is_set;
};

class html_line_plot : public html_template
{
public:
  html_line_plot();
  virtual ~html_line_plot() override;

  html_line_plot(const html_line_plot&) = delete;
  html_line_plot& operator=(const html_line_plot&) = delete;

  html_line_plot& set_legend(const std::string& value);
  html_line_plot& set_sub_title(const std::string& value);
  html_line_plot& set_title(const std::string& value);
  html_line_plot& set_title_anchor(const std::string& value);
  html_line_plot& set_values(const std::string& value);
  html_line_plot& set_width(const std::string& value);
  html_line_plot& set_x_axis(const std::string& value);
  html_line_plot& set_y_axis(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_legend;
  std::string m_sub_title;
  std::string m_title;
  bool m_title_is_set;
  std::string m_title_anchor;
  std::string m_values;
  bool m_values_is_set;
  std::string m_width;
  bool m_width_is_set;
  std::string m_x_axis;
  std::string m_y_axis;
};

class html_facet_line_plot : public html_template
{
public:
  html_facet_line_plot();
  virtual ~html_facet_line_plot() override;

  html_facet_line_plot(const html_facet_line_plot&) = delete;
  html_facet_line_plot& operator=(const html_facet_line_plot&) = delete;

  html_facet_line_plot& set_legend(const std::string& value);
  html_facet_line_plot& set_title(const std::string& value);
  html_facet_line_plot& set_values(const std::string& value);
  html_facet_line_plot& set_width(const std::string& value);
  html_facet_line_plot& set_x_axis(const std::string& value);
  html_facet_line_plot& set_y_axis(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_legend;
  std::string m_title;
  bool m_title_is_set;
  std::string m_values;
  bool m_values_is_set;
  std::string m_width;
  bool m_width_is_set;
  std::string m_x_axis;
  std::string m_y_axis;
};

class html_bar_plot : public html_template
{
public:
  html_bar_plot();
  virtual ~html_bar_plot() override;

  html_bar_plot(const html_bar_plot&) = delete;
  html_bar_plot& operator=(const html_bar_plot&) = delete;

  html_bar_plot& set_title(const std::string& value);
  html_bar_plot& set_values(const std::string& value);
  html_bar_plot& set_width(const std::string& value);
  html_bar_plot& set_x_axis(const std::string& value);
  html_bar_plot& set_y_axis(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_title;
  bool m_title_is_set;
  std::string m_values;
  bool m_values_is_set;
  std::string m_width;
  bool m_width_is_set;
  std::string m_x_axis;
  std::string m_y_axis;
};

class html_demultiplexing_head : public html_template
{
public:
  html_demultiplexing_head();
  virtual ~html_demultiplexing_head() override;

  html_demultiplexing_head(const html_demultiplexing_head&) = delete;
  html_demultiplexing_head& operator=(const html_demultiplexing_head&) = delete;

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
};

class html_demultiplexing_row : public html_template
{
public:
  html_demultiplexing_row();
  virtual ~html_demultiplexing_row() override;

  html_demultiplexing_row(const html_demultiplexing_row&) = delete;
  html_demultiplexing_row& operator=(const html_demultiplexing_row&) = delete;

  html_demultiplexing_row& set_barcode_1(const std::string& value);
  html_demultiplexing_row& set_barcode_2(const std::string& value);
  html_demultiplexing_row& set_bp(const std::string& value);
  html_demultiplexing_row& set_gc(const std::string& value);
  html_demultiplexing_row& set_length(const std::string& value);
  html_demultiplexing_row& set_n(const std::string& value);
  html_demultiplexing_row& set_name(const std::string& value);
  html_demultiplexing_row& set_pct(const std::string& value);
  html_demultiplexing_row& set_reads(const std::string& value);

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
  std::string m_barcode_1;
  bool m_barcode_1_is_set;
  std::string m_barcode_2;
  bool m_barcode_2_is_set;
  std::string m_bp;
  bool m_bp_is_set;
  std::string m_gc;
  bool m_gc_is_set;
  std::string m_length;
  bool m_length_is_set;
  std::string m_n;
  bool m_n_is_set;
  std::string m_name;
  bool m_name_is_set;
  std::string m_pct;
  bool m_pct_is_set;
  std::string m_reads;
  bool m_reads_is_set;
};

class html_demultiplexing_tail : public html_template
{
public:
  html_demultiplexing_tail();
  virtual ~html_demultiplexing_tail() override;

  html_demultiplexing_tail(const html_demultiplexing_tail&) = delete;
  html_demultiplexing_tail& operator=(const html_demultiplexing_tail&) = delete;

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
};

class html_body_end : public html_template
{
public:
  html_body_end();
  virtual ~html_body_end() override;

  html_body_end(const html_body_end&) = delete;
  html_body_end& operator=(const html_body_end&) = delete;

  virtual void write(std::ofstream& out) override;

private:
  bool m_written;
};

} // namespace adapterremoval
