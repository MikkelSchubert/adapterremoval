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
#pragma once

#include <ostream>
#include <string>
#include <vector>

namespace adapterremoval {

using string_vec = std::vector<std::string>;


class html_template
{
public:
  html_template() = default;
  virtual ~html_template() = default;
  virtual void write(std::ostream& out) = 0;

  html_template(const html_template&) = delete;
  html_template(html_template&&) = delete;
  html_template& operator=(const html_template&) = delete;
  html_template& operator=(html_template&&) = delete;
};

class html_head : public html_template
{
public:
  html_head() = default;
  ~html_head() override;

  html_head(const html_head&) = delete;
  html_head(html_head&&) = delete;
  html_head& operator=(const html_head&) = delete;
  html_head& operator=(html_head&&) = delete;

  html_head& set_name(const std::string& value);
  html_head& set_version(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_name_is_set{};
  bool m_version_is_set{};
  std::string m_name{};
  std::string m_version{};
};

class html_body_start : public html_template
{
public:
  html_body_start() = default;
  ~html_body_start() override;

  html_body_start(const html_body_start&) = delete;
  html_body_start(html_body_start&&) = delete;
  html_body_start& operator=(const html_body_start&) = delete;
  html_body_start& operator=(html_body_start&&) = delete;

  void write(std::ostream& out) override;

private:
  bool m_written{};
};

class html_summary : public html_template
{
public:
  html_summary() = default;
  ~html_summary() override;

  html_summary(const html_summary&) = delete;
  html_summary(html_summary&&) = delete;
  html_summary& operator=(const html_summary&) = delete;
  html_summary& operator=(html_summary&&) = delete;

  html_summary& set_command(const std::string& value);
  html_summary& set_date_and_time(const std::string& value);
  html_summary& set_runtime(const std::string& value);
  html_summary& set_version(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_command_is_set{};
  bool m_date_and_time_is_set{};
  bool m_runtime_is_set{};
  bool m_version_is_set{};
  std::string m_command{};
  std::string m_date_and_time{};
  std::string m_runtime{};
  std::string m_version{};
};

class html_summary_io : public html_template
{
public:
  html_summary_io() = default;
  ~html_summary_io() override;

  html_summary_io(const html_summary_io&) = delete;
  html_summary_io(html_summary_io&&) = delete;
  html_summary_io& operator=(const html_summary_io&) = delete;
  html_summary_io& operator=(html_summary_io&&) = delete;

  html_summary_io& add_columns(const std::string& value);
  html_summary_io& add_gc(const std::string& value);
  html_summary_io& add_lengths(const std::string& value);
  html_summary_io& add_n_bases(const std::string& value);
  html_summary_io& add_n_reads(const std::string& value);
  html_summary_io& add_ns(const std::string& value);
  html_summary_io& add_q20(const std::string& value);
  html_summary_io& add_q30(const std::string& value);
  html_summary_io& set_title(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_columns_is_set{};
  bool m_gc_is_set{};
  bool m_lengths_is_set{};
  bool m_n_bases_is_set{};
  bool m_n_reads_is_set{};
  bool m_ns_is_set{};
  bool m_q20_is_set{};
  bool m_q30_is_set{};
  bool m_title_is_set{};
  string_vec m_columns{};
  string_vec m_gc{};
  string_vec m_lengths{};
  string_vec m_n_bases{};
  string_vec m_n_reads{};
  string_vec m_ns{};
  string_vec m_q20{};
  string_vec m_q30{};
  std::string m_title{};
};

class html_sampling_note : public html_template
{
public:
  html_sampling_note() = default;
  ~html_sampling_note() override;

  html_sampling_note(const html_sampling_note&) = delete;
  html_sampling_note(html_sampling_note&&) = delete;
  html_sampling_note& operator=(const html_sampling_note&) = delete;
  html_sampling_note& operator=(html_sampling_note&&) = delete;

  html_sampling_note& set_label(const std::string& value);
  html_sampling_note& set_pct(const std::string& value);
  html_sampling_note& set_reads(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_label_is_set{};
  bool m_pct_is_set{};
  bool m_reads_is_set{};
  std::string m_label{};
  std::string m_pct{};
  std::string m_reads{};
};

class html_output_note : public html_template
{
public:
  html_output_note() = default;
  ~html_output_note() override;

  html_output_note(const html_output_note&) = delete;
  html_output_note(html_output_note&&) = delete;
  html_output_note& operator=(const html_output_note&) = delete;
  html_output_note& operator=(html_output_note&&) = delete;

  html_output_note& set_text(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_text_is_set{};
  std::string m_text{};
};

class html_output_footnote : public html_template
{
public:
  html_output_footnote() = default;
  ~html_output_footnote() override;

  html_output_footnote(const html_output_footnote&) = delete;
  html_output_footnote(html_output_footnote&&) = delete;
  html_output_footnote& operator=(const html_output_footnote&) = delete;
  html_output_footnote& operator=(html_output_footnote&&) = delete;

  html_output_footnote& set_symbol(const std::string& value);
  html_output_footnote& set_text(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_symbol_is_set{};
  bool m_text_is_set{};
  std::string m_symbol{};
  std::string m_text{};
};

class html_summary_trimming_head : public html_template
{
public:
  html_summary_trimming_head() = default;
  ~html_summary_trimming_head() override;

  html_summary_trimming_head(const html_summary_trimming_head&) = delete;
  html_summary_trimming_head(html_summary_trimming_head&&) = delete;
  html_summary_trimming_head& operator=(const html_summary_trimming_head&) = delete;
  html_summary_trimming_head& operator=(html_summary_trimming_head&&) = delete;

  void write(std::ostream& out) override;

private:
  bool m_written{};
};

class html_summary_trimming_row : public html_template
{
public:
  html_summary_trimming_row() = default;
  ~html_summary_trimming_row() override;

  html_summary_trimming_row(const html_summary_trimming_row&) = delete;
  html_summary_trimming_row(html_summary_trimming_row&&) = delete;
  html_summary_trimming_row& operator=(const html_summary_trimming_row&) = delete;
  html_summary_trimming_row& operator=(html_summary_trimming_row&&) = delete;

  html_summary_trimming_row& set_avg_bases(const std::string& value);
  html_summary_trimming_row& set_bases(const std::string& value);
  html_summary_trimming_row& set_label_1(const std::string& value);
  html_summary_trimming_row& set_label_2(const std::string& value);
  html_summary_trimming_row& set_pct_bases(const std::string& value);
  html_summary_trimming_row& set_pct_reads(const std::string& value);
  html_summary_trimming_row& set_reads(const std::string& value);
  html_summary_trimming_row& set_stage(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_avg_bases_is_set{};
  bool m_bases_is_set{};
  bool m_label_1_is_set{};
  bool m_label_2_is_set{};
  bool m_pct_bases_is_set{};
  bool m_pct_reads_is_set{};
  bool m_reads_is_set{};
  bool m_stage_is_set{};
  std::string m_avg_bases{};
  std::string m_bases{};
  std::string m_label_1{};
  std::string m_label_2{};
  std::string m_pct_bases{};
  std::string m_pct_reads{};
  std::string m_reads{};
  std::string m_stage{};
};

class html_summary_trimming_tail : public html_template
{
public:
  html_summary_trimming_tail() = default;
  ~html_summary_trimming_tail() override;

  html_summary_trimming_tail(const html_summary_trimming_tail&) = delete;
  html_summary_trimming_tail(html_summary_trimming_tail&&) = delete;
  html_summary_trimming_tail& operator=(const html_summary_trimming_tail&) = delete;
  html_summary_trimming_tail& operator=(html_summary_trimming_tail&&) = delete;

  html_summary_trimming_tail& set_n_enabled_filt(const std::string& value);
  html_summary_trimming_tail& set_n_enabled_proc(const std::string& value);
  html_summary_trimming_tail& set_n_total_filt(const std::string& value);
  html_summary_trimming_tail& set_n_total_proc(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_n_enabled_filt_is_set{};
  bool m_n_enabled_proc_is_set{};
  bool m_n_total_filt_is_set{};
  bool m_n_total_proc_is_set{};
  std::string m_n_enabled_filt{};
  std::string m_n_enabled_proc{};
  std::string m_n_total_filt{};
  std::string m_n_total_proc{};
};

class html_consensus_adapter_head : public html_template
{
public:
  html_consensus_adapter_head() = default;
  ~html_consensus_adapter_head() override;

  html_consensus_adapter_head(const html_consensus_adapter_head&) = delete;
  html_consensus_adapter_head(html_consensus_adapter_head&&) = delete;
  html_consensus_adapter_head& operator=(const html_consensus_adapter_head&) = delete;
  html_consensus_adapter_head& operator=(html_consensus_adapter_head&&) = delete;

  html_consensus_adapter_head& set_overlapping_pairs(const std::string& value);
  html_consensus_adapter_head& set_pairs_with_adapters(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_overlapping_pairs_is_set{};
  bool m_pairs_with_adapters_is_set{};
  std::string m_overlapping_pairs{};
  std::string m_pairs_with_adapters{};
};

class html_consensus_adapter_table : public html_template
{
public:
  html_consensus_adapter_table() = default;
  ~html_consensus_adapter_table() override;

  html_consensus_adapter_table(const html_consensus_adapter_table&) = delete;
  html_consensus_adapter_table(html_consensus_adapter_table&&) = delete;
  html_consensus_adapter_table& operator=(const html_consensus_adapter_table&) = delete;
  html_consensus_adapter_table& operator=(html_consensus_adapter_table&&) = delete;

  html_consensus_adapter_table& set_alignment_1(const std::string& value);
  html_consensus_adapter_table& set_alignment_2(const std::string& value);
  html_consensus_adapter_table& set_consensus_1(const std::string& value);
  html_consensus_adapter_table& set_consensus_2(const std::string& value);
  html_consensus_adapter_table& set_name_1(const std::string& value);
  html_consensus_adapter_table& set_name_2(const std::string& value);
  html_consensus_adapter_table& set_qualities_1(const std::string& value);
  html_consensus_adapter_table& set_qualities_2(const std::string& value);
  html_consensus_adapter_table& set_reference_1(const std::string& value);
  html_consensus_adapter_table& set_reference_2(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_alignment_1_is_set{};
  bool m_alignment_2_is_set{};
  bool m_consensus_1_is_set{};
  bool m_consensus_2_is_set{};
  bool m_name_1_is_set{};
  bool m_name_2_is_set{};
  bool m_qualities_1_is_set{};
  bool m_qualities_2_is_set{};
  bool m_reference_1_is_set{};
  bool m_reference_2_is_set{};
  std::string m_alignment_1{};
  std::string m_alignment_2{};
  std::string m_consensus_1{};
  std::string m_consensus_2{};
  std::string m_name_1{};
  std::string m_name_2{};
  std::string m_qualities_1{};
  std::string m_qualities_2{};
  std::string m_reference_1{};
  std::string m_reference_2{};
};

class html_consensus_adapter_kmer_head : public html_template
{
public:
  html_consensus_adapter_kmer_head() = default;
  ~html_consensus_adapter_kmer_head() override;

  html_consensus_adapter_kmer_head(const html_consensus_adapter_kmer_head&) = delete;
  html_consensus_adapter_kmer_head(html_consensus_adapter_kmer_head&&) = delete;
  html_consensus_adapter_kmer_head& operator=(const html_consensus_adapter_kmer_head&) = delete;
  html_consensus_adapter_kmer_head& operator=(html_consensus_adapter_kmer_head&&) = delete;

  html_consensus_adapter_kmer_head& set_kmer_length(const std::string& value);
  html_consensus_adapter_kmer_head& set_n_kmers(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_kmer_length_is_set{};
  bool m_n_kmers_is_set{};
  std::string m_kmer_length{};
  std::string m_n_kmers{};
};

class html_consensus_adapter_kmer_row : public html_template
{
public:
  html_consensus_adapter_kmer_row() = default;
  ~html_consensus_adapter_kmer_row() override;

  html_consensus_adapter_kmer_row(const html_consensus_adapter_kmer_row&) = delete;
  html_consensus_adapter_kmer_row(html_consensus_adapter_kmer_row&&) = delete;
  html_consensus_adapter_kmer_row& operator=(const html_consensus_adapter_kmer_row&) = delete;
  html_consensus_adapter_kmer_row& operator=(html_consensus_adapter_kmer_row&&) = delete;

  html_consensus_adapter_kmer_row& set_count_1(const std::string& value);
  html_consensus_adapter_kmer_row& set_count_2(const std::string& value);
  html_consensus_adapter_kmer_row& set_index(const std::string& value);
  html_consensus_adapter_kmer_row& set_kmer_1(const std::string& value);
  html_consensus_adapter_kmer_row& set_kmer_2(const std::string& value);
  html_consensus_adapter_kmer_row& set_pct_1(const std::string& value);
  html_consensus_adapter_kmer_row& set_pct_2(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_count_1_is_set{};
  bool m_count_2_is_set{};
  bool m_index_is_set{};
  bool m_kmer_1_is_set{};
  bool m_kmer_2_is_set{};
  bool m_pct_1_is_set{};
  bool m_pct_2_is_set{};
  std::string m_count_1{};
  std::string m_count_2{};
  std::string m_index{};
  std::string m_kmer_1{};
  std::string m_kmer_2{};
  std::string m_pct_1{};
  std::string m_pct_2{};
};

class html_consensus_adapter_kmer_tail : public html_template
{
public:
  html_consensus_adapter_kmer_tail() = default;
  ~html_consensus_adapter_kmer_tail() override;

  html_consensus_adapter_kmer_tail(const html_consensus_adapter_kmer_tail&) = delete;
  html_consensus_adapter_kmer_tail(html_consensus_adapter_kmer_tail&&) = delete;
  html_consensus_adapter_kmer_tail& operator=(const html_consensus_adapter_kmer_tail&) = delete;
  html_consensus_adapter_kmer_tail& operator=(html_consensus_adapter_kmer_tail&&) = delete;

  void write(std::ostream& out) override;

private:
  bool m_written{};
};

class html_h2_tag : public html_template
{
public:
  html_h2_tag() = default;
  ~html_h2_tag() override;

  html_h2_tag(const html_h2_tag&) = delete;
  html_h2_tag(html_h2_tag&&) = delete;
  html_h2_tag& operator=(const html_h2_tag&) = delete;
  html_h2_tag& operator=(html_h2_tag&&) = delete;

  html_h2_tag& set_title(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_title_is_set{};
  std::string m_title{};
};

class html_line_plot : public html_template
{
public:
  html_line_plot() = default;
  ~html_line_plot() override;

  html_line_plot(const html_line_plot&) = delete;
  html_line_plot(html_line_plot&&) = delete;
  html_line_plot& operator=(const html_line_plot&) = delete;
  html_line_plot& operator=(html_line_plot&&) = delete;

  html_line_plot& set_legend(const std::string& value);
  html_line_plot& set_sub_title(const std::string& value);
  html_line_plot& set_title(const std::string& value);
  html_line_plot& set_values(const std::string& value);
  html_line_plot& set_width(const std::string& value);
  html_line_plot& set_x_axis(const std::string& value);
  html_line_plot& set_y_axis(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_sub_title_is_set{};
  bool m_title_is_set{};
  bool m_values_is_set{};
  bool m_width_is_set{};
  std::string m_legend{"{ \"title\": \"Legend\", \"padding\": 5 }"};
  std::string m_sub_title{};
  std::string m_title{};
  std::string m_values{};
  std::string m_width{};
  std::string m_x_axis{"null"};
  std::string m_y_axis{"null"};
};

class html_facet_line_plot : public html_template
{
public:
  html_facet_line_plot() = default;
  ~html_facet_line_plot() override;

  html_facet_line_plot(const html_facet_line_plot&) = delete;
  html_facet_line_plot(html_facet_line_plot&&) = delete;
  html_facet_line_plot& operator=(const html_facet_line_plot&) = delete;
  html_facet_line_plot& operator=(html_facet_line_plot&&) = delete;

  html_facet_line_plot& set_legend(const std::string& value);
  html_facet_line_plot& set_title(const std::string& value);
  html_facet_line_plot& set_values(const std::string& value);
  html_facet_line_plot& set_width(const std::string& value);
  html_facet_line_plot& set_x_axis(const std::string& value);
  html_facet_line_plot& set_y_axis(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_title_is_set{};
  bool m_values_is_set{};
  bool m_width_is_set{};
  std::string m_legend{"{ \"title\": \"Legend\", \"padding\": 5 }"};
  std::string m_title{};
  std::string m_values{};
  std::string m_width{};
  std::string m_x_axis{"null"};
  std::string m_y_axis{"null"};
};

class html_bar_plot : public html_template
{
public:
  html_bar_plot() = default;
  ~html_bar_plot() override;

  html_bar_plot(const html_bar_plot&) = delete;
  html_bar_plot(html_bar_plot&&) = delete;
  html_bar_plot& operator=(const html_bar_plot&) = delete;
  html_bar_plot& operator=(html_bar_plot&&) = delete;

  html_bar_plot& set_title(const std::string& value);
  html_bar_plot& set_values(const std::string& value);
  html_bar_plot& set_width(const std::string& value);
  html_bar_plot& set_x_axis(const std::string& value);
  html_bar_plot& set_y_axis(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_title_is_set{};
  bool m_values_is_set{};
  bool m_width_is_set{};
  std::string m_title{};
  std::string m_values{};
  std::string m_width{};
  std::string m_x_axis{"null"};
  std::string m_y_axis{"null"};
};

class html_demultiplexing_head : public html_template
{
public:
  html_demultiplexing_head() = default;
  ~html_demultiplexing_head() override;

  html_demultiplexing_head(const html_demultiplexing_head&) = delete;
  html_demultiplexing_head(html_demultiplexing_head&&) = delete;
  html_demultiplexing_head& operator=(const html_demultiplexing_head&) = delete;
  html_demultiplexing_head& operator=(html_demultiplexing_head&&) = delete;

  void write(std::ostream& out) override;

private:
  bool m_written{};
};

class html_demultiplexing_row : public html_template
{
public:
  html_demultiplexing_row() = default;
  ~html_demultiplexing_row() override;

  html_demultiplexing_row(const html_demultiplexing_row&) = delete;
  html_demultiplexing_row(html_demultiplexing_row&&) = delete;
  html_demultiplexing_row& operator=(const html_demultiplexing_row&) = delete;
  html_demultiplexing_row& operator=(html_demultiplexing_row&&) = delete;

  html_demultiplexing_row& set_barcode_1(const std::string& value);
  html_demultiplexing_row& set_barcode_2(const std::string& value);
  html_demultiplexing_row& set_bp(const std::string& value);
  html_demultiplexing_row& set_gc(const std::string& value);
  html_demultiplexing_row& set_length(const std::string& value);
  html_demultiplexing_row& set_n(const std::string& value);
  html_demultiplexing_row& set_name(const std::string& value);
  html_demultiplexing_row& set_pct(const std::string& value);
  html_demultiplexing_row& set_reads(const std::string& value);

  void write(std::ostream& out) override;

private:
  bool m_written{};
  bool m_barcode_1_is_set{};
  bool m_barcode_2_is_set{};
  bool m_bp_is_set{};
  bool m_gc_is_set{};
  bool m_length_is_set{};
  bool m_n_is_set{};
  bool m_name_is_set{};
  bool m_pct_is_set{};
  bool m_reads_is_set{};
  std::string m_barcode_1{};
  std::string m_barcode_2{};
  std::string m_bp{};
  std::string m_gc{};
  std::string m_length{};
  std::string m_n{};
  std::string m_name{};
  std::string m_pct{};
  std::string m_reads{};
};

class html_demultiplexing_tail : public html_template
{
public:
  html_demultiplexing_tail() = default;
  ~html_demultiplexing_tail() override;

  html_demultiplexing_tail(const html_demultiplexing_tail&) = delete;
  html_demultiplexing_tail(html_demultiplexing_tail&&) = delete;
  html_demultiplexing_tail& operator=(const html_demultiplexing_tail&) = delete;
  html_demultiplexing_tail& operator=(html_demultiplexing_tail&&) = delete;

  void write(std::ostream& out) override;

private:
  bool m_written{};
};

class html_body_end : public html_template
{
public:
  html_body_end() = default;
  ~html_body_end() override;

  html_body_end(const html_body_end&) = delete;
  html_body_end(html_body_end&&) = delete;
  html_body_end& operator=(const html_body_end&) = delete;
  html_body_end& operator=(html_body_end&&) = delete;

  void write(std::ostream& out) override;

private:
  bool m_written{};
};

} // namespace adapterremoval
