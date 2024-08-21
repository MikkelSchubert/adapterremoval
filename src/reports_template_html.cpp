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
#include "reports_template_html.hpp"
#include "debug.hpp" // for AR_REQUIRE
#include <cstddef>   // for size_t

namespace adapterremoval {

size_t g_html_id = 1;

html_head::~html_head()
{
  AR_REQUIRE(m_written);
}

html_head&
html_head::set_name(const std::string& value)
{
  m_name = value;
  m_name_is_set = true;
  return *this;
}

html_head&
html_head::set_version(const std::string& value)
{
  m_version = value;
  m_version_is_set = true;
  return *this;
}

void
html_head::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_name_is_set);
  AR_REQUIRE(m_version_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "<!DOCTYPE html>\n";
  out << "<html lang='en'>\n";
  out << "\n";
  out << "<head>\n";
  out << "    <meta charset='utf-8'>\n";
  out << "    <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n";
  out << "    <title>" << m_name << " " << m_version << "</title>\n";
  out << "    <link rel='stylesheet' href='https://cdn.jsdelivr.net/npm/purecss@2.1.0/build/pure-min.css'\n";
  out << "        integrity='sha384-yHIFVG6ClnONEA5yB5DJXfW2/KC173DIQrYoZMEtBvGzmf0PKiGyNEqe9N6BNDBH' crossorigin='anonymous'>\n";
  out << "    <script src=\"https://cdn.jsdelivr.net/npm/vega@5.21.0/build/vega.min.js\"\n";
  out << "        integrity=\"sha384-s2nYi9D0FfKNopEKsfINeS1Ffhcf+5uvwIrb7Zqso2II+HPhzBTWvXClt+NdUwFc\"\n";
  out << "        crossorigin=\"anonymous\"></script>\n";
  out << "    <script src=\"https://cdn.jsdelivr.net/npm/vega-lite@5.2.0/build/vega-lite.min.js\"\n";
  out << "        integrity=\"sha384-tU6fj0fI2gxrcWwC7uBMp70QvipC9ukjcXyOs85VMmdCq33CrA7xQ3nJkJu0SmDm\"\n";
  out << "        crossorigin=\"anonymous\"></script>\n";
  out << "    <script src=\"https://cdn.jsdelivr.net/npm/vega-embed@6.20.2/build/vega-embed.min.js\"\n";
  out << "        integrity=\"sha384-oP1rwLY7weRZ5jvAVzfnJsAn+sYA69rQC4geH82Y9oMvr8ruA1oeE9Jkft2noCHR\"\n";
  out << "        crossorigin=\"anonymous\"></script>\n";
  out << "    <style type='text/css'>\n";
  out << "        body {\n";
  out << "            background-color: #E3E2DE;\n";
  out << "        }\n";
  out << "\n";
  out << "        div#layout {\n";
  out << "            max-width: 920px;\n";
  out << "            margin-left: auto;\n";
  out << "            margin-right: auto;\n";
  out << "            font-size: smaller;\n";
  out << "        }\n";
  out << "\n";
  out << "        div.title {\n";
  out << "            background-color: #8C9CC0;\n";
  out << "            margin: -10px !important;\n";
  out << "            text-align: center;\n";
  out << "            border-radius: 5px;\n";
  out << "        }\n";
  out << "\n";
  out << "        div.title>h1,\n";
  out << "        div.title>h2 {\n";
  out << "            padding: 5px;\n";
  out << "        }\n";
  out << "\n";
  out << "        h5 {\n";
  out << "            margin-bottom: 2px;\n";
  out << "            margin-left: 6px;\n";
  out << "        }\n";
  out << "\n";
  out << "        .pure-table {\n";
  out << "            margin-left: 1em;\n";
  out << "        }\n";
  out << "\n";
  out << "        .pure-table thead>tr>th {\n";
  out << "            font-weight: bold;\n";
  out << "            background-color: #C4CCDB;\n";
  out << "        }\n";
  out << "\n";
  out << "        .summary-table tr>td:first-child {\n";
  out << "            font-weight: bold;\n";
  out << "        }\n";
  out << "\n";
  out << "        .io-table tr>td:first-child {\n";
  out << "            font-weight: bold;\n";
  out << "        }\n";
  out << "\n";
  out << "        .io-table tr>td:not(:first-child) {\n";
  out << "            text-align: right;\n";
  out << "            width: 100px;\n";
  out << "        }\n";
  out << "\n";
  out << "        .trimming-table tr>td:nth-child(-n+3) {\n";
  out << "            font-weight: bold;\n";
  out << "        }\n";
  out << "\n";
  out << "        .trimming-table tr>td:nth-child(n+3) {\n";
  out << "            text-align: right;\n";
  out << "        }\n";
  out << "\n";
  out << "        .adapter-table tr>td:first-child {\n";
  out << "            font-weight: bold;\n";
  out << "        }\n";
  out << "\n";
  out << "        .kmer-table thead {\n";
  out << "            font-weight: bold;\n";
  out << "        }\n";
  out << "\n";
  out << "        .kmer-table tr>td:first-child {\n";
  out << "            font-weight: bold;\n";
  out << "        }\n";
  out << "\n";
  out << "        .kmer-table {\n";
  out << "            text-align: right;\n";
  out << "        }\n";
  out << "\n";
  out << "        .trimming-table tr>td:nth-child(4),\n";
  out << "        .trimming-table tr>td:nth-child(6),\n";
  out << "        .trimming-table tr>td:nth-child(8) {\n";
  out << "            width: 80px;\n";
  out << "        }\n";
  out << "\n";
  out << "        .trimming-table tr>td:nth-child(5),\n";
  out << "        .trimming-table tr>td:nth-child(7) {\n";
  out << "            width: 40px;\n";
  out << "        }\n";
  out << "\n";
  out << "        .trimming-table tr>td:nth-of-type(3) {\n";
  out << "            text-align: center;\n";
  out << "            width: 25px;\n";
  out << "        }\n";
  out << "\n";
  out << "        .fixed-height-table {\n";
  out << "            max-height: 300px;\n";
  out << "            overflow-y: scroll;\n";
  out << "        }\n";
  out << "\n";
  out << "        .fixed-height-table>table {\n";
  out << "            width: 100%;\n";
  out << "        }\n";
  out << "\n";
  out << "        .fixed-height-table>table>thead>tr>th {\n";
  out << "            position: sticky;\n";
  out << "            top: 0;\n";
  out << "            box-shadow: 0 2px 2px -1px rgba(0, 0, 0, 0.4);\n";
  out << "        }\n";
  out << "\n";
  out << "        .section {\n";
  out << "            background-color: #FFF;\n";
  out << "            border-radius: 5px;\n";
  out << "            margin-bottom: 10px;\n";
  out << "            padding: 10px;\n";
  out << "            padding-top: 0px;\n";
  out << "        }\n";
  out << "\n";
  out << "        .note {\n";
  out << "            color: #777;\n";
  out << "            font-size: small;\n";
  out << "        }\n";
  out << "\n";
  out << "        .epilogue {\n";
  out << "            padding-top: 10px;\n";
  out << "        }\n";
  out << "\n";
  out << "        .anchor {\n";
  out << "            color: #104fc6;\n";
  out << "            font-size: smaller;\n";
  out << "            text-decoration: none;\n";
  out << "        }\n";
  out << "\n";
  out << "        .anchor:hover {\n";
  out << "            text-decoration: underline;\n";
  out << "        }\n";
  out << "    </style>\n";
  out << "</head>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_body_start::~html_body_start()
{
  AR_REQUIRE(m_written);
}

void
html_body_start::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "<body>\n";
  out << "    <div id='layout'>\n";
  out << "        <div class=\"title\">\n";
  out << "            <h1>AdapterRemoval</h1>\n";
  out << "        </div>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_summary::~html_summary()
{
  AR_REQUIRE(m_written);
}

html_summary&
html_summary::set_command(const std::string& value)
{
  m_command = value;
  m_command_is_set = true;
  return *this;
}

html_summary&
html_summary::set_date_and_time(const std::string& value)
{
  m_date_and_time = value;
  m_date_and_time_is_set = true;
  return *this;
}

html_summary&
html_summary::set_runtime(const std::string& value)
{
  m_runtime = value;
  m_runtime_is_set = true;
  return *this;
}

html_summary&
html_summary::set_version(const std::string& value)
{
  m_version = value;
  m_version_is_set = true;
  return *this;
}

void
html_summary::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_command_is_set);
  AR_REQUIRE(m_date_and_time_is_set);
  AR_REQUIRE(m_runtime_is_set);
  AR_REQUIRE(m_version_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "        <div class=\"section\">\n";
  out << "            <div class=\"title\" id=\"summary\">\n";
  out << "                <h2>Summary <a class=\"anchor\" href=\"#summary\">#</a></h2>\n";
  out << "            </div>\n";
  out << "\n";
  out << "            <h4 id=\"summary-program\">\n";
  out << "                Program <a class=\"anchor\" href=\"#summary-program\">#</a>\n";
  out << "            </h4>\n";
  out << "            <table class=\"pure-table pure-table-striped summary-table\">\n";
  out << "                <tbody>\n";
  out << "                    <tr>\n";
  out << "                        <td>Date</td>\n";
  out << "                        <td>" << m_date_and_time << "</td>\n";
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Command</td>\n";
  out << "                        <td>" << m_command << "</td>\n";
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Version</td>\n";
  out << "                        <td>" << m_version << "</td>\n";
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Runtime</td>\n";
  out << "                        <td>" << m_runtime << "</td>\n";
  out << "                    </tr>\n";
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_summary_io::~html_summary_io()
{
  AR_REQUIRE(m_written);
}

html_summary_io&
html_summary_io::add_columns(const std::string& value)
{
  m_columns.push_back(value);
  m_columns_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::add_gc(const std::string& value)
{
  m_gc.push_back(value);
  m_gc_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::set_href(const std::string& value)
{
  m_href = value;
  m_href_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::add_lengths(const std::string& value)
{
  m_lengths.push_back(value);
  m_lengths_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::add_n_bases(const std::string& value)
{
  m_n_bases.push_back(value);
  m_n_bases_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::add_n_reads(const std::string& value)
{
  m_n_reads.push_back(value);
  m_n_reads_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::add_ns(const std::string& value)
{
  m_ns.push_back(value);
  m_ns_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::add_q20(const std::string& value)
{
  m_q20.push_back(value);
  m_q20_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::add_q30(const std::string& value)
{
  m_q30.push_back(value);
  m_q30_is_set = true;
  return *this;
}

html_summary_io&
html_summary_io::set_title(const std::string& value)
{
  m_title = value;
  m_title_is_set = true;
  return *this;
}

void
html_summary_io::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_columns_is_set);
  AR_REQUIRE(m_gc_is_set);
  AR_REQUIRE(m_href_is_set);
  AR_REQUIRE(m_lengths_is_set);
  AR_REQUIRE(m_n_bases_is_set);
  AR_REQUIRE(m_n_reads_is_set);
  AR_REQUIRE(m_ns_is_set);
  AR_REQUIRE(m_q20_is_set);
  AR_REQUIRE(m_q30_is_set);
  AR_REQUIRE(m_title_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4 id=\"" << m_href << "\">\n";
  out << "                " << m_title << " <a class=\"anchor\" href=\"#" << m_href << "\">#</a>\n";
  out << "            </h4>\n";
  out << "            <table class=\"pure-table io-table pure-table-striped\">\n";
  out << "                <thead>\n";
  out << "                    <tr>\n";
  out << "                        <th></th>\n";
  for (const auto& value : m_columns) {
    out << "                        <th>" << value << "</th>\n";
  }
  out << "                    </tr>\n";
  out << "                </thead>\n";
  out << "                <tbody>\n";
  out << "                    <tr>\n";
  out << "                        <td>#Reads</td>\n";
  for (const auto& value : m_n_reads) {
    out << "                        <td>" << value << "</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>#Bases</td>\n";
  for (const auto& value : m_n_bases) {
    out << "                        <td>" << value << "</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Read length</td>\n";
  for (const auto& value : m_lengths) {
    out << "                        <td>" << value << "</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Q20</td>\n";
  for (const auto& value : m_q20) {
    out << "                        <td>" << value << "</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Q30</td>\n";
  for (const auto& value : m_q30) {
    out << "                        <td>" << value << "</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>GC</td>\n";
  for (const auto& value : m_gc) {
    out << "                        <td>" << value << "</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>N</td>\n";
  for (const auto& value : m_ns) {
    out << "                        <td>" << value << "</td>\n";
  }
  out << "                    </tr>\n";
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_sampling_note::~html_sampling_note()
{
  AR_REQUIRE(m_written);
}

html_sampling_note&
html_sampling_note::set_label(const std::string& value)
{
  m_label = value;
  m_label_is_set = true;
  return *this;
}

html_sampling_note&
html_sampling_note::set_pct(const std::string& value)
{
  m_pct = value;
  m_pct_is_set = true;
  return *this;
}

html_sampling_note&
html_sampling_note::set_reads(const std::string& value)
{
  m_reads = value;
  m_reads_is_set = true;
  return *this;
}

void
html_sampling_note::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_label_is_set);
  AR_REQUIRE(m_pct_is_set);
  AR_REQUIRE(m_reads_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                Base composition statistics/plots are based on " << m_reads << " (" << m_pct << "%) of " << m_label << " reads sampled during\n";
  out << "                execution.\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_output_note::~html_output_note()
{
  AR_REQUIRE(m_written);
}

html_output_note&
html_output_note::set_text(const std::string& value)
{
  m_text = value;
  m_text_is_set = true;
  return *this;
}

void
html_output_note::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_text_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <p class=\"note\">" << m_text << "</p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_output_footnote::~html_output_footnote()
{
  AR_REQUIRE(m_written);
}

html_output_footnote&
html_output_footnote::set_symbol(const std::string& value)
{
  m_symbol = value;
  m_symbol_is_set = true;
  return *this;
}

html_output_footnote&
html_output_footnote::set_text(const std::string& value)
{
  m_text = value;
  m_text_is_set = true;
  return *this;
}

void
html_output_footnote::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_symbol_is_set);
  AR_REQUIRE(m_text_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <p class=\"note\"><b>" << m_symbol << "</b> " << m_text << "</p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_summary_trimming_head::~html_summary_trimming_head()
{
  AR_REQUIRE(m_written);
}

void
html_summary_trimming_head::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4 id=\"summary-process\">\n";
  out << "                Process summary <a class=\"anchor\" href=\"#summary-process\">#</a>\n";
  out << "            </h4>\n";
  out << "            <table class=\"pure-table trimming-table pure-table-striped\">\n";
  out << "                <thead>\n";
  out << "                    <tr>\n";
  out << "                        <th>Stage</th>\n";
  out << "                        <th>Step</th>\n";
  out << "                        <th></th>\n";
  out << "                        <th>Bases</th>\n";
  out << "                        <th>(%)</th>\n";
  out << "                        <th>Reads</th>\n";
  out << "                        <th>(%)</th>\n";
  out << "                        <th>Mean Bases</th>\n";
  out << "                    </tr>\n";
  out << "                </thead>\n";
  out << "                <tbody>\n";
  // clang-format on
  m_written = true;
}

html_summary_trimming_row::~html_summary_trimming_row()
{
  AR_REQUIRE(m_written);
}

html_summary_trimming_row&
html_summary_trimming_row::set_avg_bases(const std::string& value)
{
  m_avg_bases = value;
  m_avg_bases_is_set = true;
  return *this;
}

html_summary_trimming_row&
html_summary_trimming_row::set_bases(const std::string& value)
{
  m_bases = value;
  m_bases_is_set = true;
  return *this;
}

html_summary_trimming_row&
html_summary_trimming_row::set_label_1(const std::string& value)
{
  m_label_1 = value;
  m_label_1_is_set = true;
  return *this;
}

html_summary_trimming_row&
html_summary_trimming_row::set_label_2(const std::string& value)
{
  m_label_2 = value;
  m_label_2_is_set = true;
  return *this;
}

html_summary_trimming_row&
html_summary_trimming_row::set_pct_bases(const std::string& value)
{
  m_pct_bases = value;
  m_pct_bases_is_set = true;
  return *this;
}

html_summary_trimming_row&
html_summary_trimming_row::set_pct_reads(const std::string& value)
{
  m_pct_reads = value;
  m_pct_reads_is_set = true;
  return *this;
}

html_summary_trimming_row&
html_summary_trimming_row::set_reads(const std::string& value)
{
  m_reads = value;
  m_reads_is_set = true;
  return *this;
}

html_summary_trimming_row&
html_summary_trimming_row::set_stage(const std::string& value)
{
  m_stage = value;
  m_stage_is_set = true;
  return *this;
}

void
html_summary_trimming_row::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_avg_bases_is_set);
  AR_REQUIRE(m_bases_is_set);
  AR_REQUIRE(m_label_1_is_set);
  AR_REQUIRE(m_label_2_is_set);
  AR_REQUIRE(m_pct_bases_is_set);
  AR_REQUIRE(m_pct_reads_is_set);
  AR_REQUIRE(m_reads_is_set);
  AR_REQUIRE(m_stage_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "                    <tr>\n";
  out << "                        <td>" << m_stage << "</td>\n";
  out << "                        <td>" << m_label_1 << "</td>\n";
  out << "                        <td>" << m_label_2 << "</td>\n";
  out << "                        <td>" << m_bases << "</td>\n";
  out << "                        <td>" << m_pct_bases << "</td>\n";
  out << "                        <td>" << m_reads << "</td>\n";
  out << "                        <td>" << m_pct_reads << "</td>\n";
  out << "                        <td>" << m_avg_bases << "</td>\n";
  out << "                    </tr>\n";
  // clang-format on
  m_written = true;
}

html_summary_trimming_tail::~html_summary_trimming_tail()
{
  AR_REQUIRE(m_written);
}

html_summary_trimming_tail&
html_summary_trimming_tail::set_n_enabled_filt(const std::string& value)
{
  m_n_enabled_filt = value;
  m_n_enabled_filt_is_set = true;
  return *this;
}

html_summary_trimming_tail&
html_summary_trimming_tail::set_n_enabled_proc(const std::string& value)
{
  m_n_enabled_proc = value;
  m_n_enabled_proc_is_set = true;
  return *this;
}

html_summary_trimming_tail&
html_summary_trimming_tail::set_n_total_filt(const std::string& value)
{
  m_n_total_filt = value;
  m_n_total_filt_is_set = true;
  return *this;
}

html_summary_trimming_tail&
html_summary_trimming_tail::set_n_total_proc(const std::string& value)
{
  m_n_total_proc = value;
  m_n_total_proc_is_set = true;
  return *this;
}

void
html_summary_trimming_tail::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_n_enabled_filt_is_set);
  AR_REQUIRE(m_n_enabled_proc_is_set);
  AR_REQUIRE(m_n_total_filt_is_set);
  AR_REQUIRE(m_n_total_proc_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                " << m_n_enabled_proc << " of " << m_n_total_proc << " processing steps enabled, " << m_n_enabled_filt << " of " << m_n_total_filt << "\n";
  out << "                filtering steps enabled. Numbers of reads are given in terms of input reads, meaning that a merged read\n";
  out << "                counts for two. For <b>Processing</b> and <b>Filtering</b>, the numbers of bases specify how many were\n";
  out << "                lost during each step.\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_consensus_adapter_head::~html_consensus_adapter_head()
{
  AR_REQUIRE(m_written);
}

html_consensus_adapter_head&
html_consensus_adapter_head::set_overlapping_pairs(const std::string& value)
{
  m_overlapping_pairs = value;
  m_overlapping_pairs_is_set = true;
  return *this;
}

html_consensus_adapter_head&
html_consensus_adapter_head::set_pairs_with_adapters(const std::string& value)
{
  m_pairs_with_adapters = value;
  m_pairs_with_adapters_is_set = true;
  return *this;
}

void
html_consensus_adapter_head::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_overlapping_pairs_is_set);
  AR_REQUIRE(m_pairs_with_adapters_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4 id=\"analyses-consensus-adapter\">\n";
  out << "                Consensus adapter sequences <a class=\"anchor\" href=\"#analyses-consensus-adapter\">#</a>\n";
  out << "            </h4>\n";
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                Found " << m_overlapping_pairs << " overlapping pairs of which " << m_pairs_with_adapters << " contained adapter\n";
  out << "                sequence(s).\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_consensus_adapter_table::~html_consensus_adapter_table()
{
  AR_REQUIRE(m_written);
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_alignment_1(const std::string& value)
{
  m_alignment_1 = value;
  m_alignment_1_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_alignment_2(const std::string& value)
{
  m_alignment_2 = value;
  m_alignment_2_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_consensus_1(const std::string& value)
{
  m_consensus_1 = value;
  m_consensus_1_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_consensus_2(const std::string& value)
{
  m_consensus_2 = value;
  m_consensus_2_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_name_1(const std::string& value)
{
  m_name_1 = value;
  m_name_1_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_name_2(const std::string& value)
{
  m_name_2 = value;
  m_name_2_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_qualities_1(const std::string& value)
{
  m_qualities_1 = value;
  m_qualities_1_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_qualities_2(const std::string& value)
{
  m_qualities_2 = value;
  m_qualities_2_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_reference_1(const std::string& value)
{
  m_reference_1 = value;
  m_reference_1_is_set = true;
  return *this;
}

html_consensus_adapter_table&
html_consensus_adapter_table::set_reference_2(const std::string& value)
{
  m_reference_2 = value;
  m_reference_2_is_set = true;
  return *this;
}

void
html_consensus_adapter_table::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_alignment_1_is_set);
  AR_REQUIRE(m_alignment_2_is_set);
  AR_REQUIRE(m_consensus_1_is_set);
  AR_REQUIRE(m_consensus_2_is_set);
  AR_REQUIRE(m_name_1_is_set);
  AR_REQUIRE(m_name_2_is_set);
  AR_REQUIRE(m_qualities_1_is_set);
  AR_REQUIRE(m_qualities_2_is_set);
  AR_REQUIRE(m_reference_1_is_set);
  AR_REQUIRE(m_reference_2_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <table class=\"pure-table adapter-table pure-table-striped\" style=\"font-family: monospace;\">\n";
  out << "                <tbody>\n";
  out << "                    <tr>\n";
  out << "                        <td>" << m_name_1 << "</td>\n";
  out << "                        <td>" << m_reference_1 << "</td>\n";
  out << "                        <td></td>\n";
  out << "                        <td>" << m_name_2 << "</td>\n";
  out << "                        <td>" << m_reference_2 << "</td>\n";
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td></td>\n";
  out << "                        <td>" << m_alignment_1 << "</td>\n";
  out << "                        <td></td>\n";
  out << "                        <td></td>\n";
  out << "                        <td>" << m_alignment_2 << "</td>\n";
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Consensus</td>\n";
  out << "                        <td>" << m_consensus_1 << "</td>\n";
  out << "                        <td></td>\n";
  out << "                        <td>Consensus</td>\n";
  out << "                        <td>" << m_consensus_2 << "</td>\n";
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Qualities</td>\n";
  out << "                        <td>" << m_qualities_1 << "</td>\n";
  out << "                        <td></td>\n";
  out << "                        <td>Qualities</td>\n";
  out << "                        <td>" << m_qualities_2 << "</td>\n";
  out << "                    </tr>\n";
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_consensus_adapter_kmer_head::~html_consensus_adapter_kmer_head()
{
  AR_REQUIRE(m_written);
}

html_consensus_adapter_kmer_head&
html_consensus_adapter_kmer_head::set_kmer_length(const std::string& value)
{
  m_kmer_length = value;
  m_kmer_length_is_set = true;
  return *this;
}

html_consensus_adapter_kmer_head&
html_consensus_adapter_kmer_head::set_n_kmers(const std::string& value)
{
  m_n_kmers = value;
  m_n_kmers_is_set = true;
  return *this;
}

void
html_consensus_adapter_kmer_head::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_kmer_length_is_set);
  AR_REQUIRE(m_n_kmers_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4 id=\"analyses-top-kmers\">\n";
  out << "                Top " << m_n_kmers << " most common " << m_kmer_length << "-bp 5' adapter k-mers <a class=\"anchor\"\n";
  out << "                    href=\"#analyses-top-kmers\">#</a>\n";
  out << "            </h4>\n";
  out << "\n";
  out << "            <table class=\"pure-table kmer-table pure-table-striped\" style=\"font-family: monospace;\">\n";
  out << "                <thead>\n";
  out << "                    <tr>\n";
  out << "                        <td></td>\n";
  out << "                        <td>Adapter 1 k-mer</td>\n";
  out << "                        <td>Count</td>\n";
  out << "                        <td>%</td>\n";
  out << "                        <td></td>\n";
  out << "                        <td>Adapter 2 k-mer</td>\n";
  out << "                        <td>Count</td>\n";
  out << "                        <td>%</td>\n";
  out << "                    </tr>\n";
  out << "                </thead>\n";
  out << "                <tbody>\n";
  // clang-format on
  m_written = true;
}

html_consensus_adapter_kmer_row::~html_consensus_adapter_kmer_row()
{
  AR_REQUIRE(m_written);
}

html_consensus_adapter_kmer_row&
html_consensus_adapter_kmer_row::set_count_1(const std::string& value)
{
  m_count_1 = value;
  m_count_1_is_set = true;
  return *this;
}

html_consensus_adapter_kmer_row&
html_consensus_adapter_kmer_row::set_count_2(const std::string& value)
{
  m_count_2 = value;
  m_count_2_is_set = true;
  return *this;
}

html_consensus_adapter_kmer_row&
html_consensus_adapter_kmer_row::set_index(const std::string& value)
{
  m_index = value;
  m_index_is_set = true;
  return *this;
}

html_consensus_adapter_kmer_row&
html_consensus_adapter_kmer_row::set_kmer_1(const std::string& value)
{
  m_kmer_1 = value;
  m_kmer_1_is_set = true;
  return *this;
}

html_consensus_adapter_kmer_row&
html_consensus_adapter_kmer_row::set_kmer_2(const std::string& value)
{
  m_kmer_2 = value;
  m_kmer_2_is_set = true;
  return *this;
}

html_consensus_adapter_kmer_row&
html_consensus_adapter_kmer_row::set_pct_1(const std::string& value)
{
  m_pct_1 = value;
  m_pct_1_is_set = true;
  return *this;
}

html_consensus_adapter_kmer_row&
html_consensus_adapter_kmer_row::set_pct_2(const std::string& value)
{
  m_pct_2 = value;
  m_pct_2_is_set = true;
  return *this;
}

void
html_consensus_adapter_kmer_row::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_count_1_is_set);
  AR_REQUIRE(m_count_2_is_set);
  AR_REQUIRE(m_index_is_set);
  AR_REQUIRE(m_kmer_1_is_set);
  AR_REQUIRE(m_kmer_2_is_set);
  AR_REQUIRE(m_pct_1_is_set);
  AR_REQUIRE(m_pct_2_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "                    <tr>\n";
  out << "                        <td>" << m_index << "</td>\n";
  out << "                        <td>" << m_kmer_1 << "</td>\n";
  out << "                        <td>" << m_count_1 << "</td>\n";
  out << "                        <td>" << m_pct_1 << "</td>\n";
  out << "                        <td></td>\n";
  out << "                        <td>" << m_kmer_2 << "</td>\n";
  out << "                        <td>" << m_count_2 << "</td>\n";
  out << "                        <td>" << m_pct_2 << "</td>\n";
  out << "                    </tr>\n";
  // clang-format on
  m_written = true;
}

html_consensus_adapter_kmer_tail::~html_consensus_adapter_kmer_tail()
{
  AR_REQUIRE(m_written);
}

void
html_consensus_adapter_kmer_tail::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_h2_tag::~html_h2_tag()
{
  AR_REQUIRE(m_written);
}

html_h2_tag&
html_h2_tag::set_href(const std::string& value)
{
  m_href = value;
  m_href_is_set = true;
  return *this;
}

html_h2_tag&
html_h2_tag::set_title(const std::string& value)
{
  m_title = value;
  m_title_is_set = true;
  return *this;
}

void
html_h2_tag::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_href_is_set);
  AR_REQUIRE(m_title_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "        </div>\n";
  out << "\n";
  out << "        <div class=\"section\">\n";
  out << "            <div class=\"title\">\n";
  out << "                <h2>" << m_title << " <a class=\"anchor\" id=\"" << m_href << "\" href=\"#" << m_href << "\">#</a></h2>\n";
  out << "            </div>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_plot_title::~html_plot_title()
{
  AR_REQUIRE(m_written);
}

html_plot_title&
html_plot_title::set_href(const std::string& value)
{
  m_href = value;
  m_href_is_set = true;
  return *this;
}

html_plot_title&
html_plot_title::set_title(const std::string& value)
{
  m_title = value;
  m_title_is_set = true;
  return *this;
}

void
html_plot_title::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_href_is_set);
  AR_REQUIRE(m_title_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4 id=\"" << m_href << "\">\n";
  out << "                " << m_title << " <a class=\"anchor\" href=\"#" << m_href << "\">#</a>\n";
  out << "            </h4>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_plot_sub_title::~html_plot_sub_title()
{
  AR_REQUIRE(m_written);
}

html_plot_sub_title&
html_plot_sub_title::set_sub_title(const std::string& value)
{
  m_sub_title = value;
  m_sub_title_is_set = true;
  return *this;
}

void
html_plot_sub_title::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_sub_title_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <div class=\"note\">" << m_sub_title << "</div>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_frequency_plot::~html_frequency_plot()
{
  AR_REQUIRE(m_written);
}

html_frequency_plot&
html_frequency_plot::set_legend(const std::string& value)
{
  m_legend = value;
  return *this;
}

html_frequency_plot&
html_frequency_plot::set_values(const std::string& value)
{
  m_values = value;
  m_values_is_set = true;
  return *this;
}

html_frequency_plot&
html_frequency_plot::set_width(const std::string& value)
{
  m_width = value;
  m_width_is_set = true;
  return *this;
}

html_frequency_plot&
html_frequency_plot::set_x_axis(const std::string& value)
{
  m_x_axis = value;
  return *this;
}

html_frequency_plot&
html_frequency_plot::set_y_axis(const std::string& value)
{
  m_y_axis = value;
  return *this;
}

void
html_frequency_plot::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_values_is_set);
  AR_REQUIRE(m_width_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <div id=\"vis_" << id << "\"></div>\n";
  out << "\n";
  out << "            <script>\n";
  out << "                var schema = {\n";
  out << "                    \"$schema\": \"https://vega.github.io/schema/vega-lite/v5.json\",\n";
  out << "                    \"width\": " << m_width << ",\n";
  out << "                    \"height\": 250,\n";
  out << "                    \"padding\": 5,\n";
  out << "                    \"autosize\": \"pad\",\n";
  out << "                    \"encoding\": {\n";
  out << "                        \"x\": {\n";
  out << "                            \"title\": " << m_x_axis  << ",\n";
  out << "                            \"field\": \"x\",\n";
  out << "                            \"type\": \"quantitative\",\n";
  out << "                            \"axis\": { \"offset\": 5 },\n";
  out << "                        },\n";
  out << "                    },\n";
  out << "                    \"layer\": [\n";
  out << "                        {\n";
  out << "                            \"encoding\": {\n";
  out << "                                \"y\": {\n";
  out << "                                    \"title\": " << m_y_axis  << ",\n";
  out << "                                    \"field\": \"y\",\n";
  out << "                                    \"type\": \"quantitative\",\n";
  out << "                                    \"axis\": { \"offset\": 5, \"minExtent\": 40 }\n";
  out << "                                },\n";
  out << "                                \"color\": {\n";
  out << "                                    \"field\": \"group\",\n";
  out << "                                    \"sort\": null,\n";
  out << "                                    \"legend\": " << m_legend  << "\n";
  out << "                                }\n";
  out << "                            },\n";
  out << "                            \"layer\": [\n";
  out << "                                { \"mark\": \"line\" },\n";
  out << "                                {\n";
  out << "                                    \"transform\": [{ \"filter\": { \"param\": \"hover\", \"empty\": false } }],\n";
  out << "                                    \"mark\": \"point\"\n";
  out << "                                }\n";
  out << "                            ]\n";
  out << "                        },\n";
  out << "                        {\n";
  out << "                            \"transform\": [{ \"pivot\": \"group\", \"value\": \"y\", \"groupby\": [\"x\"] }],\n";
  out << "                            \"mark\": { \"type\": \"rule\", \"tooltip\": { \"content\": \"data\" } },\n";
  out << "                            \"encoding\": {\n";
  out << "                                \"opacity\": {\n";
  out << "                                    \"condition\": { \"value\": 0.3, \"param\": \"hover\", \"empty\": false },\n";
  out << "                                    \"value\": 0\n";
  out << "                                }\n";
  out << "                            },\n";
  out << "                            \"params\": [\n";
  out << "                                {\n";
  out << "                                    \"name\": \"hover\",\n";
  out << "                                    \"select\": {\n";
  out << "                                        \"type\": \"point\",\n";
  out << "                                        \"fields\": [\"x\"],\n";
  out << "                                        \"nearest\": true,\n";
  out << "                                        \"on\": \"mouseover\",\n";
  out << "                                        \"clear\": \"mouseout\"\n";
  out << "                                    }\n";
  out << "                                }\n";
  out << "                            ]\n";
  out << "                        }\n";
  out << "                    ],\n";
  out << "                    \"transform\": [\n";
  out << "                        { \"flatten\": [\"y\"] },\n";
  out << "                        { \"window\": [{ \"op\": \"row_number\", \"as\": \"x\" }], \"groupby\": [\"group\"] },\n";
  out << "                        { \"calculate\": \"datum.x - 1 + (datum.offset)\", \"as\": \"x\" },\n";
  out << "                        { \"joinaggregate\": [{ \"op\": \"sum\", \"field\": \"y\", \"as\": \"sum_y\" }], \"groupby\": [\"group\"] },\n";
  out << "                        { \"calculate\": \"datum.y / max(1, datum.sum_y)\", \"as\": \"y\" },\n";
  out << "                    ],\n";
  out << "                    \"data\": {\n";
  out << "                        \"values\": " << m_values << "\n";
  out << "                    }\n";
  out << "                };\n";
  out << "\n";
  out << "                vegaEmbed('#vis_" << id << "', schema);\n";
  out << "            </script>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_facet_line_plot::~html_facet_line_plot()
{
  AR_REQUIRE(m_written);
}

html_facet_line_plot&
html_facet_line_plot::set_legend(const std::string& value)
{
  m_legend = value;
  return *this;
}

html_facet_line_plot&
html_facet_line_plot::set_values(const std::string& value)
{
  m_values = value;
  m_values_is_set = true;
  return *this;
}

html_facet_line_plot&
html_facet_line_plot::set_width(const std::string& value)
{
  m_width = value;
  m_width_is_set = true;
  return *this;
}

html_facet_line_plot&
html_facet_line_plot::set_x_axis(const std::string& value)
{
  m_x_axis = value;
  return *this;
}

html_facet_line_plot&
html_facet_line_plot::set_y_axis(const std::string& value)
{
  m_y_axis = value;
  return *this;
}

void
html_facet_line_plot::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_values_is_set);
  AR_REQUIRE(m_width_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <div id=\"vis_" << id << "\"></div>\n";
  out << "\n";
  out << "            <script>\n";
  out << "                var schema = {\n";
  out << "                    \"$schema\": \"https://vega.github.io/schema/vega-lite/v5.json\",\n";
  out << "                    \"padding\": 5,\n";
  out << "                    \"spec\": {\n";
  out << "                        \"width\": " << m_width << ",\n";
  out << "                        \"height\": 250,\n";
  out << "                        \"encoding\": {\n";
  out << "                            \"x\": {\n";
  out << "                                \"title\": " << m_x_axis  << ",\n";
  out << "                                \"field\": \"x\",\n";
  out << "                                \"type\": \"quantitative\",\n";
  out << "                                \"axis\": { \"offset\": 5 }\n";
  out << "                            }\n";
  out << "                        },\n";
  out << "                        \"layer\": [\n";
  out << "                            {\n";
  out << "                                \"encoding\": {\n";
  out << "                                    \"y\": {\n";
  out << "                                        \"title\": " << m_y_axis  << ",\n";
  out << "                                        \"field\": \"y\",\n";
  out << "                                        \"type\": \"quantitative\",\n";
  out << "                                        \"sort\": null,\n";
  out << "                                        \"axis\": { \"offset\": 5, \"minExtent\": 40 }\n";
  out << "                                    },\n";
  out << "                                    \"color\": {\n";
  out << "                                        \"field\": \"group\",\n";
  out << "                                        \"sort\": null,\n";
  out << "                                        \"legend\": " << m_legend  << "\n";
  out << "                                    }\n";
  out << "                                },\n";
  out << "                                \"layer\": [\n";
  out << "                                    { \"mark\": \"line\" },\n";
  out << "                                    {\n";
  out << "                                        \"transform\": [{ \"filter\": { \"param\": \"hover\", \"empty\": false } }],\n";
  out << "                                        \"mark\": \"point\"\n";
  out << "                                    }\n";
  out << "                                ]\n";
  out << "                            },\n";
  out << "                            {\n";
  out << "                                \"transform\": [{ \"pivot\": \"group\", \"value\": \"y\", \"groupby\": [\"x\"] }],\n";
  out << "                                \"mark\": { \"type\": \"rule\", \"tooltip\": { \"content\": \"data\" } },\n";
  out << "                                \"encoding\": {\n";
  out << "                                    \"opacity\": {\n";
  out << "                                        \"condition\": { \"value\": 0.3, \"param\": \"hover\", \"empty\": false },\n";
  out << "                                        \"value\": 0\n";
  out << "                                    }\n";
  out << "                                },\n";
  out << "                                \"params\": [\n";
  out << "                                    {\n";
  out << "                                        \"name\": \"hover\",\n";
  out << "                                        \"select\": {\n";
  out << "                                            \"type\": \"point\",\n";
  out << "                                            \"fields\": [\"x\"],\n";
  out << "                                            \"nearest\": true,\n";
  out << "                                            \"on\": \"mouseover\",\n";
  out << "                                            \"clear\": \"mouseout\"\n";
  out << "                                        }\n";
  out << "                                    }\n";
  out << "                                ]\n";
  out << "                            }\n";
  out << "                        ]\n";
  out << "                    },\n";
  out << "                    \"facet\": { \"field\": \"read\", \"type\": \"ordinal\", \"title\": \"\", \"sort\": [] },\n";
  out << "                    \"columns\": 2,\n";
  out << "                    \"transform\": [\n";
  out << "                        { \"flatten\": [\"y\"] },\n";
  out << "                        { \"window\": [{ \"op\": \"row_number\", \"as\": \"x\" }], \"groupby\": [\"read\", \"group\"] },\n";
  out << "                        { \"calculate\": \"datum.x - 1 + (datum.offset)\", \"as\": \"x\" },\n";
  out << "                    ],\n";
  out << "                    \"data\": {\n";
  out << "                        \"values\": " << m_values << "\n";
  out << "                    }\n";
  out << "                };\n";
  out << "\n";
  out << "                vegaEmbed('#vis_" << id << "', schema);\n";
  out << "            </script>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_bar_plot::~html_bar_plot()
{
  AR_REQUIRE(m_written);
}

html_bar_plot&
html_bar_plot::set_values(const std::string& value)
{
  m_values = value;
  m_values_is_set = true;
  return *this;
}

html_bar_plot&
html_bar_plot::set_width(const std::string& value)
{
  m_width = value;
  m_width_is_set = true;
  return *this;
}

html_bar_plot&
html_bar_plot::set_x_axis(const std::string& value)
{
  m_x_axis = value;
  return *this;
}

html_bar_plot&
html_bar_plot::set_y_axis(const std::string& value)
{
  m_y_axis = value;
  return *this;
}

void
html_bar_plot::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_values_is_set);
  AR_REQUIRE(m_width_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <div id=\"vis_" << id << "\"></div>\n";
  out << "\n";
  out << "            <script>\n";
  out << "                var schema = {\n";
  out << "                    \"$schema\": \"https://vega.github.io/schema/vega-lite/v5.json\",\n";
  out << "                    \"width\": " << m_width << ",\n";
  out << "                    \"height\": 250,\n";
  out << "                    \"padding\": 5,\n";
  out << "                    \"autosize\": \"pad\",\n";
  out << "                    \"mark\": \"bar\",\n";
  out << "                    \"encoding\": {\n";
  out << "                        \"x\": {\n";
  out << "                            \"title\": " << m_x_axis  << ",\n";
  out << "                            \"field\": \"x\",\n";
  out << "                            \"type\": \"nominal\",\n";
  out << "                            \"title\": \"\",\n";
  out << "                            \"axis\": { \"labels\": false }\n";
  out << "                        },\n";
  out << "                        \"y\": {\n";
  out << "                            \"title\": " << m_y_axis  << ",\n";
  out << "                            \"field\": \"y\",\n";
  out << "                            \"type\": \"quantitative\",\n";
  out << "                            \"axis\": { \"offset\": 5, \"minExtent\": 40 }\n";
  out << "                        },\n";
  out << "                        \"tooltip\": { \"field\": \"x\" }\n";
  out << "                    },\n";
  out << "                    \"data\": {\n";
  out << "                        \"values\": " << m_values << "\n";
  out << "                    }\n";
  out << "                };\n";
  out << "\n";
  out << "                vegaEmbed('#vis_" << id << "', schema);\n";
  out << "            </script>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_demultiplexing_head::~html_demultiplexing_head()
{
  AR_REQUIRE(m_written);
}

void
html_demultiplexing_head::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <h4 id=\"demux-per-sample\">\n";
  out << "                Per sample statistics <a class=\"anchor\" href=\"#demux-per-sample\">#</a>\n";
  out << "            </h4>\n";
  out << "\n";
  out << "            <div class=\"fixed-height-table\">\n";
  out << "                <table class=\"pure-table pure-table-striped\">\n";
  out << "                    <thead>\n";
  out << "                        <tr>\n";
  out << "                            <th>#</th>\n";
  out << "                            <th>Barcode #1</th>\n";
  out << "                            <th>Barcode #2</th>\n";
  out << "                            <th>Sample</th>\n";
  out << "                            <th>%</th>\n";
  out << "                            <th>#Reads</th>\n";
  out << "                            <th>#Bases</th>\n";
  out << "                            <th>Read length</th>\n";
  out << "                            <th>GC</th>\n";
  out << "                        </tr>\n";
  out << "                    </thead>\n";
  out << "                    <tbody>\n";
  // clang-format on
  m_written = true;
}

html_demultiplexing_row::~html_demultiplexing_row()
{
  AR_REQUIRE(m_written);
}

html_demultiplexing_row&
html_demultiplexing_row::set_barcode_1(const std::string& value)
{
  m_barcode_1 = value;
  m_barcode_1_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_barcode_2(const std::string& value)
{
  m_barcode_2 = value;
  m_barcode_2_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_bp(const std::string& value)
{
  m_bp = value;
  m_bp_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_gc(const std::string& value)
{
  m_gc = value;
  m_gc_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_length(const std::string& value)
{
  m_length = value;
  m_length_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_n(const std::string& value)
{
  m_n = value;
  m_n_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_name(const std::string& value)
{
  m_name = value;
  m_name_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_pct(const std::string& value)
{
  m_pct = value;
  m_pct_is_set = true;
  return *this;
}

html_demultiplexing_row&
html_demultiplexing_row::set_reads(const std::string& value)
{
  m_reads = value;
  m_reads_is_set = true;
  return *this;
}

void
html_demultiplexing_row::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  AR_REQUIRE(m_barcode_1_is_set);
  AR_REQUIRE(m_barcode_2_is_set);
  AR_REQUIRE(m_bp_is_set);
  AR_REQUIRE(m_gc_is_set);
  AR_REQUIRE(m_length_is_set);
  AR_REQUIRE(m_n_is_set);
  AR_REQUIRE(m_name_is_set);
  AR_REQUIRE(m_pct_is_set);
  AR_REQUIRE(m_reads_is_set);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "                        <tr>\n";
  out << "                            <td>" << m_n << "</td>\n";
  out << "                            <td>" << m_barcode_1 << "</td>\n";
  out << "                            <td>" << m_barcode_2 << "</td>\n";
  out << "                            <td>" << m_name << "</td>\n";
  out << "                            <td>" << m_pct << "</td>\n";
  out << "                            <td>" << m_reads << "</td>\n";
  out << "                            <td>" << m_bp << "</td>\n";
  out << "                            <td>" << m_length << "</td>\n";
  out << "                            <td>" << m_gc << " %</td>\n";
  out << "                        </tr>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_demultiplexing_tail::~html_demultiplexing_tail()
{
  AR_REQUIRE(m_written);
}

void
html_demultiplexing_tail::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "                    </tbody>\n";
  out << "                </table>\n";
  out << "            </div>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_body_end::~html_body_end()
{
  AR_REQUIRE(m_written);
}

void
html_body_end::write(std::ostream& out)
{
  AR_REQUIRE(!m_written);
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "        </div>\n";
  out << "\n";
  out << "        <div class=\"section note epilogue\">\n";
  out << "            <p>\n";
  out << "                If you use AdapterRemoval, please cite\n";
  out << "                <a href=\"https://doi.org/10.1186/s13104-016-1900-2\">Schubert et. al. 2016</a>:\n";
  out << "            </p>\n";
  out << "\n";
  out << "            <p style=\"font-family: monospace; padding-left: 2em; padding-right: 2em; max-width: 50%;\">\n";
  out << "                Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid adapter trimming, identification, and\n";
  out << "                read merging. BMC Research Notes, 12;9(1):88.\n";
  out << "            </p>\n";
  out << "\n";
  out << "            <p>\n";
  out << "                For comments, suggestions, or feedback, please use\n";
  out << "                <a href=\"https://github.com/MikkelSchubert/adapterremoval/issues/new\">GitHub Issues</a>.\n";
  out << "            </p>\n";
  out << "            <p>\n";
  out << "                This report was generated using <a href=\"https://purecss.io/\">Pure</a> (<a\n";
  out << "                    href=\"https://github.com/pure-css/pure/blob/v2.1.0/LICENSE\">license</a>) and <a\n";
  out << "                    href=\"https://vega.github.io/\">Vega-Lite</a> (<a\n";
  out << "                    href=\"https://github.com/vega/vega-lite/blob/v5.2.0/LICENSE\">license</a>).\n";
  out << "            </p>\n";
  out << "        </div>\n";
  out << "    </div>\n";
  out << "</body>\n";
  out << "\n";
  out << "</html>";
  // clang-format on
  m_written = true;
}

} // namespace adapterremoval
