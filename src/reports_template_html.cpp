/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Stinus Lindgreen - stinus@binf.ku.dk            *
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
#include "debug.hpp"

namespace adapterremoval {

size_t g_html_id = 1;

html_template::html_template()
{
  //
}

html_template::~html_template()
{
  //
}

html_head::html_head()
  : m_written()
  , m_name()
  , m_name_is_set()
  , m_version()
  , m_version_is_set()
{
  //
}

html_head::~html_head()
{
  AR_REQUIRE(m_written, "template html_head was not written");
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
html_head::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_head already written");
  AR_REQUIRE(m_name_is_set, "html_head::name not set");
  AR_REQUIRE(m_version_is_set, "html_head::version not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "<!DOCTYPE html>\n";
  out << "<html lang='en'>\n";
  out << "\n";
  out << "<head>\n";
  out << "    <meta charset='utf-8'>\n";
  out << "    <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n";
  out << "    <title>" << m_name << " " << m_version << "</title>\n";
  out << "    <link rel='stylesheet' href='https://unpkg.com/purecss@2.1.0/build/pure-min.css'\n";
  out << "        integrity='sha384-yHIFVG6ClnONEA5yB5DJXfW2/KC173DIQrYoZMEtBvGzmf0PKiGyNEqe9N6BNDBH' crossorigin='anonymous'>\n";
  out << "    <script src=\"https://unpkg.com/vega@5.21.0/build/vega.min.js\"\n";
  out << "        integrity=\"sha384-s2nYi9D0FfKNopEKsfINeS1Ffhcf+5uvwIrb7Zqso2II+HPhzBTWvXClt+NdUwFc\"\n";
  out << "        crossorigin=\"anonymous\"></script>\n";
  out << "    <script src=\"https://unpkg.com/vega-lite@5.2.0/build/vega-lite.min.js\"\n";
  out << "        integrity=\"sha384-tU6fj0fI2gxrcWwC7uBMp70QvipC9ukjcXyOs85VMmdCq33CrA7xQ3nJkJu0SmDm\"\n";
  out << "        crossorigin=\"anonymous\"></script>\n";
  out << "    <script src=\"https://unpkg.com/vega-embed@6.20.2/build/vega-embed.min.js\"\n";
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
  out << "            width: 100px;\n";
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
  out << "        .epilogue,\n";
  out << "        .note {\n";
  out << "            color: #777;\n";
  out << "            font-size: small;\n";
  out << "            padding-top: 10px;\n";
  out << "        }\n";
  out << "\n";
  out << "    </style>\n";
  out << "</head>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_body_start::html_body_start()
  : m_written()
{
  //
}

html_body_start::~html_body_start()
{
  AR_REQUIRE(m_written, "template html_body_start was not written");
}

void
html_body_start::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_body_start already written");
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

html_summary::html_summary()
  : m_written()
  , m_command()
  , m_command_is_set()
  , m_date_and_time()
  , m_date_and_time_is_set()
  , m_runtime()
  , m_runtime_is_set()
  , m_version()
  , m_version_is_set()
{
  //
}

html_summary::~html_summary()
{
  AR_REQUIRE(m_written, "template html_summary was not written");
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
html_summary::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary already written");
  AR_REQUIRE(m_command_is_set, "html_summary::command not set");
  AR_REQUIRE(m_date_and_time_is_set, "html_summary::date_and_time not set");
  AR_REQUIRE(m_runtime_is_set, "html_summary::runtime not set");
  AR_REQUIRE(m_version_is_set, "html_summary::version not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "        <div class=\"section\">\n";
  out << "            <div class=\"title\">\n";
  out << "                <h2>Summary</h2>\n";
  out << "            </div>\n";
  out << "\n";
  out << "            <h4>Program</h4>\n";
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

html_summary_io::html_summary_io()
  : m_written()
  , m_columns()
  , m_columns_is_set()
  , m_gc()
  , m_gc_is_set()
  , m_lengths()
  , m_lengths_is_set()
  , m_n_bases()
  , m_n_bases_is_set()
  , m_n_reads()
  , m_n_reads_is_set()
  , m_ns()
  , m_ns_is_set()
  , m_q20()
  , m_q20_is_set()
  , m_q30()
  , m_q30_is_set()
  , m_title()
  , m_title_is_set()
{
  //
}

html_summary_io::~html_summary_io()
{
  AR_REQUIRE(m_written, "template html_summary_io was not written");
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
html_summary_io::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary_io already written");
  AR_REQUIRE(m_columns_is_set, "html_summary_io::columns not set");
  AR_REQUIRE(m_gc_is_set, "html_summary_io::gc not set");
  AR_REQUIRE(m_lengths_is_set, "html_summary_io::lengths not set");
  AR_REQUIRE(m_n_bases_is_set, "html_summary_io::n_bases not set");
  AR_REQUIRE(m_n_reads_is_set, "html_summary_io::n_reads not set");
  AR_REQUIRE(m_ns_is_set, "html_summary_io::ns not set");
  AR_REQUIRE(m_q20_is_set, "html_summary_io::q20 not set");
  AR_REQUIRE(m_q30_is_set, "html_summary_io::q30 not set");
  AR_REQUIRE(m_title_is_set, "html_summary_io::title not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4>" << m_title << "</h4>\n";
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
    out << "                        <td>" << value << " bp</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Q20</td>\n";
  for (const auto& value : m_q20) {
    out << "                        <td>" << value << " %</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>Q30</td>\n";
  for (const auto& value : m_q30) {
    out << "                        <td>" << value << " %</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>GC</td>\n";
  for (const auto& value : m_gc) {
    out << "                        <td>" << value << " %</td>\n";
  }
  out << "                    </tr>\n";
  out << "                    <tr>\n";
  out << "                        <td>N</td>\n";
  for (const auto& value : m_ns) {
    out << "                        <td>" << value << " %</td>\n";
  }
  out << "                    </tr>\n";
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_sampling_note::html_sampling_note()
  : m_written()
  , m_label()
  , m_label_is_set()
  , m_pct()
  , m_pct_is_set()
{
  //
}

html_sampling_note::~html_sampling_note()
{
  AR_REQUIRE(m_written, "template html_sampling_note was not written");
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

void
html_sampling_note::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_sampling_note already written");
  AR_REQUIRE(m_label_is_set, "html_sampling_note::label not set");
  AR_REQUIRE(m_pct_is_set, "html_sampling_note::pct not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                Base composition statistics/plots are based on " << m_pct << "% of " << m_label << " reads sampled during execution.\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_output_note_pe::html_output_note_pe()
  : m_written()
{
  //
}

html_output_note_pe::~html_output_note_pe()
{
  AR_REQUIRE(m_written, "template html_output_note_pe was not written");
}

void
html_output_note_pe::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_output_note_pe already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                <b>*</b> The <b>Passed</b> column includes all read types except for <b>Discarded</b> reads.\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_output_note_se::html_output_note_se()
  : m_written()
{
  //
}

html_output_note_se::~html_output_note_se()
{
  AR_REQUIRE(m_written, "template html_output_note_se was not written");
}

void
html_output_note_se::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_output_note_se already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                <b>*</b> <b>Discarded</b> reads are not included in the <b>Output</b> column.\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_summary_trimming_head::html_summary_trimming_head()
  : m_written()
{
  //
}

html_summary_trimming_head::~html_summary_trimming_head()
{
  AR_REQUIRE(m_written, "template html_summary_trimming_head was not written");
}

void
html_summary_trimming_head::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary_trimming_head already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4>Processing</h4>\n";
  out << "            <table class=\"pure-table trimming-table pure-table-striped\">\n";
  out << "                <thead>\n";
  out << "                    <tr>\n";
  out << "                        <th>Stage</th>\n";
  out << "                        <th>Trimmed</th>\n";
  out << "                        <th>X</th>\n";
  out << "                        <th>Reads</th>\n";
  out << "                        <th>Bases</th>\n";
  out << "                        <th>Bases/Read</th>\n";
  out << "                    </tr>\n";
  out << "                </thead>\n";
  out << "                <tbody>\n";
  // clang-format on
  m_written = true;
}

html_summary_trimming_row::html_summary_trimming_row()
  : m_written()
  , m_avg_bases()
  , m_avg_bases_is_set()
  , m_bases()
  , m_bases_is_set()
  , m_label_1()
  , m_label_1_is_set()
  , m_label_2()
  , m_label_2_is_set()
  , m_reads()
  , m_reads_is_set()
  , m_stage()
  , m_stage_is_set()
{
  //
}

html_summary_trimming_row::~html_summary_trimming_row()
{
  AR_REQUIRE(m_written, "template html_summary_trimming_row was not written");
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
html_summary_trimming_row::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary_trimming_row already written");
  AR_REQUIRE(m_avg_bases_is_set, "html_summary_trimming_row::avg_bases not set");
  AR_REQUIRE(m_bases_is_set, "html_summary_trimming_row::bases not set");
  AR_REQUIRE(m_label_1_is_set, "html_summary_trimming_row::label_1 not set");
  AR_REQUIRE(m_label_2_is_set, "html_summary_trimming_row::label_2 not set");
  AR_REQUIRE(m_reads_is_set, "html_summary_trimming_row::reads not set");
  AR_REQUIRE(m_stage_is_set, "html_summary_trimming_row::stage not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "                    <tr>\n";
  out << "                        <td>" << m_stage << "</td>\n";
  out << "                        <td>" << m_label_1 << "</td>\n";
  out << "                        <td>" << m_label_2 << "</td>\n";
  out << "                        <td>" << m_reads << "</td>\n";
  out << "                        <td>" << m_bases << "</td>\n";
  out << "                        <td>" << m_avg_bases << "</td>\n";
  out << "                    </tr>\n";
  // clang-format on
  m_written = true;
}

html_summary_trimming_tail::html_summary_trimming_tail()
  : m_written()
  , m_n_enabled()
  , m_n_enabled_is_set()
  , m_n_total()
  , m_n_total_is_set()
{
  //
}

html_summary_trimming_tail::~html_summary_trimming_tail()
{
  AR_REQUIRE(m_written, "template html_summary_trimming_tail was not written");
}

html_summary_trimming_tail&
html_summary_trimming_tail::set_n_enabled(const std::string& value)
{
  m_n_enabled = value;
  m_n_enabled_is_set = true;
  return *this;
}

html_summary_trimming_tail&
html_summary_trimming_tail::set_n_total(const std::string& value)
{
  m_n_total = value;
  m_n_total_is_set = true;
  return *this;
}

void
html_summary_trimming_tail::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary_trimming_tail already written");
  AR_REQUIRE(m_n_enabled_is_set, "html_summary_trimming_tail::n_enabled not set");
  AR_REQUIRE(m_n_total_is_set, "html_summary_trimming_tail::n_total not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                " << m_n_enabled << " of " << m_n_total << " trimming steps enabled. Numbers are\n";
  out << "                given in terms of input reads, also for merged reads.\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_summary_filtering_head::html_summary_filtering_head()
  : m_written()
{
  //
}

html_summary_filtering_head::~html_summary_filtering_head()
{
  AR_REQUIRE(m_written, "template html_summary_filtering_head was not written");
}

void
html_summary_filtering_head::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary_filtering_head already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "            <h4>Filtering</h4>\n";
  out << "            <table class=\"pure-table io-table pure-table-striped\">\n";
  out << "                <thead>\n";
  out << "                    <tr>\n";
  out << "                        <th></th>\n";
  out << "                        <th>Reads</th>\n";
  out << "                        <th>Bases</th>\n";
  out << "                        <th>Bases/Read</th>\n";
  out << "                    </tr>\n";
  out << "                </thead>\n";
  out << "                <tbody>\n";
  // clang-format on
  m_written = true;
}

html_summary_filtering_row::html_summary_filtering_row()
  : m_written()
  , m_avg_bases()
  , m_avg_bases_is_set()
  , m_bases()
  , m_bases_is_set()
  , m_label()
  , m_label_is_set()
  , m_reads()
  , m_reads_is_set()
{
  //
}

html_summary_filtering_row::~html_summary_filtering_row()
{
  AR_REQUIRE(m_written, "template html_summary_filtering_row was not written");
}

html_summary_filtering_row&
html_summary_filtering_row::set_avg_bases(const std::string& value)
{
  m_avg_bases = value;
  m_avg_bases_is_set = true;
  return *this;
}

html_summary_filtering_row&
html_summary_filtering_row::set_bases(const std::string& value)
{
  m_bases = value;
  m_bases_is_set = true;
  return *this;
}

html_summary_filtering_row&
html_summary_filtering_row::set_label(const std::string& value)
{
  m_label = value;
  m_label_is_set = true;
  return *this;
}

html_summary_filtering_row&
html_summary_filtering_row::set_reads(const std::string& value)
{
  m_reads = value;
  m_reads_is_set = true;
  return *this;
}

void
html_summary_filtering_row::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary_filtering_row already written");
  AR_REQUIRE(m_avg_bases_is_set, "html_summary_filtering_row::avg_bases not set");
  AR_REQUIRE(m_bases_is_set, "html_summary_filtering_row::bases not set");
  AR_REQUIRE(m_label_is_set, "html_summary_filtering_row::label not set");
  AR_REQUIRE(m_reads_is_set, "html_summary_filtering_row::reads not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "                    <tr>\n";
  out << "                        <td>" << m_label << "</td>\n";
  out << "                        <td>" << m_reads << "</td>\n";
  out << "                        <td>" << m_bases << "</td>\n";
  out << "                        <td>" << m_avg_bases << "</td>\n";
  out << "                    </tr>\n";
  // clang-format on
  m_written = true;
}

html_summary_filtering_tail::html_summary_filtering_tail()
  : m_written()
  , m_n_enabled()
  , m_n_enabled_is_set()
  , m_n_total()
  , m_n_total_is_set()
{
  //
}

html_summary_filtering_tail::~html_summary_filtering_tail()
{
  AR_REQUIRE(m_written, "template html_summary_filtering_tail was not written");
}

html_summary_filtering_tail&
html_summary_filtering_tail::set_n_enabled(const std::string& value)
{
  m_n_enabled = value;
  m_n_enabled_is_set = true;
  return *this;
}

html_summary_filtering_tail&
html_summary_filtering_tail::set_n_total(const std::string& value)
{
  m_n_total = value;
  m_n_total_is_set = true;
  return *this;
}

void
html_summary_filtering_tail::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_summary_filtering_tail already written");
  AR_REQUIRE(m_n_enabled_is_set, "html_summary_filtering_tail::n_enabled not set");
  AR_REQUIRE(m_n_total_is_set, "html_summary_filtering_tail::n_total not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "                </tbody>\n";
  out << "            </table>\n";
  out << "\n";
  out << "            <p class=\"note\">\n";
  out << "                " << m_n_enabled << " of " << m_n_total << " filtering steps enabled.\n";
  out << "            </p>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_h2_tag::html_h2_tag()
  : m_written()
  , m_title()
  , m_title_is_set()
{
  //
}

html_h2_tag::~html_h2_tag()
{
  AR_REQUIRE(m_written, "template html_h2_tag was not written");
}

html_h2_tag&
html_h2_tag::set_title(const std::string& value)
{
  m_title = value;
  m_title_is_set = true;
  return *this;
}

void
html_h2_tag::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_h2_tag already written");
  AR_REQUIRE(m_title_is_set, "html_h2_tag::title not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "        </div>\n";
  out << "\n";
  out << "        <div class=\"section\">\n";
  out << "            <div class=\"title\">\n";
  out << "                <h2>" << m_title << "</h2>\n";
  out << "            </div>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_line_plot::html_line_plot()
  : m_written()
  , m_legend("{ \"title\": \"Legend\", \"padding\": 5 }")
  , m_sub_title("\"\"")
  , m_title()
  , m_title_is_set()
  , m_title_anchor("null")
  , m_values()
  , m_values_is_set()
  , m_width()
  , m_width_is_set()
  , m_x_axis("null")
  , m_y_axis("null")
{
  //
}

html_line_plot::~html_line_plot()
{
  AR_REQUIRE(m_written, "template html_line_plot was not written");
}

html_line_plot&
html_line_plot::set_legend(const std::string& value)
{
  m_legend = value;
  return *this;
}

html_line_plot&
html_line_plot::set_sub_title(const std::string& value)
{
  m_sub_title = value;
  return *this;
}

html_line_plot&
html_line_plot::set_title(const std::string& value)
{
  m_title = value;
  m_title_is_set = true;
  return *this;
}

html_line_plot&
html_line_plot::set_title_anchor(const std::string& value)
{
  m_title_anchor = value;
  return *this;
}

html_line_plot&
html_line_plot::set_values(const std::string& value)
{
  m_values = value;
  m_values_is_set = true;
  return *this;
}

html_line_plot&
html_line_plot::set_width(const std::string& value)
{
  m_width = value;
  m_width_is_set = true;
  return *this;
}

html_line_plot&
html_line_plot::set_x_axis(const std::string& value)
{
  m_x_axis = value;
  return *this;
}

html_line_plot&
html_line_plot::set_y_axis(const std::string& value)
{
  m_y_axis = value;
  return *this;
}

void
html_line_plot::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_line_plot already written");
  AR_REQUIRE(m_title_is_set, "html_line_plot::title not set");
  AR_REQUIRE(m_values_is_set, "html_line_plot::values not set");
  AR_REQUIRE(m_width_is_set, "html_line_plot::width not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <div id=\"vis_" << id << "\"></div>\n";
  out << "\n";
  out << "            <script>\n";
  out << "                var schema = {\n";
  out << "                    \"$schema\": \"https://vega.github.io/schema/vega-lite/v5.json\",\n";
  out << "                    \"title\": {\n";
  out << "                        \"text\": " << m_title << ",\n";
  out << "                        \"subtitle\": " << m_sub_title  << ",\n";
  out << "                        \"anchor\": " << m_title_anchor  << ",\n";
  out << "                    },\n";
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

html_facet_line_plot::html_facet_line_plot()
  : m_written()
  , m_legend("{ \"title\": \"Legend\", \"padding\": 5 }")
  , m_title()
  , m_title_is_set()
  , m_values()
  , m_values_is_set()
  , m_width()
  , m_width_is_set()
  , m_x_axis("null")
  , m_y_axis("null")
{
  //
}

html_facet_line_plot::~html_facet_line_plot()
{
  AR_REQUIRE(m_written, "template html_facet_line_plot was not written");
}

html_facet_line_plot&
html_facet_line_plot::set_legend(const std::string& value)
{
  m_legend = value;
  return *this;
}

html_facet_line_plot&
html_facet_line_plot::set_title(const std::string& value)
{
  m_title = value;
  m_title_is_set = true;
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
html_facet_line_plot::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_facet_line_plot already written");
  AR_REQUIRE(m_title_is_set, "html_facet_line_plot::title not set");
  AR_REQUIRE(m_values_is_set, "html_facet_line_plot::values not set");
  AR_REQUIRE(m_width_is_set, "html_facet_line_plot::width not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <div id=\"vis_" << id << "\"></div>\n";
  out << "\n";
  out << "            <script>\n";
  out << "                var schema = {\n";
  out << "                    \"$schema\": \"https://vega.github.io/schema/vega-lite/v5.json\",\n";
  out << "                    \"title\": " << m_title << ",\n";
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

html_bar_plot::html_bar_plot()
  : m_written()
  , m_title()
  , m_title_is_set()
  , m_values()
  , m_values_is_set()
  , m_width()
  , m_width_is_set()
  , m_x_axis("null")
  , m_y_axis("null")
{
  //
}

html_bar_plot::~html_bar_plot()
{
  AR_REQUIRE(m_written, "template html_bar_plot was not written");
}

html_bar_plot&
html_bar_plot::set_title(const std::string& value)
{
  m_title = value;
  m_title_is_set = true;
  return *this;
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
html_bar_plot::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_bar_plot already written");
  AR_REQUIRE(m_title_is_set, "html_bar_plot::title not set");
  AR_REQUIRE(m_values_is_set, "html_bar_plot::values not set");
  AR_REQUIRE(m_width_is_set, "html_bar_plot::width not set");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <div id=\"vis_" << id << "\"></div>\n";
  out << "\n";
  out << "            <script>\n";
  out << "                var schema = {\n";
  out << "                    \"$schema\": \"https://vega.github.io/schema/vega-lite/v5.json\",\n";
  out << "                    \"title\": { \"text\": " << m_title << ", \"anchor\": \"start\" },\n";
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

html_demultiplexing_head::html_demultiplexing_head()
  : m_written()
{
  //
}

html_demultiplexing_head::~html_demultiplexing_head()
{
  AR_REQUIRE(m_written, "template html_demultiplexing_head was not written");
}

void
html_demultiplexing_head::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_demultiplexing_head already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "\n";
  out << "            <h5>Per sample statistics</h5>\n";
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

html_demultiplexing_row::html_demultiplexing_row()
  : m_written()
  , m_barcode_1()
  , m_barcode_1_is_set()
  , m_barcode_2()
  , m_barcode_2_is_set()
  , m_bp()
  , m_bp_is_set()
  , m_gc()
  , m_gc_is_set()
  , m_length()
  , m_length_is_set()
  , m_n()
  , m_n_is_set()
  , m_name()
  , m_name_is_set()
  , m_pct()
  , m_pct_is_set()
  , m_reads()
  , m_reads_is_set()
{
  //
}

html_demultiplexing_row::~html_demultiplexing_row()
{
  AR_REQUIRE(m_written, "template html_demultiplexing_row was not written");
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
html_demultiplexing_row::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_demultiplexing_row already written");
  AR_REQUIRE(m_barcode_1_is_set, "html_demultiplexing_row::barcode_1 not set");
  AR_REQUIRE(m_barcode_2_is_set, "html_demultiplexing_row::barcode_2 not set");
  AR_REQUIRE(m_bp_is_set, "html_demultiplexing_row::bp not set");
  AR_REQUIRE(m_gc_is_set, "html_demultiplexing_row::gc not set");
  AR_REQUIRE(m_length_is_set, "html_demultiplexing_row::length not set");
  AR_REQUIRE(m_n_is_set, "html_demultiplexing_row::n not set");
  AR_REQUIRE(m_name_is_set, "html_demultiplexing_row::name not set");
  AR_REQUIRE(m_pct_is_set, "html_demultiplexing_row::pct not set");
  AR_REQUIRE(m_reads_is_set, "html_demultiplexing_row::reads not set");
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
  out << "                            <td>" << m_length << " bp</td>\n";
  out << "                            <td>" << m_gc << " %</td>\n";
  out << "                        </tr>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_demultiplexing_tail::html_demultiplexing_tail()
  : m_written()
{
  //
}

html_demultiplexing_tail::~html_demultiplexing_tail()
{
  AR_REQUIRE(m_written, "template html_demultiplexing_tail was not written");
}

void
html_demultiplexing_tail::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_demultiplexing_tail already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "                    </tbody>\n";
  out << "                </table>\n";
  out << "            </div>\n";
  out << "\n";
  // clang-format on
  m_written = true;
}

html_body_end::html_body_end()
  : m_written()
{
  //
}

html_body_end::~html_body_end()
{
  AR_REQUIRE(m_written, "template html_body_end was not written");
}

void
html_body_end::write(std::ofstream& out)
{
  AR_REQUIRE(!m_written, "template html_body_end already written");
  // clang-format off
  auto id = g_html_id++; (void)id;
  out << "        </div>\n";
  out << "\n";
  out << "        <div class=\"section epilogue\">\n";
  out << "            <p>\n";
  out << "                If you use AdapterRemoval, please cite\n";
  out << "                <a href=\"https://doi.org/10.1186/s13104-016-1900-2\">Schubert et. al. 2016</a>:\n";
  out << "\n";
  out << "            <pre>    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid\n";
  out << "    adapter trimming, identification, and read merging.\n";
  out << "    BMC Research Notes, 12;9(1):88.</pre>\n";
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
  out << "\n";
  out << "</body>\n";
  out << "\n";
  out << "</html>\n";
  // clang-format on
  m_written = true;
}

} // namespace adapterremoval
