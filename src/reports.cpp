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
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>

#include "json.hpp"
#include "main.hpp"
#include "strutils.hpp"
#include "userconfig.hpp"

#define WITH_SECTION(writer, key)                                              \
  if (const auto section##__LINE__ = (writer).start(key))

void
write_report_meta(const userconfig& config, json_writer& writer)
{
  WITH_SECTION(writer, "meta")
  {
    writer.write("version", NAME + " " + VERSION);
    writer.write("command", config.args);
    writer.write_float("runtime", config.runtime());
  }
}

void
write_report_summary_stats(json_writer& writer,
                           const std::vector<const fastq_statistics*>& stats)
{
  size_t n_reads = 0;
  size_t n_bases = 0;
  size_t n_g = 0;
  size_t n_c = 0;
  size_t n_n = 0;
  size_t n_q20 = 0;
  size_t n_q30 = 0;

  for (const auto& it : stats) {
    n_reads += it->length_dist().sum();
    n_bases += it->length_dist().product();
    n_g += it->nucleotides_pos('G').sum();
    n_c += it->nucleotides_pos('C').sum();
    n_n += it->uncalled_pos().sum();
    n_q20 += it->quality_dist().sum(20);
    n_q30 += it->quality_dist().sum(30);
  }

  writer.write_int("reads", n_reads);
  writer.write_int("bases", n_bases);
  writer.write_float("mean_length", static_cast<double>(n_bases) / n_reads);
  writer.write_int("q20_bases", n_q20);
  writer.write_int("q30_bases", n_q30);
  writer.write_float("q20_rate", static_cast<double>(n_q20) / n_bases);
  writer.write_float("q30_rate", static_cast<double>(n_q30) / n_bases);
  writer.write_int("uncalled_bases", n_n);
  writer.write_float("uncalled_rate", static_cast<double>(n_n) / n_bases);
  writer.write_float("gc_content", static_cast<double>(n_g + n_c) / n_bases);
}

void
write_read_length(json_writer& writer,
                  const std::string& label,
                  size_t reads,
                  size_t bases)
{
  writer.write_int(label + "_reads", reads);
  writer.write_float(label + "_mean_length",
                     static_cast<double>(bases) / reads);
}

void
write_report_summary(const userconfig& config,
                     json_writer& writer,
                     const ar_statistics& stats)
{
  WITH_SECTION(writer, "summary")
  {
    WITH_SECTION(writer, "input")
    {
      std::vector<const fastq_statistics*> input = { &stats.input_1,
                                                     &stats.input_2 };
      write_report_summary_stats(writer, input);
    }

    if (config.adapters.barcode_count()) {
      WITH_SECTION(writer, "demultiplexing")
      {
        const auto& demux = stats.demultiplexing;
        const auto total = demux.total();

        writer.write_int("total_reads", total);
        writer.write_int("identified_reads",
                         total - demux.unidentified - demux.ambiguous);
        writer.write_int("unidentified_reads", demux.unidentified);
        writer.write_int("ambiguous_reads", demux.ambiguous);

        WITH_SECTION(writer, "samples")
        {
          for (size_t i = 0; i < demux.barcodes.size(); ++i) {
            writer.write_int(config.adapters.get_sample_name(i),
                             demux.barcodes.at(i));
          }
        }
      }
    }

    WITH_SECTION(writer, "output")
    {
      std::vector<const fastq_statistics*> output;
      for (const auto& it : stats.trimming) {
        output.push_back(&it.read_1);
        output.push_back(&it.read_2);
        output.push_back(&it.merged);
      }

      write_report_summary_stats(writer, output);
    }

    WITH_SECTION(writer, "discarded")
    {
      std::vector<const fastq_statistics*> discarded;
      for (const auto& it : stats.trimming) {
        discarded.push_back(&it.discarded);
      }

      write_report_summary_stats(writer, discarded);
    }
  }
}

void
write_io_section(const std::string& key,
                 json_writer& writer,
                 const fastq_statistics& stats)
{
  WITH_SECTION(writer, key)
  {
    const std::vector<const fastq_statistics*> read_1 = { &stats };

    auto total_bases = stats.uncalled_pos();
    auto total_quality = stats.uncalled_quality_pos();

    writer.write("lengths", stats.length_dist());

    WITH_SECTION(writer, "quality_curves")
    {
      for (size_t nuc_i = 0; nuc_i < 4; ++nuc_i) {
        const auto nuc = IDX_TO_ACGT(nuc_i);
        const auto& nucleotides = stats.nucleotides_pos(nuc);
        const auto& quality = stats.qualities_pos(nuc);

        total_bases += nucleotides;
        total_quality += quality;

        writer.write(std::string(1, nuc), quality / nucleotides);
      }

      writer.write("Mean", total_quality / total_bases);
    }

    WITH_SECTION(writer, "content_curves")
    {
      for (size_t nuc_i = 0; nuc_i < 4; ++nuc_i) {
        const auto nuc = IDX_TO_ACGT(nuc_i);
        const auto& bases = stats.nucleotides_pos(nuc);

        writer.write(std::string(1, nuc), bases / total_bases);
      }

      writer.write("N", stats.uncalled_quality_pos() / total_bases);
      writer.write("GC",
                   (stats.nucleotides_pos('G') + stats.nucleotides_pos('C')) /
                     total_bases);
    }
  }
}

void
write_report_input(const userconfig& config,
                   json_writer& writer,
                   const ar_statistics& stats)
{
  WITH_SECTION(writer, "input")
  {
    write_io_section("read1", writer, stats.input_1);

    if (config.paired_ended_mode) {
      write_io_section("read2", writer, stats.input_2);
    }
  }
}

void
write_report_demultiplexing(const userconfig& config,
                            json_writer& writer,
                            const ar_statistics& stats)
{
  if (config.adapters.barcode_count()) {
    WITH_SECTION(writer, "demultiplexing")
    {
      size_t assigned_reads = 0;
      const auto& demux = stats.demultiplexing;
      for (size_t it : demux.barcodes) {
        assigned_reads += it;
      }

      writer.write_int("assigned_reads", assigned_reads);
      writer.write_int("ambiguous_reads", demux.ambiguous);
      writer.write_int("unassigned_reads", demux.unidentified);

      WITH_SECTION(writer, "samples")
      {
        for (size_t i = 0; i < demux.barcodes.size(); ++i) {
          WITH_SECTION(writer, config.adapters.get_sample_name(i))
          {
            writer.write_int("reads", demux.barcodes.at(i));

            WITH_SECTION(writer, "output")
            {
              write_io_section("read1", writer, stats.trimming.at(i).read_1);

              if (config.paired_ended_mode) {
                write_io_section("read2", writer, stats.trimming.at(i).read_2);

                if (config.collapse) {
                  write_io_section(
                    "merged", writer, stats.trimming.at(i).merged);
                }
              }

              write_io_section(
                "discarded", writer, stats.trimming.at(i).discarded);
            }
          }
        }
      }
    }
  } else {
    writer.write_null("demultiplexing");
  }
}

void
write_report_trimming(const userconfig& config,
                      json_writer& writer,
                      const ar_statistics& stats)
{
  WITH_SECTION(writer, "trimming_and_filtering")
  {
#warning TODO: Trimming/Filtering statistics
  }
}

void
write_report_output(const userconfig& config,
                    json_writer& writer,
                    const ar_statistics& stats)
{
  fastq_statistics output_1;
  fastq_statistics output_2;
  fastq_statistics merged;
  fastq_statistics discarded;

  for (const auto& it : stats.trimming) {
    output_1 += it.read_1;
    output_2 += it.read_2;
    merged += it.merged;
    discarded += it.discarded;
  }

  WITH_SECTION(writer, "output")
  {
    write_io_section("read1", writer, output_1);

    if (config.paired_ended_mode) {
      write_io_section("read2", writer, output_2);

      if (config.collapse) {
        write_io_section("merged", writer, merged);
      }
    }

    write_io_section("discarded", writer, discarded);
  }
}

bool
write_report(const userconfig& config, const ar_statistics& stats)
{
  const auto filename = config.get_output_filename("--settings");

  try {
    std::ofstream output(filename, std::ofstream::out);
    if (!output.is_open()) {
      throw std::ofstream::failure(std::strerror(errno));
    }

    output.exceptions(std::ofstream::failbit | std::ofstream::badbit);

    {
      json_writer writer(output);
      write_report_meta(config, writer);
      write_report_summary(config, writer, stats);
      write_report_input(config, writer, stats);
      write_report_demultiplexing(config, writer, stats);
      write_report_trimming(config, writer, stats);
      write_report_output(config, writer, stats);
    }

    output << std::endl;
  } catch (const std::ios_base::failure& error) {
    std::cerr << "Error writing JSON report to '" << filename << "':\n"
              << cli_formatter::fmt(error.what()) << std::endl;
    return false;
  }

  return true;
}
