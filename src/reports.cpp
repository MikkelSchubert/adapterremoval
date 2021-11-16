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
#include <cstring>  // for size_t, strerror
#include <errno.h>  // for errno
#include <fstream>  // for ofstream
#include <iostream> // for ofstream, operator<<, basic_ostream, endl
#include <string>   // for operator+, string, operator<<
#include <vector>   // for vector

#include "adapterset.hpp" // for adapter_set
#include "counts.hpp"     // for counts, counts_tmpl
#include "debug.hpp"      // for AR_DEBUG_FAIL
#include "fastq.hpp"      // for fastq_pair_vec, IDX_TO_ACGT, fastq
#include "json.hpp"       // for json_writer, json_section
#include "main.hpp"       // for NAME, VERSION
#include "statistics.hpp" // for fastq_statistics, trimming_statistics, ar_...
#include "strutils.hpp"   // for cli_formatter
#include "userconfig.hpp" // for userconfig, ar_command, ar_command::demult...

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
    n_reads += it->number_of_input_reads();
    n_reads_s += it->number_of_sampled_reads();
    n_bases += it->length_dist().product();
    // The following stats are all (potentially) based on a subset of reads
    n_a += it->nucleotides_pos('A').sum();
    n_c += it->nucleotides_pos('C').sum();
    n_g += it->nucleotides_pos('G').sum();
    n_t += it->nucleotides_pos('T').sum();
    n_n += it->uncalled_pos().sum();
    n_q20 += it->quality_dist().sum(20);
    n_q30 += it->quality_dist().sum(30);
  }

  const auto n_bases_s = n_a + n_c + n_g + n_t + n_n;

  writer.write_int("reads", n_reads);
  writer.write_int("bases", n_bases);
  writer.write_float("mean_length", static_cast<double>(n_bases) / n_reads);
  writer.write_int("reads_sampled", n_reads_s);
  writer.write_float("q20_rate", static_cast<double>(n_q20) / n_bases_s);
  writer.write_float("q30_rate", static_cast<double>(n_q30) / n_bases_s);
  writer.write_float("uncalled_rate", static_cast<double>(n_n) / n_bases_s);
  writer.write_float("gc_content", static_cast<double>(n_g + n_c) / n_bases_s);
}

void
write_report_trimming(const userconfig& config,
                      json_writer& writer,
                      const trimming_statistics& totals,
                      const fastq_pair_vec& adapters)
{
  if (config.run_type == ar_command::demultiplex_sequences) {
    writer.write_null("trimming_and_filtering");
    return;
  }

  WITH_SECTION(writer, "trimming_and_filtering")
  {
    writer.write_int("adapter_trimmed_reads",
                     totals.adapter_trimmed_reads.sum());
    writer.write_int("adapter_trimmed_bases",
                     totals.adapter_trimmed_bases.sum());

    if (adapters.size() == 1) {
      writer.write("adapter_sequence_1", adapters.front().first.sequence());
      writer.write("adapter_sequence_2", adapters.front().second.sequence());

      writer.write_null("adapter_sequences");
    } else {
      writer.write_null("adapter_sequence_1");
      writer.write_null("adapter_sequence_2");

      writer.start_list("adapter_sequences");

      for (size_t i = 0; i < adapters.size(); ++i) {
        const auto _section = writer.start();

        writer.write("adapter_sequence_1", adapters.at(i).first.sequence());
        writer.write("adapter_sequence_2", adapters.at(i).second.sequence());
        writer.write_int("adapter_trimmed_reads",
                         totals.adapter_trimmed_reads.get(i));
        writer.write_int("adapter_trimmed_bases",
                         totals.adapter_trimmed_bases.get(i));
      }

      writer.end_list();
    }

    writer.write_int("overlapping_reads_merged",
                     totals.overlapping_reads_merged);
    writer.write_int("terminal_bases_trimmed", totals.terminal_bases_trimmed);
    writer.write_int("low_quality_trimmed_reads",
                     totals.low_quality_trimmed_reads);
    writer.write_int("low_quality_trimmed_bases",
                     totals.low_quality_trimmed_bases);
    writer.write_int("filtered_min_length_reads",
                     totals.filtered_min_length_reads);
    writer.write_int("filtered_min_length_bases",
                     totals.filtered_min_length_bases);
    writer.write_int("filtered_max_length_reads",
                     totals.filtered_max_length_reads);
    writer.write_int("filtered_max_length_bases",
                     totals.filtered_max_length_bases);
    writer.write_int("filtered_ambiguous_reads",
                     totals.filtered_ambiguous_reads);
    writer.write_int("filtered_ambiguous_bases",
                     totals.filtered_ambiguous_bases);
  }
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
  const bool demux_only = config.run_type == ar_command::demultiplex_sequences;

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
        writer.write_int("assigned_reads",
                         total - demux.unidentified - demux.ambiguous);
        writer.write_int("ambiguous_reads", demux.ambiguous);
        writer.write_int("unassigned_reads", demux.unidentified);

        WITH_SECTION(writer, "samples")
        {
          for (size_t i = 0; i < demux.barcodes.size(); ++i) {
            writer.write_int(config.adapters.get_sample_name(i),
                             demux.barcodes.at(i));
          }
        }
      }
    } else {
      writer.write_null("demultiplexing");
    }

    {
      trimming_statistics totals;
      for (const auto& it : stats.trimming) {
        totals += it;
      }

      write_report_trimming(
        config, writer, totals, config.adapters.get_raw_adapters());
    }

    WITH_SECTION(writer, "output")
    {
      WITH_SECTION(writer, "passed")
      {
        std::vector<const fastq_statistics*> output;
        for (const auto& it : stats.trimming) {
          output.push_back(&it.read_1);
          output.push_back(&it.read_2);
          output.push_back(&it.merged);
        }

        write_report_summary_stats(writer, output);
      }

      if (config.adapters.barcode_count()) {
        WITH_SECTION(writer, "unidentified_1")
        {
          write_report_summary_stats(
            writer, { &stats.demultiplexing.unidentified_stats_1 });
        }

        if (config.paired_ended_mode) {
          WITH_SECTION(writer, "unidentified_2")
          {
            write_report_summary_stats(
              writer, { &stats.demultiplexing.unidentified_stats_2 });
          }
        } else {
          writer.write_null("unidentified_2");
        }
      } else {
        writer.write_null("unidentified_1");
        writer.write_null("unidentified_2");
      }

      if (demux_only) {
        writer.write_null("discarded");
      } else {
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
  }
}

/** Helper struct used to simplify writing of multiple io sections. */
struct io_section
{
  io_section(read_type rtype,
             const fastq_statistics& stats,
             const output_sample_files& sample_files,
             bool multiple_files = false)
    : io_section(rtype, stats, string_vec(), multiple_files)
  {
    const auto offset = sample_files.offset(rtype);
    if (offset != output_sample_files::disabled) {
      m_filenames.push_back(sample_files.filenames.at(offset));
    }
  }

  io_section(read_type rtype,
             const fastq_statistics& stats,
             const string_vec& filenames,
             bool multiple_files = false)
    : m_read_type(rtype)
    , m_stats(stats)
    , m_filenames(filenames)
    , m_multiple_files(multiple_files)
  {}

  const char* name() const
  {
    switch (m_read_type) {
      case read_type::mate_1:
        return "read1";
      case read_type::mate_2:
        return "read2";
      case read_type::merged:
        return "merged";
      case read_type::singleton_1:
      case read_type::singleton_2:
        return "singleton";
      case read_type::discarded_1:
      case read_type::discarded_2:
        return "discarded";
      case read_type::unidentified_1:
        return "unidentified_1";
      case read_type::unidentified_2:
        return "unidentified_2";
      default:
        AR_DEBUG_FAIL("unknown read type");
    }
  }

  void write_to_if(json_writer& writer, bool enabled = true)
  {
    if (!enabled) {
      writer.write_null(name());
      return;
    }

    WITH_SECTION(writer, name())
    {
      if (m_filenames.empty()) {
        writer.write_null("filenames");
      } else {
        writer.write("filenames", m_filenames);
      }

      writer.write_int("reads", m_stats.number_of_input_reads());
      writer.write_int("reads_sampled", m_stats.number_of_sampled_reads());
      writer.write("lengths", m_stats.length_dist());

      auto total_bases = m_stats.uncalled_pos();
      auto total_quality = m_stats.uncalled_quality_pos();

      WITH_SECTION(writer, "quality_curves")
      {
        for (size_t nuc_i = 0; nuc_i < 4; ++nuc_i) {
          const auto nuc = IDX_TO_ACGT(nuc_i);
          const auto& nucleotides = m_stats.nucleotides_pos(nuc);
          const auto& quality = m_stats.qualities_pos(nuc);

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
          const auto& bases = m_stats.nucleotides_pos(nuc);

          writer.write(std::string(1, nuc), bases / total_bases);
        }

        writer.write("N", m_stats.uncalled_quality_pos() / total_bases);
        writer.write(
          "GC",
          (m_stats.nucleotides_pos('G') + m_stats.nucleotides_pos('C')) /
            total_bases);
      }
    }
  }

private:
  const read_type m_read_type;
  const fastq_statistics& m_stats;
  string_vec m_filenames;
  const bool m_multiple_files;
};

void
write_report_input(const userconfig& config,
                   json_writer& writer,
                   const ar_statistics& stats)
{
  WITH_SECTION(writer, "input")
  {
    const auto mate_2_filenames =
      config.interleaved_input ? config.input_files_1 : config.input_files_2;

    io_section(read_type::mate_1, stats.input_1, config.input_files_1, true)
      .write_to_if(writer);
    io_section(read_type::mate_2, stats.input_2, mate_2_filenames, true)
      .write_to_if(writer, config.paired_ended_mode);
  }
}

void
write_report_demultiplexing(const userconfig& config,
                            json_writer& writer,
                            const ar_statistics& sample_stats)
{
  const bool demux_only = config.run_type == ar_command::demultiplex_sequences;
  const auto out_files = config.get_output_filenames();

  if (config.adapters.barcode_count()) {
    WITH_SECTION(writer, "demultiplexing")
    {
      size_t assigned_reads = 0;
      const auto& demux = sample_stats.demultiplexing;
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
            const auto& stats = sample_stats.trimming.at(i);
            const auto& files = out_files.samples.at(i);

            writer.write_int("reads", demux.barcodes.at(i));

            write_report_trimming(
              config, writer, stats, config.adapters.get_adapter_set(i));

            WITH_SECTION(writer, "output")
            {
              io_section(read_type::mate_1, stats.read_1, files)
                .write_to_if(writer);

              io_section(read_type::mate_2, stats.read_2, files)
                .write_to_if(writer, config.paired_ended_mode);

              io_section(read_type::merged, stats.merged, files)
                .write_to_if(writer, config.merge && !demux_only);

              io_section(read_type::discarded_1, stats.discarded, files)
                .write_to_if(writer, !demux_only);
            }
          }
        }
      }
    }
  } else {
    writer.write_null("demultiplexing");
  }
}

string_vec
collect_files(const output_files& files, read_type rtype)
{
  string_vec filenames;

  for (const auto& sample_files : files.samples) {
    const auto offset = sample_files.offset(rtype);
    if (offset != output_sample_files::disabled) {
      filenames.push_back(sample_files.filenames.at(offset));
    }
  }

  return filenames;
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

  const auto out_files = config.get_output_filenames();
  const auto mate_1_files = collect_files(out_files, read_type::mate_1);
  const auto mate_2_files = collect_files(out_files, read_type::mate_2);
  const auto merged_files = collect_files(out_files, read_type::merged);
  const auto discarded_files = collect_files(out_files, read_type::discarded_1);

  const bool demux_only = config.run_type == ar_command::demultiplex_sequences;

  WITH_SECTION(writer, "output")
  {
    io_section(read_type::mate_1, output_1, mate_1_files)
      .write_to_if(writer, true);
    io_section(read_type::mate_2, output_2, mate_2_files)
      .write_to_if(writer, config.paired_ended_mode);
    io_section(read_type::merged, merged, merged_files)
      .write_to_if(writer, config.merge && !demux_only);

    io_section(read_type::unidentified_1,
               stats.demultiplexing.unidentified_stats_1,
               { out_files.unidentified_1 })
      .write_to_if(writer, config.adapters.barcode_count());
    io_section(read_type::unidentified_2,
               stats.demultiplexing.unidentified_stats_2,
               { out_files.unidentified_2 })
      .write_to_if(writer,
                   config.adapters.barcode_count() && config.paired_ended_mode);

    io_section(read_type::discarded_1, discarded, discarded_files)
      .write_to_if(writer, !demux_only);
  }
}

bool
write_report(const userconfig& config,
             const ar_statistics& stats,
             const std::string& filename)
{
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
