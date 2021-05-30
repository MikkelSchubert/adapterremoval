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

#include "demultiplexing.hpp"
#include "fastq.hpp"
#include "main.hpp"
#include "strutils.hpp"
#include "trimming.hpp"
#include "userconfig.hpp"

std::ostream&
operator<<(std::ostream& stream, const fastq::ntrimmed& ntrim)
{
  stream << ntrim.first;
  if (ntrim.first != ntrim.second) {
    stream << " " << ntrim.second;
  }

  return stream;
}

void
write_trimming_settings(const userconfig& config, std::ostream& output, int nth)
{
  output << NAME << " " << VERSION << "\nTrimming of ";

  if (config.adapters.barcode_count()) {
    if (config.adapters.get_barcodes().front().second.length()) {
      output << "double-indexed ";
    } else {
      output << "single-indexed ";
    }
  }

  if (config.paired_ended_mode) {
    if (config.interleaved_input) {
      output << "interleaved ";
    }

    output << "paired-end reads\n";
  } else {
    output << "single-end reads\n";
  }

  if (config.adapters.barcode_count()) {
    output << "\n\n[Demultiplexing]"
           << "\nMaximum mismatches (total): " << config.barcode_mm;

    if (config.paired_ended_mode) {
      output << "\nMaximum mate 1 mismatches: " << config.barcode_mm_r1;
      output << "\nMaximum mate 2 mismatches: " << config.barcode_mm_r2;
    }

    output << "\n\n\n[Demultiplexing samples]"
           << "\nName\tBarcode_1\tBarcode_2\n";

    const fastq_pair_vec barcodes = config.adapters.get_barcodes();
    for (size_t idx = 0; idx < barcodes.size(); ++idx) {
      output << config.adapters.get_sample_name(idx);
      if (static_cast<int>(idx) == nth) {
        output << "*";
      }

      const fastq_pair& current = barcodes.at(idx);
      output << "\t" << current.first.sequence();

      if (current.second.length()) {
        output << "\t" << current.second.sequence() << "\n";
      } else {
        output << "\t*\n";
      }
    }
  }

  output << "\n\n[Adapter sequences]";
  if (nth == -1) {
    const fastq_pair_vec adapters = config.adapters.get_raw_adapters();
    size_t adapter_id = 0;
    for (auto it = adapters.cbegin(); it != adapters.cend();
         ++it, ++adapter_id) {
      output << "\nAdapter1[" << adapter_id + 1
             << "]: " << it->first.sequence();

      fastq adapter_2 = it->second;
      adapter_2.reverse_complement();

      output << "\nAdapter2[" << adapter_id + 1 << "]: " << adapter_2.sequence()
             << "\n";
    }
  } else {
    const string_pair_vec adapters =
      config.adapters.get_pretty_adapter_set(nth);
    size_t adapter_id = 0;
    for (auto it = adapters.cbegin(); it != adapters.cend();
         ++it, ++adapter_id) {
      output << "\nAdapter1[" << adapter_id + 1 << "]: " << it->first;
      output << "\nAdapter2[" << adapter_id + 1 << "]: " << it->second << "\n";
    }
  }

  output << "\n\n[Adapter trimming]";
  if (config.max_threads > 1 || config.deterministic) {
    output << "\nRNG seed: NA";
  } else {
    output << "\nRNG seed: " << config.seed;
  }

  output << "\nAlignment shift value: " << config.shift
         << "\nGlobal mismatch threshold: " << config.mismatch_threshold
         << "\nQuality format (input): " << config.quality_input_fmt->name()
         << "\nQuality score max (input): "
         << config.quality_input_fmt->max_score()
         << "\nQuality format (output): " << config.quality_output_fmt->name()
         << "\nQuality score max (output): "
         << config.quality_output_fmt->max_score()
         << "\nMate-number separator (input): '" << config.mate_separator << "'"
         << "\nTrimming 5p: " << config.trim_fixed_5p
         << "\nTrimming 3p: " << config.trim_fixed_3p
         << "\nTrimming Ns: " << ((config.trim_ambiguous_bases) ? "Yes" : "No")
         << "\nTrimming Phred scores <= " << config.low_quality_score << ": "
         << (config.trim_by_quality ? "Yes" : "No")
         << "\nTrimming using sliding windows: ";

  if (config.trim_window_length >= 1) {
    output << static_cast<size_t>(config.trim_window_length);
  } else if (config.trim_window_length >= 0) {
    output << config.trim_window_length;
  } else {
    output << "No";
  }

  output << "\nMinimum genomic length: " << config.min_genomic_length
         << "\nMaximum genomic length: " << config.max_genomic_length
         << "\nCollapse overlapping reads: "
         << ((config.collapse) ? "Yes" : "No") << "\nDeterministic collapse: "
         << (config.deterministic ? "Yes" : "No") << "\nConservative collapse: "
         << (config.collapse_conservatively ? "Yes" : "No")
         << "\nMinimum overlap (in case of collapse): "
         << config.min_alignment_length;

  if (!config.paired_ended_mode) {
    output << "\nMinimum adapter overlap: " << config.min_adapter_overlap;
  }
}

void
write_trimming_statistics(const userconfig& config,
                          std::ostream& settings,
                          const statistics& stats)
{
  const std::string reads_type =
    (config.paired_ended_mode ? "read pairs: " : "reads: ");
  settings << "\n\n\n[Trimming statistics]"
           << "\nTotal number of " << reads_type << stats.records
           << "\nNumber of unaligned " << reads_type << stats.unaligned_reads
           << "\nNumber of well aligned " << reads_type
           << stats.well_aligned_reads
           << "\nNumber of discarded mate 1 reads: " << stats.discard1
           << "\nNumber of singleton mate 1 reads: " << stats.keep1;

  if (config.paired_ended_mode) {
    settings << "\nNumber of discarded mate 2 reads: " << stats.discard2
             << "\nNumber of singleton mate 2 reads: " << stats.keep2;
  }

  for (size_t adapter_id = 0;
       adapter_id < stats.number_of_reads_with_adapter.size();
       ++adapter_id) {
    const size_t count = stats.number_of_reads_with_adapter.at(adapter_id);
    // Value between 0 and stats.records for SE, and 0 and 2*stats.records
    // for N PE pairs. For PE reads, mate 1 and mate 2 reads being of
    // unequal length can cause uneven numbers.
    settings << "\nNumber of reads with adapters[" << adapter_id + 1
             << "]: " << count;
  }

  if (config.collapse) {
    settings << "\nNumber of collapsed pairs: " << stats.number_of_collapsed;
  }

  settings << "\nNumber of retained reads: " << stats.total_number_of_good_reads
           << "\nNumber of retained nucleotides: "
           << stats.total_number_of_nucleotides
           << "\nAverage length of retained reads: "
           << (stats.total_number_of_good_reads
                 ? (static_cast<double>(stats.total_number_of_nucleotides) /
                    stats.total_number_of_good_reads)
                 : 0);

  settings << "\n\n\n[Length distribution]"
           << "\nLength\tMate1\t";
  if (config.paired_ended_mode) {
    settings << "Mate2\tSingleton\t";
  }

  if (config.collapse) {
    settings << "Collapsed\t";
  }

  settings << "Discarded\tAll\n";

  for (size_t length = 0; length < stats.read_lengths.size(); ++length) {
    const std::vector<size_t>& lengths = stats.read_lengths.at(length);
    const size_t total = std::accumulate(lengths.begin(), lengths.end(), 0);

    settings << length << '\t'
             << lengths.at(static_cast<size_t>(read_type::mate_1));

    if (config.paired_ended_mode) {
      settings << '\t' << lengths.at(static_cast<size_t>(read_type::mate_2))
               << '\t' << lengths.at(static_cast<size_t>(read_type::singleton));
    }

    if (config.collapse) {
      settings << '\t' << lengths.at(static_cast<size_t>(read_type::collapsed));
    }

    settings << '\t' << lengths.at(static_cast<size_t>(read_type::discarded))
             << '\t' << total << '\n';
  }

  settings.flush();
}

void
write_demultiplex_statistics(const userconfig& config,
                             std::ofstream& output,
                             const demux_statistics& stats)
{
  const size_t total = stats.total();

  output.precision(3);
  output << std::fixed << std::setw(3) << "\n\n[Demultiplexing statistics]"
         << "\nName\tBarcode_1\tBarcode_2\tHits\tFraction\n"
         << "unidentified\tNA\tNA\t" << stats.unidentified << "\t"
         << stats.unidentified / static_cast<double>(total) << "\n"
         << "ambiguous\tNA\tNA\t" << stats.ambiguous << "\t"
         << stats.ambiguous / static_cast<double>(total) << "\n";

  const fastq_pair_vec barcodes = config.adapters.get_barcodes();
  for (size_t nth = 0; nth < barcodes.size(); ++nth) {
    const fastq_pair& current = barcodes.at(nth);

    output << config.adapters.get_sample_name(nth) << "\t"
           << current.first.sequence() << "\t";
    if (current.second.length()) {
      output << current.second.sequence() << "\t";
    } else {
      output << "*\t";
    }

    output << stats.barcodes.at(nth) << "\t"
           << stats.barcodes.at(nth) / static_cast<double>(total) << "\n";
  }

  output << "*\t*\t*\t" << total << "\t" << 1.0 << std::endl;
}

void
write_demultiplex_settings(const userconfig& config,
                           std::ostream& output,
                           int nth)
{
  output << NAME << " " << VERSION << "\nDemultiplexing of ";

  if (config.adapters.barcode_count()) {
    if (config.adapters.get_barcodes().front().second.length()) {
      output << "double-indexed ";
    } else {
      output << "single-indexed ";
    }
  }

  if (config.paired_ended_mode) {
    if (config.interleaved_input) {
      output << "interleaved ";
    }

    output << "paired-end reads";
  } else {
    output << "single-end reads";
  }

  output << "\n\n\n[Demultiplexing]"
         << "\nMaximum mismatches (total): " << config.barcode_mm;

  if (config.paired_ended_mode) {
    output << "\nMaximum mate 1 mismatches: " << config.barcode_mm_r1;
    output << "\nMaximum mate 2 mismatches: " << config.barcode_mm_r2;
  }

  output << "\n\n\n[Demultiplexing samples]"
         << "\nName\tBarcode_1\tBarcode_2\n";

  const fastq_pair_vec barcodes = config.adapters.get_barcodes();
  for (size_t idx = 0; idx < barcodes.size(); ++idx) {
    output << config.adapters.get_sample_name(idx);
    if (static_cast<size_t>(nth) == idx) {
      output << "*";
    }

    const fastq_pair& current = barcodes.at(idx);
    output << "\t" << current.first.sequence();

    if (current.second.length()) {
      output << "\t" << current.second.sequence() << "\n";
    } else {
      output << "\t*\n";
    }
  }

  output << "\n\n[Adapter sequences]";
  if (nth == -1) {
    const fastq_pair_vec adapters = config.adapters.get_raw_adapters();
    size_t adapter_id = 0;
    for (fastq_pair_vec::const_iterator it = adapters.begin();
         it != adapters.end();
         ++it, ++adapter_id) {
      output << "\nAdapter1[" << adapter_id + 1
             << "]: " << it->first.sequence();

      fastq adapter_2 = it->second;
      adapter_2.reverse_complement();
      output << "\nAdapter2[" << adapter_id + 1 << "]: " << adapter_2.sequence()
             << "\n";
    }

  } else {
    const string_pair_vec adapters =
      config.adapters.get_pretty_adapter_set(nth);
    size_t adapter_id = 0;
    for (string_pair_vec::const_iterator it = adapters.begin();
         it != adapters.end();
         ++it, ++adapter_id) {
      output << "\nAdapter1[" << adapter_id + 1 << "]: " << it->first;
      output << "\nAdapter2[" << adapter_id + 1 << "]: " << it->second << "\n";
    }
  }
}
