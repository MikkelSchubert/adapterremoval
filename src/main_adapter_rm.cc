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
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cerrno>
#include <cstring>
#include <stdexcept>

#include "main.h"
#include "main_adapter_rm.h"
#include "fastq.h"
#include "alignment.h"
#include "userconfig.h"
#include "timer.h"


std::string describe_phred_format(quality_format fmt)
{
    switch (fmt) {
        case phred_33: return "Phred+33";
        case phred_64: return "Phred+64";
        case solexa: return "Solexa";
        default: throw std::invalid_argument("invalid quality score format");
    }

}


std::ostream& write_settings(const userconfig& config,
                             std::ostream& settings)
{
    settings << "Running " << NAME << " " << VERSION << " using the following options:\n";
    settings << "RNG seed: " << config.seed << "\n";

    if (config.paired_ended_mode) {
        settings << "Paired end mode\n";
    } else {
        settings << "Single end mode\n";
    }

    size_t adapter_id = 0;
    for (fastq_pair_vec::const_iterator it = config.adapters.begin(); it != config.adapters.end(); ++it, adapter_id++) {
        settings << "PCR1[" << adapter_id << "]: " << it->first.sequence() << "\n";
        if (config.paired_ended_mode) {
            settings << "PCR2[" << adapter_id << "]: " << it->second.sequence() << "\n";
        }
    }

    if (config.trim_barcodes_mode) {
        adapter_id = 0;
        for (fastq_pair_vec::const_iterator it = config.barcodes.begin(); it != config.barcodes.end(); ++it, adapter_id++) {
            settings << "Mate 1 5' barcode[" << adapter_id << "]: " << it->first.sequence() << "\n";
        }
    }

    settings << "Alignment shift value: " << config.shift << "\n";
    settings << "Global mismatch threshold: " << config.mismatch_threshold << "\n";
    settings << "Quality format (input): " << describe_phred_format(config.quality_input_fmt) << "\n";
    settings << "Quality format (output): " << describe_phred_format(config.quality_output_fmt) << "\n";
    settings << "Trimming Ns: " << ((config.trim_ambiguous_bases) ? "Yes" : "No") << "\n";

    settings << "Trimming Phred scores <= " << config.low_quality_score << ": " << (config.trim_by_quality ? "yes" : "no") << "\n";
    settings << "Minimum genomic length: " << config.min_genomic_length << "\n";
    settings << "Collapse overlapping reads: " << ((config.collapse) ? "Yes" : "No") << "\n";
    settings << "Minimum overlap (in case of collapse): " << config.min_alignment_length << "\n";

    settings.flush();

    return settings;
}


std::ostream& write_statistics(const userconfig& config, std::ostream& settings, const statistics& stats)
{
    const std::string reads_type = (config.paired_ended_mode ? "read pairs: " : "reads: ");

    settings << "\n";
    settings << "Total number of " << reads_type << stats.records << "\n";
    settings << "Number of unaligned " << reads_type << stats.unaligned_reads << "\n";
    settings << "Number of well aligned " << reads_type << stats.well_aligned_reads << "\n";
    settings << "Number of inadequate alignments: " << stats.poorly_aligned_reads << "\n";
    settings << "Number of discarded mate 1 reads: " << stats.discard1 << "\n";
    settings << "Number of singleton mate 1 reads: " << stats.keep1 << "\n";

    if (config.paired_ended_mode) {
        settings << "Number of discarded mate 2 reads: " << stats.discard2 << "\n";
        settings << "Number of singleton mate 2 reads: " << stats.keep2 << "\n";
    }

    settings << "\n";
    if (config.trim_barcodes_mode) {
        for (size_t barcode_id = 0; barcode_id < stats.number_of_barcodes_trimmed.size(); ++barcode_id) {
            const size_t count = stats.number_of_barcodes_trimmed.at(barcode_id);
            settings << "Number of reads with barcode[" << barcode_id << "]: " << count << "\n";
        }
    }

    for (size_t adapter_id = 0; adapter_id < stats.number_of_reads_with_adapter.size(); ++adapter_id) {
        const size_t count = stats.number_of_reads_with_adapter.at(adapter_id);
        settings << "Number of reads with adapters[" << adapter_id << "]: " << count << "\n";
    }

    if (config.collapse) {
        settings << "Number of full-length collapsed pairs: " << stats.number_of_full_length_collapsed << "\n";
        settings << "Number of truncated collapsed pairs: " << stats.number_of_truncated_collapsed << "\n";
    }

    settings << "Number of retained reads: " << stats.total_number_of_good_reads << "\n";
    settings << "Number of retained nucleotides: " << stats.total_number_of_nucleotides << "\n";
    settings << "Average read length of trimmed reads: "
             << (stats.total_number_of_good_reads ? ( static_cast<double>(stats.total_number_of_nucleotides) / stats.total_number_of_good_reads) : 0)
             << "\n";

    settings.flush();

    return settings;
}


void process_collapsed_read(const userconfig& config, statistics& stats,
                            fastq& collapsed_read,
                            std::ostream& io_discarded,
                            std::ostream& io_collapsed,
                            std::ostream& io_collapsed_truncated)
{
    const fastq::ntrimmed trimmed = config.trim_sequence_by_quality_if_enabled(collapsed_read);

    // If trimmed, the external coordinates are no longer reliable
    // for determining the size of the original template.
    const bool was_trimmed = trimmed.first || trimmed.second;
    if (was_trimmed) {
        collapsed_read.add_prefix_to_header("MT_");
        stats.number_of_truncated_collapsed++;
    } else {
        collapsed_read.add_prefix_to_header("M_");
        stats.number_of_full_length_collapsed++;
    }

    if (config.is_acceptable_read(collapsed_read)) {
        stats.total_number_of_nucleotides += collapsed_read.length();
        stats.total_number_of_good_reads++;
        collapsed_read.write((was_trimmed ? io_collapsed_truncated : io_collapsed), config.quality_output_fmt);
    } else {
        stats.discard1++;
        stats.discard2++;
        collapsed_read.write(io_discarded, config.quality_output_fmt);
    }
}


bool process_single_ended_reads(const userconfig& config, statistics& stats)
{
    std::auto_ptr<std::istream> io_input;
    std::auto_ptr<std::ostream> io_output;
    std::auto_ptr<std::ostream> io_discarded;
    std::auto_ptr<std::ostream> io_collapsed;
    std::auto_ptr<std::ostream> io_collapsed_truncated;

    try {
        io_input = config.open_ifstream(config.input_file_1);

        io_discarded = config.open_with_default_filename("--discarded", ".discarded");
        io_output = config.open_with_default_filename("--output1", ".truncated");

        if (config.collapse) {
            io_collapsed = config.open_with_default_filename("--outputcollapsed", ".collapsed");
            io_collapsed_truncated = config.open_with_default_filename("--outputcollapsedtruncated", ".collapsed.truncated");
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n    " << error.what() << std::endl;
        return false;
    }

    try {
        fastq read;
        timer progress("reads", config.quiet);
        for ( ; read.read(*io_input, config.quality_input_fmt); ++stats.records) {
            progress.increment();
            config.trim_barcodes_if_enabled(read, stats);

            const alignment_info alignment = align_single_ended_sequence(read, config.adapters, config.shift);
            const userconfig::alignment_type aln_type = config.evaluate_alignment(alignment);

            if (aln_type == userconfig::valid_alignment) {
                truncate_single_ended_sequence(alignment, read);
                stats.number_of_reads_with_adapter.at(alignment.adapter_id)++;
                stats.well_aligned_reads++;

                if (config.is_alignment_collapsible(alignment)) {
                    process_collapsed_read(config, stats, read, *io_discarded, *io_collapsed, *io_collapsed_truncated);
                    continue;
                }
            } else if (aln_type == userconfig::poor_alignment) {
                stats.poorly_aligned_reads++;
            } else {
                stats.unaligned_reads++;
            }

            config.trim_sequence_by_quality_if_enabled(read);
            if (config.is_acceptable_read(read)) {
                stats.keep1++;
                stats.total_number_of_good_reads++;
                stats.total_number_of_nucleotides += read.length();

                read.write(*io_output, config.quality_output_fmt);
            } else {
                stats.discard1++;

                read.write(*io_discarded, config.quality_output_fmt);
            }
        }

        progress.finalize();
    } catch (const fastq_error& error) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << error.what() << std::endl;
        return false;
    } catch (const std::ios_base::failure&) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << std::strerror(errno) << std::endl;
        return false;
    }

    io_output->flush();
    io_discarded->flush();
    if (config.collapse) {
        io_collapsed->flush();
        io_collapsed_truncated->flush();
    }

    return true;
}


bool process_paired_ended_reads(const userconfig& config, statistics& stats)
{
    std::auto_ptr<std::istream> io_input_1;
    std::auto_ptr<std::istream> io_input_2;

    std::auto_ptr<std::ostream> io_output_1;
    std::auto_ptr<std::ostream> io_output_2;
    std::auto_ptr<std::ostream> io_singleton;
    std::auto_ptr<std::ostream> io_collapsed;
    std::auto_ptr<std::ostream> io_collapsed_truncated;
    std::auto_ptr<std::ostream> io_discarded;

    try {
        io_input_1 = config.open_ifstream(config.input_file_1);
        io_input_2 = config.open_ifstream(config.input_file_2);

        io_discarded = config.open_with_default_filename("--discarded", ".discarded");

        io_output_1 = config.open_with_default_filename("--output1", ".pair1.truncated");
        io_output_2 = config.open_with_default_filename("--output2", ".pair2.truncated");

        io_singleton = config.open_with_default_filename("--singleton", ".singleton.truncated");

        if (config.collapse) {
            io_collapsed = config.open_with_default_filename("--outputcollapsed", ".collapsed");
            io_collapsed_truncated = config.open_with_default_filename("--outputcollapsedtruncated", ".collapsed.truncated");
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n    " << error.what() << std::endl;
        return false;
    }

    fastq read1;
    fastq read2;
    try {
        timer progress("pairs", config.quiet);
        for (; ; ++stats.records) {
            const bool read_file_1_ok = read1.read(*io_input_1, config.quality_input_fmt);
            const bool read_file_2_ok = read2.read(*io_input_2, config.quality_input_fmt);

            if (read_file_1_ok != read_file_2_ok) {
                throw fastq_error("files contain unequal number of records");
            } else if (!read_file_1_ok) {
                break;
            }

            progress.increment();

            // Throws if read-names or mate numbering does not match
            fastq::validate_paired_reads(read1, read2);

            config.trim_barcodes_if_enabled(read1, stats);

            // Reverse complement to match the orientation of read1
            read2.reverse_complement();

            const alignment_info alignment = align_paired_ended_sequences(read1, read2, config.adapters, config.shift);
            const userconfig::alignment_type aln_type = config.evaluate_alignment(alignment);
            if (aln_type == userconfig::valid_alignment) {
                stats.well_aligned_reads++;
                const size_t n_adapters = truncate_paired_ended_sequences(alignment, read1, read2);
                stats.number_of_reads_with_adapter.at(alignment.adapter_id) += n_adapters;

                if (config.is_alignment_collapsible(alignment)) {
                    fastq collapsed_read = collapse_paired_ended_sequences(alignment, read1, read2);
                    process_collapsed_read(config, stats, collapsed_read, *io_discarded, *io_collapsed, *io_collapsed_truncated);
                    continue;
                }
            } else if (aln_type == userconfig::poor_alignment) {
                stats.poorly_aligned_reads++;
            } else {
                stats.unaligned_reads++;
            }

            // Reads were not aligned or collapsing is not enabled
            // Undo reverse complementation (post truncation of adapters)
            read2.reverse_complement();

            // Are the reads good enough? Not too many Ns?
            config.trim_sequence_by_quality_if_enabled(read1);
            config.trim_sequence_by_quality_if_enabled(read2);
            const bool read_1_acceptable = config.is_acceptable_read(read1);
            const bool read_2_acceptable = config.is_acceptable_read(read2);

            stats.total_number_of_nucleotides += read_1_acceptable ? read1.length() : 0u;
            stats.total_number_of_nucleotides += read_1_acceptable ? read2.length() : 0u;
            stats.total_number_of_good_reads += read_1_acceptable;
            stats.total_number_of_good_reads += read_2_acceptable;

            if (read_1_acceptable && read_2_acceptable) {
                read1.write(*io_output_1, config.quality_output_fmt);
                read2.write(*io_output_2, config.quality_output_fmt);
            } else {
                // Keep one or none of the reads ...
                stats.keep1 += read_1_acceptable;
                stats.keep2 += read_2_acceptable;
                stats.discard1 += !read_1_acceptable;
                stats.discard2 += !read_2_acceptable;
                read1.write((read_1_acceptable ? *io_singleton : *io_discarded), config.quality_output_fmt);
                read2.write((read_2_acceptable ? *io_singleton : *io_discarded), config.quality_output_fmt);
            }
        }

        progress.finalize();
    } catch (const fastq_error& error) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << error.what() << std::endl;
        return false;
    } catch (const std::ios_base::failure&) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << std::strerror(errno) << std::endl;
        return false;
    }

    io_output_1->flush();
    io_output_2->flush();
    io_singleton->flush();
    io_discarded->flush();
    if (config.collapse) {
        io_collapsed->flush();
        io_collapsed_truncated->flush();
    }

    return true;
}



int remove_adapter_sequences(const userconfig& config)
{
    std::auto_ptr<std::ostream> settings;
    try {
        settings = config.open_with_default_filename("--settings",
                                                     ".settings",
                                                     false);
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening settings file; aborting:\n    "
                  << error.what() << std::endl;
        return 1;
    }

    if (!write_settings(config, *settings)) {
        std::cerr << "Error writing settings file; aborting!" << std::endl;
        return 1;
    }

    statistics stats = config.create_stats();
    if (config.paired_ended_mode) {
        if (!process_paired_ended_reads(config, stats)) {
            return 1;
        }
    } else if (!process_single_ended_reads(config, stats)) {
        return 1;
    }

    if (!write_statistics(config, *settings, stats)) {
        std::cerr << "Error writing statistics to settings file!" << std::endl;
        return 1;
    }

    settings->flush();

    return 0;
}
