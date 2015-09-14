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
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "alignment.h"
#include "debug.h"
#include "demultiplex.h"
#include "fastq.h"
#include "fastq_io.h"
#include "main.h"
#include "strutils.h"
#include "userconfig.h"


typedef std::auto_ptr<fastq_output_chunk> output_chunk_ptr;


void add_chunk(chunk_vec& chunks, size_t target, std::auto_ptr<fastq_output_chunk> chunk)
{
    try {
        if (chunk.get()) {
            chunks.push_back(chunk_pair(target, chunk.release()));
        }
    } catch (...) {
        for (chunk_vec::iterator it = chunks.begin(); it != chunks.end(); ++it) {
            delete it->second;
        }

        throw;
    }
}



void write_settings(const userconfig& config, std::ostream& output, int nth)
{
    output << NAME << " " << VERSION
             << "\nTrimming of ";

    if (config.adapters.barcode_count()) {
        if (config.adapters.get_barcodes().front().second.length()) {
            output << "double-indexed ";
        } else {
            output << "single-indexed ";
        }
    }

    if (config.paired_ended_mode) {
        output << "paired-end reads";
    } else {
        output << "single-end reads\n";
    }

    if (config.adapters.barcode_count()) {
        output << "\n\n\n[Demultiplexing]"
               << "\nMaximum mismatches (total): " << config.barcode_mm;
        if (config.paired_ended_mode) {
            output << "\nMaximum mate 1 mismatches: " << config.barcode_mm_r1;
            output << "\nMaximum mate 2 mismatches: " << config.barcode_mm_r2;
        }

        output << "\n\n\n[Demultiplexing samples]"
               << "\nName\tBarcode_1\tBarcode_2";

        const fastq_pair_vec barcodes = config.adapters.get_barcodes();
        for (size_t idx = 0; idx < barcodes.size(); ++idx) {
            output << "\n" << config.adapters.get_sample_name(idx);
            if (static_cast<int>(idx) == nth) {
                output << "*";
            }

            const fastq_pair& current = barcodes.at(idx);
            output << "\t" << current.first.sequence();

            if (current.second.length()) {
                output << "\t" << current.second.sequence();
            } else {
                output << "\t*";
            }
        }
    }

    output << "\n\n\n[Adapter sequences]";
    if (nth == -1) {
        const fastq_pair_vec adapters = config.adapters.get_raw_adapters();
        size_t adapter_id = 0;
        for (fastq_pair_vec::const_iterator it = adapters.begin(); it != adapters.end(); ++it, ++adapter_id) {
            output << "\nAdapter1[" << adapter_id << "]: " << it->first.sequence();
            if (config.paired_ended_mode) {
                output << "\nAdapter2[" << adapter_id << "]: " << it->second.sequence() << "\n";
            }
        }
    } else {
        const string_pair_vec adapters = config.adapters.get_pretty_adapter_set(nth);
        size_t adapter_id = 0;
        for (string_pair_vec::const_iterator it = adapters.begin(); it != adapters.end(); ++it, ++adapter_id) {
            output << "\nAdapter1[" << adapter_id << "]: " << it->first;
            if (config.paired_ended_mode) {
                output << "\nAdapter2[" << adapter_id << "]: " << it->second << "\n";
            }
        }
    }

    output << "\n\n[Adapter trimming]"
           << "\nRNG seed: " << config.seed
           << "\nAlignment shift value: " << config.shift
           << "\nGlobal mismatch threshold: " << config.mismatch_threshold
           << "\nQuality format (input): " << config.quality_input_fmt->name()
           << "\nQuality score max (input): " << config.quality_input_fmt->max_score()
           << "\nQuality format (output): " << config.quality_output_fmt->name()
           << "\nQuality score max (output): " << config.quality_output_fmt->max_score()
           << "\nTrimming Ns: " << ((config.trim_ambiguous_bases) ? "Yes" : "No")
           << "\nTrimming Phred scores <= " << config.low_quality_score
           << ": " << (config.trim_by_quality ? "yes" : "no")
           << "\nMinimum genomic length: " << config.min_genomic_length
           << "\nMaximum genomic length: " << config.max_genomic_length
           << "\nCollapse overlapping reads: " << ((config.collapse) ? "Yes" : "No")
           << "\nMinimum overlap (in case of collapse): " << config.min_alignment_length;
}



void write_trimming_settings(const userconfig& config,
                             const statistics& stats,
                             size_t nth,
                             std::ostream& settings)
{
    write_settings(config, settings, nth);

    const std::string reads_type = (config.paired_ended_mode ? "read pairs: " : "reads: ");
    settings << "\n\n\n[Trimming statistics]"
             << "\nTotal number of " << reads_type << stats.records
             << "\nNumber of unaligned " << reads_type << stats.unaligned_reads
             << "\nNumber of well aligned " << reads_type << stats.well_aligned_reads
             << "\nNumber of inadequate alignments: " << stats.poorly_aligned_reads
             << "\nNumber of discarded mate 1 reads: " << stats.discard1
             << "\nNumber of singleton mate 1 reads: " << stats.keep1;

    if (config.paired_ended_mode) {
        settings << "\nNumber of discarded mate 2 reads: " << stats.discard2
                 << "\nNumber of singleton mate 2 reads: " << stats.keep2;
    }

    for (size_t adapter_id = 0; adapter_id < stats.number_of_reads_with_adapter.size(); ++adapter_id) {
        const size_t count = stats.number_of_reads_with_adapter.at(adapter_id);
        settings << "\nNumber of reads with adapters[" << adapter_id << "]: " << count;
    }

    if (config.collapse) {
        settings << "\nNumber of full-length collapsed pairs: " << stats.number_of_full_length_collapsed
                 << "\nNumber of truncated collapsed pairs: " << stats.number_of_truncated_collapsed;
    }

    settings << "\nNumber of retained reads: " << stats.total_number_of_good_reads
             << "\nNumber of retained nucleotides: " << stats.total_number_of_nucleotides
             << "\nAverage read length of trimmed reads: "
             << (stats.total_number_of_good_reads ? ( static_cast<double>(stats.total_number_of_nucleotides) / stats.total_number_of_good_reads) : 0);

    settings << "\n\n\n[Length distribution]"
             << "\nLength\tMate1\t";
    if (config.paired_ended_mode) {
        settings << "Mate2\tSingleton\t";
    }

    if (config.collapse) {
        settings << "Collapsed\tCollapsedTruncated\t";
    }

    settings << "Discarded\tAll\n";

    for (size_t length = 0; length < stats.read_lengths.size(); ++length) {
        const std::vector<size_t>& lengths = stats.read_lengths.at(length);
        const size_t total = std::accumulate(lengths.begin(), lengths.end(), 0);

        settings << length << '\t' << lengths.at(rt_mate_1);

        if (config.paired_ended_mode) {
            settings << '\t' << lengths.at(rt_mate_2)
                     << '\t' << lengths.at(rt_singleton);
        }

        if (config.collapse) {
            settings << '\t' << lengths.at(rt_collapsed)
                     << '\t' << lengths.at(rt_collapsed_truncated);
        }

        settings << '\t' << lengths.at(rt_discarded)
                 << '\t' << total << '\n';
    }

    settings.flush();
}



bool write_demux_settings(const userconfig& config,
                          const demultiplex_reads* step)
{
    if (!step) {
        // Demultiplexing not enabled; nothing to do
        return true;
    }

    const demux_statistics stats = step->statistics();
    const std::string filename = config.get_output_filename("demux_stats");

    try {
        std::ofstream output(filename.c_str(), std::ofstream::out);
        if (!output.is_open()) {
            std::string message = std::string("Failed to open file '") + filename + "': ";
            throw std::ofstream::failure(message + std::strerror(errno));
        }

        output.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        write_settings(config, output, -1);

        const size_t total = stats.total();

        output.precision(3);
        output << std::fixed << std::setw(3)
               << "\n\n\n[Demultiplexing statistics]"
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
                   << stats.barcodes.at(nth) / static_cast<double>(total)
                   << "\n";
        }

        output << "*\t*\t*\t" << total << "\t" << 1.0 << std::endl;
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error writing demultiplexing statistics; aborting:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;
        return false;
    }

    return true;
}


void process_collapsed_read(const userconfig& config, statistics& stats,
                            fastq& collapsed_read,
                            fastq_output_chunk& out_collapsed,
                            fastq_output_chunk& out_collapsed_truncated,
                            fastq_output_chunk& out_discarded)
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

    const size_t read_count = config.paired_ended_mode ? 2 : 1;
    if (config.is_acceptable_read(collapsed_read)) {
        stats.total_number_of_nucleotides += collapsed_read.length();
        stats.total_number_of_good_reads++;
        stats.inc_length_count(was_trimmed ? rt_collapsed_truncated : rt_collapsed,
                               collapsed_read.length());

        if (was_trimmed) {
            out_collapsed_truncated.add(*config.quality_output_fmt, collapsed_read, read_count);
        } else {
            out_collapsed.add(*config.quality_output_fmt, collapsed_read, read_count);
        }
    } else {
        stats.discard1++;
        stats.discard2++;
        stats.inc_length_count(rt_discarded, collapsed_read.length());
        out_discarded.add(*config.quality_output_fmt, collapsed_read, read_count);
    }
}


class reads_processor : public analytical_step
{
public:
    reads_processor(const userconfig& config, size_t nth)
      : analytical_step(analytical_step::unordered)
      , m_config(config)
      , m_adapters(config.adapters.get_adapter_set(nth))
      , m_stats(config)
      , m_nth(nth)
    {

    }

    statistics* get_final_statistics() {
        return m_stats.finalize();
    }

protected:
    class stats_sink : public statistics_sink<statistics>
    {
    public:
        stats_sink(const userconfig& config)
          : m_config(config)
        {
        }

    protected:
        virtual statistics* new_sink() const {
            return m_config.create_stats().release();
        }

        const userconfig& m_config;
    };

    const userconfig& m_config;
    const fastq_pair_vec m_adapters;
    stats_sink m_stats;
    const size_t m_nth;
};



class se_reads_processor : public reads_processor
{
public:
    se_reads_processor(const userconfig& config, size_t nth = 0)
      : reads_processor(config, nth)
    {
    }

    chunk_vec process(analytical_chunk* chunk)
    {
        std::auto_ptr<fastq_read_chunk> read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));

        std::auto_ptr<statistics> stats(m_stats.get_sink());

        const fastq_encoding& encoding = *m_config.quality_output_fmt;
        output_chunk_ptr out_mate_1(new fastq_output_chunk(read_chunk->eof));
        output_chunk_ptr out_collapsed;
        output_chunk_ptr out_collapsed_truncated;
        output_chunk_ptr out_discarded(new fastq_output_chunk(read_chunk->eof));

        if (m_config.collapse) {
            out_collapsed.reset(new fastq_output_chunk(read_chunk->eof));
            out_collapsed_truncated.reset(new fastq_output_chunk(read_chunk->eof));
        }

        for (fastq_vec::iterator it = read_chunk->reads_1.begin(); it != read_chunk->reads_1.end(); ++it) {
            fastq& read = *it;

            const alignment_info alignment = align_single_ended_sequence(read, m_adapters, m_config.shift);
            const userconfig::alignment_type aln_type = m_config.evaluate_alignment(alignment);

            if (aln_type == userconfig::valid_alignment) {
                truncate_single_ended_sequence(alignment, read);
                stats->number_of_reads_with_adapter.at(alignment.adapter_id)++;
                stats->well_aligned_reads++;

                if (m_config.is_alignment_collapsible(alignment)) {
                    process_collapsed_read(m_config, *stats, read,
                                           *out_collapsed,
                                           *out_collapsed_truncated,
                                           *out_discarded);
                    continue;
                }
            } else if (aln_type == userconfig::poor_alignment) {
                stats->poorly_aligned_reads++;
            } else {
                stats->unaligned_reads++;
            }

            m_config.trim_sequence_by_quality_if_enabled(read);
            if (m_config.is_acceptable_read(read)) {
                stats->keep1++;
                stats->total_number_of_good_reads++;
                stats->total_number_of_nucleotides += read.length();

                out_mate_1->add(encoding, read);
                stats->inc_length_count(rt_mate_1, read.length());
            } else {
                stats->discard1++;
                stats->inc_length_count(rt_discarded, read.length());

                out_discarded->add(encoding, read);
            }
        }

        stats->records += read_chunk->reads_1.size();
        m_stats.return_sink(stats.release());

        chunk_vec chunks;
        const size_t offset = m_nth * ai_analyses_offset;
        add_chunk(chunks, offset + ai_write_mate_1, out_mate_1);
        add_chunk(chunks, offset + ai_write_collapsed, out_collapsed);
        add_chunk(chunks, offset + ai_write_collapsed_truncated, out_collapsed_truncated);
        add_chunk(chunks, offset + ai_write_discarded, out_discarded);

        return chunks;
    }
};


class pe_reads_processor : public reads_processor
{
public:
    pe_reads_processor(const userconfig& config, size_t nth)
      : reads_processor(config, nth)
    {
    }

    chunk_vec process(analytical_chunk* chunk)
    {
        std::auto_ptr<fastq_read_chunk> read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));

        std::auto_ptr<statistics> stats(m_stats.get_sink());

        const fastq_encoding& encoding = *m_config.quality_output_fmt;
        output_chunk_ptr out_mate_1(new fastq_output_chunk(read_chunk->eof));
        output_chunk_ptr out_mate_2(new fastq_output_chunk(read_chunk->eof));
        output_chunk_ptr out_singleton(new fastq_output_chunk(read_chunk->eof));
        output_chunk_ptr out_collapsed;
        output_chunk_ptr out_collapsed_truncated;
        output_chunk_ptr out_discarded(new fastq_output_chunk(read_chunk->eof));

        if (m_config.collapse) {
            out_collapsed.reset(new fastq_output_chunk(read_chunk->eof));
            out_collapsed_truncated.reset(new fastq_output_chunk(read_chunk->eof));
        }

        AR_DEBUG_ASSERT(read_chunk->reads_1.size() == read_chunk->reads_2.size());

        fastq_vec::iterator it_1 = read_chunk->reads_1.begin();
        fastq_vec::iterator it_2 = read_chunk->reads_2.begin();
        while (it_1 != read_chunk->reads_1.end()) {
            fastq read1 = *it_1++;
            fastq read2 = *it_2++;

            // Throws if read-names or mate numbering does not match
            fastq::validate_paired_reads(read1, read2);

            // Reverse complement to match the orientation of read1
            read2.reverse_complement();

            const alignment_info alignment = align_paired_ended_sequences(read1, read2, m_adapters, m_config.shift);
            const userconfig::alignment_type aln_type = m_config.evaluate_alignment(alignment);
            if (aln_type == userconfig::valid_alignment) {
                stats->well_aligned_reads++;
                const size_t n_adapters = truncate_paired_ended_sequences(alignment, read1, read2);
                stats->number_of_reads_with_adapter.at(alignment.adapter_id) += n_adapters;

                if (m_config.is_alignment_collapsible(alignment)) {
                    fastq collapsed_read = collapse_paired_ended_sequences(alignment, read1, read2);
                    process_collapsed_read(m_config, *stats, collapsed_read,
                                           *out_collapsed,
                                           *out_collapsed_truncated,
                                           *out_discarded);
                    continue;
                }
            } else if (aln_type == userconfig::poor_alignment) {
                stats->poorly_aligned_reads++;
            } else {
                stats->unaligned_reads++;
            }

            // Reads were not aligned or collapsing is not enabled
            // Undo reverse complementation (post truncation of adapters)
            read2.reverse_complement();

            // Are the reads good enough? Not too many Ns?
            m_config.trim_sequence_by_quality_if_enabled(read1);
            m_config.trim_sequence_by_quality_if_enabled(read2);
            const bool read_1_acceptable = m_config.is_acceptable_read(read1);
            const bool read_2_acceptable = m_config.is_acceptable_read(read2);

            stats->total_number_of_nucleotides += read_1_acceptable ? read1.length() : 0u;
            stats->total_number_of_nucleotides += read_1_acceptable ? read2.length() : 0u;
            stats->total_number_of_good_reads += read_1_acceptable;
            stats->total_number_of_good_reads += read_2_acceptable;

            if (read_1_acceptable && read_2_acceptable) {
                out_mate_1->add(encoding, read1);
                out_mate_2->add(encoding, read2);

                stats->inc_length_count(rt_mate_1, read1.length());
                stats->inc_length_count(rt_mate_2, read2.length());
            } else {
                // Keep one or none of the reads ...
                stats->keep1 += read_1_acceptable;
                stats->keep2 += read_2_acceptable;
                stats->discard1 += !read_1_acceptable;
                stats->discard2 += !read_2_acceptable;
                stats->inc_length_count(read_1_acceptable ? rt_mate_1 : rt_discarded, read1.length());
                stats->inc_length_count(read_2_acceptable ? rt_mate_2 : rt_discarded, read2.length());

                if (read_1_acceptable) {
                    out_singleton->add(encoding, read1);
                } else {
                    out_discarded->add(encoding, read1);
                }

                if (read_2_acceptable) {
                    out_singleton->add(encoding, read2);
                } else {
                    out_discarded->add(encoding, read2);
                }
            }
        }

        stats->records += read_chunk->reads_1.size();
        m_stats.return_sink(stats.release());

        chunk_vec chunks;
        const size_t offset = m_nth * ai_analyses_offset;
        add_chunk(chunks, offset + ai_write_mate_1, out_mate_1);
        add_chunk(chunks, offset + ai_write_mate_2, out_mate_2);
        add_chunk(chunks, offset + ai_write_singleton, out_singleton);
        add_chunk(chunks, offset + ai_write_collapsed, out_collapsed);
        add_chunk(chunks, offset + ai_write_collapsed_truncated, out_collapsed_truncated);
        add_chunk(chunks, offset + ai_write_discarded, out_discarded);

        return chunks;
    }
};


bool write_settings(const userconfig& config, const std::vector<reads_processor*>& processors)
{
    for (size_t nth = 0; nth < processors.size(); ++nth) {
        const std::string filename = config.get_output_filename("--settings", nth);

        const std::auto_ptr<statistics> stats(processors.at(nth)->get_final_statistics());

        try {
            std::ofstream output(filename.c_str(), std::ofstream::out);

            if (!output.is_open()) {
                std::string message = std::string("Failed to open file '") + filename + "': ";
                throw std::ofstream::failure(message + std::strerror(errno));
            }

            output.exceptions(std::ofstream::failbit | std::ofstream::badbit);
            write_trimming_settings(config, *stats, nth, output);
        } catch (const std::ios_base::failure& error) {
            std::cerr << "IO error writing settings file; aborting:\n"
                      << cli_formatter::fmt(error.what()) << std::endl;
            return false;
        }
    }

    return true;
}


void add_write_step(const userconfig& config, scheduler& sch, size_t offset,
                    analytical_step* step)
{
#ifdef AR_GZIP_SUPPORT
    if (config.gzip) {
        sch.add_step(offset + ai_zip_offset, step);
        sch.add_step(offset, new gzip_paired_fastq(config, offset + ai_zip_offset));
    } else
#endif

#ifdef AR_BZIP2_SUPPORT
    if (config.bzip2) {
        sch.add_step(offset + ai_zip_offset, step);
        sch.add_step(offset, new bzip2_paired_fastq(config, offset + ai_zip_offset));
    } else
#endif
    {
        sch.add_step(offset, step);
    }
}


int remove_adapter_sequences_se(const userconfig& config)
{
    std::cout << "Trimming single ended reads ..." << std::endl;

    scheduler sch;
    std::vector<reads_processor*> processors;
    demultiplex_reads* demultiplexer = NULL;

    try {
        if (config.adapters.barcode_count()) {
            // Step 1: Read input file
            sch.add_step(ai_read_fastq, new read_single_fastq(config.quality_input_fmt.get(),
                                                              config.input_file_1,
                                                              ai_demultiplex));

            // Step 2: Parse and demultiplex reads based on single or double indices
            sch.add_step(ai_demultiplex, demultiplexer = new demultiplex_se_reads(&config));

            add_write_step(config, sch, ai_write_unidentified_1,
                           new write_paired_fastq(config.get_output_filename("demux_unknown")));
        } else {
            sch.add_step(ai_read_fastq, new read_single_fastq(config.quality_input_fmt.get(),
                                                              config.input_file_1,
                                                              ai_analyses_offset));
        }

        // Step 3 - N: Trim and write demultiplexed readss
        for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
            const size_t offset = nth * ai_analyses_offset;

            processors.push_back(new se_reads_processor(config, nth));
            sch.add_step(offset + ai_trim_se, processors.back());

            add_write_step(config, sch, offset + ai_write_mate_1,
                           new write_paired_fastq(config.get_output_filename("--output1", nth)));
            add_write_step(config, sch, offset + ai_write_discarded,
                         new write_paired_fastq(config.get_output_filename("--discarded", nth)));

            if (config.collapse) {
                add_write_step(config, sch, offset + ai_write_collapsed,
                               new write_paired_fastq(config.get_output_filename("--outputcollapsed", nth)));
                add_write_step(config, sch, offset + ai_write_collapsed_truncated,
                               new write_paired_fastq(config.get_output_filename("--outputcollapsedtruncated", nth)));
            }
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;
        return 1;
    }

    if (!sch.run(config.max_threads, config.seed)) {
        return 1;
    } else if (!write_settings(config, processors)) {
        return 1;
    } else if (!write_demux_settings(config, demultiplexer)) {
        return 1;
    }

    return 0;
}


int remove_adapter_sequences_pe(const userconfig& config)
{
    std::cout << "Trimming paired end reads ..." << std::endl;

    scheduler sch;
    std::vector<reads_processor*> processors;
    demultiplex_reads* demultiplexer = NULL;

    try {
        if (config.adapters.barcode_count()) {
            // Step 1: Read input file
            sch.add_step(ai_read_fastq, new read_paired_fastq(config.quality_input_fmt.get(),
                                                              config.input_file_1,
                                                              config.input_file_2,
                                                              ai_demultiplex));

            // Step 2: Parse and demultiplex reads based on single or double indices
            sch.add_step(ai_demultiplex, demultiplexer = new demultiplex_pe_reads(&config));

            add_write_step(config, sch, ai_write_unidentified_1,
                           new write_paired_fastq(config.get_output_filename("demux_unknown", 1)));
            add_write_step(config, sch, ai_write_unidentified_2,
                           new write_paired_fastq(config.get_output_filename("demux_unknown", 2)));
        } else {
            sch.add_step(ai_read_fastq, new read_paired_fastq(config.quality_input_fmt.get(),
                                                              config.input_file_1,
                                                              config.input_file_2,
                                                              ai_analyses_offset));
        }

        // Step 3 - N: Trim and write demultiplexed reads
        for (size_t nth = 0; nth < config.adapters.adapter_set_count(); ++nth) {
            const size_t offset = nth * ai_analyses_offset;

            processors.push_back(new pe_reads_processor(config, nth));
            sch.add_step(offset + ai_trim_pe, processors.back());

            add_write_step(config, sch, offset + ai_write_mate_1,
                           new write_paired_fastq(config.get_output_filename("--output1", nth)));
            add_write_step(config, sch, offset + ai_write_mate_2,
                           new write_paired_fastq(config.get_output_filename("--output2", nth)));
            add_write_step(config, sch, offset + ai_write_discarded,
                           new write_paired_fastq(config.get_output_filename("--discarded", nth)));
            add_write_step(config, sch, offset + ai_write_singleton,
                           new write_paired_fastq(config.get_output_filename("--singleton", nth)));

            if (config.collapse) {
                add_write_step(config, sch, offset + ai_write_collapsed,
                               new write_paired_fastq(config.get_output_filename("--outputcollapsed", nth)));
                add_write_step(config, sch, offset + ai_write_collapsed_truncated,
                               new write_paired_fastq(config.get_output_filename("--outputcollapsedtruncated", nth)));
            }
        }
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;
        return 1;
    }

    if (!sch.run(config.max_threads, config.seed)) {
        return 1;
    } else if (!write_settings(config, processors)) {
        return 1;
    } else if (!write_demux_settings(config, demultiplexer)) {
        return 1;
    }

    return 0;
}


int remove_adapter_sequences(const userconfig& config)
{
    if (config.paired_ended_mode) {
        return remove_adapter_sequences_pe(config);
    } else {
        return remove_adapter_sequences_se(config);
    }

    return 0;
}
