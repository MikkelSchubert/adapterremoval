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
#include <numeric>

#include "main.h"
#include "main_adapter_rm.h"
#include "fastq.h"
#include "fastq_io.h"
#include "alignment.h"
#include "userconfig.h"


typedef std::auto_ptr<fastq_output_chunk> output_chunk_ptr;


void add_chunk(chunk_list& chunks, size_t target, std::auto_ptr<fastq_output_chunk> chunk)
{
    try {
        if (chunk.get()) {
            chunks.push_back(chunk_pair(target, chunk.release()));
        }
    } catch (...) {
        for (chunk_list::iterator it = chunks.begin(); it != chunks.end(); ++it) {
            delete it->second;
        }

        throw;
    }
}


std::string describe_phred_format(const fastq::quality_format fmt)
{
    switch (fmt) {
        case fastq::phred_33: return "Phred+33";
        case fastq::phred_64: return "Phred+64";
        case fastq::solexa: return "Solexa";
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
    for (fastq_pair_vec::const_iterator it = config.adapters.begin(); it != config.adapters.end(); ++it, ++adapter_id) {
        settings << "Adapter1[" << adapter_id << "]: " << it->first.sequence() << "\n";
        if (config.paired_ended_mode) {
            fastq adapter = it->second;
            adapter.reverse_complement();

            settings << "Adapter2[" << adapter_id << "]: " << adapter.sequence() << "\n";
        }
    }

    adapter_id = 0;
    for (fastq_pair_vec::const_iterator it = config.barcodes.begin(); it != config.barcodes.end(); ++it, ++adapter_id) {
        settings << "Mate 1 5' barcode[" << adapter_id << "]: " << it->first.sequence() << "\n";
    }

    settings << "Alignment shift value: " << config.shift << "\n";
    settings << "Global mismatch threshold: " << config.mismatch_threshold << "\n";
    settings << "Quality format (input): " << describe_phred_format(config.quality_input_fmt) << "\n";
    settings << "Quality format (output): " << describe_phred_format(config.quality_output_fmt) << "\n";
    settings << "Trimming Ns: " << ((config.trim_ambiguous_bases) ? "Yes" : "No") << "\n";

    settings << "Trimming Phred scores <= " << config.low_quality_score << ": " << (config.trim_by_quality ? "yes" : "no") << "\n";
    settings << "Minimum genomic length: " << config.min_genomic_length << "\n";
    settings << "Maximum genomic length: " << config.max_genomic_length << "\n";
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
    for (size_t barcode_id = 0; barcode_id < stats.number_of_barcodes_trimmed.size(); ++barcode_id) {
        const size_t count = stats.number_of_barcodes_trimmed.at(barcode_id);
        settings << "Number of reads with barcode[" << barcode_id << "]: " << count << "\n";
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
             << "\n\n";

    const std::string prefix = "Length distribution: ";
    settings << prefix << "Length\tMate1\tMate2\tSingleton\tCollapsed\tCollapsedTruncated\tDiscarded\tAll\n";
    for (size_t length = 0; length < stats.read_lengths.size(); ++length) {
        const std::vector<size_t>& lengths = stats.read_lengths.at(length);
        const size_t total = std::accumulate(lengths.begin(), lengths.end(), 0);

        settings << prefix << length
                 << '\t' << lengths.at(rt_mate_1)
                 << '\t' << lengths.at(rt_mate_2)
                 << '\t' << lengths.at(rt_singleton)
                 << '\t' << lengths.at(rt_collapsed)
                 << '\t' << lengths.at(rt_collapsed_truncated)
                 << '\t' << lengths.at(rt_discarded)
                 << '\t' << total << '\n';
    }

    settings.flush();

    return settings;
}


void process_collapsed_read(const userconfig& config, statistics& stats,
                            fastq& collapsed_read,
                            string_vec& out_collapsed,
                            string_vec& out_collapsed_truncated,
                            string_vec& out_discarded)
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
        stats.inc_length_count(was_trimmed ? rt_collapsed_truncated : rt_collapsed,
                               collapsed_read.length());

        const std::string& line = collapsed_read.to_str(config.quality_output_fmt);
        if (was_trimmed) {
            out_collapsed_truncated.push_back(line);
        } else {
            out_collapsed.push_back(line);
        }
    } else {
        stats.discard1++;
        stats.discard2++;
        stats.inc_length_count(rt_discarded, collapsed_read.length());
        out_discarded.push_back(collapsed_read.to_str(config.quality_output_fmt));
    }
}



class reads_processor : public analytical_step
{
public:
    reads_processor(const userconfig& config)
      : analytical_step(analytical_step::unordered)
      , m_config(config)
      , m_stats(config)
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
    stats_sink m_stats;
};



class se_reads_processor : public reads_processor
{
public:
    se_reads_processor(const userconfig& config)
      : reads_processor(config)
    {
    }

    chunk_list process(analytical_chunk* chunk)
    {
        std::auto_ptr<fastq_file_chunk> file_chunk(dynamic_cast<fastq_file_chunk*>(chunk));
        string_vec_citer file_1_it = file_chunk->mates.at(rt_mate_1).begin();
        const string_vec_citer file_1_end = file_chunk->mates.at(rt_mate_1).end();

        std::auto_ptr<statistics> stats(m_stats.get_sink());

        output_chunk_ptr out_mate_1(new fastq_output_chunk(file_chunk->eof));
        output_chunk_ptr out_collapsed;
        output_chunk_ptr out_collapsed_truncated;
        output_chunk_ptr out_discarded(new fastq_output_chunk(file_chunk->eof));

        if (m_config.collapse) {
            out_collapsed.reset(new fastq_output_chunk(file_chunk->eof));
            out_collapsed_truncated.reset(new fastq_output_chunk(file_chunk->eof));
        }

        try {
            for (fastq read; read.read(file_1_it, file_1_end, m_config.quality_input_fmt); file_chunk->offset += 4) {
                m_config.trim_barcodes_if_enabled(read, *stats);

                const alignment_info alignment = align_single_ended_sequence(read, m_config.adapters, m_config.shift);
                const userconfig::alignment_type aln_type = m_config.evaluate_alignment(alignment);

                if (aln_type == userconfig::valid_alignment) {
                    truncate_single_ended_sequence(alignment, read);
                    stats->number_of_reads_with_adapter.at(alignment.adapter_id)++;
                    stats->well_aligned_reads++;

                    if (m_config.is_alignment_collapsible(alignment)) {
                        process_collapsed_read(m_config, *stats, read,
                                               out_collapsed->reads,
                                               out_collapsed_truncated->reads,
                                               out_discarded->reads);
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

                    out_mate_1->reads.push_back(read.to_str(m_config.quality_output_fmt));
                    stats->inc_length_count(rt_mate_1, read.length());
                } else {
                    stats->discard1++;
                    stats->inc_length_count(rt_discarded, read.length());

                    out_discarded->reads.push_back(read.to_str(m_config.quality_output_fmt));
                }
            }
        } catch (const fastq_error& error) {
            print_locker lock;
            std::cerr << "Error reading FASTQ record at line " << file_chunk->offset << "; aborting:\n    " << error.what() << std::endl;
            throw thread_abort();
        }

        stats->records += file_chunk->mates.at(rt_mate_1).size() / 4;
        m_stats.return_sink(stats.release());

        chunk_list chunks;
        const bool compress = m_config.gzip || m_config.bzip2;
        add_chunk(chunks, compress ? ai_zip_mate_1 : ai_write_mate_1, out_mate_1);
        add_chunk(chunks, compress ? ai_zip_collapsed : ai_write_collapsed, out_collapsed);
        add_chunk(chunks, compress ? ai_zip_collapsed_truncated : ai_write_collapsed_truncated, out_collapsed_truncated);
        add_chunk(chunks, compress ? ai_zip_discarded : ai_write_discarded, out_discarded);

        return chunks;
    }
};


class pe_reads_processor : public reads_processor
{
public:
    pe_reads_processor(const userconfig& config)
      : reads_processor(config)
    {
    }

    chunk_list process(analytical_chunk* chunk)
    {
        std::auto_ptr<fastq_file_chunk> file_chunk(dynamic_cast<fastq_file_chunk*>(chunk));

        string_vec_citer file_1_it = file_chunk->mates.at(rt_mate_1).begin();
        const string_vec_citer file_1_end = file_chunk->mates.at(rt_mate_1).end();
        string_vec_citer file_2_it = file_chunk->mates.at(rt_mate_2).begin();
        const string_vec_citer file_2_end = file_chunk->mates.at(rt_mate_2).end();

        std::auto_ptr<statistics> stats(m_stats.get_sink());

        output_chunk_ptr out_mate_1(new fastq_output_chunk(file_chunk->eof));
        output_chunk_ptr out_mate_2(new fastq_output_chunk(file_chunk->eof));
        output_chunk_ptr out_singleton(new fastq_output_chunk(file_chunk->eof));
        output_chunk_ptr out_collapsed;
        output_chunk_ptr out_collapsed_truncated;
        output_chunk_ptr out_discarded(new fastq_output_chunk(file_chunk->eof));

        if (m_config.collapse) {
            out_collapsed.reset(new fastq_output_chunk(file_chunk->eof));
            out_collapsed_truncated.reset(new fastq_output_chunk(file_chunk->eof));
        }

        fastq read1;
        fastq read2;
        try {
            for (; ; file_chunk->offset += 4) {
                const bool read_file_1_ok = read1.read(file_1_it, file_1_end, m_config.quality_input_fmt);
                const bool read_file_2_ok = read2.read(file_2_it, file_2_end, m_config.quality_input_fmt);

                if (read_file_1_ok != read_file_2_ok) {
                    throw fastq_error("files contain unequal number of records");
                } else if (!read_file_1_ok) {
                    break;
                }

                // Throws if read-names or mate numbering does not match
                fastq::validate_paired_reads(read1, read2);

                m_config.trim_barcodes_if_enabled(read1, *stats);

                // Reverse complement to match the orientation of read1
                read2.reverse_complement();

                const alignment_info alignment = align_paired_ended_sequences(read1, read2, m_config.adapters, m_config.shift);
                const userconfig::alignment_type aln_type = m_config.evaluate_alignment(alignment);
                if (aln_type == userconfig::valid_alignment) {
                    stats->well_aligned_reads++;
                    const size_t n_adapters = truncate_paired_ended_sequences(alignment, read1, read2);
                    stats->number_of_reads_with_adapter.at(alignment.adapter_id) += n_adapters;

                    if (m_config.is_alignment_collapsible(alignment)) {
                        fastq collapsed_read = collapse_paired_ended_sequences(alignment, read1, read2);
                        process_collapsed_read(m_config, *stats, collapsed_read,
                                               out_collapsed->reads,
                                               out_collapsed_truncated->reads,
                                               out_discarded->reads);
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
                    out_mate_1->reads.push_back(read1.to_str(m_config.quality_output_fmt));
                    out_mate_2->reads.push_back(read2.to_str(m_config.quality_output_fmt));

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
                        out_mate_1->reads.push_back(read1.to_str(m_config.quality_output_fmt));
                    } else {
                        out_singleton->reads.push_back(read1.to_str(m_config.quality_output_fmt));
                    }

                    if (read_2_acceptable) {
                        out_mate_2->reads.push_back(read2.to_str(m_config.quality_output_fmt));
                    } else {
                        out_singleton->reads.push_back(read2.to_str(m_config.quality_output_fmt));
                    }
                }
            }
        } catch (const fastq_error& error) {
            print_locker lock;
            std::cerr << "Error reading FASTQ record at line " << file_chunk->offset << "; aborting:\n    " << error.what() << std::endl;
            throw thread_abort();
        }

        stats->records += file_chunk->mates.at(rt_mate_1).size() / 4;
        m_stats.return_sink(stats.release());

        chunk_list chunks;
        const bool compress = m_config.gzip || m_config.bzip2;
        add_chunk(chunks, compress ? ai_zip_mate_1 : ai_write_mate_1, out_mate_1);
        add_chunk(chunks, compress ? ai_zip_mate_2 : ai_write_mate_2, out_mate_2);
        add_chunk(chunks, compress ? ai_zip_singleton : ai_write_singleton, out_singleton);
        add_chunk(chunks, compress ? ai_zip_collapsed : ai_write_collapsed, out_collapsed);
        add_chunk(chunks, compress ? ai_zip_collapsed_truncated : ai_write_collapsed_truncated, out_collapsed_truncated);
        add_chunk(chunks, compress ? ai_zip_discarded : ai_write_discarded, out_discarded);

        return chunks;
    }
};


int remove_adapter_sequences(const userconfig& config)
{
    if (!config.quiet) {
        if (config.paired_ended_mode) {
            std::cerr << "Trimming pair ended reads ..." << std::endl;
        } else {
            std::cerr << "Trimming single ended reads ..." << std::endl;
        }
    }

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

    scheduler sch;
    reads_processor* processor = NULL;
    try {
        if (config.paired_ended_mode) {
            sch.add_step(ai_read_mate_1, new read_paired_fastq(config, rt_mate_1, ai_read_mate_2));
            sch.add_step(ai_read_mate_2, new read_paired_fastq(config, rt_mate_2, ai_trim_pe));
            sch.add_step(ai_trim_pe, processor = new pe_reads_processor(config));
            sch.add_step(ai_write_mate_1, new write_paired_fastq(config, rt_mate_1));
            sch.add_step(ai_write_mate_2, new write_paired_fastq(config, rt_mate_2));
            sch.add_step(ai_write_singleton, new write_paired_fastq(config, rt_singleton));
        } else {
            sch.add_step(ai_read_mate_1, new read_paired_fastq(config, rt_mate_1, ai_trim_se));
            sch.add_step(ai_trim_se, processor = new se_reads_processor(config));
            sch.add_step(ai_write_mate_1, new write_paired_fastq(config, rt_mate_1));
        }

        if (config.collapse) {
            sch.add_step(ai_write_collapsed, new write_paired_fastq(config, rt_collapsed));
            sch.add_step(ai_write_collapsed_truncated, new write_paired_fastq(config, rt_collapsed_truncated));
        }

        if (config.gzip) {
            sch.add_step(ai_zip_mate_1, new gzip_paired_fastq(config, ai_write_mate_1));
            sch.add_step(ai_zip_mate_2, new gzip_paired_fastq(config, ai_write_mate_2));
            sch.add_step(ai_zip_singleton, new gzip_paired_fastq(config, ai_write_singleton));
            sch.add_step(ai_zip_collapsed, new gzip_paired_fastq(config, ai_write_collapsed));
            sch.add_step(ai_zip_collapsed_truncated, new gzip_paired_fastq(config, ai_write_collapsed_truncated));
            sch.add_step(ai_zip_discarded, new gzip_paired_fastq(config, ai_write_discarded));
#ifdef AR_BZIP2_SUPPORT
        } else if (config.bzip2) {
            sch.add_step(ai_zip_mate_1, new bzip2_paired_fastq(config, ai_write_mate_1));
            sch.add_step(ai_zip_mate_2, new bzip2_paired_fastq(config, ai_write_mate_2));
            sch.add_step(ai_zip_singleton, new bzip2_paired_fastq(config, ai_write_singleton));
            sch.add_step(ai_zip_collapsed, new bzip2_paired_fastq(config, ai_write_collapsed));
            sch.add_step(ai_zip_collapsed_truncated, new bzip2_paired_fastq(config, ai_write_collapsed_truncated));
            sch.add_step(ai_zip_discarded, new bzip2_paired_fastq(config, ai_write_discarded));
#endif
        }

        // Progress reporting enabled for final writer
        sch.add_step(ai_write_discarded, new write_paired_fastq(config, rt_discarded));
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n    " << error.what() << std::endl;
        return 1;
    }

    if (!sch.run(config.max_threads, config.seed)) {
        return 1;
    }

    std::auto_ptr<statistics> stats(processor->get_final_statistics());
    if (!write_statistics(config, *settings, *stats)) {
        std::cerr << "Error writing statistics to settings file!" << std::endl;
        return 1;
    }

    settings->flush();

    return 0;
}
