#include <algorithm>
#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <vector>
#include <queue>
#include <iomanip>


#include "main_adapter_id.h"
#include "userconfig.h"
#include "alignment.h"
#include "timer.h"


const size_t KMER_LENGTH = 12;
const size_t N_KMERS = 2 << (2 * KMER_LENGTH);
const size_t TOP_N_KMERS = 5;


inline size_t ACGT_TO_IDX(char nt)
{
    return (nt >> 1) & 0x3;
}


inline size_t kmer_to_size_t(const std::string& kmer)
{
    size_t index = 0;
    for (size_t i = 0; i < kmer.length(); ++i) {
        index = (index << 2) | ACGT_TO_IDX(kmer.at(i));
    }

    return index;
}


std::string size_t_to_kmer(size_t kmer)
{
    std::string kmer_s(KMER_LENGTH, 'N');
    for (size_t i = 1; i <= KMER_LENGTH; ++i) {
        kmer_s.at(KMER_LENGTH - i) = "ACTG"[kmer & 0x3];
        kmer = kmer >> 2;
    }

    return kmer_s;
}


struct char_counts
{
    char_counts()
      : counts(4, 0)
    {
    }

    std::vector<size_t> counts;
};


typedef std::vector<char_counts> char_count_vec;
typedef std::vector<unsigned> kmer_map;
typedef std::pair<size_t, unsigned> kmer_count;

struct cmp_kmer_count
{
    bool operator()(const kmer_count& a, const kmer_count& b) const
    {
        return (a.second > b.second);
    }
};

typedef std::priority_queue<kmer_count, std::vector<kmer_count>, cmp_kmer_count> kmer_queue;
typedef std::vector<kmer_count> kmer_vector;


void print_most_common_kmers(const kmer_map& kmers)
{
    size_t total = 0;
    kmer_queue queue;
    for (size_t i = 0; i < kmers.size(); ++i) {
        kmer_count value(i, kmers.at(i));
        total += value.second;

        if (queue.size() >= TOP_N_KMERS) {
            if (queue.top().second < value.second) {
                queue.pop();
                queue.push(value);
            }
        } else if (value.second) {
            queue.push(value);
        }
    }

    kmer_vector top_n_kmers;
    while (!queue.empty()) {
        top_n_kmers.push_back(queue.top());
        queue.pop();
    }

    std::cout.precision(2);
    std::cout << std::fixed;
    std::cout << "    Top 5 most common " << KMER_LENGTH << "-bp 5'-kmers:\n";

    std::reverse(top_n_kmers.begin(), top_n_kmers.end());
    for (size_t i = 0; i < top_n_kmers.size(); ++i) {
        const kmer_count& count = top_n_kmers.at(i);
        std::string kmer_s = size_t_to_kmer(count.first);

        std::cout << std::string(12, ' ');
        std::cout << i + 1<< ": " << kmer_s
                  << " = " << std::right << std::setw(5) << (100.0 * count.second) / total
                  << "% (" << count.second << ")"
                  << "\n";
    }
}


void print_consensus_adapter(const char_count_vec& counts,
                             const kmer_map& kmers,
                             const std::string& name)
{
    const char* NTs = "ACGT";
    std::stringstream sequence;
    std::stringstream qualities;

    for(char_count_vec::const_iterator it = counts.begin(); it != counts.end(); ++it) {
        char best_nt = 'N';
        size_t best_count = 0;
        // Always assume one non-consensus observation; this is more reasonable
        // than allowing an error-rate of 0, especially for few observations.
        size_t total = 1;
        for (size_t nt_i = 0; nt_i < 4; ++nt_i) {
            const size_t cur_count = it->counts.at(ACGT_TO_IDX(NTs[nt_i]));
            total += cur_count;

            if (cur_count > best_count) {
                best_nt = NTs[nt_i];
                best_count = cur_count;
            }
        }

        sequence << best_nt;
        qualities << fastq::p_to_phred_33(1.0 - best_count / static_cast<double>(total));
    }

    const std::string consensus = sequence.str();
    std::cout << "  " << name << ":  " << consensus << "\n"
              << "               " << qualities.str() << "\n\n";

    print_most_common_kmers(kmers);
}


void process_adapter(const std::string& sequence, char_count_vec& counts, kmer_map& kmers)
{
    if (counts.size() < sequence.length()) {
        counts.resize(sequence.length());
    }

    for (size_t i = 0; i < std::min(counts.size(), sequence.length()); ++i) {
        counts.at(i).counts.at(ACGT_TO_IDX(sequence.at(i)))++;
    }

    if (sequence.length() >= KMER_LENGTH) {
        const std::string kmer = sequence.substr(0, KMER_LENGTH);
        if (!std::count(kmer.begin(), kmer.end(), 'N')) {
            kmers.at(kmer_to_size_t(kmer)) += 1;
        }
    }
}


int identify_adapter_sequences(const userconfig& config)
{
    std::auto_ptr<std::istream> io_input_1;
    std::auto_ptr<std::istream> io_input_2;

    try {
        io_input_1 = config.open_ifstream(config.input_file_1);
        io_input_2 = config.open_ifstream(config.input_file_2);
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n    " << error.what() << std::endl;
        return 1;
    }

    const fastq empty_adapter("dummy", "", "");
    fastq_pair_vec adapters;
    adapters.push_back(fastq_pair(empty_adapter, empty_adapter));

    std::auto_ptr<statistics> stats_ptr = config.create_stats();
    statistics& stats = *stats_ptr;
    fastq read1;
    fastq read2;

    char_count_vec pcr1_counts;
    char_count_vec pcr2_counts;
    kmer_map pcr1_kmers(N_KMERS, 0);
    kmer_map pcr2_kmers(N_KMERS, 0);

    std::cout << "Attemping to identify adapter sequences ..." << std::endl;

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

            const alignment_info alignment = align_paired_ended_sequences(read1, read2, adapters, config.shift);
            const userconfig::alignment_type aln_type = config.evaluate_alignment(alignment);
            if (aln_type == userconfig::valid_alignment) {
                stats.well_aligned_reads++;
                if (!config.is_alignment_collapsible(alignment)) {
                    continue;
                }

                if (extract_adapter_sequences(alignment, read1, read2)) {
                    stats.number_of_reads_with_adapter.at(0)++;

                    process_adapter(read1.sequence(), pcr1_counts, pcr1_kmers);

                    read2.reverse_complement();
                    process_adapter(read2.sequence(), pcr2_counts, pcr2_kmers);
                }
            } else if (aln_type == userconfig::poor_alignment) {
                stats.poorly_aligned_reads++;
            } else {
                stats.unaligned_reads++;
            }
        }

        progress.finalize();
    } catch (const fastq_error& error) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << error.what() << std::endl;
        return 1;
    } catch (const std::ios_base::failure&) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << std::strerror(errno) << std::endl;
        return 1;
    }

    std::cout << "   Found " << stats.well_aligned_reads << " overlapping pairs ...\n"
              << "   Of which " << stats.number_of_reads_with_adapter.at(0) << " contained adapter sequence(s) ...\n\n"
              << "Printing adapter sequences, including poly-A tails:"
              << std::endl;

    print_consensus_adapter(pcr1_counts, pcr1_kmers, "--adapter1");
    std::cout << "\n\n";
    print_consensus_adapter(pcr2_counts, pcr2_kmers, "--adapter2");

    return 0;
}
