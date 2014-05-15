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
#include "statistics.h"


const size_t KMER_LENGTH = 13;
const size_t N_KMERS = 2 << (2 * KMER_LENGTH);
const size_t TOP_N_KMERS = 5;


inline size_t ACGT_TO_IDX(char nt)
{
    return (nt >> 1) & 0x3;
}


size_t kmer_to_size_t(const std::string& kmer)
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


void most_common_kmers(const kmer_map& kmers, bool reverse, size_t rjust)
{
    size_t total = 0;
    kmer_queue queue;
    for (size_t i = 0; i < kmers.size(); ++i) {
        kmer_count value(i, kmers.at(i));
        total += value.second;

        #warning remove
//        if (value.second) std::cout << size_t_to_kmer(i) << " " << value.second << "\n";

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

    std::cerr.precision(2);
    std::cerr << std::fixed;
    std::reverse(top_n_kmers.begin(), top_n_kmers.end());
    for (size_t i = 0; i < top_n_kmers.size(); ++i) {
        const kmer_count& count = top_n_kmers.at(i);
        std::string kmer_s = size_t_to_kmer(count.first);
        if (reverse) {
            std::reverse(kmer_s.begin(), kmer_s.end());
        }

        std::cerr << std::string(rjust, ' ');
        std::cerr << i + 1<< ": " << kmer_s
                  << " = " << std::right << std::setw(5) << (100.0 * count.second) / total
                  << "% (" << count.second << ")"
                  << "\n";
    }

    std::cerr.flush();
}



void print_consensus_adapter(const char_count_vec& counts,
                             const kmer_map& kmers,
                             const std::string& name,
                             const std::string& expected,
                             bool rjust = false)
{
    const char* NTs = "ACGT";
    std::stringstream sequence;
    std::stringstream qualities;

    for(char_count_vec::const_iterator it = counts.begin(); it != counts.end(); ++it) {
        char best_nt = 'N';
        size_t best_count = 0;
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
        qualities << p_to_phred_33(1.0 - best_count / static_cast<double>(total));
    }

    const std::string consensus = sequence.str();
    std::cerr << name << ":     " << expected << "\n";
    std::cerr << "Consensus:  " << consensus << "\n";
    std::cerr << "            " << qualities.str() << "\n";

    std::cerr << "            ";

    for (size_t i = 0; i < std::min(consensus.length(), expected.length()); ++i) {
        if (consensus.at(i) != expected.at(i)) {
            std::cerr << "X";
        } else {
            std::cerr << "-";
        }
    }

    std::cerr << "\n\n";
    most_common_kmers(kmers, rjust, rjust ? consensus.size() - KMER_LENGTH + 9 : 8);
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
    std::ifstream io_input_1;
    std::ifstream io_input_2;

    try {
        config.open_ifstream(io_input_1, config.input_file_1);
        config.open_ifstream(io_input_2, config.input_file_2);
    } catch (const std::ios_base::failure& error) {
        std::cerr << "IO error opening file; aborting:\n    " << error.what() << std::endl;
        return 1;
    }

    const fastq empty_adapter("PCR1", "", "");
    statistics stats;
    fastq read1;
    fastq read2;

    char_count_vec pcr1_counts;
    char_count_vec pcr2_counts;
    kmer_map pcr1_kmers(N_KMERS, 0);
    kmer_map pcr2_kmers(N_KMERS, 0);


    try {
        for (; ; ++stats.records) {
            const bool read_file_1_ok = read1.read(io_input_1, config.quality_input_fmt);
            const bool read_file_2_ok = read2.read(io_input_2, config.quality_input_fmt);

            if (read_file_1_ok != read_file_2_ok) {
                throw fastq_error("files contain unequal number of records");
            } else if (!read_file_1_ok) {
                break;
            }

            if (config.trim_barcode_if_enabled(read1)) {
                stats.number_of_barcodes_trimmed++;
            }

            // Reverse complement to match the orientation of read1
            read2.reverse_complement();

            const alignment_info alignment = align_paired_ended_sequences(read1, read2, empty_adapter, empty_adapter, config.shift);
            const userconfig::alignment_type aln_type = config.evaluate_alignment(alignment);
            if (aln_type == userconfig::valid_alignment) {
                stats.well_aligned_reads++;
                if (extract_adapter_sequences(alignment, read1, read2)) {
                    stats.number_of_reads_with_adapter++;

                    process_adapter(read1.sequence(), pcr1_counts, pcr1_kmers);

                    const std::string seq2 = std::string(read2.sequence().rbegin(), read2.sequence().rend());
                    process_adapter(seq2, pcr2_counts, pcr2_kmers);
                }
            } else if (aln_type == userconfig::poor_alignment) {
                stats.poorly_aligned_reads++;
            } else {
                stats.unaligned_reads++;
            }
        }
    } catch (const fastq_error& error) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << error.what() << std::endl;
        return 1;
    } catch (const std::ios_base::failure&) {
        std::cerr << "Error reading FASTQ record (" << stats.records << "); aborting:\n    " << std::strerror(errno) << std::endl;
        return 1;
    }

    std::cout << "Processed " << stats.records << " read pairs ...\n"
              << "   Found " << stats.well_aligned_reads << " overlapping pairs ...\n"
              << "   Of which " << stats.number_of_reads_with_adapter << " contained adapter sequence(s) ...\n"
              << std::endl;

    std::string padding = std::string(std::max<int>(0, pcr1_counts.size() - config.PCR1.length()),  'N');
    print_consensus_adapter(pcr1_counts, pcr1_kmers, "--pcr1", config.PCR1 + padding);
    std::cerr << "\n";

    std::reverse(pcr2_counts.begin(), pcr2_counts.end());
    padding = std::string(std::max<int>(0, pcr2_counts.size() - config.PCR2.length()),  'N');
    print_consensus_adapter(pcr2_counts, pcr2_kmers, "--pcr2", padding + config.PCR2, true);

    return 0;
}
