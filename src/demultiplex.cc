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

#include "debug.h"
#include "demultiplex.h"
#include "commontypes.h"
#include "fastq_io.h"
#include "userconfig.h"
#include "strutils.h"


typedef std::vector<unsigned> int_vec;
typedef demux_node_vec::iterator node_vec_iter;


///////////////////////////////////////////////////////////////////////////////

/**
 * Struct representing node in quad-tree; children are referenced using the
 * corresponding indice in the vector representing the tree; -1 is used to
 * represent unasigned children.
 */
struct demultiplexer_node
{
public:
    demultiplexer_node()
        : children()
        , barcodes()
    {
        children[0] = -1;
        children[1] = -1;
        children[2] = -1;
        children[3] = -1;
    }

    bool has_children() const
    {
        for (size_t nt_idx = 0; nt_idx < 4; ++nt_idx) {
            if (children[nt_idx] != -1) {
                 return true;
            }
        }

        return false;
    }

    int children[4];
    int_vec barcodes;
};


///////////////////////////////////////////////////////////////////////////////

/** Fills a quadtree with nodes representing keys of len 'max_key_len' */
void fill_demux_tree(demux_node_vec& tree, const size_t max_key_len)
{
    const size_t n_keys = (1 - (1 << (2 * (max_key_len + 1)))) / (1 - 4);
    const size_t n_internal = n_keys - (1 << (2 * max_key_len));

    tree.resize(n_keys);

    for (size_t index = 0; index < n_keys; ++index) {
        demultiplexer_node& node = tree.at(index);

        for (size_t nt_idx = 0; nt_idx < 4; ++nt_idx) {
            node.children[nt_idx] = (index < n_internal) ? 4 * index + nt_idx + 1 : -1;
        }
    }
}


/** Adds a nucleotide sequence with a given ID to a quadtree. */
void add_sequence_to_tree(demux_node_vec& tree,
                          const std::string& sequence,
                          const size_t barcode_id)
{
    demultiplexer_node* parent = &tree.at(0);
    for (std::string::const_iterator it = sequence.begin(); it != sequence.end(); ++it) {
        const int child = parent->children[ACGT_TO_IDX(*it)];
        AR_DEBUG_ASSERT(child != -1);

        parent = &tree.at(child);
    }

    parent->barcodes.push_back(barcode_id);
}


/**
 * Converts a full quad-tree to a sparse quad-tree, by removing any nodes
 * which does not contain any barcodes and which does not have any children.
 * Once these have been removed, the resulting list is compressed such that
 * only populated parts of the tree are included in the vector.
 */
void prune_tree(demux_node_vec& tree)
{
    if (tree.empty()) {
        return;
    }

    for (size_t idx = tree.size() - 1; idx; --idx) {
        const demultiplexer_node& node = tree.at(0);

        if (node.barcodes.empty() && !node.has_children()) {
            tree.at((idx - 1) / 4).children[(idx - 1) % 4] = -1;
        }
    }

    size_t last_index = 0;
    std::vector<int> new_indices(tree.size(), -1);
    for (size_t idx = 0; idx < tree.size(); ++idx) {
        const demultiplexer_node& node = tree.at(0);

        if (!node.barcodes.empty() || node.has_children()) {
            if (idx != last_index) {
                tree.at(last_index) = node;
            }

            new_indices.at(idx) = last_index;
            last_index++;
        }
    }

    tree.resize(last_index);

    for (size_t idx = 0; idx < tree.size(); ++idx) {
        demultiplexer_node& node = tree.at(0);

        for (size_t nt_idx = 0; nt_idx < 4; ++nt_idx) {
            if (node.children[nt_idx] != -1) {
                node.children[nt_idx] = new_indices.at(node.children[nt_idx]);
            }
        }
    }
}


/**
 * Builds a sparse quad tree using the first sequence in a set of unique
 * barcodes pairs; duplicate pairs will negatively impact the identification of
 * these, since all hits will be considered ambiguous.
 */
demux_node_vec build_demux_tree(const fastq_pair_vec& barcodes)
{
    demux_node_vec tree;

    if (barcodes.empty()) {
        return tree;
    }

    const size_t max_key_1_len = barcodes.front().first.length();
    const size_t max_key_2_len = barcodes.front().second.length();
    for (fastq_pair_vec::const_iterator it = barcodes.begin(); it != barcodes.end(); ++it) {
        AR_DEBUG_ASSERT(it->first.length() == max_key_1_len);
        AR_DEBUG_ASSERT(it->second.length() == max_key_2_len);
    }

    fill_demux_tree(tree, max_key_1_len);

    for (fastq_pair_vec::const_iterator it = barcodes.begin(); it != barcodes.end(); ++it) {
        add_sequence_to_tree(tree, it->first.sequence(), it - barcodes.begin());
    }

    prune_tree(tree);

    return tree;
}


///////////////////////////////////////////////////////////////////////////////

typedef std::vector<std::pair<int, size_t> > candidate_vec;


void rec_lookup_sequence_no_mm(candidate_vec& candidates,
                               const demux_node_vec& tree,
                               const std::string& seq,
                               size_t seq_pos = 0,
                               int parent = 0,
                               size_t mismatches = 0)
{
    const demultiplexer_node& node = tree.at(parent);
    for (int_vec::const_iterator it = node.barcodes.begin(); it != node.barcodes.end(); ++it) {
        candidates.push_back(std::pair<int, size_t>(*it, mismatches));
    }

    if (seq_pos < seq.length()) {
        const size_t nt_idx = ACGT_TO_IDX(seq.at(seq_pos));

        if (node.children[nt_idx] != -1) {
            rec_lookup_sequence_no_mm(candidates, tree, seq, seq_pos + 1,
                                      node.children[nt_idx], mismatches);
        }
    }
}


void rec_lookup_sequence(candidate_vec& candidates,
                         const demux_node_vec& tree,
                         const std::string& seq,
                         size_t max_mismatches,
                         size_t seq_pos = 0,
                         int parent = 0,
                         size_t mismatches = 0)
{
    const demultiplexer_node& node = tree.at(parent);
    for (int_vec::const_iterator it = node.barcodes.begin(); it != node.barcodes.end(); ++it) {
        candidates.push_back(std::pair<int, size_t>(*it, mismatches));
    }

    if (seq_pos < seq.length()) {
        const size_t current_nt = ACGT_TO_IDX(seq.at(seq_pos));

        for (size_t nt_idx = 0; nt_idx < 4; ++nt_idx) {
            if (node.children[nt_idx] != -1) {
                if (nt_idx == current_nt) {
                    rec_lookup_sequence(candidates, tree, seq, max_mismatches, seq_pos + 1, node.children[nt_idx], mismatches);
                } else if (mismatches + 1 < max_mismatches) {
                    rec_lookup_sequence(candidates, tree, seq, max_mismatches, seq_pos + 1, node.children[nt_idx], mismatches + 1);
                } else if (mismatches + 1 == max_mismatches) {
                    rec_lookup_sequence_no_mm(candidates, tree, seq, seq_pos + 1, node.children[nt_idx], mismatches + 1);
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////


demultiplex_reads::demultiplex_reads(const userconfig* config)
    : analytical_step(analytical_step::ordered)
    , m_barcodes(config->adapters.get_barcodes())
    , m_tree(build_demux_tree(m_barcodes))
    , m_max_mismatches(config->barcode_mm)
    , m_max_mismatches_r1(std::min<size_t>(config->barcode_mm, config->barcode_mm_r1))
    , m_max_mismatches_r2(std::min<size_t>(config->barcode_mm, config->barcode_mm_r2))
    , m_config(config)
    , m_cache(m_barcodes.size(), NULL)
    , m_unidentified_1(new fastq_output_chunk())
    , m_unidentified_2(new fastq_output_chunk())
    , m_statistics(m_barcodes.size())
{
    AR_DEBUG_ASSERT(!m_barcodes.empty());

    for (demultiplexed_cache::iterator it = m_cache.begin(); it != m_cache.end(); ++it) {
        *it = new fastq_read_chunk();
    }
}


demultiplex_reads::~demultiplex_reads()
{
    for (demultiplexed_cache::iterator it = m_cache.begin(); it != m_cache.end(); ++it) {
        delete *it;
    }

    delete m_unidentified_1;
    delete m_unidentified_2;
}



size_t count_mismatches(const std::string& barcode,
                        const std::string& sequence,
                        const size_t max_mismatches)
{
    std::string::const_iterator b_iter = barcode.begin();
    std::string::const_iterator r2_iter = sequence.begin();

    size_t mismatches = 0;
    size_t n_check = std::min(barcode.length(), sequence.length());
    if (n_check < barcode.length()) {
        // Missing bases are considered mismatching
        mismatches += barcode.length() - n_check;
    }

    while (n_check-- && mismatches <= max_mismatches) {
        if (*b_iter++ != *r2_iter++) {
            ++mismatches;
        }
    }

    return mismatches;
}


/**
 * Returns the best matching barcode (pair) for sequences read_r1 and read_r2
 *
 */
int demultiplex_reads::select_barcode(const fastq& read_r1, const fastq& read_r2)
{
    candidate_vec candidates;
    if (m_max_mismatches_r1) {
        rec_lookup_sequence(candidates, m_tree, read_r1.sequence(), m_max_mismatches_r1);
    } else {
        rec_lookup_sequence_no_mm(candidates, m_tree, read_r1.sequence());
    }

    int best_barcode = -1;
    size_t min_mismatches = m_max_mismatches + 1;
    for (candidate_vec::iterator it = candidates.begin(); it != candidates.end(); ++it) {
        if (m_config->paired_ended_mode) {
            const std::string& barcode = m_barcodes.at(it->first).second.sequence();
            const size_t max_mismatches_r2 = std::min(m_max_mismatches - it->second,
                                                      m_max_mismatches_r2);

            const size_t mismatches = count_mismatches(barcode,
                                                       read_r2.sequence(),
                                                       max_mismatches_r2);

            if (mismatches > max_mismatches_r2) {
                continue;
            }

            it->second += mismatches;
        }

        if (it->second < min_mismatches) {
            best_barcode = it->first;
            min_mismatches = it->second;
        } else if (it->second == min_mismatches) {
            // Ambiguous results; multiple best matches
            best_barcode = -1;
        }
    }

    if (best_barcode >= 0) {
        return best_barcode;
    } else if (min_mismatches == m_max_mismatches + 1) {
        // No viable candidates
        return -1;
    } else {
        // Ambiguous results
        return -2;
    }
}


chunk_vec demultiplex_reads::flush_cache(bool eof)
{
    chunk_vec output;

    if (eof || m_unidentified_1->count >= FASTQ_CHUNK_SIZE) {
        output.push_back(chunk_pair(ai_write_unidentified_1, m_unidentified_1));
        m_unidentified_1->eof = eof;
        m_unidentified_1 = new fastq_output_chunk();
    }

    if (m_config->paired_ended_mode && (eof || m_unidentified_2->count >= FASTQ_CHUNK_SIZE)) {
        output.push_back(chunk_pair(ai_write_unidentified_2, m_unidentified_2));
        m_unidentified_2->eof = eof;
        m_unidentified_2 = new fastq_output_chunk();
    }

    for (size_t nth = 0; nth < m_cache.size(); ++nth) {
        fastq_read_chunk* chunk = m_cache.at(nth);
        if (eof || chunk->reads_1.size() >= FASTQ_CHUNK_SIZE) {
            chunk->eof = eof;

            const size_t step_id = (nth + 1) * ai_analyses_offset;
            output.push_back(chunk_pair(step_id, chunk));
            m_cache.at(nth) = new fastq_read_chunk();
        }
    }

    return output;
}


demux_statistics demultiplex_reads::statistics() const
{
    return m_statistics;
}


///////////////////////////////////////////////////////////////////////////////

demultiplex_se_reads::demultiplex_se_reads(const userconfig* config)
    : demultiplex_reads(config)
{
}


chunk_vec demultiplex_se_reads::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_read_chunk> read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));

    fastq empty_read;
    for (fastq_vec::iterator it = read_chunk->reads_1.begin(); it != read_chunk->reads_1.end(); ++it) {
        const int best_barcode = select_barcode(*it, empty_read);

        if (best_barcode < 0) {
            m_unidentified_1->add(*m_config->quality_output_fmt, *it);

            if (best_barcode == -1) {
                m_statistics.unidentified += 1;
            } else {
                m_statistics.ambiguous += 1;
            }
        } else {
            fastq_read_chunk* dst = m_cache.at(best_barcode);
            dst->reads_1.push_back(*it);
            dst->reads_1.back().truncate(m_barcodes.at(best_barcode).first.length());

           m_statistics.barcodes.at(best_barcode) += 1;
        }
    }

    return flush_cache(read_chunk->eof);
}


///////////////////////////////////////////////////////////////////////////////

demultiplex_pe_reads::demultiplex_pe_reads(const userconfig* config)
    : demultiplex_reads(config)
{
}


chunk_vec demultiplex_pe_reads::process(analytical_chunk* chunk)
{
    std::auto_ptr<fastq_read_chunk> read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
    AR_DEBUG_ASSERT(read_chunk->reads_1.size() == read_chunk->reads_2.size());

    fastq_vec::iterator it_1 = read_chunk->reads_1.begin();
    fastq_vec::iterator it_2 = read_chunk->reads_2.begin();
    for (; it_1 != read_chunk->reads_1.end(); ++it_1, ++it_2) {
        const int best_barcode = select_barcode(*it_1, *it_2);

        if (best_barcode < 0) {
            m_unidentified_1->add(*m_config->quality_output_fmt, *it_1);
            m_unidentified_2->add(*m_config->quality_output_fmt, *it_2);

            if (best_barcode == -1) {
                m_statistics.unidentified += 1;
            } else {
                m_statistics.ambiguous += 1;
            }
        } else {
            fastq_read_chunk* dst = m_cache.at(best_barcode);

            it_1->truncate(m_barcodes.at(best_barcode).first.length());
            dst->reads_1.push_back(*it_1);
            it_2->truncate(m_barcodes.at(best_barcode).second.length());
            dst->reads_2.push_back(*it_2);

            m_statistics.barcodes.at(best_barcode) += 1;
        }
    }

    return flush_cache(read_chunk->eof);
}
