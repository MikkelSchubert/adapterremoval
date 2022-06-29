/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include "debug.hpp" // for AR_REQUIRE
#include <algorithm> // for min, max, sort
#include <utility>   // for pair

#include "barcode_table.hpp"

namespace adapterremoval {

typedef std::pair<std::string, size_t> barcode_pair;
typedef std::vector<barcode_pair> barcode_vec;

struct next_subsequence
{
  explicit next_subsequence(const char* seq_,
                            const size_t max_local_mismatches_)
    : seq(seq_)
    , max_local_mismatches(max_local_mismatches_)
  {
  }

  const char* seq;
  const size_t max_local_mismatches;
};

///////////////////////////////////////////////////////////////////////////////
// barcode_error

barcode_error::barcode_error(const std::string& message)
  : std::exception()
  , m_message(message)
{
}

barcode_error::barcode_error(const barcode_error& error)
  : std::exception()
  , m_message(error.m_message)
{
}

barcode_error::~barcode_error() {}

const char*
barcode_error::what() const noexcept
{
  return m_message.c_str();
}

///////////////////////////////////////////////////////////////////////////////

demultiplexer_node::demultiplexer_node()
  : children()
  , value(barcode_table::no_match)
{
  children.fill(barcode_table::no_match);
}

barcode_table::candidate::candidate(int barcode_, size_t mismatches_)
  : barcode(barcode_)
  , mismatches(mismatches_)
{
}

///////////////////////////////////////////////////////////////////////////////

/**
 * Returns a lexicographically sorted list of merged barcodes, each paired with
 * the 0-based index of corresponding barcode in the source vector.
 */
barcode_vec
sort_barcodes(const fastq_pair_vec& barcodes)
{
  AR_REQUIRE(!barcodes.empty());
  barcode_vec sorted_barcodes;

  const size_t max_key_1_len = barcodes.front().first.length();
  const size_t max_key_2_len = barcodes.front().second.length();
  for (auto it = barcodes.begin(); it != barcodes.end(); ++it) {
    if (it->first.length() != max_key_1_len) {
      throw barcode_error("mate 1 barcodes do not have the same length");
    } else if (it->second.length() != max_key_2_len) {
      throw barcode_error("mate 2 barcodes do not have the same length");
    }

    std::string barcode;
    barcode.reserve(max_key_1_len + max_key_2_len);
    barcode.append(it->first.sequence());
    barcode.append(it->second.sequence());

    sorted_barcodes.push_back(barcode_pair(barcode, it - barcodes.begin()));
  }

  std::sort(sorted_barcodes.begin(), sorted_barcodes.end());

  return sorted_barcodes;
}

/** Adds a nucleotide sequence with a given ID to a quad-tree. */
void
add_sequence_to_tree(demux_node_vec& tree,
                     const std::string& sequence,
                     const size_t barcode_id)
{
  size_t node_idx = 0;
  bool added_last_node = false;
  for (auto nuc : sequence) {
    auto& node = tree.at(node_idx);
    // Indicate when PE barcodes can be unambigiously identified from SE
    // reads
    node.value = (node.value == barcode_table::no_match)
                   ? barcode_id
                   : barcode_table::ambigious;

    const auto nuc_idx = ACGT::to_index(nuc);
    auto child = node.children[nuc_idx];

    added_last_node = (child == barcode_table::no_match);
    if (added_last_node) {
      // New nodes are added to the end of the list; as barcodes are
      // added in lexicographic order, this helps ensure that a set of
      // similar barcodes will be placed in mostly contiguous runs
      // of the vector representation.
      child = node.children[nuc_idx] = tree.size();
      tree.push_back(demultiplexer_node());
    }

    node_idx = child;
  }

  if (!added_last_node) {
    throw barcode_error(std::string("duplicate barcode (pair): ") + sequence);
  }

  tree.at(node_idx).value = barcode_id;
}

/**
 * Builds a sparse quad tree using the first sequence in a set of unique
 * barcodes pairs; duplicate pairs will negatively impact the identification of
 * these, since all hits will be considered ambiguous.
 */
demux_node_vec
build_demux_tree(const fastq_pair_vec& barcodes)
{
  // Step 1: Construct list of merged, sorted barcodes barcodes;
  //         this allows construction of the sparse tree in one pass.
  const barcode_vec sorted_barcodes = sort_barcodes(barcodes);

  // Step 2: Create empty tree containing just the root node; creating
  //         the root here simplifies the 'add_sequence_to_tree' function.
  demux_node_vec tree;
  tree.push_back(demultiplexer_node());

  // Step 3: Add each barcode to the tree, in sorted order
  for (auto& pair : sorted_barcodes) {
    add_sequence_to_tree(tree, pair.first, pair.second);
  }

  return tree;
}

///////////////////////////////////////////////////////////////////////////////

const int barcode_table::no_match;
const int barcode_table::ambigious;

barcode_table::barcode_table(const fastq_pair_vec& barcodes,
                             size_t mismatches,
                             size_t mm_r1,
                             size_t mm_r2)
  : m_nodes()
  , m_max_mismatches()
  , m_max_mismatches_r1()
  , m_max_mismatches_r2()
  , m_barcode_1_len()
  , m_barcode_2_len()
{
  m_max_mismatches = std::min<size_t>(mismatches, mm_r1 + mm_r2);
  m_max_mismatches_r1 = std::min<size_t>(m_max_mismatches, mm_r1);
  m_max_mismatches_r2 = std::min<size_t>(m_max_mismatches, mm_r2);

  if (!barcodes.empty()) {
    m_nodes = build_demux_tree(barcodes);
    m_barcode_1_len = barcodes.front().first.length();
    m_barcode_2_len = barcodes.front().second.length();
  }
}

int
barcode_table::identify(const fastq& read_r1) const
{
  if (read_r1.length() < m_barcode_1_len) {
    return barcode_table::no_match;
  }

  const std::string barcode = read_r1.sequence().substr(0, m_barcode_1_len);
  auto match = lookup(barcode.c_str(), 0, 0, nullptr);
  if (match.barcode == no_match && m_max_mismatches) {
    match = lookup_with_mm(
      barcode.c_str(), 0, m_max_mismatches, m_max_mismatches_r1, nullptr);
  }

  return match.barcode;
}

int
barcode_table::identify(const fastq& read_r1, const fastq& read_r2) const
{
  if (read_r1.length() < m_barcode_1_len ||
      read_r2.length() < m_barcode_2_len) {
    return no_match;
  }

  const auto barcode_1 = read_r1.sequence().substr(0, m_barcode_1_len);
  const auto barcode_2 = read_r2.sequence().substr(0, m_barcode_2_len);
  const auto combined_barcode = barcode_1 + barcode_2;

  auto match = lookup(combined_barcode.c_str(), 0, 0, nullptr);
  if (match.barcode == no_match && m_max_mismatches) {
    const next_subsequence next(barcode_2.c_str(), m_max_mismatches_r2);

    if (m_max_mismatches_r1) {
      match = lookup_with_mm(
        barcode_1.c_str(), 0, m_max_mismatches, m_max_mismatches_r1, &next);
    } else {
      match = lookup(barcode_1.c_str(), 0, m_max_mismatches, &next);
    }
  }

  return match.barcode;
}

barcode_table::candidate
barcode_table::lookup(const char* seq,
                      int parent,
                      const size_t max_global_mismatches,
                      const next_subsequence* next) const

{
  for (; *seq && parent != no_match; ++seq) {
    if (*seq == 'N') {
      return candidate(no_match);
    }

    const auto encoded_nuc = ACGT::to_index(*seq);
    const auto& node = m_nodes.at(parent);

    parent = node.children.at(encoded_nuc);
  }

  if (parent == no_match) {
    return candidate(no_match);
  } else if (next) {
    const auto max_local_mismatches =
      std::min(max_global_mismatches, next->max_local_mismatches);

    if (max_local_mismatches) {
      return lookup_with_mm(next->seq,
                            parent,
                            max_global_mismatches,
                            max_local_mismatches,
                            nullptr);
    } else {
      return lookup(next->seq, parent, max_global_mismatches, nullptr);
    }
  } else {
    const auto& node = m_nodes.at(parent);
    return candidate(node.value, m_max_mismatches - max_global_mismatches);
  }
}

barcode_table::candidate
barcode_table::lookup_with_mm(const char* seq,
                              int parent,
                              const size_t max_global_mismatches,
                              const size_t max_local_mismatches,
                              const next_subsequence* next) const

{
  const demultiplexer_node& node = m_nodes.at(parent);
  const auto nucleotide = *seq;

  if (nucleotide) {
    candidate best_candidate;

    for (size_t encoded_i = 0; encoded_i < 4; ++encoded_i) {
      const auto child = node.children.at(encoded_i);

      if (child != -1) {
        candidate current_candidate;

        if (nucleotide == ACGT::to_value(encoded_i)) {
          current_candidate = lookup_with_mm(
            seq + 1, child, max_global_mismatches, max_local_mismatches, next);
        } else if (max_local_mismatches == 1) {
          current_candidate =
            lookup(seq + 1, child, max_global_mismatches - 1, next);
        } else if (max_local_mismatches) {
          current_candidate = lookup_with_mm(seq + 1,
                                             child,
                                             max_global_mismatches - 1,
                                             max_local_mismatches - 1,
                                             next);
        }

        if (current_candidate.mismatches < best_candidate.mismatches) {
          best_candidate = current_candidate;
        } else if (current_candidate.mismatches == best_candidate.mismatches) {
          if (current_candidate.barcode != no_match) {
            best_candidate.barcode = ambigious;
          }
        }
      }
    }

    return best_candidate;
  } else if (next) {
    const size_t next_max_local_mismatches =
      std::min(max_global_mismatches, next->max_local_mismatches);

    if (next_max_local_mismatches) {
      return lookup_with_mm(next->seq,
                            parent,
                            max_global_mismatches,
                            next_max_local_mismatches,
                            nullptr);
    } else {
      return lookup(next->seq, parent, max_global_mismatches, nullptr);
    }
  } else {
    return candidate(node.value, m_max_mismatches - max_global_mismatches);
  }
}

} // namespace adapterremoval
