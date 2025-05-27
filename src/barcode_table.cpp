// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "barcode_table.hpp" // declarations
#include "debug.hpp"         // for AR_REQUIRE
#include "errors.hpp"        // for parsing_error
#include "fastq.hpp"         // for fastq
#include "fastq_enc.hpp"     // for ACGT
#include "sequence.hpp"      // for sequence_pair
#include "sequence_sets.hpp" // for sample_set
#include "strutils.hpp"      // for stringify
#include <algorithm>         // for min, max, sort
#include <cstddef>           // for size_t
#include <cstdint>           // for int32_t, uint32_t
#include <ostream>           // for ostream
#include <string>            // for std::string
#include <utility>           // for pair

namespace adapterremoval {

using barcode_pair = std::pair<sequence_pair, barcode_key>;
using barcode_vec = std::vector<barcode_pair>;

struct barcode_match
{
  barcode_match() = default;

  explicit barcode_match(const barcode_key& key_, uint32_t mismatches_ = -1)
    : key(key_)
    , mismatches(mismatches_)
  {
  }

  barcode_key key{};
  uint32_t mismatches = -1;
};

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

barcode_node::barcode_node()
  : children()
  , key()
{
  children.fill(barcode_key::unidentified);
}

///////////////////////////////////////////////////////////////////////////////

namespace {

/** Adds a nucleotide sequence with a given ID to a quad-tree. */
void
add_sequence_to_tree(std::vector<barcode_node>& tree,
                     const std::string& sequence,
                     const barcode_key key)
{
  size_t node_idx = 0;
  bool added_last_node = false;
  for (auto nuc : sequence) {
    auto& node = tree.at(node_idx);
    // Indicate when PE barcodes can be unambiguously identified from SE reads
    if (node.key.sample == barcode_key::unidentified) {
      node.key.sample = key.sample;
      node.key.barcode = key.barcode;
    } else if (node.key.sample == key.sample) {
      node.key.barcode = barcode_key::ambiguous;
    } else {
      node.key.sample = barcode_key::ambiguous;
      node.key.barcode = barcode_key::ambiguous;
    }

    const auto nuc_idx = ACGT::to_index(nuc);
    auto child = node.children[nuc_idx];

    added_last_node = (child == barcode_key::unidentified);
    if (added_last_node) {
      // New nodes are added to the end of the list; as barcodes are added in
      // lexicographic order, this helps ensure that a set of similar barcodes
      // will be placed in mostly contiguous runs of the vector representation.
      child = node.children[nuc_idx] = tree.size();
      tree.emplace_back();
    }

    node_idx = child;
  }

  if (!added_last_node) {
    throw parsing_error(std::string("duplicate barcode(s): ") + sequence);
  }

  tree.at(node_idx).key = key;
}

barcode_vec
build_barcode_vec(const sample_set& samples)
{
  int32_t sample_i = 0;
  barcode_vec barcodes;
  for (const auto& sample : samples) {
    AR_REQUIRE(!sample.name().empty());

    int32_t barcode_i = 0;
    for (const auto& sequences : sample) {
      barcodes.emplace_back(
        sequence_pair{ sequences.barcode_1, sequences.barcode_2 },
        barcode_key{ sample_i, barcode_i++ });
    }

    sample_i++;
  }

  std::sort(barcodes.begin(), barcodes.end());

  return barcodes;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////

barcode_table::barcode_table(const sample_set& samples,
                             size_t max_mm,
                             size_t max_mm_r1,
                             size_t max_mm_r2)
{
  AR_REQUIRE(samples.size());
  m_max_mismatches = std::min<size_t>(max_mm, max_mm_r1 + max_mm_r2);
  m_max_mismatches_r1 = std::min<size_t>(m_max_mismatches, max_mm_r1);
  m_max_mismatches_r2 = std::min<size_t>(m_max_mismatches, max_mm_r2);

  // Flatten and lexiographically sort barcodes to simplify tree building
  auto barcodes = build_barcode_vec(samples);

  // Create empty tree containing just the root node; creating the root here
  // simplifies the 'add_sequence_to_tree' function.
  m_nodes.emplace_back();

  const auto& front = barcodes.front().first;
  m_barcode_1_len = front.first.length();
  m_barcode_2_len = front.second.length();

  // Step 3: Add each barcode to the tree, in sorted order
  for (const auto& [sequences, key] : barcodes) {
    AR_REQUIRE(m_barcode_1_len && sequences.first.length() == m_barcode_1_len &&
               sequences.second.length() == m_barcode_2_len);

    std::string barcode;
    barcode.reserve(m_barcode_1_len + m_barcode_2_len);
    barcode.append(sequences.first.as_string());
    barcode.append(sequences.second.as_string());

    add_sequence_to_tree(m_nodes, barcode, key);
  }
}

barcode_key
barcode_table::identify(const fastq& read_r1) const
{
  if (read_r1.length() < m_barcode_1_len) {
    return barcode_key{};
  }

  const std::string barcode = read_r1.sequence().substr(0, m_barcode_1_len);
  auto match = lookup(barcode.c_str(), 0, 0, nullptr);
  if (match.key.sample == barcode_key::unidentified && m_max_mismatches_r1) {
    match = lookup_with_mm(barcode.c_str(),
                           0,
                           m_max_mismatches_r1,
                           m_max_mismatches_r1,
                           nullptr);
  }

  return match.key;
}

barcode_key
barcode_table::identify(const fastq& read_r1, const fastq& read_r2) const
{
  if (read_r1.length() < m_barcode_1_len ||
      read_r2.length() < m_barcode_2_len) {
    return barcode_key{};
  }

  const auto barcode_1 = read_r1.sequence().substr(0, m_barcode_1_len);
  const auto barcode_2 = read_r2.sequence().substr(0, m_barcode_2_len);

  auto match = lookup(barcode_1.c_str(), 0, 0, barcode_2.c_str());
  if (match.key.sample == barcode_key::unidentified && m_max_mismatches) {
    if (m_max_mismatches_r1) {
      match = lookup_with_mm(barcode_1.c_str(),
                             0,
                             m_max_mismatches,
                             m_max_mismatches_r1,
                             barcode_2.c_str());
    } else {
      match = lookup(barcode_1.c_str(), 0, m_max_mismatches, barcode_2.c_str());
    }
  }

  return match.key;
}

barcode_match
barcode_table::lookup(const char* seq,
                      int32_t parent,
                      const size_t max_global_mismatches,
                      const char* next) const
{
  for (; *seq && parent != barcode_key::unidentified; ++seq) {
    if (*seq == 'N') {
      return {};
    }

    const auto encoded_nuc = ACGT::to_index(*seq);
    const auto& node = m_nodes.at(parent);

    parent = node.children.at(encoded_nuc);
  }

  if (parent == barcode_key::unidentified) {
    return {};
  } else if (next) {
    const auto max_local_mismatches =
      std::min(max_global_mismatches, m_max_mismatches_r2);

    if (max_local_mismatches) {
      return lookup_with_mm(next,
                            parent,
                            max_global_mismatches,
                            max_local_mismatches,
                            nullptr);
    } else {
      return lookup(next, parent, max_global_mismatches, nullptr);
    }
  } else {
    const auto& node = m_nodes.at(parent);
    return barcode_match(node.key, m_max_mismatches - max_global_mismatches);
  }
}

barcode_match
barcode_table::lookup_with_mm(const char* seq,
                              int32_t parent,
                              const size_t max_global_mismatches,
                              size_t max_local_mismatches,
                              const char* next) const

{
  const barcode_node& node = m_nodes.at(parent);
  const auto nucleotide = *seq;

  if (nucleotide) {
    barcode_match best_candidate;

    for (size_t encoded_i = 0; encoded_i < 4; ++encoded_i) {
      const auto child = node.children.at(encoded_i);

      if (child != barcode_key::unidentified) {
        barcode_match current_candidate;

        // Nucleotide may be 'N', so test has to be done by converting the index
        if (nucleotide == ACGT::to_value(encoded_i)) {
          current_candidate = lookup_with_mm(seq + 1,
                                             child,
                                             max_global_mismatches,
                                             max_local_mismatches,
                                             next);
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
        } else if (current_candidate.mismatches == best_candidate.mismatches &&
                   current_candidate.key.sample != barcode_key::unidentified) {
          if (current_candidate.key.sample == best_candidate.key.sample) {
            best_candidate.key.barcode = barcode_key::ambiguous;
          } else {
            best_candidate.key.sample = barcode_key::ambiguous;
            best_candidate.key.barcode = barcode_key::ambiguous;
          }
        }
      }
    }

    return best_candidate;
  } else if (next) {
    max_local_mismatches = std::min(max_global_mismatches, m_max_mismatches_r2);

    if (max_local_mismatches) {
      return lookup_with_mm(next,
                            parent,
                            max_global_mismatches,
                            max_local_mismatches,
                            nullptr);
    } else {
      return lookup(next, parent, max_global_mismatches, nullptr);
    }
  } else {
    return barcode_match(node.key, m_max_mismatches - max_global_mismatches);
  }
}

std::ostream&
operator<<(std::ostream& os, const barcode_key& value)
{
  const auto to_str_ = [](int32_t idx) -> std::string {
    switch (idx) {
      case barcode_key::unidentified:
        return "unidentified";
      case barcode_key::ambiguous:
        return "ambiguous";
      default:
        return stringify(idx);
    }
  };

  return os << "barcode_key{sample=" << to_str_(value.sample)
            << ", barcode=" << to_str_(value.barcode) << "}";
}

} // namespace adapterremoval
