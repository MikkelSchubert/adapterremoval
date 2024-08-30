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
#pragma once

#include "barcode_table.hpp"
#include "sequence.hpp"     // for for dna_sequence
#include <cstddef>          // for size_t
#include <initializer_list> // for initializer_list
#include <string>           // for string
#include <string_view>      // for string_view
#include <vector>           // for vector

namespace adapterremoval {

using string_view_pair = std::pair<std::string_view, std::string_view>;

/** Contains SAM/BAM read-group information */
class read_group
{
public:
  read_group();

  /**
   * Parses a read-group string in the form "ID:1\tSM:sample" (optionally
   * including a leading "@RG\t"). Throws std::invalid_argument if the value
   * is invalid.
   */
  explicit read_group(std::string_view value);

  /** Returns the read-group ID for use in per-read 'RG' tags */
  [[nodiscard]] std::string_view id() const { return m_id; }

  /** Returns the full @RG header, not including a trailing new-line */
  [[nodiscard]] std::string_view header() const { return m_header; }

  /** Adds/replaces the barcode (ID) tag */
  void set_id(std::string_view id) { update_tag("ID", id); }

  /** Adds/replaces the sample (SM) tag */
  void set_sample(std::string_view name) { update_tag("SM", name); }

  /** Adds/replaces the barcode (BC) tag */
  void set_barcodes(std::string_view value) { update_tag("BC", value); }

private:
  /** Updates or adds the specified tag; sets `m_id` if key is `ID` */
  void update_tag(std::string_view key, std::string_view value);

  //! The full read_group header, including leading `@RG\t`
  std::string m_header{};
  //! Value mapping reads (via `RG:Z:${ID}`) to the @RG header
  std::string m_id{};
};

/**
 * Class for loading/handling adapter adapter sequences.
 *
 * Adapter sequences are found in one of two orientations:
 *  - Read orientation, corresponding to the sequence in input fastq reads
 *  - Alignment orientation, corresponding to the orientation used during
 *    sequence alignment. For the mate 1 adapter, this is read orientation,
 *    but for the mate 2 adapter this is the reverse complement.
 */
class adapter_set
{
public:
  /** Initialize empty adapter list. */
  adapter_set() = default;

  /** Initializes with adapters in read orientation */
  adapter_set(std::initializer_list<string_view_pair> args);

  /** Adds a pair of adapters to the set in read orientation */
  void add(dna_sequence adapter1, dna_sequence adapter2);

  /** Adds a pair of adapters to the set in read orientation */
  void add(std::string adapter1, std::string adapter2);

  /** Generate new adapter set with these barcodes (in read orientation) */
  [[nodiscard]] adapter_set add_barcodes(const dna_sequence& barcode1,
                                         const dna_sequence& barcode2) const;

  /**
   * Loads adapters in read orientation from a TSV file, throwing on failure.
   * Two adapter sequences are expected if 'paired_end_mode' is set.
   */
  void load(const std::string& filename, bool paired_end_mode);

  /** Returns the number of adapters/adapter pairs added/loaded */
  [[nodiscard]] size_t size() const { return m_adapters.size(); }

  /** Iterator over adapter sequences in alignment orientation */
  [[nodiscard]] auto begin() const { return m_adapters.begin(); }

  /** Terminal iterator over adapter sequences in alignment orientation */
  [[nodiscard]] auto end() const { return m_adapters.end(); }

  [[nodiscard]] const auto& at(size_t n) const { return m_adapters.at(n); }

  /** Returns the adapters in read orientation */
  [[nodiscard]] sequence_pair_vec to_read_orientation() const;

private:
  //! Adapter sequences in alignment orientation
  sequence_pair_vec m_adapters{};
};

/** Represents a demultiplexing sample with one or more barcodes */
class sample
{
public:
  explicit sample(std::string name,
                  dna_sequence barcode1,
                  dna_sequence barcode2)
    : m_name(std::move(name))
    , m_barcodes()
  {
    add(std::move(barcode1), std::move(barcode2));
  };

  explicit sample(std::string name, std::string barcode1, std::string barcode2)
    : sample(name, dna_sequence{ barcode1 }, dna_sequence{ barcode2 }){};

  /** Adds a pair of barcodes in read orientation */
  void add(dna_sequence adapter1, dna_sequence adapter2);

  /** Adds barcodes in read orientation */
  void add(std::string barcode1, std::string barcode2);

  /** Returns the unique name of this sample */
  [[nodiscard]] const auto& name() const { return m_name; }

  /** Returns the number of barcode sequences loaded */
  [[nodiscard]] size_t size() const { return m_barcodes.size(); }

  /** Iterator over adapter sequences in alignment orientation */
  [[nodiscard]] auto begin() const { return m_barcodes.begin(); }

  /** Terminal iterator over adapter sequences in alignment orientation */
  [[nodiscard]] auto end() const { return m_barcodes.end(); }

  /** Returns the nth barcode / pair of barcodes */
  [[nodiscard]] const auto& at(size_t n) const { return m_barcodes.at(n); }

private:
  //! Unique name associated with this sample
  std::string m_name;
  //! Barcodes identifying this sample
  sequence_pair_vec m_barcodes;
};

/**
 * Class for handling samples for  demultiplexing. The class further checks for
 * the correctness of these sequences, and detects duplicate barcode sequences /
 * pairs of sequences.
 */
class barcode_set
{
public:
  barcode_set() = default;

  /**
   * In SE mode the first barcode must be unique, in PE mode the individual
   * sequences may be non-unique as long as the combination is unique
   */
  void set_paired_end_mode(bool b) { m_paired_end_mode = b; }

  /** Enable or disable support for multiple barcodes for the same sample */
  void set_allow_multiple_barcodes(bool b) { m_allow_multiple_barcodes = b; }

  /** Adds a unnamed sample with empty barcodes */
  void add_default_sample();

  /**
   * Adds a pair of barcodes for the given sample; intended for small / test
   * tables as the table is fully validated after each addition
   **/
  void add(std::string name, std::string barcode1, std::string barcode2);

  /**
   * Adds reverse (complemented) barcodes for barcode for each sample; if the
   * reverse (complemented) barcode is already present for a sample, then it is
   * ignored. Duplicates across samples trigger an exception.
   */
  void add_reversed_barcodes();

  /** Loads barcodes from a TSV file, throwing on error */
  void load(const std::string& filename);

  /** Returns the number of (demultiplexing) samples */
  [[nodiscard]] size_t size() const { return m_samples.size(); }

  /** Iterator over (demultiplexing) samples */
  [[nodiscard]] auto begin() const { return m_samples.begin(); }

  /** Terminal iterator over (demultiplexing) samples */
  [[nodiscard]] auto end() const { return m_samples.end(); }

  /** Returns the nth (demultiplexing) sample */
  [[nodiscard]] const auto& at(size_t n) const { return m_samples.at(n); }

private:
  //! Demultiplexing samples. Names and barcode pairs are both unique
  std::vector<sample> m_samples{};
  //! Whether running in paired or single end mode; is used to determine whether
  //! or not b can be uniquely identified from the barcodes provided
  bool m_paired_end_mode = false;
  //! cates whether or not multiple barcodes/barcode pairs can bused b
  //! identify the same sample.
  bool m_allow_multiple_barcodes = false;
};

} // namespace adapterremoval
