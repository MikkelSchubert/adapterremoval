// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp"  // for barcode_orientation, ...
#include "read_group.hpp"   // for read_group
#include "sequence.hpp"     // for for dna_sequence
#include <cstddef>          // for size_t
#include <initializer_list> // for initializer_list
#include <iosfwd>           // for ostream
#include <string>           // for string
#include <string_view>      // for string_view
#include <vector>           // for vector

namespace adapterremoval {

class line_reader_base;

/** Maps a name to the corresponding barcode_table_orientation enum */
barcode_table_orientation
parse_table_orientation(std::string_view value);

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
  void add(dna_sequence adapter1, const dna_sequence& adapter2);

  /** Generate new adapter set with these barcodes (in read orientation) */
  [[nodiscard]] adapter_set add_barcodes(const dna_sequence& barcode1,
                                         const dna_sequence& barcode2) const;

  /**
   * Loads adapters in read orientation, clearing existing adapters. Two adapter
   * sequences are expected if 'paired_end_mode' is set.
   */
  void load(const std::string& filename, bool paired_end_mode);
  /**
   * Loads adapters in read orientation, clearing existing adapters. Two adapter
   * sequences are expected if 'paired_end_mode' is set.
   */
  void load(line_reader_base& reader, bool paired_end_mode);

  /** Returns the number of adapters/adapter pairs added/loaded */
  [[nodiscard]] size_t size() const { return m_adapters.size(); }

  /** Returns true if the adapter set is empty */
  [[nodiscard]] bool empty() const { return m_adapters.empty(); }

  /** Iterator over adapter sequences in alignment orientation */
  [[nodiscard]] auto begin() const { return m_adapters.begin(); }

  /** Terminal iterator over adapter sequences in alignment orientation */
  [[nodiscard]] auto end() const { return m_adapters.end(); }

  /** Returns the nth adapter sequences */
  [[nodiscard]] const auto& at(size_t n) const { return m_adapters.at(n); }

  /** Returns the adapters in read orientation */
  [[nodiscard]] sequence_pair_vec to_read_orientation() const;

  /** Returns true if the adapters (including ordering) are identical  */
  [[nodiscard]] bool operator==(const adapter_set& other) const;

  /** Creates debug representation of an adapter set */
  friend std::ostream& operator<<(std::ostream& os, const adapter_set& value);

private:
  //! Adapter sequences in alignment orientation
  sequence_pair_vec m_adapters{};
};

/** Represents sequences used for identifying/processing a sample */
struct sample_sequences
{
  sample_sequences() = default;

  sample_sequences(dna_sequence barcode_1,
                   dna_sequence barcode_2,
                   barcode_orientation orientation);

  //! Whether read groups are specified for this set of sequences
  bool has_read_group{};
  //! Optional read-group for this sample/barcode combination
  read_group read_group_{};
  //! Barcode expected to be found in mate 1 reads, if any (read orientation)
  dna_sequence barcode_1{};
  //! Barcode expected to be found in mate 2 reads, if any (read orientation)
  dna_sequence barcode_2{};
  //! User specified orientation of the barcodes
  barcode_orientation orientation = barcode_orientation::unspecified;

  /** Returns the sample adapters. Aborts if the adapters are uninitialized */
  [[nodiscard]] const adapter_set& adapters() const;

  /** Replaces the current adapters and clears the uninitialized flag, if set */
  void set_adapters(adapter_set as);

  /** Mark adapters as uninitialized; see `adapters()` */
  void flag_uninitialized_adapters();

  /** Returns true if all members are identical */
  [[nodiscard]] bool operator==(const sample_sequences& other) const;

  /** Creates debug representation of sample sequences */
  friend std::ostream& operator<<(std::ostream& os,
                                  const sample_sequences& value);

private:
  //! Adapter set with the above barcodes added
  adapter_set m_adapters{};
  //! The adapters have yet to be determined, and accessing those are prohibited
  bool m_uninitialized_adapters = false;
};

/** Represents a (demultiplexing) sample with one or more barcodes */
class sample
{
public:
  /** Creates basic unnamed sample without barcodes */
  sample();

  /** Creates named sample with the specified barcodes */
  sample(std::string name,
         dna_sequence barcode1,
         dna_sequence barcode2,
         barcode_orientation orientation);

  /** Adds a pair of barcodes in read orientation */
  void add_barcodes(dna_sequence barcode1,
                    dna_sequence barcode2,
                    barcode_orientation orientation);

  /**
   * Assigns adapter sequences for each pair of barcodes. If the adapters have
   * been flagged as uninitialized, then this clears that flag
   */
  void set_adapters(const adapter_set& adapters);

  /** Assigns read groups for each pair of barcodes */
  void set_read_group(const read_group& read_group_);

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

  /**
   * Mark adapters as uninitialized, i.e. before automatic selection is done.
   * When set, adapters may not be accessed until set_adapters` has been called
   */
  void flag_uninitialized_adapters();

  /** Returns true if name and barcodes are identical */
  [[nodiscard]] bool operator==(const sample& other) const;

  /** Creates debug representation of a sample */
  friend std::ostream& operator<<(std::ostream& os, const sample& value);

private:
  //! Unique name associated with this sample
  std::string m_name{};
  //! Barcodes identifying this sample
  std::vector<sample_sequences> m_barcodes{};
};

/** Configuration for loading of barcode tables */
class barcode_config
{
public:
  constexpr barcode_config() = default;

  /**
   * If PE mode is enabled, barcode 1 and 2 together must be unique, otherwise
   * barcode 1 sequences must be unique to allow samples to be identified.
   */
  constexpr auto& paired_end_mode(bool value = true) noexcept
  {
    m_paired_end_mode = value;
    return *this;
  }

  /** Specifies if barcodes are expected in one or both orientations */
  constexpr auto& orientation(barcode_table_orientation value) noexcept
  {
    m_orientation = value;
    return *this;
  }

  /** Enable or disable support for multiple barcodes for the same sample */
  constexpr auto& allow_multiple_barcodes(bool value = true) noexcept
  {
    m_allow_multiple_barcodes = value;
    return *this;
  }

private:
  friend class sample_set;

  //! Whether running in paired or single end mode; is used to determine whether
  //! or not samples can be uniquely identified from the barcodes provided
  bool m_paired_end_mode = false;
  //! Indicates if multiple barcodes/barcode pairs are allowed per sample
  bool m_allow_multiple_barcodes = false;
  //! Indicates the orientation of barcodes in the table
  barcode_table_orientation m_orientation =
    barcode_table_orientation::unspecified;
};

/**
 * Class for handling user-specified samples for demultiplexing, in addition to
 * an 'unidentified' sample representing reads that could not be assigned to a
 * sample. In non-demultiplexing mode, the set contains a single, unnamed sample
 * with an optional read-group, and no barcode sequences.
 */
class sample_set
{
public:
  /** Creates sample set with single unnamed sample with empty barcodes */
  sample_set();
  /** Creates sample set from  lines representing a barcode table */
  sample_set(std::initializer_list<std::string_view> lines,
             barcode_config config = {});

  /**
   * Sets adapter sequences for all samples, generating unique sequences for
   * each set of adapters for each set of barcodes. If the adapters have been
   * flagged as uninitialized, then calling `set_adapters` clears that flag
   */
  void set_adapters(adapter_set adapters);

  /** Parses read group string and updates existing samples */
  void set_read_group(std::string_view value);

  /**
   * Overrides table with samples specified in a whitespace separated table,
   * containing a name column, and one or two barcode columns. Samples are
   * updated with the current read group and adapters set.
   */
  void load(const std::string& filename, const barcode_config& config);
  /** See `load(const std::string&, const barcode_config&)` */
  void load(line_reader_base& reader, const barcode_config& config);

  /** Returns the number of (demultiplexing) samples */
  [[nodiscard]] size_t size() const { return m_samples.size(); }

  /** Iterator over (demultiplexing) samples */
  [[nodiscard]] auto begin() const { return m_samples.begin(); }

  /** Terminal iterator over (demultiplexing) samples */
  [[nodiscard]] auto end() const { return m_samples.end(); }

  /** Returns the nth (demultiplexing) sample */
  [[nodiscard]] const auto& at(size_t n) const { return m_samples.at(n); }

  /** Returns the original, user-supplied adapter sequences */
  [[nodiscard]] const adapter_set& adapters() const;

  /** Returns the original, user-supplied adapter sequences */
  [[nodiscard]] const read_group& readgroup() const { return m_read_group; }

  /** Returns special sample representing uidentified reads */
  [[nodiscard]] const auto& unidentified() const { return m_unidentified; }

  /** Returns the vector of samples */
  [[nodiscard]] const auto& samples() const { return m_samples; }

  /**
   * Mark adapters as uninitialized, i.e. before automatic selection is done.
   * When set, adapters may not be accessed until set_adapters` has been called
   */
  void flag_uninitialized_adapters();

  /** Access partially initialized adapters; for adapter selection only */
  [[nodiscard]] const adapter_set& uninitialized_adapters() const;

  /** Creates debug representation of a sample set */
  friend std::ostream& operator<<(std::ostream& os, const sample_set& value);

private:
  /** Sets read-group for unidentified reads */
  void set_unidentified_read_group(read_group tmpl);

  //! Demultiplexing samples. Names and barcode pairs are both unique
  std::vector<sample> m_samples{};
  //! Special sample representing unidentified samples;
  sample m_unidentified{};
  //! User-supplied read group used to generate per-sample read-groups
  read_group m_read_group{};
  //! User-supplied adapter sequences used to generate per-barcode adapters
  adapter_set m_adapters{};
  //! The adapters have yet to be determined, and read-access is prohibited
  bool m_uninitialized_adapters = false;
};

////////////////////////////////////////////////////////////////////////////////
// Stream operators for debugging output

std::ostream&
operator<<(std::ostream& os, const barcode_orientation& value);

std::ostream&
operator<<(std::ostream& os, const barcode_table_orientation& value);

} // namespace adapterremoval
