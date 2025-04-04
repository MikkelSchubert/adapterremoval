// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

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
  void set_id(std::string_view id);

  /** Adds/replaces the sample (SM) tag */
  void set_sample(std::string_view name) { update_tag("SM", name); }

  /** Adds/replaces the barcode (BC) tag */
  void set_barcodes(std::string_view value) { update_tag("BC", value); }

  /** Adds/replaces the comment (CO) tag */
  void set_comment(std::string_view value) { update_tag("CO", value); }

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

/** Represents sequences used for identifying/processing a sample */
struct sample_sequences
{
  sample_sequences() = default;

  sample_sequences(dna_sequence barcode1, dna_sequence barcode2)
    : barcode_1(std::move(barcode1))
    , barcode_2(std::move(barcode2))
  {
  }

  //! Whether read groups are specified for this set of sequences
  bool has_read_group{};
  //! Read-group for this sample/barcode combination
  read_group info{};
  //! Barcode expected to be found in mate 1 reads, if any (read orientation)
  dna_sequence barcode_1{};
  //! Barcode expected to be found in mate 2 reads, if any (read orientation)
  dna_sequence barcode_2{};
  //! Adapter set with the above barcodes added
  adapter_set adapters{};
};

/** Represents a demultiplexing sample with one or more barcodes */
class sample
{
public:
  sample() { add(dna_sequence{}, dna_sequence{}); }

  explicit sample(std::string name,
                  dna_sequence barcode1,
                  dna_sequence barcode2)
    : m_name(std::move(name))
  {
    add(std::move(barcode1), std::move(barcode2));
  };

  explicit sample(std::string name, std::string barcode1, std::string barcode2)
    : sample(name, dna_sequence{ barcode1 }, dna_sequence{ barcode2 }) {};

  /** Adds a pair of barcodes in read orientation */
  void add(dna_sequence barcode1, dna_sequence barcode2);

  /** Adds barcodes in read orientation */
  void add(std::string barcode1, std::string barcode2);

  /** Assigns adapter sequences for each pair of barcodes */
  void set_adapters(const adapter_set& adapters);

  /** Assigns read groups for each pair of barcodes */
  void set_read_group(const read_group& info);

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
  std::string m_name{};
  //! Barcodes identifying this sample
  std::vector<sample_sequences> m_barcodes{};
};

/** Configuration for loading of barcode tables */
class barcode_config
{
public:
  barcode_config() = default;

  /**
   * If PE mode is enabled, barcode 1 and 2 together must be unique, otherwise
   * barcode 1 sequences alone must be unique to allow unambiguous
   * identification of samples
   */
  auto& paired_end_mode(bool value)
  {
    m_paired_end_mode = value;
    return *this;
  }

  /** Specifies if barcodes are expected in one or both orientations */
  auto& unidirectional_barcodes(bool value)
  {
    m_unidirectional_barcodes = value;
    return *this;
  }

  /** Enable or disable support for multiple barcodes for the same sample */
  auto& allow_multiple_barcodes(bool value)
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
  //! Indicates if barcode pairs can be annealed in both orientations
  bool m_unidirectional_barcodes = true;
};

/**
 * Class for handling samples for demultiplexing. The class further checks for
 * the correctness of these sequences, and detects duplicate barcode sequences /
 * pairs of sequences.
 */
class sample_set
{
public:
  /** Creates barcode set with single unnamed sample with empty barcodes */
  sample_set();
  /** Creates barcode set from set of samples. Allows multiple barcodes */
  sample_set(std::initializer_list<sample> args);

  /** Sets adapter sequences for all samples */
  void set_adapters(adapter_set adapters);

  /** Sets read group for samples using information parsed using `read_group` */
  void set_read_group(std::string_view value);

  /** Clears existing samples and loads barcodes from a TSV file */
  void load(const std::string& filename, const barcode_config& config);

  /** Convenience function to get sequences for sample / barcode pair */
  [[nodiscard]] const auto& get_sequences(const size_t sample,
                                          const size_t barcodes) const
  {
    return m_samples.at(sample).at(barcodes);
  }

  /** Returns the number of (demultiplexing) samples */
  [[nodiscard]] size_t size() const { return m_samples.size(); }

  /** Iterator over (demultiplexing) samples */
  [[nodiscard]] auto begin() const { return m_samples.begin(); }

  /** Terminal iterator over (demultiplexing) samples */
  [[nodiscard]] auto end() const { return m_samples.end(); }

  /** Returns the nth (demultiplexing) sample */
  [[nodiscard]] const auto& at(size_t n) const { return m_samples.at(n); }

  /** Returns the original, user-supplied adapter sequences */
  [[nodiscard]] const adapter_set& adapters() const { return m_adapters; }

  /** Returns special sample representing uidentified reads */
  [[nodiscard]] const auto& unidentified() const { return m_unidentified; }

private:
  /** Sets read-group for unidentified reads */
  void set_unidentified_read_group(read_group tmpl);

  /** Adds the reverse complement of barcodes for all samples, if missing */
  void add_reversed_barcodes(const barcode_config& config);

  //! Demultiplexing samples. Names and barcode pairs are both unique
  std::vector<sample> m_samples{};
  //! Special sample representing unidentified samples;
  sample m_unidentified{};
  //! User-supplied read group used to generate per-sample read-groups
  read_group m_read_group{};
  //! User-supplied adapter sequences used to generate per-barcode adapters
  adapter_set m_adapters{};
};

} // namespace adapterremoval
