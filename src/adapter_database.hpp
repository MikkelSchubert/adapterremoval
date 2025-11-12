// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp"   // for read_mate
#include "sequence.hpp"      // for dna_sequences
#include "sequence_sets.hpp" // for adapter_set
#include <initializer_list>  // for initializer_list
#include <iosfwd>            // for ostream
#include <string_view>       // for string_view
#include <vector>            // for vector

namespace adapterremoval {

/** Class representing recommended published sequences for trimming */
class known_adapters
{
public:
  known_adapters(std::string name,
                 dna_sequence adapter_1,
                 dna_sequence adapter_2 = {},
                 bool user_provided = false);

  known_adapters(std::string name,
                 std::initializer_list<std::string_view> adapter_1,
                 std::initializer_list<std::string_view> adapter_2 = {},
                 bool user_provided = false);

  /** Returns the canonnical company / technology for these sequences */
  [[nodiscard]] const std::string& name() const { return m_name; }

  /** Adapter sequences used to trim mate 1 reads */
  [[nodiscard]] const sequence_vec& adapter_1() const { return m_adapter_1; }

  /** Adapter sequences used to trim mate 2 reads; defaults to adapter_1 */
  [[nodiscard]] const sequence_vec& adapter_2() const
  {
    return m_adapter_2.empty() ? m_adapter_1 : m_adapter_2;
  }

  [[nodiscard]] bool has_adapter_2() const { return !m_adapter_2.empty(); }

  /** Was the adapter(s) provided by the user? */
  [[nodiscard]] bool user_provided() const { return m_user_provided; }

private:
  //! The name of sequencing companies / technologies that use these adapters
  std::string m_name;
  //! Adapter 1 sequences
  sequence_vec m_adapter_1;
  //! Optional adapter 2 sequences; if unspecified, adapter 1 sequences are used
  sequence_vec m_adapter_2;
  //! Was this adapter set user-provided; used to favor informative records
  bool m_user_provided;
};

/** A lookup match, either exact or closest. May be empty for no match */
struct identified_adapter
{
  //! The canonical name of the adapter sequence
  std::string name{};
  //! The published adapter sequence matching the detected adapter sequence
  dna_sequence sequence{};
  //! Is this adapter expected to be found in mate 1 or mate 2 reads?
  read_mate mate = read_mate::_1;

  /** Returns true if the two objects have identical properties */
  [[nodiscard]] bool operator==(const identified_adapter& other) const noexcept;

  /** Stream operator for debugging output */
  friend std::ostream& operator<<(std::ostream& os,
                                  const identified_adapter& match);
};

/** A pair of identified adapters; sequences may be empty for no matches */
using identified_adapter_pair =
  std::pair<identified_adapter, identified_adapter>;

class adapter_database
{
public:
  /** Create empty adapter database */
  adapter_database() = default;

  /** Adds a set of (user-provided) adapters to the database */
  void add(const adapter_set& adapters);
  /** Adds all known adapters to the database */
  void add_known();

  /** Returns the number of known / user-provided adapter sets */
  [[nodiscard]] size_t size() const { return m_adapters.size(); }

  /** Returns start iterator for known / user provided adapters */
  [[nodiscard]] auto begin() const { return m_adapters.begin(); }

  /** Returns end iterator for known / user provided adapters */
  [[nodiscard]] auto end() const { return m_adapters.end(); }

  /** Returns the nth known / user provided adapters */
  [[nodiscard]] const known_adapters& at(size_t i) const
  {
    return m_adapters.at(i);
  }

  /**
   * Returns the closest sequence/pair of sequences matching the arguments.
   * Attempted to return sequences combinations that make sense, favoring
   * sequences in the expected reads and from the same sources.
   */
  [[nodiscard]] identified_adapter_pair identify_closest(
    const dna_sequence& seq_1,
    const dna_sequence& seq_2) const;

  /** Like `identify_closest`, but assume/require that exact matches exist */
  [[nodiscard]] identified_adapter_pair identify_exact(
    const dna_sequence& seq_1,
    const dna_sequence& seq_2) const;

  /** Specifies the output format when exporting the database */
  enum class export_fmt
  {
    /** Tab separated values */
    tsv,
    /** Pretty printed JSON list of records */
    json
  };

  /** Return table of known adapters as TSV or JSON */
  [[nodiscard]] static std::string export_known(export_fmt format);

private:
  //! List of known sequences used for adapter trimming
  std::vector<known_adapters> m_adapters{};
};

/** Writes "read 1" or "read 2" */
std::ostream&
operator<<(std::ostream& os, read_mate mate);

} // namespace adapterremoval
