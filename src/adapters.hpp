// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
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
  known_adapters(std::vector<std::string> sources,
                 std::initializer_list<std::string_view> adapter_1,
                 std::initializer_list<std::string_view> adapter_2 = {},
                 bool user_provided = false);

  /** Returns the canonnical company / technology for these sequences */
  [[nodiscard]] std::string_view source() const { return m_sources.front(); }

  /** Returns sources using the following adapter trimming sequences  */
  [[nodiscard]] const std::vector<std::string>& sources() const
  {
    return m_sources;
  }

  /** Adapter sequences used to trim mate 1 reads */
  [[nodiscard]] const std::vector<dna_sequence>& adapter_1() const
  {
    return m_adapter_1;
  }

  /** Adapter sequences used to trim mate 2 reads; defaults to adapter_1 */
  [[nodiscard]] const std::vector<dna_sequence>& adapter_2() const
  {
    return m_adapter_2.empty() ? m_adapter_1 : m_adapter_2;
  }

  /** Was the adapter(s) provided by the user */
  [[nodiscard]] bool user_provided() const { return m_user_provided; }

private:
  //! The name of sequencing companies / technologies that use these adapters
  const std::vector<std::string> m_sources;
  //! Adapter 1 sequences
  const std::vector<dna_sequence> m_adapter_1;
  //! Optional adapter 2 sequences; if unspecified, adapter 1 sequences are used
  const std::vector<dna_sequence> m_adapter_2;
  //! Was this adapter set user-provided; used to favor informative records
  const bool m_user_provided;
};

struct identified_adapter
{
  identified_adapter() = default;

  identified_adapter(std::string_view source_,
                     dna_sequence sequence_,
                     read_mate mate_)
    : source(source_)
    , sequence(std::move(sequence_))
    , mate(mate_)
  {
  }

  std::string source{};
  dna_sequence sequence{};
  read_mate mate = read_mate::_1;
};

using identified_adapter_pair =
  std::pair<identified_adapter, identified_adapter>;

class adapter_database
{
public:
  /** Create database containing only predefined adapter sequences */
  adapter_database();
  /** Create database with pre-defined and user-defined adapter sequences */
  explicit adapter_database(const adapter_set& user_adapters);

  [[nodiscard]] size_t size() const { return m_adapters.size(); }

  [[nodiscard]] auto begin() const { return m_adapters.begin(); }

  [[nodiscard]] auto end() const { return m_adapters.end(); }

  [[nodiscard]] const known_adapters& at(size_t i) const
  {
    return m_adapters.at(i);
  }

  /**
   * Returns the closest sequence/pair of sequences matching the arguments. It
   * is attempted to return sequences combinations that make sense, favoring
   * sequences in the expected reads and from the same sources.
   */
  [[nodiscard]] identified_adapter_pair identify(
    const dna_sequence& seq_1,
    const dna_sequence& seq_2) const;

private:
  //! List of known sequences used for adapter trimming
  std::vector<known_adapters> m_adapters;
};

/** Writes "read 1" or "read 2" */
std::ostream&
operator<<(std::ostream& os, read_mate mate);

} // namespace adapterremoval
