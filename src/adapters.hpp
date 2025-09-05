/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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

  [[nodiscard]] std::pair<identified_adapter, identified_adapter> identify(
    const dna_sequence& seq1,
    const dna_sequence& seq2) const;

private:
  //! List of known sequences used for adapter trimming
  std::vector<known_adapters> m_adapters;
};

/** Writes "read 1" or "read 2" */
std::ostream&
operator<<(std::ostream& os, read_mate mate);

} // namespace adapterremoval
