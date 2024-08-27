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

#include "fastq_enc.hpp" // for FASTQ_ENCODING_33, MATE_SEPARATOR
#include <array>         // for array
#include <string>        // for string
#include <string_view>   // for string_view
#include <vector>        // for vector

namespace adapterremoval {

/** Lightweight represention of a DNA sequence limited to bases ACGTN */
class dna_sequence
{
public:
  /** Creates empty sequence */
  dna_sequence() = default;

  /** Creates validated DNA sequence using the supplied base encoding rules */
  explicit dna_sequence(std::string seq,
                        const fastq_encoding& encoding = FASTQ_ENCODING_33)
    : m_sequence(std::move(seq))
  {
    encoding.process_nucleotides(m_sequence);
  }

  /** Creates validated DNA sequence using the supplied base encoding rules */
  explicit dna_sequence(std::string_view seq,
                        const fastq_encoding& encoding = FASTQ_ENCODING_33)
    : dna_sequence(std::string{ seq }, encoding)
  {
  }

  /** Creates validated DNA sequence using the supplied base encoding rules */
  explicit dna_sequence(const char* seq,
                        const fastq_encoding& encoding = FASTQ_ENCODING_33)
    : dna_sequence(std::string{ seq }, encoding)
  {
  }

  /** Returns true IFF the sequences are identical. **/
  bool operator==(const dna_sequence& other) const noexcept
  {
    return m_sequence == other.m_sequence;
  }

  /** Returns true IFF this sequence is less than the other sequence. **/
  bool operator<(const dna_sequence& other) const noexcept
  {
    return m_sequence < other.m_sequence;
  }

  dna_sequence operator+(const dna_sequence& other) const
  {
    dna_sequence result;
    result.m_sequence.reserve(length() + other.length());
    result.m_sequence += m_sequence;
    result.m_sequence += other.m_sequence;

    return result;
  }

  operator std::string_view() const noexcept { return m_sequence; }

  /** Returns the length of the sequence. */
  [[nodiscard]] size_t length() const noexcept { return m_sequence.length(); }

  [[nodiscard]] auto begin() const noexcept { return m_sequence.begin(); }

  [[nodiscard]] auto end() const noexcept { return m_sequence.end(); }

  /** Returns the reverse complements of the sequence */
  [[nodiscard]] dna_sequence reverse_complement() const
  {
    // Lookup table for complementary bases based only on the last 4 bits
    static constexpr std::array<char, 16> complements{ '-', 'T', '-', 'G',
                                                       'A', '-', '-', 'C',
                                                       '-', '-', '-', '-',
                                                       '-', '-', 'N', '-' };

    dna_sequence rc;
    rc.m_sequence.reserve(length());
    for (auto it = m_sequence.rbegin(); it != m_sequence.rend(); ++it) {
      rc.m_sequence.push_back(complements[*it & 0xf]);
    }

    return rc;
  }

private:
  //! Nucleotide sequence; contains only uppercase letters "ACGTN"
  std::string m_sequence{};
};

using sequence_pair = std::pair<dna_sequence, dna_sequence>;
using sequence_pair_vec = std::vector<sequence_pair>;

} // namespace adapterremoval
