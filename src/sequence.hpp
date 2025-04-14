// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include <cstddef>     // for size_t
#include <ostream>     // for ostream
#include <string>      // for string
#include <string_view> // for string_view
#include <utility>     // for pair
#include <vector>      // for vector

namespace adapterremoval {

/** Lightweight representation of a DNA sequence limited to bases ACGTN */
class dna_sequence
{
public:
  /** Creates empty sequence */
  dna_sequence() = default;

  /** Creates validated DNA sequence using the supplied base encoding rules */
  explicit dna_sequence(std::string seq);

  /** Creates validated DNA sequence using the supplied base encoding rules */
  explicit dna_sequence(std::string_view seq)
    : dna_sequence(std::string{ seq })
  {
  }

  /** Creates validated DNA sequence using the supplied base encoding rules */
  explicit dna_sequence(const char* seq)
    : dna_sequence(std::string{ seq })
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

  dna_sequence operator+(const dna_sequence& other) const;

  /** Returns true if the sequence is empty. */
  [[nodiscard]] bool empty() const noexcept { return m_sequence.empty(); }

  /** Returns the length of the sequence. */
  [[nodiscard]] size_t length() const noexcept { return m_sequence.length(); }

  [[nodiscard]] auto begin() const noexcept { return m_sequence.begin(); }

  [[nodiscard]] auto end() const noexcept { return m_sequence.end(); }

  /** Returns the reverse complements of the sequence */
  [[nodiscard]] dna_sequence reverse_complement() const;

  [[nodiscard]] const std::string& as_string() const { return m_sequence; }

private:
  //! Nucleotide sequence; contains only uppercase letters "ACGTN"
  std::string m_sequence{};
};

using sequence_pair = std::pair<dna_sequence, dna_sequence>;
using sequence_pair_vec = std::vector<sequence_pair>;

/** String literal; mostly for testing */
inline dna_sequence
operator""_dna(const char* seq, size_t length)
{
  return dna_sequence{ std::string_view{ seq, length } };
}

/** Stream operator for debugging output */
std::ostream&
operator<<(std::ostream& os, const dna_sequence& value);

} // namespace adapterremoval
