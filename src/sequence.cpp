// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#include "sequence.hpp"  // declarations
#include "fastq_enc.hpp" // for FASTQ_ENCODING_33
#include "strutils.hpp"  // for log_escape
#include <array>         // for array

namespace adapterremoval {

dna_sequence::dna_sequence(std::string seq)
  : m_sequence(std::move(seq))
{
  FASTQ_ENCODING_SAM.process_nucleotides(m_sequence);
}

dna_sequence
dna_sequence::operator+(const dna_sequence& other) const
{
  dna_sequence result;
  result.m_sequence.reserve(length() + other.length());
  result.m_sequence += m_sequence;
  result.m_sequence += other.m_sequence;

  return result;
}

dna_sequence
dna_sequence::reverse_complement() const
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

std::ostream&
operator<<(std::ostream& os, const dna_sequence& value)
{
  return os << "dna_sequence{" << log_escape(value.as_string()) << "}";
}

} // namespace adapterremoval
