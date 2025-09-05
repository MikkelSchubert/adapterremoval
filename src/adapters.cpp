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
#include "adapters.hpp"      // declarations
#include "commontypes.hpp"   // for read_mate
#include "debug.hpp"         // for AR_FAIL, AR_REQUIRE
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for adapter_set
#include "utilities.hpp"     // for underlying_value
#include <algorithm>         // for find
#include <sstream>           // for ostringstream

namespace adapterremoval {

#if 0
// TODO
// Review: https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/FrontPages/AdapterSeq.htm
#endif

namespace {

const std::array<known_adapters, 9> ADAPTER_DATABASE = {
  // Sourced from
  // https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
  // https://web.archive.org/web/20240819030322/https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
  known_adapters{
    {
      "Illumina TruSeq",
      "Illumina TruSeq DNA CD",
      "Illumina TruSeq DNA HT",
      "Illumina TruSeq DNA LT",
      "Illumina TruSeq single index",
    },
    { "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" },
    { "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" },
  },
  // https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/FrontPages/AdapterSeq.htm
  known_adapters{
    {
      "Illumina AmpliSeq",
      "Illumina DNA Prep (M) Tagmentation",
      "Illumina DNA Prep with Enrichment (S) Tagmentation",
      "Illumina RNA Prep with Enrichment, Ligation",
      "Illumina Stranded mRNA Prep, Ligation",
      "Illumina Stranded Total RNA Prep, Ligation with Ribo-Zero Plus",
      "Illumina Nextera DNA Flex",
      "Illumina Nextera DNA",
      "Illumina Nextera Enrichment",
      "Illumina Nextera Flex for Enrichment",
      "Illumina Nextera Rapid Capture Enrichment",
      "Illumina Nextera XT",
      "Illumina TruSight Enrichment",
      "Illumina TruSight HLA",
      "Illumina TruSight Rapid Capture Enrichment",
    },
    { "CTGTCTCTTATACACATCT" },
  },
  known_adapters{
    { "Illumina DNA PCR-Free Prep, Tagmentation" },
    {
      "CTGTCTCTTATACACATCT",
      "ATGTGTATAAGAGACA",
    },
  },
  known_adapters{
    {
      "Illumina ScriptSeq",
      "Illumina TruSeq DNA Methylation",
    },
    { "AGATCGGAAGAGCACACGTCTGAAC" },
    { "AGATCGGAAGAGCGTCGTGTAGGGA" },
  },
  known_adapters{
    { "Illumina TruSeq Small RNA" },
    { "TGGAATTCTCGGGTGCCAAGG" },
  },
  // Example data = PRJEB8559 / ERR760535
  known_adapters{
    { "Illumina Nextera Mate Pair" },
    {
      "CTGTCTCTTATACACATCT",
      "AGATGTGTATAAGAGACAG",
    },
  },
  // Sourced from
  // https://en.mgitech.cn/Download/download_file/id/71
  // https://web.archive.org/web/20240528231846/https://en.mgitech.cn/Download/download_file/id/71
  known_adapters{
    {
      "MGI Tech BGISEQ",
      "MGI Tech DNBSEQ",
      "MGI Tech MGISEQ",
    },
    { "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA" },
    { "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG" },
  },
  // Sourced from
  // https://www.qiagen.com/us/resources/faq?id=f12b85b4-df4f-43b5-9e82-a4fd0ddbdcc0
  // https://web.archive.org/web/20240908090002/https://www.qiagen.com/us/resources/faq?id=f12b85b4-df4f-43b5-9e82-a4fd0ddbdcc0
  // Example data = SRR8557389
  known_adapters{
    { "QIAGEN QIAseq miRNA" },
    { "AACTGTAGGCACCATCAAT" },
  },
  // Example data = PRJNA162397 / SRR491337
  known_adapters{
    { "Life Technologies SOLiD small RNA" },
    { "CGCCTTGGCCGTACAGCAG" },
  },
};

//! Simplifies checking boolean conditions that may be false or true
const std::array<bool, 2> FALSE_AND_TRUE{ false, true };

} // namespace

known_adapters::known_adapters(
  std::vector<std::string> sources,
  std::initializer_list<std::string_view> adapter_1,
  std::initializer_list<std::string_view> adapter_2,
  bool user_provided)
  : m_sources(std::move(sources))
  , m_adapter_1(adapter_1.begin(), adapter_1.end())
  , m_adapter_2(adapter_2.begin(), adapter_2.end())
  , m_user_provided(user_provided)
{
  AR_REQUIRE(m_sources.size());
}

adapter_database::adapter_database()
  : adapter_database::adapter_database(adapter_set{})
{
}

adapter_database::adapter_database(const adapter_set& user_adapters)
  : m_adapters(ADAPTER_DATABASE.begin(), ADAPTER_DATABASE.end())
{
  size_t nth = 1;
  for (const auto& it : user_adapters) {
    std::ostringstream ss;
    ss << "User adapters #" << nth++;

    m_adapters.push_back(known_adapters({ ss.str() },
                                        { it.first.as_string() },
                                        { it.second.as_string() },
                                        true));
  }
}

namespace {

bool
includes(const std::vector<dna_sequence>& items, const dna_sequence& seq)
{
  return std::find(items.begin(), items.end(), seq) != items.end();
}

bool
includes(const known_adapters& items, const dna_sequence& seq, read_mate mate)
{
  switch (mate) {
    case read_mate::_1:
      return includes(items.adapter_1(), seq);
    case read_mate::_2:
      return includes(items.adapter_2(), seq);
    default:
      AR_FAIL("invalid sequence_orientation");
  }
}

read_mate
swap_mate(read_mate mate)
{
  switch (mate) {
    case read_mate::_1:
      return read_mate::_2;
    case read_mate::_2:
      return read_mate::_1;
    default:
      AR_FAIL("invalid sequence_orientation");
  }
}

identified_adapter
select_se(const std::vector<known_adapters>& candidates,
          const dna_sequence& seq,
          read_mate mate)
{
  // prefer adapters in the expected reads over swapped adapters
  for (bool swapped : FALSE_AND_TRUE) {
    auto swapped_mate = swapped ? mate : swap_mate(mate);

    // prefer known adapters over user-provided adapters, for the extra info
    for (bool user_provided : FALSE_AND_TRUE) {
      for (const auto& it : candidates) {
        if (it.user_provided() == user_provided) {
          if (includes(it, seq, swapped_mate)) {
            return { it.source(), seq, swapped_mate };
          }
        }
      }
    }
  }

  AR_FAIL("Attempted to retrieve unknown adapter sequence");
}

} // namespace

std::pair<identified_adapter, identified_adapter>
adapter_database::identify(const dna_sequence& seq1,
                           const dna_sequence& seq2) const
{
  if (seq1.empty() && seq2.empty()) {
    return {};
  } else if (seq2.empty()) {
    return { select_se(m_adapters, seq1, read_mate::_1), {} };
  } else if (seq1.empty()) {
    return { {}, select_se(m_adapters, seq2, read_mate::_2) };
  }

  // prefer adapters in the expected reads over swapped adapters
  for (bool flipped : FALSE_AND_TRUE) {
    auto flipped_seq1 = flipped ? seq2 : seq1;
    auto flipped_seq2 = flipped ? seq1 : seq2;
    auto orientation = flipped ? read_mate::_2 : read_mate::_1;

    // prefer known adapters over user-provided adapters, for the extra info
    for (bool user_provided : FALSE_AND_TRUE) {
      for (const auto& it : m_adapters) {
        if (it.user_provided() == user_provided) {
          if (includes(it.adapter_1(), flipped_seq1) &&
              includes(it.adapter_2(), flipped_seq2)) {

            return {
              identified_adapter{ it.source(), flipped_seq1, orientation },
              identified_adapter{ it.source(),
                                  flipped_seq2,
                                  swap_mate(orientation) },
            };
          }
        }
      }
    }
  }

  // fall back to arbitrary picks among all candidates
  return { select_se(m_adapters, seq1, read_mate::_1),
           select_se(m_adapters, seq2, read_mate::_2) };
}

std::ostream&
operator<<(std::ostream& os, const read_mate mate)
{
  switch (mate) {
    case read_mate::_1:
    case read_mate::_2:
      return os << "read " << underlying_value(mate);
    default:
      AR_FAIL("invalid sequence_orientation");
  }
}

} // namespace adapterremoval

#if 0

// https://www.mdpi.com/2311-553X/5/4/49
TGGAATTCTCGGGTGCCAAGG PerkinElmer, NEXTflex Small RNA
TGGAATTCTCGGGTGCCAAGG SeqMatic, TailorMix miRNA
TGGAATTCTCGGGTGCCAAGG TriLink, CleanTag Small RNA
AGATCGGAAGAGCACACGTCT New England Biolabs, NEBNext Multiplex Small RNA (SRR10713837)
S Lexogen, Small RNA-Seq Library Prep Kit
GATCGGAAGAGCACACGTCTG Diagenode CATS small RNA-seq Kit *

// Include, but comment out with reasson
Takara 	SMARTer smRNA-Seq Kit for Illumina 	Version 	5'-AAAAAAAAAA-3'

    // https://github.com/joey0214/adapter4srna/blob/master/README.md
    // https://www.biorxiv.org/content/biorxiv/early/2019/07/14/702456/DC3/embed/media-3.pdf?download=true

#endif
