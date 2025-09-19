// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>

#include "adapters.hpp"      // declarations
#include "commontypes.hpp"   // for read_mate
#include "debug.hpp"         // for AR_FAIL, AR_REQUIRE
#include "sequence.hpp"      // for dna_sequence
#include "sequence_sets.hpp" // for adapter_set
#include "strutils.hpp"
#include "utilities.hpp" // for underlying_value
#include <algorithm>     // for find
#include <limits>
#include <sstream> // for ostringstream

namespace adapterremoval {

namespace {

const std::array ADAPTER_DATABASE = {
  // Sourced from
  // https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
  // https://web.archive.org/web/20240819030322/https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
  known_adapters{
    {
      "Illumina TruSeq",
      "Illumina TruSeq DNA and RNA CD",
      "Illumina TruSeq DNA and RNA UD",
      "Illumina TruSeq DNA HT",
      "Illumina TruSeq DNA LT",
      "Illumina TruSeq single index",
      "Illumina TruSight Oncology 500",
      "Illumina TruSight Oncology ctDNA",
      "Illumina TruSight RNA Pan-Cancer Panel",
      "Illumina TruSight Tumor 170",
    },
    { "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" },
    { "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" },
  },
  // https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/FrontPages/AdapterSeq.htm
  // https://web.archive.org/web/20250920194612/https://support-docs.illumina.com/SHARE/AdapterSequences/1000000002694_21_illumina_adapter_sequences.pdf
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
  known_adapters{
    {
      "Illumina Stranded mRNA",
      "Illumina Stranded Total RNA",
    },
    {
      "ACTGTCTCTTATACACATCT",
    },
  },
  known_adapters{
    { "Illumina TruSeq Ribo Profile" },
    {
      "AGATCGGAAGAGCACACGTCT",
    },
  },
  // Sourced from
  // https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/miseq-sample-sheet-quick-ref-guide-15028392-j.pdf
  known_adapters{
    { "Illumina MiSeq" },
    { "TGGAATTCTCGGGTGCCAAGGC" },
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
  // Sourced from
  // https://www.mdpi.com/2311-553X/5/4/49
  known_adapters{
    { "Diagenode CATS small RNA-seq Kit" },
    { "GATCGGAAGAGCACACGTCTG" },
  },
  // Sourced from
  // https://faqs.lexogen.com/faq/what-is-the-adapter-sequence-i-need-to-use-for-t-2
  // https://web.archive.org/web/20250719105852/https://faqs.lexogen.com/faq/what-is-the-adapter-sequence-i-need-to-use-for-t-2
  known_adapters{
    { "Lexogene Small RNA-Seq" },
    { "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" },
  }
#ifdef TAKARA_SMARTER
  // Not included due to high probability of false positives
  known_adapters{
    { "Takara SMARTer smRNA-Seq" },
    { "AAAAAAAAAA" },
  },
#endif
};

//! Simplifies checking boolean conditions that may be false or true
const std::array FALSE_AND_TRUE{ false, true };

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

std::pair<size_t, dna_sequence>
find_best_match(const std::vector<dna_sequence>& items, const dna_sequence& seq)
{
  std::pair<size_t, dna_sequence> match{ std::numeric_limits<size_t>::max(),
                                         {} };

  const std::string_view seq1_view = seq.as_string();
  for (const auto& it : items) {
    const auto length = std::min(it.length(), seq.length());
    const std::string_view seq2_view = seq.as_string();

    const auto dist =
      levenshtein(seq1_view.substr(0, length), seq2_view.substr(0, length));
    if (dist < match.first) {
      match.first = dist;
      match.second = it;
    }
  }

  return match;
}

using best_match = std::pair<size_t, identified_adapter>;

best_match
find_best_match(const known_adapters& items,
                const dna_sequence& seq,
                read_mate mate)
{
  switch (mate) {
    case read_mate::_1: {
      auto [diff, match] = find_best_match(items.adapter_1(), seq);

      return { diff, { items.source(), std::move(match), read_mate::_1 } };
    }
    case read_mate::_2: {
      auto [diff, match] = find_best_match(items.adapter_2(), seq);

      return { diff, { items.source(), std::move(match), read_mate::_2 } };
    }
    default:
      AR_FAIL("invalid sequence_orientation");
  }
}

best_match
select_se(const std::vector<known_adapters>& candidates,
          const dna_sequence& seq,
          read_mate mate)
{
  std::array<read_mate, 2> mates{ read_mate::_1, read_mate::_2 };
  if (mates.front() != mate) {
    std::swap(mates.front(), mates.back());
  }

  best_match best{ std::numeric_limits<size_t>::max(), {} };

  // prefer adapters in the expected reads over swapped adapters
  for (const auto swapped_mate : mates) {
    // prefer known adapters over user-provided adapters, for the extra info
    for (bool user_provided : FALSE_AND_TRUE) {
      for (const auto& it : candidates) {
        if (it.user_provided() == user_provided) {
          auto current = find_best_match(it, seq, swapped_mate);

          if (current.first == 0) {
            return current;
          } else if (current.first < best.first) {
            std::swap(best, current);
          }
        }
      }
    }
  }

  return best;
}

} // namespace

std::pair<identified_adapter, identified_adapter>
adapter_database::identify(const dna_sequence& seq_1,
                           const dna_sequence& seq_2) const
{
  if (seq_1.empty() && seq_2.empty()) {
    return {};
  } else if (seq_2.empty()) {
    return { select_se(m_adapters, seq_1, read_mate::_1).second, {} };
  } else if (seq_1.empty()) {
    return { {}, select_se(m_adapters, seq_2, read_mate::_2).second };
  }

  best_match best_1{ std::numeric_limits<size_t>::max(), {} };
  best_match best_2{ std::numeric_limits<size_t>::max(), {} };

  // prefer adapters in the expected reads over swapped adapters
  for (const bool swapped : FALSE_AND_TRUE) {
    const auto& swapped_seq1 = swapped ? seq_2 : seq_1;
    const auto& swapped_seq2 = swapped ? seq_1 : seq_2;
    const auto mate_1 = swapped ? read_mate::_2 : read_mate::_1;
    const auto mate_2 = swapped ? read_mate::_1 : read_mate::_2;

    // prefer known adapters over user-provided adapters, for the extra info
    for (const bool user_provided : FALSE_AND_TRUE) {
      for (const auto& it : m_adapters) {
        if (it.user_provided() == user_provided) {
          auto match_1 = find_best_match(it.adapter_1(), swapped_seq1);
          auto match_2 = find_best_match(it.adapter_2(), swapped_seq2);

          if (match_1.first == 0 && match_2.first == 0) {
            return {
              identified_adapter{ it.source(), swapped_seq1, mate_1 },
              identified_adapter{ it.source(), swapped_seq2, mate_2 },
            };
          } else if (match_1.first + match_2.first <
                     best_1.first + best_2.first) {
            best_1 = {
              match_1.first,
              identified_adapter{ it.source(), swapped_seq1, mate_1 }
            };
            best_2 = {
              match_2.first,
              identified_adapter{ it.source(), swapped_seq2, mate_2 }
            };
          }
        }
      }
    }
  }

  // fall back to arbitrary picks among all candidates
  auto match_1 = select_se(m_adapters, seq_1, read_mate::_1);
  auto match_2 = select_se(m_adapters, seq_2, read_mate::_2);

  if (match_1.first + match_2.first < best_1.first + best_2.first) {
    return { match_1.second, match_2.second };
  } else {
    return { best_1.second, best_2.second };
  }
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
