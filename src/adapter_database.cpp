// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "adapter_database.hpp" // declarations
#include "alignment.hpp"        // for simple_alignment
#include "commontypes.hpp"      // for read_mate
#include "debug.hpp"            // for AR_FAIL, AR_REQUIRE
#include "json.hpp"             // for json
#include "sequence.hpp"         // for dna_sequence
#include "sequence_sets.hpp"    // for adapter_set
#include "strutils.hpp"         // for log_escape
#include "utilities.hpp"        // for underlying_value
#include <algorithm>            // for swap, none_of
#include <sstream>              // for ostringstream
#include <string_view>          // for string_view
#include <utility>              // for move, pair
#include <vector>               // for vector

namespace adapterremoval {

namespace {

const std::array ADAPTER_DATABASE = {
  // Sourced from
  // https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
  // https://web.archive.org/web/20240819030322/https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
  known_adapters{
    "Illumina TruSeq",
    { "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" },
    { "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" },
  },
  // https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/FrontPages/AdapterSeq.htm
  // https://web.archive.org/web/20250920194612/https://support-docs.illumina.com/SHARE/AdapterSequences/1000000002694_21_illumina_adapter_sequences.pdf
  known_adapters{
    "Illumina AmpliSeq",
    { "CTGTCTCTTATACACATCT" },
  },
  known_adapters{
    "Illumina DNA PCR-Free Prep, Tagmentation",
    {
      "CTGTCTCTTATACACATCT",
      "ATGTGTATAAGAGACA",
    },
  },
  known_adapters{
    "Illumina ScriptSeq",
    { "AGATCGGAAGAGCACACGTCTGAAC" },
    { "AGATCGGAAGAGCGTCGTGTAGGGA" },
  },
  known_adapters{
    "Illumina TruSeq Small RNA",
    { "TGGAATTCTCGGGTGCCAAGG" },
  },
  // Example data = PRJEB8559 / ERR760535
  known_adapters{
    "Illumina Nextera Mate Pair",
    {
      "CTGTCTCTTATACACATCT",
      "AGATGTGTATAAGAGACAG",
    },
  },
  known_adapters{
    "Illumina Stranded mRNA",
    {
      "ACTGTCTCTTATACACATCT",
    },
  },
  known_adapters{
    "Illumina TruSeq Ribo Profile",
    {
      "AGATCGGAAGAGCACACGTCT",
    },
  },
  // Sourced from
  // https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/miseq-sample-sheet-quick-ref-guide-15028392-j.pdf
  known_adapters{
    "Illumina MiSeq",
    { "TGGAATTCTCGGGTGCCAAGGC" },
  },
  // Sourced from
  // https://en.mgitech.cn/Download/download_file/id/71
  // https://web.archive.org/web/20240528231846/https://en.mgitech.cn/Download/download_file/id/71
  known_adapters{
    "MGI Tech DNBSEQ",
    { "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA" },
    { "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG" },
  },
  // Sourced from
  // https://www.qiagen.com/us/resources/faq?id=f12b85b4-df4f-43b5-9e82-a4fd0ddbdcc0
  // https://web.archive.org/web/20240908090002/https://www.qiagen.com/us/resources/faq?id=f12b85b4-df4f-43b5-9e82-a4fd0ddbdcc0
  // Example data = SRR8557389
  known_adapters{
    "QIAGEN QIAseq miRNA",
    { "AACTGTAGGCACCATCAAT" },
  },
  // Example data = PRJNA162397 / SRR491337
  known_adapters{
    "Life Technologies SOLiD small RNA",
    { "CGCCTTGGCCGTACAGCAG" },
  },
  // Sourced from
  // https://www.mdpi.com/2311-553X/5/4/49
  known_adapters{
    "Diagenode CATS small RNA-seq Kit",
    { "GATCGGAAGAGCACACGTCTG" },
  },
  // Sourced from
  // https://faqs.lexogen.com/faq/what-is-the-adapter-sequence-i-need-to-use-for-t-2
  // https://web.archive.org/web/20250719105852/https://faqs.lexogen.com/faq/what-is-the-adapter-sequence-i-need-to-use-for-t-2
  known_adapters{
    "Lexogen Small RNA-Seq",
    { "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" },
  }
#ifdef TAKARA_SMARTER
  // Not included due to high probability of false positives
  known_adapters{
    "Takara SMARTer smRNA-Seq",
    { "AAAAAAAAAA" },
  },
#endif
};

//! Simplifies checking boolean conditions that may be false or true
const std::array FALSE_AND_TRUE{ false, true };

} // namespace

////////////////////////////////////////////////////////////////////////////////
// known_adapters

known_adapters::known_adapters(std::string name,
                               dna_sequence adapter_1,
                               dna_sequence adapter_2,
                               bool user_provided)
  : m_name(std::move(name))
  , m_adapter_1()
  , m_adapter_2()
  , m_user_provided(user_provided)
{
  AR_REQUIRE(!m_name.empty());
  m_adapter_1.emplace_back(std::move(adapter_1));

  if (!adapter_2.empty()) {
    m_adapter_2.emplace_back(std::move(adapter_2));
  }
}

known_adapters::known_adapters(
  std::string name,
  std::initializer_list<std::string_view> adapter_1,
  std::initializer_list<std::string_view> adapter_2,
  bool user_provided)
  : m_name(std::move(name))
  , m_adapter_1(adapter_1.begin(), adapter_1.end())
  , m_adapter_2(adapter_2.begin(), adapter_2.end())
  , m_user_provided(user_provided)
{
  AR_REQUIRE(!m_name.empty());

  // read 1 adapters are required, but read 2 adapters may be missing / implicit
  AR_REQUIRE(!m_adapter_1.empty());

  // if any source or sequence is empty, then I probably made a mistake
  const auto empty = [](const std::string_view& it) { return it.empty(); };
  AR_REQUIRE(std::none_of(adapter_1.begin(), adapter_1.end(), empty));
  AR_REQUIRE(std::none_of(adapter_2.begin(), adapter_2.end(), empty));
}

////////////////////////////////////////////////////////////////////////////////
// identified_adapter

bool
identified_adapter::operator==(const identified_adapter& other) const noexcept
{
  return (this->name == other.name) && (this->sequence == other.sequence) &&
         (this->mate == other.mate);
}

std::ostream&
operator<<(std::ostream& os, const identified_adapter& match)
{
  return os << "identified_adapter{source=" << log_escape(match.name)
            << ", sequence=" << match.sequence << ", mate=" << match.mate
            << "}";
}

////////////////////////////////////////////////////////////////////////////////
// adapter_database

namespace {

bool
includes(const sequence_vec& seqs, const dna_sequence& seq)
{
  return std::find(seqs.begin(), seqs.end(), seq) != seqs.end();
}

bool
is_better_candidate(const alignment_info& match_1,
                    const size_t length_1,
                    const alignment_info& match_2,
                    const size_t length_2)
{
  return match_1.is_better_than(match_2) ||
         (match_1 == match_2 && length_1 < length_2);
}

std::pair<alignment_info, dna_sequence>
best_match(const sequence_vec& candidates, const dna_sequence& seq)
{
  alignment_info best_aln;
  dna_sequence best_match;

  if (!seq.empty()) {
    for (const auto& it : candidates) {
      auto alignment = simple_alignment(it, seq);
      if (is_better_candidate(alignment,
                              it.length(),
                              best_aln,
                              best_match.length())) {
        best_aln = alignment;
        best_match = it;
      }
    }
  }

  return std::pair{ best_aln, best_match };
}

identified_adapter
possible_match(const std::string& source,
               const dna_sequence& seq,
               read_mate mate)
{
  return seq.empty() ? identified_adapter{}
                     : identified_adapter{ source, seq, mate };
}

std::pair<identified_adapter, identified_adapter>
exact_match_pe(const std::vector<known_adapters>& candidates,
               const dna_sequence& seq_1,
               const dna_sequence& seq_2)
{
  // prefer adapters in the expected reads over swapped adapters
  for (const bool swapped : FALSE_AND_TRUE) {
    // prefer known adapters over user-provided adapters, for the extra info
    for (const bool user_provided : FALSE_AND_TRUE) {
      for (const auto& it : candidates) {
        if (it.user_provided() == user_provided) {
          const auto& adapters_1 = swapped ? it.adapter_2() : it.adapter_1();
          const auto& adapters_2 = swapped ? it.adapter_1() : it.adapter_2();

          const auto includes_1 = seq_1.empty() || includes(adapters_1, seq_1);
          const auto includes_2 = seq_2.empty() || includes(adapters_2, seq_2);

          if (includes_1 && includes_2) {
            const auto mate_1 = swapped ? read_mate::_2 : read_mate::_1;
            const auto mate_2 = swapped ? read_mate::_1 : read_mate::_2;

            return std::pair{
              possible_match(it.name(), seq_1, mate_1),
              possible_match(it.name(), seq_2, mate_2),
            };
          }
        }
      }
    }
  }

  return {};
}

std::pair<alignment_info, std::pair<identified_adapter, identified_adapter>>
best_match_pe(const std::vector<known_adapters>& candidates,
              const dna_sequence& seq_1,
              const dna_sequence& seq_2)
{
  // missing/empty sequences are allowed, to simplify doing lookups
  if (seq_1.empty() && seq_2.empty()) {
    return {};
  }

  identified_adapter best_seq_1;
  identified_adapter best_seq_2;
  alignment_info best_aln;
  size_t best_len = 0;

  // prefer adapters in the expected reads over swapped adapters
  for (const bool swapped : FALSE_AND_TRUE) {
    // prefer known adapters over user-provided adapters, for the extra info
    for (const bool user_provided : FALSE_AND_TRUE) {
      for (const auto& it : candidates) {
        const auto& adapters_1 = swapped ? it.adapter_2() : it.adapter_1();
        const auto& adapters_2 = swapped ? it.adapter_1() : it.adapter_2();

        if (it.user_provided() == user_provided) {
          auto [aln_1, match_1] = best_match(adapters_1, seq_1);
          auto [aln_2, match_2] = best_match(adapters_2, seq_2);
          auto merged_aln = aln_1.add(aln_2);
          auto merged_len = match_1.length() + match_2.length();

          if (is_better_candidate(merged_aln, merged_len, best_aln, best_len)) {
            const auto mate_1 = swapped ? read_mate::_2 : read_mate::_1;
            const auto mate_2 = swapped ? read_mate::_1 : read_mate::_2;

            best_seq_1 = possible_match(it.name(), match_1, mate_1);
            best_seq_2 = possible_match(it.name(), match_2, mate_2);
            best_aln = merged_aln;
            best_len = merged_len;
          }
        }
      }
    }
  }

  return std::pair{ best_aln, std::pair{ best_seq_1, best_seq_2 } };
}

} // namespace

void
adapter_database::add(const adapter_set& adapters)
{
  size_t nth = 1;
  for (const auto& it : adapters) {
    std::ostringstream ss;
    ss << "User adapters #" << nth++;

    m_adapters.emplace_back(ss.str(),
                            it.first,
                            // convert alignment orientation to read orientation
                            it.second.reverse_complement(),
                            true);
  }
}

void
adapter_database::add_known()
{
  m_adapters.insert(m_adapters.end(),
                    ADAPTER_DATABASE.begin(),
                    ADAPTER_DATABASE.end());
}

std::pair<identified_adapter, identified_adapter>
adapter_database::identify_closest(const dna_sequence& seq_1,
                                   const dna_sequence& seq_2) const
{
  // best match with seq_1 and seq_2 belonging to the same source
  auto [best_aln, best_seqs] = best_match_pe(m_adapters, seq_1, seq_2);
  auto [best_seq_1, best_seq_2] = best_seqs;
  auto best_len = best_seq_1.sequence.length() + best_seq_2.sequence.length();

  // fall back to arbitrary picks among all sources if we have both sequences
  if (!seq_1.empty() && !seq_2.empty()) {
    auto [aln_1, seqs_1] = best_match_pe(m_adapters, seq_1, {});
    auto [aln_2, seqs_2] = best_match_pe(m_adapters, {}, seq_2);
    auto merged_aln = aln_1.add(aln_2);
    auto merged_len =
      seqs_1.first.sequence.length() + seqs_2.second.sequence.length();

    if (is_better_candidate(merged_aln, merged_len, best_aln, best_len)) {
      return std::pair{ seqs_1.first, seqs_2.second };
    }
  }
  return std::pair{ best_seq_1, best_seq_2 };
}

std::pair<identified_adapter, identified_adapter>
adapter_database::identify_exact(const dna_sequence& seq_1,
                                 const dna_sequence& seq_2) const
{
  std::pair<identified_adapter, identified_adapter> match;

  // missing/empty sequences are allowed, to simplify doing lookups
  if (seq_1.empty() && seq_2.empty()) {
    return match;
  }

  match = exact_match_pe(m_adapters, seq_1, seq_2);
  if (!match.first.sequence.empty() || !match.second.sequence.empty()) {
    return match;
  }

  // fall back to arbitrary picks among all candidates
  if (!seq_1.empty() && !seq_2.empty()) {
    match.first = exact_match_pe(m_adapters, seq_1, {}).first;
    match.second = exact_match_pe(m_adapters, {}, seq_2).second;

    if (!match.first.sequence.empty() && !match.second.sequence.empty()) {
      return match;
    }
  }

  // This function is intended to be used to pick among known sequences
  AR_FAIL("attempted to look up unknown adapter sequences");
}

std::string
adapter_database::export_known(export_fmt format)
{
  auto to_vec = [&](const sequence_vec& seqs) {
    string_vec sequences;
    for (const auto& seq : seqs) {
      sequences.emplace_back(seq.as_string());
    }

    return sequences;
  };

  auto adapters = ADAPTER_DATABASE;
  std::sort(adapters.begin(), adapters.end(), [](const auto& a, const auto& b) {
    return a.name() < b.name();
  });

  if (format == export_fmt::tsv) {
    std::ostringstream os;

    os << "Name\tAdapter1\tAdapter2";
    for (const auto& it : adapters) {
      os << "\n" << it.name() << "\t" << join_text(to_vec(it.adapter_1()), ",");

      if (it.has_adapter_2()) {
        os << "\t" << join_text(to_vec(it.adapter_2()), ",");
      } else {
        os << "\t.";
      }
    }

    return os.str();
  } else if (format == export_fmt::json) {
    json_list values;
    for (const auto& it : adapters) {
      auto entry = values.dict();
      entry->str("name", it.name());
      entry->str_vec("adapter1", to_vec(it.adapter_1()));

      if (it.has_adapter_2()) {
        entry->str_vec("adapter2", to_vec(it.adapter_2()));
      } else {
        entry->null("adapter2");
      }
    }

    return values.to_string();
  } else {
    AR_FAIL("invalid dump_format");
  }
}

std::ostream&
operator<<(std::ostream& os, const read_mate mate)
{
  switch (mate) {
    case read_mate::_1:
    case read_mate::_2:
      return os << "read " << underlying_value(mate);
    default:                                   // GCOVR_EXCL_LINE
      AR_FAIL("invalid sequence_orientation"); // GCOVR_EXCL_LINE
  }
}

} // namespace adapterremoval
