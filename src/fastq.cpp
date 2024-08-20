/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "fastq.hpp"
#include "debug.hpp"      // for AR_REQUIRE, AR_FAIL
#include "errors.hpp"     // for fastq_error
#include "linereader.hpp" // for line_reader_base
#include <algorithm>      // for reverse, count, max, min
#include <cmath>          // for log10, pow
#include <numeric>        // for accumulate
#include <sstream>        // for ostringstream
#include <string_view>    // for string_view

namespace adapterremoval {

namespace {

std::vector<double>
init_phred_to_p_values()
{
  std::vector<double> result;
  for (size_t i = PHRED_SCORE_MIN; i <= PHRED_SCORE_MAX; ++i) {
    result.push_back(std::pow(10.0, static_cast<double>(i) / -10.0));
  }

  return result;
}

const std::vector<double> g_phred_to_p = init_phred_to_p_values();

enum class read_mate
{
  unknown,
  mate_1,
  mate_2,
};

struct mate_info
{
  std::string_view desc() const
  {
    switch (mate) {
      case read_mate::unknown:
        return "unknown";
      case read_mate::mate_1:
        return "mate 1";
      case read_mate::mate_2:
        return "mate 2";
      default:
        AR_FAIL("Invalid mate in mate_info::desc");
    }
  }

  //! Read name without mate number or meta-data
  std::string_view name{};
  //! Which mate in a pair, if identified
  read_mate mate = read_mate::unknown;
  //! Position of the separator character in the header (if any)
  size_t sep_pos = std::string::npos;
};

mate_info
get_mate_info(const fastq& read, char mate_separator)
{
  const std::string_view header = read.header();

  size_t pos = header.find_first_of(' ');
  if (pos == std::string::npos) {
    pos = header.length();
  }

  mate_info info;
  if (pos >= 2 && header.at(pos - 2) == mate_separator) {
    const char digit = header.at(pos - 1);

    if (digit == '1') {
      info.mate = read_mate::mate_1;
      pos -= 2;
      info.sep_pos = pos;
    } else if (digit == '2') {
      info.mate = read_mate::mate_2;
      pos -= 2;
      info.sep_pos = pos;
    }
  }

  info.name = header.substr(0, pos);
  return info;
}

size_t
count_poly_x_tail(const std::string& m_sequence,
                  const char nucleotide,
                  const size_t min_length)
{
  // Maximum number of sequential mismatches
  const size_t max_seq_mismatches = 2;
  // Number of called bases required per mismatch (via fastp)
  const size_t min_bases_per_mismatch = 8;

  //! Number of bases in the alignment to trim, excluding leading mismatches
  size_t n_trim = 0;
  //! Number of bases in the alignment
  size_t n_bases = 0;
  //! Number of uncalled bases (Ns) in the alignment
  size_t n_uncalled = 0;
  //! Number of mismatches in the alignment
  size_t n_mismatches = 0;
  //! Current number of sequential mismatches in the alignment
  size_t n_seq_mismatches = 0;

  for (auto it = m_sequence.rbegin(); it != m_sequence.rend(); ++it) {
    n_bases++;

    if (*it == nucleotide) {
      n_trim = n_bases;
      n_seq_mismatches = 0;
    } else if (*it == 'N') {
      n_uncalled++;
      // Trailing Ns are allowed only after a match
      if (!n_seq_mismatches) {
        n_trim = n_bases;
      }
    } else {
      n_mismatches++;
      n_seq_mismatches++;
      if (n_seq_mismatches > max_seq_mismatches ||
          n_mismatches > std::max(min_length, n_bases - n_uncalled) /
                           min_bases_per_mismatch) {
        // The final mismatch is not counted as part of the alignment
        n_bases--;
        break;
      }
    }
  }

  if (n_bases - n_uncalled >= min_length) {
    return n_trim;
  }

  return 0;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// fastq

fastq::fastq()
  : m_header()
  , m_sequence()
  , m_qualities()
{
}

fastq::fastq(std::string header,
             std::string sequence,
             std::string qualities,
             const fastq_encoding& encoding)
  : m_header(std::move(header))
  , m_sequence(std::move(sequence))
  , m_qualities(std::move(qualities))
{
  if (m_qualities.length() != m_sequence.length()) {
    throw fastq_error(
      "invalid FASTQ record; sequence/quality length does not match");
  }

  post_process(encoding);
}

fastq::fastq(std::string header, std::string sequence)
  : m_header(std::move(header))
  , m_sequence(std::move(sequence))
  , m_qualities(std::string(m_sequence.length(), '!'))
{
  post_process(FASTQ_ENCODING_33);
}

bool
fastq::operator==(const fastq& other) const
{
  return (m_header == other.m_header) && (m_sequence == other.m_sequence) &&
         (m_qualities == other.m_qualities);
}

std::string_view
fastq::name(const char mate_separator) const
{
  std::string_view header = m_header;
  const size_t pos = header.find_first_of(' ');
  if (pos != std::string::npos) {
    header = header.substr(0, pos);
  }

  if (mate_separator && header.size() > 1 &&
      (header.back() == '1' || header.back() == '2') &&
      header.at(header.length() - 2) == mate_separator) {
    header = header.substr(0, header.length() - 2);
  }

  return header;
}

size_t
fastq::count_ns() const
{
  return static_cast<size_t>(
    std::count(m_sequence.begin(), m_sequence.end(), 'N'));
}

namespace {

/**
 * Calculate the absolute sequence complexity score, under the assumption that
 * the sequence does not contain Ns. Should the sequence contain Ns, then this
 * algorithm would overestimate the sequence complexity, and therefore returns
 * -1 to indicate failure.
 */
int
fast_calculate_complexity(const std::string& sequence)
{
  // The last base is not checked in the loop below
  if (sequence.back() == 'N') {
    return -1;
  }

  const size_t length = sequence.length() - 1;
  size_t i = 0;
  int score = 0;

  // Fixed block sizes allows gcc/clang to optimize the loop
  const size_t BLOCK_SIZE = 16;
  for (; i + BLOCK_SIZE < length; i += BLOCK_SIZE) {
    for (size_t j = 0; j < BLOCK_SIZE; ++j, ++i) {
      if (sequence[i] != sequence[i + 1]) {
        score++;
      }

      if (sequence[i] == 'N') {
        return -1;
      }
    }
  }

  for (; i < length; ++i) {
    if (sequence[i] != sequence[i + 1]) {
      score++;
    }

    if (sequence[i] == 'N') {
      return -1;
    }
  }

  return score;
}

} // namespace

double
fastq::complexity() const
{
  if (m_sequence.length() < 2) {
    return 0.0;
  }

  // Try to use unrolled/vectorized algorithm
  int score = fast_calculate_complexity(m_sequence);

  if (score < 0) {
    // If the sequence contains Ns then use the slower calculation, that does
    // not treat Ns as distinct bases and thereby does not inflate the score
    char prev = 'N';
    for (const auto nuc : m_sequence) {
      if (nuc != 'N' && nuc != prev) {
        prev = nuc;
        score++;
      }
    }
  }

  return std::max(0.0, score / static_cast<double>(m_sequence.length() - 1));
}

double
fastq::mean_quality() const
{
  int64_t total = 0;
  AR_REQUIRE(!m_qualities.empty());
  for (const auto c : m_qualities) {
    total += c - PHRED_OFFSET_MIN;
  }

  return static_cast<double>(total) / static_cast<double>(m_qualities.length());
}

fastq::ntrimmed
fastq::trim_trailing_bases(const bool trim_ns,
                           char low_quality,
                           const bool preserve5p)
{
  low_quality += PHRED_OFFSET_MIN;
  auto is_quality_base = [&](size_t i) {
    return m_qualities.at(i) > low_quality &&
           (!trim_ns || m_sequence.at(i) != 'N');
  };

  size_t right_exclusive = 0;
  for (size_t i = m_sequence.length(); i; --i) {
    if (is_quality_base(i - 1)) {
      right_exclusive = i;
      break;
    }
  }

  size_t left_inclusive = 0;
  for (size_t i = 0; !preserve5p && i < right_exclusive; ++i) {
    if (is_quality_base(i)) {
      left_inclusive = i;
      break;
    }
  }

  return trim_sequence_and_qualities(left_inclusive, right_exclusive);
}

//! Calculates the size of the sliding window for quality trimming given a
//! read length and a user-defined window-size (fraction or whole number).
size_t
calculate_winlen(const size_t read_length, double window_size)
{
  if (window_size < 1.0) {
    window_size = window_size * static_cast<double>(read_length);
  }

  const auto winlen = static_cast<size_t>(window_size);
  if (winlen == 0 || winlen > read_length) {
    return read_length;
  }

  return winlen;
}

fastq::ntrimmed
fastq::trim_windowed_bases(const bool trim_ns,
                           char low_quality,
                           const double window_size,
                           const bool preserve5p)
{
  AR_REQUIRE(window_size >= 0.0);
  if (m_sequence.empty()) {
    return {};
  }

  low_quality += PHRED_OFFSET_MIN;
  auto is_quality_base = [&](size_t i) {
    return m_qualities.at(i) > low_quality &&
           (!trim_ns || m_sequence.at(i) != 'N');
  };

  const size_t winlen = calculate_winlen(length(), window_size);
  long running_sum =
    std::accumulate(m_qualities.begin(), m_qualities.begin() + winlen, 0);

  size_t left_inclusive = std::string::npos;
  size_t right_exclusive = std::string::npos;
  for (size_t offset = 0; offset + winlen <= length(); ++offset) {
    const long running_avg = running_sum / static_cast<long>(winlen);

    // We trim away low quality bases and Ns from the start of reads,
    // **before** we consider windows.
    if (left_inclusive == std::string::npos && is_quality_base(offset) &&
        running_avg > low_quality) {
      left_inclusive = offset;
    }

    if (left_inclusive != std::string::npos &&
        (running_avg <= low_quality || offset + winlen == length())) {
      right_exclusive = offset;
      while (right_exclusive < length() && is_quality_base(right_exclusive)) {
        right_exclusive++;
      }

      break;
    }

    running_sum -= m_qualities.at(offset);
    if (offset + winlen < length()) {
      running_sum += m_qualities.at(offset + winlen);
    }
  }

  if (left_inclusive == std::string::npos) {
    // No starting window found. Trim all bases starting from start.
    return trim_sequence_and_qualities(length(), length());
  } else if (preserve5p) {
    left_inclusive = 0;
  }

  AR_REQUIRE(right_exclusive != std::string::npos);
  return trim_sequence_and_qualities(left_inclusive, right_exclusive);
}

fastq::ntrimmed
fastq::mott_trimming(const double error_limit, const bool preserve5p)
{
  AR_REQUIRE(error_limit >= 0 && error_limit <= 1);

  size_t left_inclusive_temp = 0;
  size_t left_inclusive = 0;
  size_t right_exclusive = 0;

  double error_sum = 0.0;
  double error_sum_max = 0.0;

  for (size_t i = 0; i < length(); i++) {
    char phred = m_qualities.at(i) - PHRED_OFFSET_MIN;

    // Reduce weighting of very low-quality bases (inspired by seqtk) and
    // normalize Ns. The latter is not expected to matter for most data, but may
    // be relevant for some old/weird data and masked FASTQ reads.
    if (phred < 3 || m_sequence.at(i) == 'N') {
      phred = 3;
    }

    error_sum += error_limit - g_phred_to_p.at(phred);

    if (error_sum < 0.0) {
      // End of current segment (if any)
      left_inclusive_temp = i + 1;
      error_sum = 0;
    } else if (error_sum > error_sum_max) {
      // Extend best segment, possibly replacing the previous candidate
      left_inclusive = left_inclusive_temp;
      right_exclusive = i + 1;
      error_sum_max = error_sum;
    }
  }

  return trim_sequence_and_qualities(preserve5p ? 0 : left_inclusive,
                                     right_exclusive);
}

std::pair<char, size_t>
fastq::poly_x_trimming(const std::string& nucleotides, size_t min_length)
{
  size_t best_count = 0;
  char best_nucleotide = 'N';
  if (m_sequence.length() >= min_length && !nucleotides.empty()) {
    // Looping over all nucleotides ended up faster than a single pass algorithm
    for (const auto nucleotide : nucleotides) {
      const auto count = count_poly_x_tail(m_sequence, nucleotide, min_length);
      if (count > best_count) {
        best_nucleotide = nucleotide;
        best_count = count;
      }
    }

    truncate(0, length() - best_count);
  }

  return { best_nucleotide, best_count };
}

void
fastq::truncate(size_t pos, size_t len)
{
  AR_REQUIRE(pos == 0 || pos <= length());

  if (pos) {
    m_sequence.erase(0, pos);
    m_qualities.erase(0, pos);
  }

  if (len < length()) {
    m_sequence.erase(len);
    m_qualities.erase(len);
  }
}

void
fastq::reverse_complement()
{
  std::reverse(m_sequence.begin(), m_sequence.end());
  std::reverse(m_qualities.begin(), m_qualities.end());

  // Lookup table for complementary bases based only on the last 4 bits
  static const char complements[] = "-T-GA--C------N-";
  for (auto& nuc : m_sequence) {
    nuc = complements[nuc & 0xf];
  }
}

void
fastq::add_prefix_to_name(const std::string& prefix)
{
  if (!prefix.empty()) {
    m_header.insert(0, prefix);
  }
}

bool
fastq::read(line_reader_base& reader, const fastq_encoding& encoding)
{
  if (read_unsafe(reader)) {
    post_process(encoding);
    return true;
  }

  return false;
}

bool
fastq::read_unsafe(line_reader_base& reader)
{
  do {
    if (!reader.getline(m_header)) {
      // End of file; terminate gracefully
      return false;
    }
  } while (m_header.empty());

  if (m_header.size() < 2 || m_header.at(0) != '@') {
    throw fastq_error("Malformed or empty FASTQ header");
  } else {
    // FIXME: Erasing the '@' and then re-adding it later is pointless work
    m_header.erase(0, 1);
  }

  if (!reader.getline(m_sequence)) {
    throw fastq_error("partial FASTQ record; cut off after header");
  } else if (m_sequence.empty()) {
    throw fastq_error("sequence is empty");
  }

  // Most of the time this will only be '+' and not require an allocation
  std::string line;
  if (!reader.getline(line)) {
    throw fastq_error("partial FASTQ record; cut off after sequence");
  } else if (line.empty() || line.at(0) != '+') {
    throw fastq_error("FASTQ record lacks separator character (+)");
  }

  if (!reader.getline(m_qualities)) {
    throw fastq_error("partial FASTQ record; cut off after separator");
  } else if (m_qualities.length() != m_sequence.length()) {
    throw fastq_error("sequence/quality lengths do not match");
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
// Public helper functions

char
fastq::p_to_phred_33(double p)
{
  // Lowest possible error rate representable is '~' (~5e-10)
  const auto min_p = std::max(5e-10, p);
  const auto raw_score = static_cast<int>(-10.0 * std::log10(min_p));

  return std::min<int>(PHRED_OFFSET_MAX, PHRED_OFFSET_MIN + raw_score);
}

char
fastq::guess_mate_separator(const std::vector<fastq>& reads_1,
                            const std::vector<fastq>& reads_2)
{
  AR_REQUIRE(reads_1.size() == reads_2.size());

  // Commonly used characters
  const std::string candidates = "/.:";

  for (auto candidate : candidates) {
    auto it_1 = reads_1.begin();
    auto it_2 = reads_2.begin();

    bool any_failures = false;
    while (it_1 != reads_1.end()) {
      const auto info1 = get_mate_info(*it_1++, candidate);
      const auto info2 = get_mate_info(*it_2++, candidate);

      if (info1.name != info2.name) {
        any_failures = true;
        break;
      }

      const auto mate_1 = info1.mate;
      const auto mate_2 = info2.mate;

      if (mate_1 != read_mate::unknown || mate_2 != read_mate::unknown) {
        if (mate_1 == mate_2) {
          // This could be valid data that just happens to include a known
          // mate separator in the name. But this could also happen if the
          // same reads are used for both mate 1 and mate 2, so we cannot
          // safely guess.
          return 0;
        } else if (mate_1 != read_mate::mate_1 || mate_2 != read_mate::mate_2) {
          // The mate separator seems to be correct, but the mate information
          // does not match: One mate is missing information or the order is
          // wrong. Return the identified separator and raise an error later.
          return candidate;
        }
      }
    }

    if (!any_failures) {
      return candidate;
    }
  }

  return 0;
}

void
fastq::normalize_paired_reads(fastq& mate1, fastq& mate2, char mate_separator)
{
  if (mate1.length() == 0 || mate2.length() == 0) {
    throw fastq_error("Pair contains empty reads");
  }

  const auto info1 = get_mate_info(mate1, mate_separator);
  const auto info2 = get_mate_info(mate2, mate_separator);

  if (info1.name != info2.name) {
    std::ostringstream error;
    error << "Pair contains reads with mismatching names:\n"
          << " - '" << info1.name << "'\n"
          << " - '" << info2.name << "'";

    if (info1.mate == read_mate::unknown || info2.mate == read_mate::unknown) {
      error << "\n\nNote that AdapterRemoval by determines the mate "
               "numbers as the digit found at the end of the read name, "
               "if this is preceded by";

      if (mate_separator) {
        error << "the character '" << mate_separator << "'";
      } else {
        error << "a character such as '/'";
      }

      error << "; if these data makes use of a different character to "
               "separate the mate number from the read name, then you "
               "will need to set the --mate-separator command-line "
               "option to the appropriate character.";
    }

    throw fastq_error(error.str());
  }

  if (info1.mate != read_mate::unknown || info2.mate != read_mate::unknown) {
    if (info1.mate != read_mate::mate_1 || info2.mate != read_mate::mate_2) {
      std::ostringstream error;
      error << "Inconsistent mate numbering; please verify data:\n"
            << "\nRead 1 identified as " << info1.desc() << ": " << mate1.name()
            << "\nRead 2 identified as " << info2.desc() << ": "
            << mate2.name();

      throw fastq_error(error.str());
    }

    AR_REQUIRE(info1.sep_pos == info2.sep_pos);
    mate1.m_header.at(info1.sep_pos) = MATE_SEPARATOR;
    mate2.m_header.at(info2.sep_pos) = MATE_SEPARATOR;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Private helper functions

void
fastq::post_process(const fastq_encoding& encoding)
{
  encoding.process_nucleotides(m_sequence);
  encoding.process_qualities(m_qualities);
}

fastq::ntrimmed
fastq::trim_sequence_and_qualities(const size_t left_inclusive,
                                   const size_t right_exclusive)
{
  const ntrimmed summary(left_inclusive, length() - right_exclusive);
  truncate(left_inclusive, right_exclusive - left_inclusive);

  return summary;
}

} // namespace adapterremoval
