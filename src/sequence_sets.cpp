// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "sequence_sets.hpp" // declarations
#include "debug.hpp"         // for AR_REQUIRE
#include "errors.hpp"        // for fastq_error
#include "linereader.hpp"    // for line_reader
#include "sequence.hpp"      // for dna_sequence
#include "strutils.hpp"      // for log_escape
#include "table_reader.hpp"  // for table_reader
#include <algorithm>         // for max, sort, find
#include <cstddef>           // for size_t
#include <initializer_list>  // for initializer_list
#include <sstream>           // for operator<<, basic_ostream, ostringstream
#include <stdexcept>         // for invalid_argument
#include <string>            // for to_string
#include <string_view>       // for string_view
#include <utility>           // for pair
#include <vector>            // for vector, vector<>::const_iterator

namespace adapterremoval {

namespace {

void
validate_sample_name(std::string_view name)
{
  AR_REQUIRE(!name.empty());
  if (name == "unidentified") {
    throw parsing_error("The sample name 'unidentified' is a reserved "
                        "name, and cannot be used!");
  }

  for (const char c : name) {
    if (!is_ascii_letter_or_digit(c) && (c != '_')) {
      std::ostringstream error;
      error << "The sample name " << log_escape(name)
            << " is not a valid sample name; only letters ('a' to 'z' and 'A' "
               "to 'Z'), numbers (0 to 9) and underscores (_) are allowed.";

      throw parsing_error(error.str());
    }
  }
}

dna_sequence
parse_barcode_sequence(const std::string_view seq, const int mate)
{
  for (const auto nuc : seq) {
    switch (to_upper(nuc)) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
        break;

      default:
        std::ostringstream error;
        error << "Unsupported character found in mate " << mate << " barcode "
              << "sequence " << log_escape(seq) << ". Only bases A, C, G, "
              << "and T, are supported; please fix before continuing";

        throw parsing_error(error.str());
    }
  }

  return dna_sequence{ seq };
}

void
validate_barcode_length(const std::string_view seq,
                        const size_t expected_length,
                        const int mate)
{
  if (seq.length() != expected_length) {
    std::ostringstream error;
    error << "Inconsistent mate " << mate << " barcode lengths found: Last "
          << "barcode was " << expected_length << " base-pairs long, but "
          << "barcode " << log_escape(seq) << " is " << seq.length() << " "
          << "base-pairs long. Variable length barcodes are not supported";

    throw parsing_error(error.str());
  }
}

/**
 * Checks for multiple barcodes, while allowing that both the forward and
 * reverse barcode has been specified for a sample
 */
void
disallow_multiple_barcodes(const sample& s)
{
  if (s.size() > 1) {
    if (s.size() == 2) {
      const auto& ss_0 = s.at(0);
      const auto& ss_1 = s.at(1);

      // It is allowed to explicitly specify both forward and reverse barcodes
      if (ss_0.barcode_1 == ss_1.barcode_2 &&
          ss_0.barcode_2 == ss_1.barcode_1 &&
          ss_0.orientation != ss_1.orientation) {
        // At this point samples should not have been configured further
        AR_REQUIRE(ss_0.orientation != barcode_orientation::unspecified);
        AR_REQUIRE(ss_1.orientation != barcode_orientation::unspecified);
        AR_REQUIRE(ss_0.has_read_group == ss_1.has_read_group);
        AR_REQUIRE(ss_0.adapters() == ss_1.adapters());
        AR_REQUIRE(ss_0.read_group_ == ss_1.read_group_);
        return;
      }
    }

    std::ostringstream error;
    error << "Duplicate sample name " << log_escape(s.name())
          << "; multiple barcodes per samples is not enabled. Either ensure "
             "that all sample names are unique or use --multiple-barcodes to "
             "map multiple barcodes to a single sample";

    throw parsing_error(error.str());
  }
}

/** Check that sample names are valid and not overlapping (case-insensitive) */
void
check_sample_names(const std::vector<sample>& samples)
{
  std::vector<std::pair<std::string, std::string>> names;
  for (const auto& sample : samples) {
    const auto& name = sample.name();

    validate_sample_name(name);
    names.emplace_back(to_lower(name), name);
  }

  std::sort(names.begin(), names.end());
  for (size_t i = 1; i < names.size(); ++i) {
    const auto& [key_0, name_0] = names.at(i - 1);
    const auto& [key_1, name_1] = names.at(i);

    if ((key_0 == key_1) && (name_0 != name_1)) {
      std::ostringstream error;
      error << "Samples with names " << log_escape(name_0) << " and "
            << log_escape(name_1) << " differ only by case. Either use the "
            << "exact same name for both, if they the same sample, or give "
               "them distinct names";

      throw parsing_error(error.str());
    }
  }
}

/** Checks for basic barcode validity, but not for overlapping barcodes */
void
check_barcode_sequences(const std::vector<sample>& samples,
                        bool allow_multiple_barcodes)
{
  if (samples.empty()) {
    throw parsing_error("No samples/barcodes found in table");
  }

  auto mate_1_len = static_cast<size_t>(-1);
  auto mate_2_len = static_cast<size_t>(-1);

  for (const auto& sample : samples) {
    if (!allow_multiple_barcodes) {
      disallow_multiple_barcodes(sample);
    }

    for (const auto& it : sample) {
      if (mate_1_len == static_cast<size_t>(-1)) {
        mate_1_len = it.barcode_1.length();
        mate_2_len = it.barcode_2.length();
      }

      validate_barcode_length(it.barcode_1.as_string(), mate_1_len, 1);
      validate_barcode_length(it.barcode_2.as_string(), mate_2_len, 2);
    }
  }
}

/** Helper class used to validate barcodes */
struct barcode_key
{
  std::string_view barcode_1;
  std::string_view barcode_2;
  barcode_orientation orientation;
  std::string_view sample;

  [[nodiscard]] std::string describe() const
  {
    std::ostringstream out;
    out << log_escape(this->sample) << " (" << this->barcode_1;

    if (!this->barcode_2.empty()) {
      out << "-" << this->barcode_2;
    }

    switch (this->orientation) {
      case barcode_orientation::unspecified:
        out << ")";
        break;
      case barcode_orientation::forward:
        out << "; forward)";
        break;
      case barcode_orientation::reverse:
        out << "; reverse)";
        break;
      default:                                  // GCOVR_EXCL_LINE
        AR_FAIL("invalid barcode_orientation"); // GCOVR_EXCL_LINE
    }

    return out.str();
  }

  bool operator<(const barcode_key& other) const
  {
    // Sorted by barcodes first to allow easy duplicate checks
    if (this->barcode_1 != other.barcode_1) {
      return this->barcode_1 < other.barcode_1;
    } else if (this->barcode_2 != other.barcode_2) {
      return this->barcode_2 < other.barcode_2;
    } else if (this->sample != other.sample) {
      // Sorted by sample name after barcodes for more sensible looking output
      // when overlapping sequences are printed
      return this->sample < other.sample;
    } else {
      return this->orientation < other.orientation;
    }
  }
};

/** Checks for overlapping barcodes */
void
check_barcode_overlap(const std::vector<sample>& samples, bool paired_end)
{
  std::vector<barcode_key> sequences;
  for (const auto& sample : samples) {
    for (const auto& it : sample) {
      sequences.emplace_back(barcode_key{ it.barcode_1.as_string(),
                                          it.barcode_2.as_string(),
                                          it.orientation,
                                          sample.name() });
    }
  }

  std::sort(sequences.begin(), sequences.end());
  for (size_t i = 1; i < sequences.size(); ++i) {
    const auto& it_0 = sequences.at(i - 1);
    const auto& it_1 = sequences.at(i);

    if (it_0.barcode_1 == it_1.barcode_1 &&
        (!paired_end || it_0.barcode_2 == it_1.barcode_2)) {
      std::ostringstream error;
      error << "Sample " << it_0.describe() << " and sample " << it_1.describe()
            << " have overlapping barcodes";

      if (it_0.barcode_2 != it_1.barcode_2) {
        error << ". Note that AdapterRemoval cannot distinguish these barcodes "
              << "in single-end mode, even though the second barcodes differ";
      }

      error << ". Please remove any duplicate entries from the barcode table "
               "before continuing";

      throw parsing_error(error.str());
    }
  }
}

/** Returns the reverse of a user specified orientation */
barcode_orientation
reverse_orientation(barcode_table_orientation orientation)
{
  switch (orientation) {
    case barcode_table_orientation::forward:
      return barcode_orientation::reverse;
    case barcode_table_orientation::reverse:
      return barcode_orientation::forward;
    // Cannot be reversed
    case barcode_table_orientation::unspecified:
    case barcode_table_orientation::explicit_:
      return barcode_orientation::unspecified;
    default:                                        // GCOVR_EXCL_LINE
      AR_FAIL("invalid barcode_table_orientation"); // GCOVR_EXCL_LINE
  }
}

/** Build list of samples from individual barcodes ordered by sample */
void
append_sample(std::vector<sample>& samples,
              const std::string& name,
              const dna_sequence& barcode_1,
              const dna_sequence& barcode_2,
              barcode_orientation orientation)
{

  if (samples.empty() || samples.back().name() != name) {
    samples.emplace_back(name, barcode_1, barcode_2, orientation);
  } else {
    samples.back().add_barcodes(barcode_1, barcode_2, orientation);
  }
}

/** Updates a vector of barcodes to include the reverse sequences  */
void
create_reversed_barcodes(std::vector<sample>& samples,
                         barcode_table_orientation orientation)
{
  const auto rev_orientation = reverse_orientation(orientation);
  if (rev_orientation == barcode_orientation::unspecified) {
    return;
  }

  // This is inefficient, but easier than trying to sort the list afterwards
  std::vector<sample> rsamples;
  for (const auto& sample : samples) {
    for (const auto& seqs : sample) {
      const auto& name = sample.name();
      const auto& barcode_1 = seqs.barcode_1;
      const auto& barcode_2 = seqs.barcode_2;

      append_sample(rsamples, name, barcode_1, barcode_2, seqs.orientation);
      // NOLINTNEXTLINE(readability-suspicious-call-argument)
      append_sample(rsamples, name, barcode_2, barcode_1, rev_orientation);
    }
  }

  std::swap(samples, rsamples);
}

/** Returns expected number of columns for different types of barcode tables */
std::pair<int, int>
barcode_table_columns(barcode_table_orientation orientation)
{
  auto min_columns = 2; // Name and one barcode
  auto max_columns = 3; // Name and two barcodes
  switch (orientation) {
    case barcode_table_orientation::unspecified:
      break;
    case barcode_table_orientation::forward:
    case barcode_table_orientation::reverse:
      min_columns = 3;
      break;
    case barcode_table_orientation::explicit_:
      min_columns = 4; // Name, two barcodes, and orientation
      max_columns = 4;
      break;
    default:                                  // GCOVR_EXCL_LINE
      AR_FAIL("invalid barcode_orientation"); // GCOVR_EXCL_LINE
  }

  return { min_columns, max_columns };
}

/** Parse barcode orientation from a table file */
barcode_orientation
parse_barcode_orientation(std::string_view name, std::string value)
{
  value = to_lower(value);
  if (value == "forward" || value == "fwd" || value == "+") {
    return barcode_orientation::forward;
  } else if (value == "reverse" || value == "rev" || value == "-") {
    return barcode_orientation::reverse;
  }

  throw parsing_error("Invalid barcode orientation for sample " +
                      log_escape(name) + ": " + log_escape(value));
}

} // namespace

///////////////////////////////////////////////////////////////////////////////

barcode_table_orientation
parse_table_orientation(std::string_view value)
{
  value = trim_ascii_whitespace(value);
  auto value_l = to_lower(std::string{ value });

  if (value_l == "forward") {
    return barcode_table_orientation::forward;
  } else if (value_l == "reverse") {
    return barcode_table_orientation::reverse;
  } else if (value_l == "explicit") {
    return barcode_table_orientation::explicit_;
  } else if (value_l == "unspecified") {
    return barcode_table_orientation::unspecified;
  } else {
    throw std::invalid_argument("invalid barcode table orientation ");
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'adapter_set' class

adapter_set::adapter_set(std::initializer_list<string_view_pair> args)
{
  for (const auto& [first, second] : args) {
    add(dna_sequence{ first }, dna_sequence{ second });
  }
}

void
adapter_set::add(dna_sequence adapter1, const dna_sequence& adapter2)
{
  m_adapters.emplace_back(std::move(adapter1), adapter2.reverse_complement());
}

adapter_set
adapter_set::add_barcodes(const dna_sequence& barcode1,
                          const dna_sequence& barcode2) const
{
  adapter_set adapters;
  const auto barcode2rc = barcode2.reverse_complement();
  for (const auto& [first, second] : m_adapters) {
    // Add sequences directly in alignment orientation
    adapters.m_adapters.emplace_back(barcode2rc + first, second + barcode1);
  }

  return adapters;
}

void
adapter_set::load(std::string filename, bool paired_end_mode)
{
  line_reader reader(std::move(filename));
  load(reader, paired_end_mode);
}

void
adapter_set::load(line_reader_base& reader, bool paired_end_mode)
{
  const auto table = table_reader()
                       .with_comment_char('#')
                       .with_min_columns(1 + paired_end_mode)
                       .with_max_columns(2)
                       .parse(reader);

  sequence_pair_vec adapters;
  for (const auto& row : table) {
    dna_sequence adapter_1{ row.at(0) };
    dna_sequence adapter_2;
    if (row.size() > 1) {
      adapter_2 = dna_sequence{ row.at(1) };
    }

    // Convert from read to alignment orientation
    adapters.emplace_back(std::move(adapter_1), adapter_2.reverse_complement());
  }

  if (adapters.empty()) {
    throw parsing_error("No adapter sequences in table");
  }

  std::swap(m_adapters, adapters);
}

sequence_pair_vec
adapter_set::to_read_orientation() const
{
  sequence_pair_vec adapters;
  for (const auto& [first, second] : m_adapters) {
    adapters.emplace_back(first, second.reverse_complement());
  }

  return adapters;
}

bool
adapter_set::operator==(const adapter_set& other) const
{
  return m_adapters == other.m_adapters;
}

std::ostream&
operator<<(std::ostream& os, const adapter_set& value)
{
  os << "adapter_set{[";

  bool is_first = true;
  for (const auto& [first, second] : value) {
    if (!is_first) {
      os << ", ";
    }

    is_first = false;
    os << "pair{first=" << first << ", second=" << second << "}";
  }

  return os << "]}";
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'sample_sequences' class

sample_sequences::sample_sequences(dna_sequence barcode_1_,
                                   dna_sequence barcode_2_,
                                   barcode_orientation orientation_)
  : barcode_1(std::move(barcode_1_))
  , barcode_2(std::move(barcode_2_))
  , orientation(orientation_)
{
}

const adapter_set&
sample_sequences::adapters() const
{
  AR_REQUIRE(!m_uninitialized_adapters);
  return m_adapters;
}

void
sample_sequences::set_adapters(adapter_set as)
{
  m_adapters = std::move(as);
  m_uninitialized_adapters = false;
}

void
sample_sequences::flag_uninitialized_adapters()
{
  m_uninitialized_adapters = true;
}

bool
sample_sequences::operator==(const sample_sequences& other) const
{
  return this->has_read_group == other.has_read_group &&
         this->read_group_ == other.read_group_ &&
         this->barcode_1 == other.barcode_1 &&
         this->barcode_2 == other.barcode_2 &&
         this->orientation == other.orientation &&
         this->m_adapters == other.m_adapters &&
         this->m_uninitialized_adapters == other.m_uninitialized_adapters;
}

std::ostream&
operator<<(std::ostream& os, const sample_sequences& value)
{
  return os << "sample_sequences{has_read_group="
            << (value.has_read_group ? "true" : "false")
            << ", read_group=" << value.read_group_
            << ", barcode_1=" << value.barcode_1
            << ", barcode_2=" << value.barcode_2
            << ", orientation=" << value.orientation
            << ", adapters=" << value.m_adapters << ", uninitialized_adapters="
            << (value.m_uninitialized_adapters ? "true" : "false") << "}";
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'sample' class

sample::sample()
  : sample(std::string{},
           dna_sequence{},
           dna_sequence{},
           barcode_orientation::unspecified)
{
}

sample::sample(std::string name,
               dna_sequence barcode1,
               dna_sequence barcode2,
               barcode_orientation orientation)
  : m_name(std::move(name))
{
  add_barcodes(std::move(barcode1), std::move(barcode2), orientation);
}

void
sample::add_barcodes(dna_sequence barcode1,
                     dna_sequence barcode2,
                     barcode_orientation orientation)
{
  AR_REQUIRE(barcode2.empty() || !barcode1.empty());
  m_barcodes.emplace_back(std::move(barcode1),
                          std::move(barcode2),
                          orientation);
}

void
sample::set_adapters(const adapter_set& adapters)
{
  for (auto& it : m_barcodes) {
    it.set_adapters(adapters.add_barcodes(it.barcode_1, it.barcode_2));
  }
}

void
sample::set_read_group(const read_group& read_group_)
{
  for (auto it = m_barcodes.begin(); it != m_barcodes.end(); ++it) {
    it->read_group_ = read_group_;
    it->has_read_group = true;

    if (!m_name.empty()) {
      it->read_group_.set_sample(m_name);

      if (m_barcodes.size() > 1) {
        std::string id = m_name;
        id.push_back('.');
        id.append(std::to_string((it - m_barcodes.begin()) + 1));

        it->read_group_.set_id(id);
      } else {
        it->read_group_.set_id(m_name);
      }
    }

    if (it->barcode_1.length() || it->barcode_2.length()) {
      std::string barcodes;
      barcodes.append(it->barcode_1.as_string());
      if (it->barcode_2.length()) {
        barcodes.push_back('-');
        barcodes.append(it->barcode_2.as_string());
      }

      it->read_group_.set_barcodes(barcodes);
    }

    it->read_group_.set_orientation(it->orientation);
  }
}

void
sample::flag_uninitialized_adapters()
{
  for (auto& it : m_barcodes) {
    it.flag_uninitialized_adapters();
  }
}

bool
sample::operator==(const sample& other) const
{
  return m_name == other.m_name && m_barcodes == other.m_barcodes;
}

std::ostream&
operator<<(std::ostream& os, const sample& value)
{
  return os << "sample{name=" << log_escape(value.name()) << ", barcodes=["
            << join_text(value, ", ") << "]}";
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'sample_set' class

sample_set::sample_set()
  : m_samples{ sample{} }
  , m_unidentified({},
                   dna_sequence{},
                   dna_sequence{},
                   barcode_orientation::unspecified)
{
  set_unidentified_read_group(m_read_group);
}

sample_set::sample_set(std::initializer_list<std::string_view> lines,
                       barcode_config config)
{
  vec_reader reader(lines);
  load(reader, config);
}

void
sample_set::set_adapters(adapter_set adapters)
{
  m_adapters = std::move(adapters);
  for (auto& sample : m_samples) {
    sample.set_adapters(m_adapters);
  }

  m_unidentified.set_adapters(m_adapters);
  m_uninitialized_adapters = false;
}

void
sample_set::set_read_group(std::string_view value)
{
  m_read_group = read_group(value);
  for (auto& sample : m_samples) {
    sample.set_read_group(m_read_group);
  }

  set_unidentified_read_group(m_read_group);
}

void
sample_set::load(std::string filename, const barcode_config& config)
{
  line_reader reader(std::move(filename));
  load(reader, config);
}

void
sample_set::load(line_reader_base& reader, const barcode_config& config)
{
  auto [min_columns, max_columns] = barcode_table_columns(config.m_orientation);

  auto barcodes = table_reader()
                    .with_comment_char('#')
                    .with_min_columns(min_columns)
                    .with_max_columns(max_columns)
                    .parse(reader);

  // Sort by sample name to simplify sample set construction
  std::sort(barcodes.begin(), barcodes.end(), [](const auto& a, const auto& b) {
    return a.at(0) < b.at(0);
  });

  std::vector<sample> samples;
  for (const auto& row : barcodes) {
    const auto& name = row.at(0);
    const auto barcode_1 = parse_barcode_sequence(row.at(1), 1);
    dna_sequence barcode_2;

    if (row.size() > 2) {
      barcode_2 = parse_barcode_sequence(row.at(2), 2);
    }

    auto orientation = barcode_orientation::unspecified;
    if (config.m_orientation == barcode_table_orientation::explicit_) {
      orientation = parse_barcode_orientation(name, row.at(3));
    } else {
      orientation = static_cast<barcode_orientation>(config.m_orientation);
    }

    append_sample(samples, name, barcode_1, barcode_2, orientation);
  }

  // Check sample names for overlap (case-insensitive) and disallowed characters
  check_sample_names(samples);

  // Basic properties are checked first, before (potentially) reversing barcodes
  check_barcode_sequences(samples, config.m_allow_multiple_barcodes);

  // Create reversed barcodes if enabled
  create_reversed_barcodes(samples, config.m_orientation);

  // Check for overlap between user barcodes, generated barcodes, or both
  check_barcode_overlap(samples, config.m_paired_end_mode);

  // Update read-group information and generate adapters based on all barcodes
  for (auto& s : samples) {
    s.set_read_group(m_read_group);
    if (m_uninitialized_adapters) {
      s.flag_uninitialized_adapters();
    } else {
      s.set_adapters(m_adapters);
    }
  }

  std::swap(m_samples, samples);
}

const adapter_set&
sample_set::adapters() const
{
  AR_REQUIRE(!m_uninitialized_adapters);
  return m_adapters;
}

void
sample_set::flag_uninitialized_adapters()
{
  m_uninitialized_adapters = true;
  for (auto& it : m_samples) {
    it.flag_uninitialized_adapters();
  }
}

const adapter_set&
sample_set::uninitialized_adapters() const
{
  AR_REQUIRE(m_uninitialized_adapters);
  return m_adapters;
}

void
sample_set::set_unidentified_read_group(read_group tmpl)
{
  // Unidentified reads lack a SM tag, so add a description instead
  tmpl.set_sample("");
  tmpl.set_description("unidentified");
  m_unidentified.set_read_group(tmpl);
}

std::ostream&
operator<<(std::ostream& os, const sample_set& value)
{
  return os << "sample_set{samples=[" << join_text(value, ", ") << "]"
            << ", unidentified=" << value.unidentified()
            << ", read_group=" << value.readgroup()
            << ", adapters=" << value.adapters() << "}";
}

////////////////////////////////////////////////////////////////////////////////

std::ostream&
operator<<(std::ostream& os, const barcode_orientation& value)
{
  switch (value) {
    case barcode_orientation::unspecified:
      return os << "barcode_orientation::unspecified";
    case barcode_orientation::forward:
      return os << "barcode_orientation::forward";
    case barcode_orientation::reverse:
      return os << "barcode_orientation::reverse";
    default:                                 // GCOVR_EXCL_LINE
      return os << "barcode_orientation{?}"; // GCOVR_EXCL_LINE
  }
}

std::ostream&
operator<<(std::ostream& os, const barcode_table_orientation& value)
{
  switch (value) {
    case barcode_table_orientation::unspecified:
      return os << "barcode_table_orientation::unspecified";
    case barcode_table_orientation::forward:
      return os << "barcode_table_orientation::forward";
    case barcode_table_orientation::reverse:
      return os << "barcode_table_orientation::reverse";
    case barcode_table_orientation::explicit_:
      return os << "barcode_table_orientation::explicit_";
    default:                                       // GCOVR_EXCL_LINE
      return os << "barcode_table_orientation{?}"; // GCOVR_EXCL_LINE
  }
}

} // namespace adapterremoval
