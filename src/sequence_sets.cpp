// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "sequence_sets.hpp" // declarations
#include "debug.hpp"         // for AR_REQUIRE
#include "errors.hpp"        // for fastq_error
#include "linereader.hpp"    // for line_reader
#include "sequence.hpp"      // for dna_sequence
#include "strutils.hpp"      // for string_vec, indent_lines
#include "table_reader.hpp"  // for table_reader
#include <algorithm>         // for max, sort, find
#include <sstream>           // for operator<<, basic_ostream, ostringstream
#include <stdexcept>         // for invalid_argument
#include <string_view>       // for string_view
#include <utility>           // for pair
#include <vector>            // for vector, vector<>::const_iterator

namespace adapterremoval {

namespace {

void
validate_sample_name(std::string_view name)
{
  if (name.empty()) {
    throw parsing_error("Sample names must not be empty");
  } else if (name == "unidentified") {
    throw parsing_error("The sample name 'unidentified' is a reserved "
                        "name, and cannot be used!");
  }

  for (const char c : name) {
    if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
          (c >= '0' && c <= '9') || (c == '_'))) {
      std::ostringstream error;
      error << "The sample name " << log_escape(name)
            << " is not a valid sample name; only letters ('a' to 'z' and 'A' "
               "to 'Z'), numbers (0 to 9) and underscores (_) are allowed.";

      throw parsing_error(error.str());
    }
  }
}

void
validate_barcode_sequence(const dna_sequence& sequence,
                          const size_t expected_length,
                          const int mate)
{
  std::string_view seq = sequence;

  if (seq.find('N') != std::string::npos) {
    std::ostringstream error;
    error << "Degenerate base (N) found in mate " << mate << " barcode "
          << "sequence " << log_escape(seq) << ". Degenerate bases are not "
          << "supported for demultiplexing; please remove before continuing!";

    throw parsing_error(error.str());
  }

  if (sequence.length() != expected_length) {
    std::ostringstream error;
    error << "Inconsistent mate " << mate << " barcode lengths found: Last "
          << "barcode was " << expected_length << " base-pairs long, but "
          << "barcode " << log_escape(seq) << " is " << seq.length() << " "
          << "base-pairs long. Variable length barcodes are not supported";

    throw parsing_error(error.str());
  }
}

bool
check_barcodes_sequences(const std::vector<sample>& samples,
                         const std::string& source,
                         bool paired_end = false)
{
  if (samples.empty()) {
    throw parsing_error("No samples/barcodes provided");
  }

  auto mate_1_len = static_cast<size_t>(-1);
  auto mate_2_len = static_cast<size_t>(-1);

  std::vector<std::pair<std::string_view, std::string_view>> sequences;
  for (const auto& sample : samples) {
    validate_sample_name(sample.name());

    for (const auto& it : sample) {
      if (mate_1_len == static_cast<size_t>(-1)) {
        mate_1_len = it.barcode_1.length();
        mate_2_len = it.barcode_2.length();

        if (!mate_1_len) {
          throw parsing_error("Empty barcode 1 sequence for " + sample.name());
        }
      }

      validate_barcode_sequence(it.barcode_1, mate_1_len, 1);
      validate_barcode_sequence(it.barcode_2, mate_2_len, 2);

      if (paired_end) {
        sequences.emplace_back(it.barcode_1, it.barcode_2);
      } else {
        sequences.emplace_back(it.barcode_1, dna_sequence{});
      }
    }
  }

  std::sort(sequences.begin(), sequences.end());
  for (size_t i = 1; i < sequences.size(); ++i) {
    if (sequences.at(i - 1) == sequences.at(i)) {
      const auto& it = sequences.at(i);

      if (paired_end) {
        std::ostringstream error;
        error << "Duplicate barcode pairs found in " << log_escape(source)
              << " with barcodes " << log_escape(it.first) << " and "
              << log_escape(it.second) << ". please verify correctness of "
              << "the barcode table and remove any duplicate entries!";

        throw parsing_error(error.str());
      } else {
        std::ostringstream error;
        error << "Duplicate mate 1 barcodes found in " << log_escape(source)
              << ": " << log_escape(it.first) << ". Even if these are "
              << "associated with different mate 2 barcodes, it is not "
                 "possible to distinguish between these in single-end mode!";

        throw parsing_error(error.str());
      }
    }
  }

  return true;
}

constexpr bool
is_ascii_alpha(const char c)
{
  return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

constexpr bool
is_ascii_alphanum(const char c)
{
  return is_ascii_alpha(c) || (c >= '0' && c <= '9');
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Implementations for `read_group`

read_group::read_group()
  : m_header{ "@RG\tID:1" }
  , m_id{ "1" }
{
}

read_group::read_group(std::string_view value)
  : m_header{ "@RG" }
{
  using invalid = std::invalid_argument;

  // It's not unreasonable to except users to try to specify a full @RG line
  if (starts_with(value, "@RG\t")) {
    value = value.substr(4);
  }

  if (!value.empty()) {
    for (const auto& field : split_text(value, '\t')) {
      for (const auto c : field) {
        if (c < ' ' || c > '~') {
          throw invalid("only characters in the range ' ' to '~' are allowed");
        }
      }

      if (field.size() < 4) {
        throw invalid("tags must be at least 4 characters long");
      } else if (!is_ascii_alpha(field.at(0))) {
        throw invalid("first character in tag name must be a letter");
      } else if (!is_ascii_alphanum(field.at(1))) {
        throw invalid(
          "second character in tag name must be a letter or number");
      } else if (field.at(2) != ':') {
        throw invalid("third character in tag must be a colon");
      } else if (starts_with(field, "ID:")) {
        if (!m_id.empty()) {
          throw invalid("multiple ID tags found");
        }

        m_id = field.substr(3);
      }
    }
  }

  // Generic read-group key if the user did not supply one
  if (m_id.empty()) {
    m_id = "1";
    m_header += "\tID:1";
  }

  if (!value.empty()) {
    m_header += "\t";
    m_header += value;
  }
}

void
read_group::set_id(std::string_view id)
{
  AR_REQUIRE(!id.empty());
  update_tag("ID", id);
  m_id = id;
}

void
read_group::update_tag(std::string_view key, std::string_view value)
{
  AR_REQUIRE(!key.empty());
  std::string cache;
  cache.push_back('\t');
  cache.append(key);
  cache.push_back(':');

  auto index = m_header.find(cache);
  if (index != std::string::npos) {
    if (value.empty()) {
      cache = m_header.substr(0, index);
    } else {
      cache = m_header.substr(0, index + cache.size());
      cache.append(value);
    }

    index = m_header.find('\t', index + 1);
    if (index != std::string::npos) {
      cache.append(m_header.substr(index));
    }

    m_header = cache;
  } else if (!value.empty()) {
    m_header.append(cache);
    m_header.append(value);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'adapters' class

adapter_set::adapter_set(std::initializer_list<string_view_pair> args)
{
  for (const auto& pair : args) {
    add(std::string{ pair.first }, std::string{ pair.second });
  }
}

void
adapter_set::add(dna_sequence adapter1, dna_sequence adapter2)
{
  m_adapters.emplace_back(adapter1, adapter2.reverse_complement());
}

void
adapter_set::add(std::string adapter1, const std::string adapter2)
{
  add(dna_sequence{ adapter1 }, dna_sequence{ adapter2 });
}

adapter_set
adapter_set::add_barcodes(const dna_sequence& barcode1,
                          const dna_sequence& barcode2) const
{
  adapter_set adapters;
  for (const auto& it : m_adapters) {
    // Add sequences directly in alignment orientation
    adapters.m_adapters.emplace_back(barcode2.reverse_complement() + it.first,
                                     it.second + barcode1);
  }

  return adapters;
}

void
adapter_set::load(const std::string& filename, bool paired_end)
{
  line_reader reader(filename);
  const auto adapters = table_reader()
                          .with_comment_char('#')
                          .with_min_columns(1 + paired_end)
                          .with_max_columns(2)
                          .parse(reader);

  for (const auto& row : adapters) {
    dna_sequence adapter_1{ row.at(0) };
    dna_sequence adapter_2;
    if (row.size() > 1) {
      adapter_2 = dna_sequence{ row.at(1) };
    }

    // Convert from read to alignment orientation
    m_adapters.emplace_back(adapter_1, adapter_2.reverse_complement());
  }
}

sequence_pair_vec
adapter_set::to_read_orientation() const
{
  sequence_pair_vec adapters;
  for (const auto& it : m_adapters) {
    adapters.emplace_back(it.first, it.second.reverse_complement());
  }

  return adapters;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'sample' class

void
sample::add(dna_sequence barcode1, dna_sequence barcode2)
{
  m_barcodes.emplace_back(std::move(barcode1), std::move(barcode2));
}

void
sample::add(std::string barcode1, std::string barcode2)
{
  add(dna_sequence(barcode1), dna_sequence(barcode2));
}

void
sample::set_adapters(const adapter_set& adapters)
{
  for (auto& it : m_barcodes) {
    it.adapters = adapters.add_barcodes(it.barcode_1, it.barcode_2);
  }
}

void
sample::set_read_group(const read_group& info)
{
  for (auto it = m_barcodes.begin(); it != m_barcodes.end(); ++it) {
    it->info = info;
    it->has_read_group = true;

    if (!m_name.empty()) {
      it->info.set_sample(m_name);

      if (m_barcodes.size() > 1) {
        std::string id = m_name;
        id.push_back('.');
        id.append(std::to_string((it - m_barcodes.begin()) + 1));

        it->info.set_id(id);
      } else {
        it->info.set_id(m_name);
      }
    }

    if (it->barcode_1.length() || it->barcode_2.length()) {
      std::string barcodes;
      barcodes.append(it->barcode_1);
      if (it->barcode_2.length()) {
        barcodes.push_back('-');
        barcodes.append(it->barcode_2);
      }

      it->info.set_barcodes(barcodes);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'sample_set' class

sample_set::sample_set()
  : m_samples{ sample{} }
  , m_unidentified("", dna_sequence{}, dna_sequence{})
{
  set_unidentified_read_group(m_read_group);
}

sample_set::sample_set(std::initializer_list<sample> args)
  : m_samples(args)
  , m_unidentified("", dna_sequence{}, dna_sequence{})
{
  set_unidentified_read_group(m_read_group);

  std::sort(m_samples.begin(),
            m_samples.end(),
            [](const auto& a, const auto& b) { return a.name() < b.name(); });

  std::string_view name;
  for (const auto& sample : m_samples) {
    validate_sample_name(sample.name());
    if (sample.name() == name) {
      throw parsing_error("duplicate sample name: " + sample.name());
    }

    name = sample.name();
  }

  check_barcodes_sequences(m_samples, "initializer_list", true);
}

void
sample_set::set_adapters(adapter_set adapters)
{
  m_adapters = std::move(adapters);
  for (auto& sample : m_samples) {
    sample.set_adapters(m_adapters);
  }

  m_unidentified.set_adapters(m_adapters);
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
sample_set::load(const std::string& filename, const barcode_config& config)
{
  line_reader reader(filename);
  auto barcodes = table_reader()
                    .with_comment_char('#')
                    .with_min_columns(2 + !config.m_unidirectional_barcodes)
                    .with_max_columns(3)
                    .parse(reader);

  std::sort(barcodes.begin(), barcodes.end(), [](const auto& a, const auto& b) {
    return a.at(0) < b.at(0);
  });

  std::vector<sample> samples{};
  for (const auto& row : barcodes) {
    const auto& name = row.at(0);
    dna_sequence barcode_1{ row.at(1) };
    dna_sequence barcode_2;

    if (row.size() > 2) {
      barcode_2 = dna_sequence{ row.at(2) };
    }

    if (samples.empty() || samples.back().name() != name) {
      samples.emplace_back(name, barcode_1, barcode_2);
    } else if (config.m_allow_multiple_barcodes) {
      samples.back().add(barcode_1, barcode_2);
    } else {
      std::ostringstream error;
      error << "Duplicate sample name " << log_escape(name)
            << "; multiple barcodes per samples is not enabled. Either ensure "
               "that all sample names are unique or use --multiple-barcodes to "
               "map multiple barcodes to a single sample";

      throw parsing_error(error.str());
    }
  }

  // Check before adding reversed barcodes, to prevent misleading error messages
  check_barcodes_sequences(samples, filename, config.m_paired_end_mode);

  std::swap(m_samples, samples);
  if (!config.m_unidirectional_barcodes) {
    add_reversed_barcodes(config);
  }

  for (auto& s : m_samples) {
    s.set_read_group(m_read_group);
    s.set_adapters(m_adapters);
  }
}

void
sample_set::set_unidentified_read_group(read_group tmpl)
{
  // Unidentified reads lack a SM tag, so add a comment instead
  tmpl.set_comment("unidentified");
  m_unidentified.set_read_group(tmpl);
}

void
sample_set::add_reversed_barcodes(const barcode_config& config)
{
  for (auto& sample : m_samples) {
    // Original number of barcodes
    const size_t count = sample.size();

    for (size_t i = 0; i < count; ++i) {
      const auto& sequences = sample.at(i);
      auto barcode_1 = sequences.barcode_2.reverse_complement();
      auto barcode_2 = sequences.barcode_1.reverse_complement();
      AR_REQUIRE(barcode_1.length() && barcode_2.length());

      bool reverse_found = false;
      for (const auto& it : sample) {
        if (it.barcode_1 == barcode_1 && it.barcode_2 == barcode_2) {
          reverse_found = true;
          break;
        }
      }

      if (!reverse_found) {
        sample.add(std::move(barcode_1), std::move(barcode_2));
      }
    }
  }

  check_barcodes_sequences(
    m_samples, "reversed barcodes", config.m_paired_end_mode);
}

} // namespace adapterremoval
