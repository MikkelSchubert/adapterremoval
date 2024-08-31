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
#include "sequence_sets.hpp"
#include "commontypes.hpp" // for fastq_vec
#include "debug.hpp"       // for AR_REQUIRE
#include "errors.hpp"      // for fastq_error
#include "linereader.hpp"  // for line_reader
#include "sequence.hpp"    // for dna_sequence
#include "strutils.hpp"    // for string_vec, indent_lines
#include <algorithm>       // for max, sort, find
#include <sstream>         // for operator<<, basic_ostream, ostringstream
#include <string_view>     // for string_view
#include <utility>         // for pair
#include <vector>          // for vector, vector<>::const_iterator

namespace adapterremoval {

namespace {

std::string
trim_comments(std::string line)
{
  const size_t index = line.find('#');
  if (index != std::string::npos) {
    line.resize(index);
  }

  return line;
}

struct table_row
{
  std::string name{};
  dna_sequence sequence_1{};
  dna_sequence sequence_2{};
};

std::vector<table_row>
read_table(const std::string& filename,
           const bool paired_required = false,
           const bool row_names = false)
{
  size_t last_row_size = 0;
  size_t line_num = 1;
  const size_t min_col = 1 + paired_required + row_names;
  const size_t max_col = 2 + row_names;
  std::vector<table_row> table;

  try {
    line_reader adapter_file(filename);

    for (std::string line; adapter_file.getline(line); ++line_num) {
      string_vec row;
      std::string field;

      std::istringstream instream(trim_comments(line));
      while (instream >> field) {
        row.push_back(field);
        field.clear();
      }

      if (row.empty()) {
        // Ignore empty lines, e.g. those containing only comments
        continue;
      } else if (row.size() < min_col) {
        std::ostringstream message;
        message << "Expected at least " << min_col << " columns in " << "table "
                << log_escape(filename) << " at line " << line_num
                << ", but only found " << row.size() << " column(s)!";

        throw parsing_error(message.str());
      } else if (row.size() > max_col) {
        std::ostringstream message;
        message << "Expected at most " << min_col << " columns in " << "table "
                << log_escape(filename) << " at line " << line_num
                << ", but found " << row.size() << " column(s)!";

        throw parsing_error(message.str());
      } else if (last_row_size && last_row_size != row.size()) {
        std::ostringstream message;
        message << "Error reading " << log_escape(filename) << " at line "
                << line_num << "; rows contain unequal number of "
                << "columns; last row contained " << last_row_size
                << " column(s) but current row contains " << row.size()
                << " column(s)!";

        throw parsing_error(message.str());
      } else {
        last_row_size = row.size();
      }

      table_row record;
      if (row_names) {
        record.name = row.at(0);
      }

      record.sequence_1 = dna_sequence{ row.at(row_names) };
      if (row.size() > row_names + 1U) {
        record.sequence_2 = dna_sequence{ row.at(row_names + 1U) };
      }

      table.push_back(std::move(record));
    }
  } catch (const std::ios_base::failure& error) {
    throw parsing_error(error.what());
  } catch (const fastq_error& error) {
    std::ostringstream message;
    message << "invalid DNA sequence on line " << line_num << ": "
            << error.what();

    throw parsing_error(message.str());
  }

  return table;
}

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
    error << "Inconsistent mate " << mate << "barcode lengths found: Last "
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
  auto mate_1_len = static_cast<size_t>(-1);
  auto mate_2_len = static_cast<size_t>(-1);

  std::vector<std::pair<std::string_view, std::string_view>> sequences;
  for (const auto& it : samples) {
    validate_sample_name(it.name());

    for (const auto& it : it) {
      if (mate_1_len == static_cast<size_t>(-1)) {
        mate_1_len = it.barcode_1.length();
        mate_2_len = it.barcode_2.length();
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

/** Unescapes the escape sequences "\\" and "\t" in a read-group string */
std::string
unescape_read_group(std::string_view value)
{
  std::string result;

  bool in_escape = false;
  for (auto c : value) {
    if (in_escape) {
      if (c == '\\') {
        result.push_back('\\');
      } else if (c == 't') {
        result.push_back('\t');
      } else {
        throw std::invalid_argument("invalid escape sequence " +
                                    std::string("'\\") + c + '\'');
      }

      in_escape = false;
    } else if (c == '\\') {
      in_escape = true;
    } else {
      result.push_back(c);
    }
  }

  if (in_escape) {
    throw std::invalid_argument("incomplete escape sequence at end of string");
  }

  return result;
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
  : m_header{ "@RG\tID:1\tPG:adapterremoval" }
  , m_id{ "1" }
{
}

read_group::read_group(std::string_view value_)
  : m_header{ "@RG" }
{
  using invalid = std::invalid_argument;
  auto value = unescape_read_group(value_);

  // It's not unreasonable to except users to try to specify a full @RG line
  if (starts_with(value, "@RG\t")) {
    value = value.substr(4);
  }

  bool has_pg_tag = false;
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
      } else if (starts_with(field, "PG:")) {
        // Allow the user to override the PG field, just in case
        has_pg_tag = true;
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

  if (!has_pg_tag) {
    m_header.append("\tPG:adapterremoval");
  }
}

void
read_group::update_tag(std::string_view key, std::string_view value)
{
  AR_REQUIRE(!key.empty() && !value.empty());
  std::string cache;
  cache.push_back('\t');
  cache.append(key);
  cache.push_back(':');

  auto index = m_header.find(cache);
  if (index != std::string::npos) {
    cache = m_header.substr(0, index + cache.size());
    cache.append(value);

    index = m_header.find('\t', index + 1);
    if (index != std::string::npos) {
      cache.append(m_header.substr(index));
    }

    m_header = cache;
  } else {
    m_header.append(cache);
    m_header.append(value);
  }

  if (key == "ID") {
    m_id = value;
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
  auto table = read_table(filename, paired_end);

  for (auto& row : table) {
    // Convert from read to alignment orientation
    m_adapters.emplace_back(row.sequence_1,
                            row.sequence_2.reverse_complement());
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
  , m_unidentified("unidentified", dna_sequence{}, dna_sequence{})
{
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

  m_unidentified.set_read_group(m_read_group);
}

void
sample_set::load(const std::string& filename, const barcode_config& config)
{
  auto barcodes = read_table(filename, !config.m_unidirectional_barcodes, true);
  std::sort(barcodes.begin(), barcodes.end(), [](const auto& a, const auto& b) {
    return a.name < b.name;
  });

  m_samples.clear();
  for (const auto& row : barcodes) {
    if (m_samples.empty() || m_samples.back().name() != row.name) {
      sample s{ row.name, row.sequence_1, row.sequence_2 };
      s.set_adapters(m_adapters);
      s.set_read_group(m_read_group);

      m_samples.push_back(std::move(s));
    } else if (config.m_allow_multiple_barcodes) {
      m_samples.back().add(row.sequence_1, row.sequence_2);
    } else {
      std::ostringstream error;
      error << "Duplicate sample name " << log_escape(row.name)
            << "; combining different barcodes for one sample is not "
               "supported. Please ensure that all sample names are unique!";

      throw parsing_error(error.str());
    }
  }

  // Check before adding reversed barcodes, to prevent misleading error messages
  check_barcodes_sequences(m_samples, filename, config.m_paired_end_mode);

  if (!config.m_unidirectional_barcodes) {
    add_reversed_barcodes(config);
  }
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
