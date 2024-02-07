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
#include "adapterset.hpp"
#include "debug.hpp"      // for AR_REQUIRE
#include "fastq_enc.hpp"  // for fastq_error
#include "linereader.hpp" // for line_reader
#include "logging.hpp"    // for error, log_stream
#include "strutils.hpp"   // for string_vec, indent_lines
#include <algorithm>      // for max, sort
#include <sstream>        // for operator<<, basic_ostream, ostringstream
#include <utility>        // for pair
#include <vector>         // for vector, vector<>::const_iterator

namespace adapterremoval {

using named_fastq_row = std::pair<std::string, fastq_vec>;
using fastq_table = std::vector<named_fastq_row>;
using fastq_table_citer = fastq_table::const_iterator;

bool
print_parse_error(const std::ostringstream& message)
{
  log::error() << "Error reading table:\n" << indent_lines(message.str());

  return false;
}

std::string
trim_comments(std::string line)
{
  const size_t index = line.find('#');
  if (index != std::string::npos) {
    line.resize(index);
  }

  return line;
}

bool
read_table(const std::string& filename,
           fastq_table& dst,
           size_t min_col,
           size_t max_col,
           bool row_names = false)
{
  AR_REQUIRE(min_col <= max_col);
  AR_REQUIRE(min_col >= 1);

  size_t last_row_size = 0;
  size_t line_num = 1;
  size_t col_num = 1;
  try {
    line_reader adapter_file(filename);

    for (std::string line; adapter_file.getline(line); ++line_num) {
      fastq_vec row;
      std::string name;
      std::string field;
      std::istringstream instream(trim_comments(line));

      for (col_num = 1; instream >> field; ++col_num) {
        if (col_num == 1 && row_names) {
          name = field;
        } else {
          row.emplace_back("sequence", field);
        }
      }

      if (name.empty() && row.empty()) {
        // Ignore empty lines, e.g. those containing only comments
        continue;
      } else if (row.size() < min_col) {
        std::ostringstream message;
        message << "Expected at least " << min_col << " columns in "
                << "table '" << filename << "' at line " << line_num
                << ", but only found " << row.size() << " column(s)!";

        return print_parse_error(message);
      } else if (row.size() > max_col) {
        std::ostringstream message;
        message << "Expected at most " << min_col << " columns in "
                << "table '" << filename << "' at line " << line_num
                << ", but found " << row.size() << " column(s)!";

        return print_parse_error(message);
      } else if (last_row_size && last_row_size != row.size()) {
        std::ostringstream message;
        message << "Error reading '" << filename << "' at line " << line_num
                << "; rows contain unequal number of "
                << "columns; last row contained " << last_row_size
                << " column(s) but current row contains " << row.size()
                << " column(s)!";

        return print_parse_error(message);
      } else {
        last_row_size = row.size();
      }

      dst.emplace_back(name, row);
    }
  } catch (const std::ios_base::failure& error) {
    std::ostringstream message;
    message << "IO error reading '" << filename << "' at line " << line_num
            << ": " << error.what();

    return print_parse_error(message);
  } catch (const fastq_error& error) {
    std::ostringstream message;
    message << "Failed to parse sequence in '" << filename << "' at line "
            << line_num << ", column " << col_num << ": " << error.what();

    return print_parse_error(message);
  }

  return true;
}

bool
check_barcodes_sequences(const fastq_pair_vec& barcodes,
                         const std::string& filename,
                         bool paired_end = false)
{
  if (barcodes.empty()) {
    return true;
  }

  const size_t mate_1_len = barcodes.front().first.length();
  const size_t mate_2_len = barcodes.front().second.length();

  string_pair_vec sequences;
  for (auto it = barcodes.begin(); it != barcodes.end(); ++it) {
    const std::string& mate_1 = it->first.sequence();
    const std::string& mate_2 = it->second.sequence();

    if (mate_1.find('N') != std::string::npos) {
      std::ostringstream error;
      error << "Degenerate base (N) found in mate 1 barcode sequence '"
            << mate_1 << "'. Degenerate bases are not supported for "
            << "demultiplexing; please remove before continuing!";

      return print_parse_error(error);
    } else if (mate_2.find('N') != std::string::npos) {
      std::ostringstream error;
      error << "Degenerate base (N) found in mate 2 barcode sequence '"
            << mate_2 << "'. Degenerate bases are not supported for "
            << "demultiplexing; please remove before continuing!";

      return print_parse_error(error);
    } else if (mate_1.length() != mate_1_len) {
      std::ostringstream error;
      error << "Inconsistent mate 1 barcode lengths found; last barcode "
               "was "
            << mate_1_len << " basepairs long, but barcode "
            << (it - barcodes.begin()) + 1 << " mate 1 sequence is "
            << mate_1.length()
            << " basepairs long! Variable length "
               "barcodes are not supported!";

      return print_parse_error(error);
    } else if (mate_2.length() != mate_2_len) {
      std::ostringstream error;
      error << "Inconsistent mate 2 barcode lengths found; last barcode "
               "was "
            << mate_2_len << " basepairs long, but barcode "
            << (it - barcodes.begin()) + 1 << " mate 2 sequence is "
            << mate_2.length()
            << " basepairs long! Variable length "
               "barcodes are not supported!";

      return print_parse_error(error);
    }

    sequences.emplace_back(it->first.sequence(), it->second.sequence());
  }

  std::sort(sequences.begin(), sequences.end());
  string_pair_vec::const_iterator prev = sequences.begin();
  string_pair_vec::const_iterator curr = prev + 1;

  for (; curr != sequences.end(); ++prev, ++curr) {
    if (prev->first == curr->first) {
      if (paired_end) {
        if (prev->second == curr->second) {
          std::ostringstream error;
          error << "Duplicate barcode pairs found in '" << filename
                << "' with barcodes " << prev->first << " and " << prev->second
                << ". please verify "
                   "correctness of the barcode table and remove any "
                   "duplicate entries!";

          return print_parse_error(error);
        }
      } else {
        std::ostringstream error;
        error << "Duplicate mate 1 barcodes found in '" << filename
              << "': " << prev->first
              << ". Even if these "
                 "are associated with different mate 2 barcodes, it "
                 "is not possible to distinguish between these in "
                 "single-end mode!";

        return print_parse_error(error);
      }
    }
  }

  return true;
}

bool
valid_sample_name(const std::string& name)
{
  for (const char c : name) {
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
        (c >= '0' && c <= '9') || (c == '_')) {
      continue;
    }

    std::ostringstream error;
    error << "The sample name '" << name
          << "' is not a valid sample "
             "name; only letters ('a' to 'z' and 'A' to 'Z'), numbers (0 "
             "to 9) and underscores (_) are allowed.";

    return print_parse_error(error);
  }

  if (name == "unidentified") {
    std::ostringstream error;
    error << "The sample name '" << name
          << "' is a reserved name, and "
             "cannot be used!";

    return print_parse_error(error);
  }

  return true;
}

bool
check_sample_names(const string_vec& names)
{
  if (names.empty()) {
    return true;
  }

  for (const auto& name : names) {
    if (!valid_sample_name(name)) {
      return false;
    }
  }

  string_vec sorted_names = names;
  std::sort(sorted_names.begin(), sorted_names.end());

  string_vec::const_iterator prev = sorted_names.begin();
  string_vec::const_iterator curr = prev + 1;
  for (; curr != sorted_names.end(); ++prev, ++curr) {
    if (*prev == *curr) {
      std::ostringstream error;
      error << "Duplicate sample name '" << *prev
            << "'; combining "
               "different barcodes for one sample is not supported. "
               "Please ensure that all sample names are unique!";

      return print_parse_error(error);
    }
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'adapters' class

adapter_set::adapter_set()
  : m_samples()
  , m_barcodes()
  , m_adapters()
{
  // Default name if no barcodes are used
  m_samples.emplace_back("main");
}

void
adapter_set::add_adapters(const std::string& adapter1,
                          const std::string& adapter2)
{
  fastq adapter1_fq("adapter1", adapter1);
  fastq adapter2_fq("adapter2", adapter2);
  adapter2_fq.reverse_complement();

  m_adapters.emplace_back(adapter1_fq, adapter2_fq);
}

bool
adapter_set::load_adapters(const std::string& filename, bool paired_end)
{
  fastq_table raw_adapters;
  if (!read_table(filename, raw_adapters, paired_end ? 2 : 1, 2)) {
    return false;
  }

  for (const auto& row : raw_adapters) {
    fastq adapter_5p = row.second.at(0);
    fastq adapter_3p;

    if (row.second.size() > 1) {
      adapter_3p = row.second.at(1);
      adapter_3p.reverse_complement();
    }

    m_adapters.emplace_back(adapter_5p, adapter_3p);
  }

  return true;
}

bool
adapter_set::load_barcodes(const std::string& filename, bool paired_end)
{
  fastq_table raw_barcodes;
  if (!read_table(filename, raw_barcodes, 1, 2, true)) {
    return false;
  }

  m_samples.clear();
  m_barcodes.clear();

  for (const auto& row : raw_barcodes) {
    fastq barcode_5p = row.second.at(0);
    fastq barcode_3p;

    if (row.second.size() > 1) {
      barcode_3p = row.second.at(1);
    }

    m_samples.push_back(row.first);
    m_barcodes.emplace_back(barcode_5p, barcode_3p);
  }

  return check_barcodes_sequences(m_barcodes, filename, paired_end) &&
         check_sample_names(m_samples);
}

size_t
adapter_set::adapter_count() const
{
  return m_adapters.size();
}

size_t
adapter_set::adapter_set_count() const
{
  if (m_barcodes.empty()) {
    return 1;
  } else {
    return barcode_count();
  }
}

size_t
adapter_set::barcode_count() const
{
  return m_barcodes.size();
}

fastq_pair_vec
adapter_set::get_adapter_set(size_t nth) const
{
  AR_REQUIRE(nth == 0 || !m_barcodes.empty(),
             "tried to get non-existing adapter set");

  if (m_barcodes.empty()) {
    return m_adapters;
  }

  fastq_pair barcodes = m_barcodes.at(nth);
  barcodes.second.reverse_complement();

  fastq_pair_vec adapters;
  for (const auto& adapter_pair : m_adapters) {
    fastq adapter_1("adapter_1",
                    barcodes.second.sequence() + adapter_pair.first.sequence());
    fastq adapter_2("adapter_2",
                    adapter_pair.second.sequence() + barcodes.first.sequence());

    adapters.emplace_back(adapter_1, adapter_2);
  }

  return adapters;
}

const fastq_pair_vec&
adapter_set::get_raw_adapters() const
{
  return m_adapters;
}

const fastq_pair_vec&
adapter_set::get_barcodes() const
{
  return m_barcodes;
}

const std::string&
adapter_set::get_sample_name(size_t nth) const
{
  return m_samples.at(nth);
}

} // namespace adapterremoval
