/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <algorithm>
#include <iostream>
#include <sstream>

#include "adapterset.h"
#include "debug.h"
#include "linereader.h"
#include "strutils.h"

namespace ar
{

typedef std::pair<std::string, fastq_vec> named_fastq_row;
typedef std::vector<named_fastq_row> fastq_table;
typedef fastq_table::const_iterator fastq_table_citer;


bool print_parse_error(const std::stringstream& message)
{
    std::cerr << "ERROR READING TABLE:\n"
              << cli_formatter::fmt(message.str()) << std::endl;

    return false;
}


std::string trim_comments(std::string line)
{
    const size_t index = line.find('#');
    if (index != std::string::npos) {
        line.resize(index);
    }

    return line;
}


bool read_table(const std::string& filename, fastq_table& dst,
                size_t min_col, size_t max_col,
                bool row_names = false)
{
    AR_DEBUG_ASSERT(min_col <= max_col);
    AR_DEBUG_ASSERT(min_col >= 1);

    size_t last_row_size = 0;
    size_t line_num = 1;
    try {
        line_reader adapter_file(filename);

        for (std::string line; adapter_file.getline(line); ++line_num) {
            fastq_vec row;
            std::string name;
            std::string field;
            std::stringstream instream(trim_comments(line));

            for (size_t index = 1; instream >> field; ++index) {
                try {
                    if (index == 1 && row_names) {
                        name = field;
                    } else {
                        row.push_back(fastq("sequence", field));
                    }
                } catch (const fastq_error& error) {
                    std::stringstream message;
                    message << "Failed to parse sequence in '" << filename
                            << "' at line " << line_num << ", column " << index
                            << ": " << error.what();

                    return print_parse_error(message);
                }
            }

            if (name.empty() && row.empty()) {
                // Ignore empty lines, e.g. those containing only comments
                continue;
            } else if (row.size() < min_col) {
                std::stringstream message;
                message << "Expected at least " << min_col << " columns in "
                        << "table '" << filename << "' at line " << line_num
                        << ", but only found " << row.size() << " column(s)!";

                return print_parse_error(message);
            } else if (row.size() > max_col) {
                std::stringstream message;
                message << "Expected at most " << min_col << " columns in "
                        << "table '" << filename << "' at line " << line_num
                        << ", but found " << row.size() << " column(s)!";

                return print_parse_error(message);
            } else if (last_row_size && last_row_size != row.size()) {
                std::stringstream message;
                message << "Error reading '" << filename << "' at line "
                          << line_num << "; rows contain unequal number of "
                          << "columns; last row contained " << last_row_size
                          << " column(s) but current row contains "
                          << row.size() << " column(s)!";

                return print_parse_error(message);
            } else {
                last_row_size = row.size();
            }

            dst.push_back(named_fastq_row(name, row));
        }
    } catch (const std::ios_base::failure& error) {
        std::stringstream message;
        message << "IO error reading '" << filename << "' at line "
                << line_num << ": " << error.what();

        return print_parse_error(message);
    }

    return true;
}


bool check_barcodes_sequences(const fastq_pair_vec& barcodes,
                              const std::string& filename,
                              bool paired_end = false)
{
    if (barcodes.empty()) {
        return true;
    }

    const size_t mate_1_len = barcodes.front().first.length();
    const size_t mate_2_len = barcodes.front().second.length();

    string_pair_vec sequences;
    for (fastq_pair_vec::const_iterator it = barcodes.begin(); it != barcodes.end(); ++it) {
        const std::string& mate_1 = it->first.sequence();
        const std::string& mate_2 = it->second.sequence();

        if (mate_1.find('N') != std::string::npos) {
            std::stringstream error;
            error << "Degenerate base (N) found in mate 1 barcode sequence '"
                  << mate_1 << "'. Degenerate bases are not supported for "
                  << "demultiplexing; please remove before continuing!";

            return print_parse_error(error);
        } else if (mate_2.find('N') != std::string::npos) {
            std::stringstream error;
            error << "Degenerate base (N) found in mate 2 barcode sequence '"
                  << mate_2 << "'. Degenerate bases are not supported for "
                  << "demultiplexing; please remove before continuing!";

            return print_parse_error(error);
        } else if (mate_1.length() != mate_1_len) {
            std::stringstream error;
            error << "Inconsistent mate 1 barcode lengths found; last barcode "
                     "was " << mate_1_len << " basepairs long, but barcode "
                  << (it - barcodes.begin()) + 1 << " mate 1 sequence is "
                  << mate_1.length() << " basepairs long! Variable length "
                     "barcodes are not supported!";

            return print_parse_error(error);
        } else if (mate_2.length() != mate_2_len) {
            std::stringstream error;
            error << "Inconsistent mate 2 barcode lengths found; last barcode "
                     "was " << mate_2_len << " basepairs long, but barcode "
                  << (it - barcodes.begin()) + 1 << " mate 2 sequence is "
                  << mate_2.length() << " basepairs long! Variable length "
                     "barcodes are not supported!";

            return print_parse_error(error);
        }

        sequences.push_back(string_pair(it->first.sequence(),
                                        it->second.sequence()));
    }

    std::sort(sequences.begin(), sequences.end());
    string_pair_vec::const_iterator prev = sequences.begin();
    string_pair_vec::const_iterator curr = prev + 1;

    for (; curr != sequences.end(); ++prev, ++curr) {
        if (prev->first == curr->first) {
            if (paired_end) {
                if (prev->second == curr->second) {
                    std::stringstream error;
                    error << "Duplicate barcode pairs found in '"
                          << filename << "' with barcodes "<< prev->first
                          << " and " << prev->second << ". please verify "
                             "correctness of the barcode table and remove any "
                             "duplicate entries!";

                    return print_parse_error(error);
                }
            } else {
                std::stringstream error;
                error << "Duplicate mate 1 barcodes found in '"
                      << filename << "': "<< prev->first << ". Even if these "
                         "are assosiated with different mate 2 barcodes, it "
                         "is not possible to distinguish between these in "
                         "single-end mode!";

                return print_parse_error(error);
            }
        }
    }

    return true;
}


bool valid_sample_name(const std::string& name)
{
    std::string::const_iterator it = name.begin();
    for (; it != name.end(); ++it) {
        const char c = *it;

        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || (c == '_')) {
            continue;
        }

        std::stringstream error;
        error << "The sample name '" << name << "' is not a valid sample "
                 "name; only letters ('a' to 'z' and 'A' to 'Z'), numbers (0 "
                 "to 9) and underscores (_) are allowed.";

        return print_parse_error(error);
    }

    if (name == "unidentified") {
        std::stringstream error;
        error << "The sample name '" << name << "' is a reserved name, and "
                 "cannot be used!";

        return print_parse_error(error);
    }

    return true;
}


bool check_sample_names(const string_vec& names)
{
    if (names.empty()) {
        return true;
    }

    for (string_vec_citer it = names.begin(); it != names.end(); ++it) {
        if (!valid_sample_name(*it)) {
            return false;
        }
    }

    string_vec sorted_names = names;
    std::sort(sorted_names.begin(), sorted_names.end());

    string_vec::const_iterator prev = sorted_names.begin();
    string_vec::const_iterator curr = prev + 1;
    for (; curr != sorted_names.end(); ++prev, ++curr) {
        if (*prev == *curr) {
            std::stringstream error;
            error << "Duplicate sample name '" << *prev << "'; combining "
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
    m_samples.push_back("main");
}


void adapter_set::add_adapters(const std::string& adapter1,
                               const std::string& adapter2,
                               bool adapter2_read_orientation)
{
    fastq adapter1_fq("adapter1", adapter1);
    fastq adapter2_fq("adapter2", adapter2);

    if (adapter2_read_orientation) {
        adapter2_fq.reverse_complement();
    }

    m_adapters.push_back(fastq_pair(adapter1_fq, adapter2_fq));
}


bool adapter_set::load_adapters(const std::string& filename, bool paired_end)
{
    fastq_table raw_adapters;
    if (!read_table(filename, raw_adapters, paired_end ? 2 : 1, 2)) {
        return false;
    }

    for (fastq_table_citer it = raw_adapters.begin(); it != raw_adapters.end(); ++it) {
        fastq adapter_5p = it->second.at(0);
        fastq adapter_3p;

        if (it->second.size() > 1) {
            adapter_3p = it->second.at(1);
        }

        m_adapters.push_back(fastq_pair(adapter_5p, adapter_3p));
    }

    return true;
}


bool adapter_set::load_barcodes(const std::string& filename, bool paired_end)
{
    fastq_table raw_barcodes;
    if (!read_table(filename, raw_barcodes, 1, 2, true)) {
        return false;
    }

    m_samples.clear();
    m_barcodes.clear();

    for (fastq_table_citer it = raw_barcodes.begin(); it != raw_barcodes.end(); ++it) {
        fastq barcode_5p = it->second.at(0);
        fastq barcode_3p;

        if (it->second.size() > 1) {
            barcode_3p = it->second.at(1);
        }

        m_samples.push_back(it->first);
        m_barcodes.push_back(fastq_pair(barcode_5p, barcode_3p));
    }

    return check_barcodes_sequences(m_barcodes, filename, paired_end)
        && check_sample_names(m_samples);
}


size_t adapter_set::adapter_count() const
{
    return m_adapters.size();
}


size_t adapter_set::adapter_set_count() const
{
    if (m_barcodes.empty()) {
        return 1;
    } else {
        return barcode_count();
    }
}


size_t adapter_set::barcode_count() const
{
    return m_barcodes.size();
}


fastq_pair_vec adapter_set::get_adapter_set(size_t nth) const
{
    if (m_barcodes.empty()) {
        AR_DEBUG_ASSERT(nth == 0);
        return m_adapters;
    }

    fastq_pair barcodes = m_barcodes.at(nth);
    barcodes.second.reverse_complement();

    fastq_pair_vec adapters;
    for (fastq_pair_vec::const_iterator it = m_adapters.begin(); it != m_adapters.end(); ++it) {
        fastq adapter_1("adapter_1", barcodes.second.sequence() + it->first.sequence());
        fastq adapter_2("adapter_2", it->second.sequence() + barcodes.first.sequence());

        adapters.push_back(fastq_pair(adapter_1, adapter_2));
    }

    return adapters;
}


string_pair_vec adapter_set::get_pretty_adapter_set(size_t nth) const
{
    fastq_pair barcodes;
    if (m_barcodes.empty()) {
        AR_DEBUG_ASSERT(nth == 0);
    } else {
        barcodes = m_barcodes.at(nth);
    }

    string_pair_vec adapters;
    const fastq_pair_vec adapter_pairs = get_adapter_set(nth);
    for (fastq_pair_vec::const_iterator it = adapter_pairs.begin(); it != adapter_pairs.end(); ++it) {
        fastq adapter_1 = it->first;
        fastq adapter_2 = it->second;
        adapter_2.reverse_complement();

        std::string seq_1 = adapter_1.sequence();
        std::string seq_2 = adapter_2.sequence();

        if (barcodes.first.length()) {
            seq_2.insert(barcodes.first.length(), 1, '_');
        }

        if (barcodes.second.length()) {
            seq_1.insert(barcodes.second.length(), 1, '_');
        }

        adapters.push_back(string_pair(seq_1, seq_2));
    }

    return adapters;
}


const fastq_pair_vec& adapter_set::get_raw_adapters() const
{
    return m_adapters;
}


const fastq_pair_vec& adapter_set::get_barcodes() const
{
    return m_barcodes;
}

const std::string& adapter_set::get_sample_name(size_t nth) const
{
    return m_samples.at(nth);
}

} // namespace ar
