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
#ifndef AR_ADAPTERS_H
#define AR_ADAPTERS_H

#include "commontypes.h"
#include "fastq.h"

namespace ar
{

/**
 * Class for reading sets of adapters and barcodes, and for generating
 * per-barcode sets of adapter sequences as needed. The class further checks
 * for the correctness of these sequences, and detects duplicate barcode
 * sequences / pairs of sequences.
 */
class adapter_set
{
public:
    /** Initialize empty adapter list. */
    adapter_set();

    /**
     * Adds a pair of adapters to the set; by default it is assumed that
     * the adapter 2 sequence is in read orientation (e.g. can be found as is
     * in the raw mate 2 reads).
     */
    void add_adapters(const std::string& adapter1,
                      const std::string& adapter2,
                      bool adapter2_read_orientation = true);

    /**
     * Loads barcodes from a table, returning true on success. The value of
     * 'paired_end_mode' is used to set the expected number of values.
     */
    bool load_adapters(const std::string& filename, bool paired_end_mode);

    /**
     * Loads barcodes from a table, returning true on success. The value of
     * 'paired_end_mode' to correctly identify duplicate sequences.
     */
    bool load_barcodes(const std::string& filename, bool paired_end_mode);

    /** Returns the number of adapters per set. */
    size_t adapter_count() const;

    /** Returns the number of adapter sets; namely 1 or barcode_count() */
    size_t adapter_set_count() const;

    /** Returns the number of barcodes. */
    size_t barcode_count() const;

    /**
     * Returns the nth set of adapters; when barcodes are specified, the
     * raw adapters are merged with the 'nth' barcodes. If no barcodes are
     * set, the get_adapters(0) == get_raw_adapters().
     */
    fastq_pair_vec get_adapter_set(size_t nth) const;

    /** Returns get_adapter_set(nth) formatted for printing. */
    string_pair_vec get_pretty_adapter_set(size_t nth) const;

    /** Returns adapter sequences absent of any barcodes. */
    const fastq_pair_vec& get_raw_adapters() const;

    /** Returns the (pairs of) barcodes. */
    const fastq_pair_vec& get_barcodes() const;

    /** Returns the name associated with the nth set of barcodes. */
    const std::string& get_sample_name(size_t nth) const;

private:
    //! Names associated with barcodes
    string_vec m_samples;
    //! User-supplied barcodes
    fastq_pair_vec m_barcodes;
    //! User-supplied adapter sequences, without barcodes added
    fastq_pair_vec m_adapters;
};

} // namespace ar

#endif
