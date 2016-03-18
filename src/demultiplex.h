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
#ifndef DEMULTIPLEX_H
#define DEMULTIPLEX_H

#include "fastq.h"
#include "scheduler.h"
#include "statistics.h"
#include "fastq_io.h"

namespace ar
{

class userconfig;
struct demultiplexer_node;

typedef std::vector<demultiplexer_node> demux_node_vec;



/**
 * Baseclass for demultiplexing of reads; responsible for building the quadtree
 * representing the set of adapter sequences, and for maintaining the cache of
 * demultiplexed reads.
 */
class demultiplex_reads : public analytical_step
{
public:
    /** Setup demultiplexer; keeps pointer to config object. */
    demultiplex_reads(const userconfig* config);

    /** Frees any unflushed caches. */
    virtual ~demultiplex_reads();

    /** Returns a statistics object summarizing the results up till now. */
    demux_statistics statistics() const;

protected:
    /**
     * Returns the id of the best matching barcode(s), or -1 if no matches were
     * found or if no single best match was found.
     */
    int select_barcode(const fastq& read_r1, const fastq& read_r2);

    //! List of barcode (pairs) supplied by caller
    const fastq_pair_vec& m_barcodes;
    //! Quadtree representing all mate 1 adapters; for search with n mismatches
    const demux_node_vec m_tree;
    //! Maximum number of mismatches allowed between the mate 1 and mate 2 read
    const size_t m_max_mismatches;
    //! Maximum number of mismatches allowed for the mate 1 read
    const size_t m_max_mismatches_r1;
    //! Maximum number of mismatches allowed for the mate 2 read
    const size_t m_max_mismatches_r2;
    //! Pointer to user settings used for output format for unidentified reads
    const userconfig* m_config;

    //! Returns a chunk-list with any set of reads exceeding the max cache size
    //! If 'eof' is true, all chunks are returned, and the 'eof' values in the
    //! chunks are set to true.
    chunk_vec flush_cache(bool eof = false);

    typedef std::vector<read_chunk_ptr> demultiplexed_cache;

    //! Cache of demultiplex reads; used to reduce the number of output chunks
    //! generated from each processed chunk, which would otherwise increase
    //! linearly with the number of barcodes.
    demultiplexed_cache m_cache;
    //! Cache of unidentified mate 1 reads
    output_chunk_ptr m_unidentified_1;
    //! Cache of unidentified mate 2 reads
    output_chunk_ptr m_unidentified_2;

    //! Sink for demultiplexing statistics; used by subclasses.
    demux_statistics m_statistics;

private:
    //! Not implemented
    demultiplex_reads(const demultiplex_reads&);
    //! Not implemented
    demultiplex_reads& operator=(const demultiplex_reads&);
};


/** Demultiplexer for single-end reads. */
class demultiplex_se_reads : public demultiplex_reads
{
public:
    /** See demultiplex_reads::demultiplex_reads. */
    demultiplex_se_reads(const userconfig* config);

    /**
     * Processes a read chunk, and forwards chunks to downstream steps, with
     * the IDs corresponding to ai_analyses_offset * (nth + 1) for the nth
     * barcode (pair). Unidentified reads are sent to ai_write_unidentified_1.
     */
    chunk_vec process(analytical_chunk* chunk);
};


/** Demultiplexer for paired-end reads. */
class demultiplex_pe_reads : public demultiplex_reads
{
public:
    /** See demultiplex_reads::demultiplex_reads. */
    demultiplex_pe_reads(const userconfig* config);

    /**
     * Processes a read chunk, and forwards chunks to downstream steps, with
     * the IDs corresponding to ai_analyses_offset * (nth + 1) for the nth
     * barcode (pair). Unidentified reads are sent to ai_write_unidentified_1
     * and ai_write_unidentified_2.
     */
    chunk_vec process(analytical_chunk* chunk);
};

} // namespace ar

#endif
