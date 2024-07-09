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
#pragma once

#include "fastq.hpp"      // for fastq_pair_vec
#include "fastq_io.hpp"   // for output_chunk_ptr
#include "scheduler.hpp"  // for chunk_vec, chunk_ptr, threadstate, analyt...
#include "statistics.hpp" // for trimming_statistics, trim_stats_ptr
#include <cstddef>        // for size_t
#include <random>         // for mt19937
#include <vector>         // for vector

namespace adapterremoval {

enum class read_type : size_t;
class output_sample_files;
class userconfig;

/** Helper class used to generate per file-type chunks for processed reads . */
class trimmed_reads
{
public:
  trimmed_reads(const output_sample_files& map, const bool eof);

  /**
   * Adds a read of the given type.
   *
   * @param read A read to be distributed in the pipeline; may be modified.
   * @param type The read type to store the read as.
   */
  void add(const fastq& read, const read_type type);

  /** Returns a chunk for each generated type of processed reads. */
  chunk_vec finalize();

private:
  const output_sample_files& m_map;

  //! A set output chunks being created; typically fewer than read_type::max.
  std::vector<output_chunk_ptr> m_chunks;
};

class reads_processor : public analytical_step
{
public:
  reads_processor(const userconfig& config,
                  const output_sample_files& output,
                  const size_t nth,
                  trim_stats_ptr sink);

  void finalize() override;

protected:
  const userconfig& m_config;
  const fastq_pair_vec m_adapters;
  const output_sample_files& m_output;
  const size_t m_nth;

  threadstate<trimming_statistics> m_stats;
  trim_stats_ptr m_stats_sink;
};

class se_reads_processor : public reads_processor
{
public:
  se_reads_processor(const userconfig& config,
                     const output_sample_files& output,
                     const size_t nth,
                     trim_stats_ptr sink);

  chunk_vec process(chunk_ptr chunk) override;
};

class pe_reads_processor : public reads_processor
{
public:
  pe_reads_processor(const userconfig& config,
                     const output_sample_files& output,
                     const size_t nth,
                     trim_stats_ptr sink);

  chunk_vec process(chunk_ptr chunk) override;

private:
  threadstate<std::mt19937> m_rngs;
};

} // namespace adapterremoval
