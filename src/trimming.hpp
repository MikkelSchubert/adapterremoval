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
#ifndef TRIMMING_H
#define TRIMMING_H

#include <random>

#include "fastq.hpp"
#include "scheduler.hpp"
#include "statistics.hpp"

namespace ar {

class userconfig;

typedef std::unique_ptr<statistics> statistics_ptr;

class reads_processor : public analytical_step
{
public:
  reads_processor(const userconfig& config, size_t nth);

  statistics_ptr get_final_statistics();

protected:
  class stats_sink : public statistics_sink<statistics>
  {
  public:
    stats_sink(const userconfig& config);

  protected:
    virtual pointer new_sink() const;

    virtual void reduce(pointer& dst, const pointer& src) const;

    const userconfig& m_config;
  };

  const userconfig& m_config;
  const fastq_pair_vec m_adapters;
  stats_sink m_stats;
  const size_t m_nth;
};

class se_reads_processor : public reads_processor
{
public:
  se_reads_processor(const userconfig& config, size_t nth = 0);

  chunk_vec process(analytical_chunk* chunk);
};

/** Class for building RNGs on demand. */
class rng_sink : public statistics_sink<std::mt19937>
{
public:
  rng_sink(unsigned seed);

protected:
  virtual pointer new_sink() const;

  virtual void reduce(pointer&, const pointer&) const;

  //! Copy construction not supported
  rng_sink(const rng_sink&) = delete;
  //! Assignment not supported
  rng_sink& operator=(const rng_sink&) = delete;

private:
  mutable std::mt19937 m_seed;
};

class pe_reads_processor : public reads_processor
{
public:
  pe_reads_processor(const userconfig& config, size_t nth);

  chunk_vec process(analytical_chunk* chunk);

private:
  rng_sink m_rngs;
};

}

#endif
