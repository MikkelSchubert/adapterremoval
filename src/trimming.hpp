// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "scheduler.hpp"     // for chunk_vec, chunk_ptr, threadstate, analyt...
#include "sequence_sets.hpp" // for adapter_set
#include "statistics.hpp"    // for trimming_statistics, trim_stats_ptr
#include <cstddef>           // for size_t
#include <vector>            // for vector

namespace adapterremoval {

enum class read_type : size_t;
class sample_output_files;
class userconfig;

class reads_processor : public analytical_step
{
public:
  reads_processor(const userconfig& config,
                  const sample_output_files& output,
                  size_t nth,
                  trim_stats_ptr sink);

  void finalize() override;

protected:
  const userconfig& m_config;
  const sample_output_files& m_output;
  const size_t m_sample;

  threadstate<trimming_statistics> m_stats{};
  trim_stats_ptr m_stats_sink;
};

class se_reads_processor : public reads_processor
{
public:
  se_reads_processor(const userconfig& config,
                     const sample_output_files& output,
                     size_t nth,
                     trim_stats_ptr sink);

  chunk_vec process(chunk_ptr chunk) override;
};

class pe_reads_processor : public reads_processor
{
public:
  pe_reads_processor(const userconfig& config,
                     const sample_output_files& output,
                     size_t nth,
                     trim_stats_ptr sink);

  chunk_vec process(chunk_ptr chunk) override;
};

} // namespace adapterremoval
