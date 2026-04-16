// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "scheduler.hpp"  // for chunk_vec, chunk_ptr, analytical_step
#include "statistics.hpp" // for trim_stats_ptr, trimming_statistics
#include "threading.hpp"  // for threadstate
#include <cstddef>        // for size_t

namespace adapterremoval {

class sample_output_files;
class userconfig;
enum class read_type;

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

  chunk_vec process(chunk_ptr data) override;
};

class pe_reads_processor : public reads_processor
{
public:
  pe_reads_processor(const userconfig& config,
                     const sample_output_files& output,
                     size_t nth,
                     trim_stats_ptr sink);

  chunk_vec process(chunk_ptr data) override;
};

/** Tracks merged reads, to account for trimming affecting the overlap */
class merged_reads
{
public:
  /** Track trimming for the specified overlapping reads */
  merged_reads(const fastq& read1, const fastq& read2, int offset);

  /** Increments bases trimmed and returns true if both reads were affected by
   * the current operation, in which case statistics should count both */
  bool increment(size_t trim5p, size_t trim3p);

private:
  //! The remaining number of bases that uniquely belong to read 1
  int64_t m_unique_1 = 0;
  //! The remaining number of bases that uniquely belong to read 2
  int64_t m_unique_2 = 0;
};

} // namespace adapterremoval
