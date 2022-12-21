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
#include <algorithm> // for min
#include <memory>    // for make_unique
#include <stdexcept> // for invalid_argument
#include <utility>   // for pair, move

#include "adapterset.hpp" // for adapter_set
#include "alignment.hpp"  // for alignment_info, sequence_merger, align_pai...
#include "counts.hpp"     // for counts
#include "debug.hpp"      // for AR_REQUIRE
#include "fastq_io.hpp"   // for output_chunk_ptr, fastq_read_chunk, fastq_...
#include "statistics.hpp" // for trimming_statistics, fastq_statistics
#include "trimming.hpp"
#include "userconfig.hpp" // for userconfig, output_sample_files, output_sa...

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// Helper functions

namespace {

/** Tracks overlapping reads, to account for trimming affecting the overlap */
class merged_reads
{
public:
  merged_reads(const fastq& read1, const fastq& read2, size_t insert_size)
    : m_len_1(read1.length())
    , m_len_2(read2.length())
    , m_overlap(m_len_1 + m_len_2 - insert_size)
    , m_trimmed_5p()
    , m_trimmed_3p()
  {
  }

  /** Increments bases trimmed and returns true if overlap was trimmed */
  bool increment(const size_t trim5p, const size_t trim3p)
  {
    if (m_trimmed_5p + m_trimmed_3p < m_len_1 + m_len_2 - m_overlap) {
      m_trimmed_5p += trim5p;
      m_trimmed_3p += trim3p;

      return (trim5p && trim3p) ||
             ((trim5p || trim3p) && ((m_trimmed_5p > m_len_1 - m_overlap) ||
                                     (m_trimmed_3p > m_len_2 - m_overlap)));
    } else {
      return false;
    }
  }

private:
  const size_t m_len_1;
  const size_t m_len_2;
  const size_t m_overlap;

  size_t m_trimmed_5p;
  size_t m_trimmed_3p;
};

/** Trims poly-X tails from sequence prior to adapter trimming **/
void
pre_trim_poly_x_tail(const userconfig& config,
                     trimming_statistics& stats,
                     fastq& read)
{
  if (config.pre_trim_poly_x.size()) {
    const auto result = read.poly_x_trimming(config.pre_trim_poly_x,
                                             config.trim_poly_x_threshold);

    if (result.second) {
      stats.poly_x_pre_trimmed_reads.inc(result.first);
      stats.poly_x_pre_trimmed_bases.inc(result.first, result.second);
    }
  }
}

/** Trims poly-X tails from sequence after adapter trimming **/
void
post_trim_poly_x_tail(const userconfig& config,
                      trimming_statistics& stats,
                      fastq& read)
{
  if (config.post_trim_poly_x.size()) {
    const auto result = read.poly_x_trimming(config.post_trim_poly_x,
                                             config.trim_poly_x_threshold);

    if (result.second) {
      stats.poly_x_post_trimmed_reads.inc(result.first);
      stats.poly_x_post_trimmed_bases.inc(result.first, result.second);
    }
  }
}

/** Trims fixed bases from read termini and returns the number trimmed. **/
void
trim_read_termini(reads_and_bases& stats,
                  fastq& read,
                  size_t trim_5p,
                  size_t trim_3p)
{
  const auto length = read.length();
  if ((trim_5p || trim_3p) && length) {
    if (trim_5p + trim_3p < length) {
      read.truncate(trim_5p, length - trim_5p - trim_3p);
    } else {
      read.truncate(0, 0);
    }

    stats.inc(length - read.length());
  }
}

/** Trims fixed number of 5'/3' bases prior to adapter trimming */
void
pre_trim_read_termini(const userconfig& config,
                      trimming_statistics& stats,
                      fastq& read,
                      read_type type)
{
  size_t trim_5p = 0;
  size_t trim_3p = 0;

  if (type == read_type::mate_1) {
    trim_5p = config.pre_trim_fixed_5p.first;
    trim_3p = config.pre_trim_fixed_3p.first;
  } else if (type == read_type::mate_2) {
    trim_5p = config.pre_trim_fixed_5p.second;
    trim_3p = config.pre_trim_fixed_3p.second;
  } else {
    AR_FAIL("invalid read type in pre_trim_read_termini");
  }

  trim_read_termini(stats.terminal_pre_trimmed, read, trim_5p, trim_3p);
}

/** Trims fixed number of 5'/3' bases after adapter trimming */
void
post_trim_read_termini(const userconfig& config,
                       trimming_statistics& stats,
                       fastq& read,
                       read_type type,
                       merged_reads* mstats = nullptr)
{
  size_t trim_5p = 0;
  size_t trim_3p = 0;

  switch (type) {
    case read_type::mate_1:
      trim_5p = config.post_trim_fixed_5p.first;
      trim_3p = config.post_trim_fixed_3p.first;
      break;

    case read_type::mate_2:
      trim_5p = config.post_trim_fixed_5p.second;
      trim_3p = config.post_trim_fixed_3p.second;
      break;

    case read_type::merged:
      AR_REQUIRE(config.paired_ended_mode);
      trim_5p = config.post_trim_fixed_5p.first;
      trim_3p = config.post_trim_fixed_5p.second;
      break;

    default:
      AR_FAIL("invalid read type in post_trim_read_termini");
  }

  trim_read_termini(stats.terminal_post_trimmed, read, trim_5p, trim_3p);
  if (mstats && mstats->increment(trim_5p, trim_3p)) {
    stats.terminal_post_trimmed.inc_reads();
  }
}

/** Quality trims a read after adapter trimming */
void
post_trim_read_by_quality(const userconfig& config,
                          trimming_statistics& stats,
                          fastq& read,
                          merged_reads* mstats = nullptr)
{
  fastq::ntrimmed trimmed;
  switch (config.trim) {
    case trimming_strategy::none:
      break;

    case trimming_strategy::mott:
      trimmed = read.mott_trimming(config.trim_mott_rate, config.preserve5p);
      break;

    case trimming_strategy::window:
      trimmed = read.trim_windowed_bases(config.trim_ambiguous_bases,
                                         config.trim_quality_score,
                                         config.trim_window_length,
                                         config.preserve5p);
      break;

    case trimming_strategy::per_base:
      trimmed = read.trim_trailing_bases(
        config.trim_ambiguous_bases,
        config.trim_low_quality_bases ? config.trim_quality_score : -1,
        config.preserve5p);
      break;

    default:
      AR_FAIL("not implemented");
  }

  if (trimmed.first || trimmed.second) {
    stats.low_quality_trimmed.inc(trimmed.first + trimmed.second);

    if (mstats && mstats->increment(trimmed.first, trimmed.second)) {
      stats.low_quality_trimmed.inc_reads();
    }
  }
}

bool
is_acceptable_read(const userconfig& config,
                   trimming_statistics& stats,
                   const fastq& seq)
{
  const auto length = seq.length();

  if (length < config.min_genomic_length) {
    stats.filtered_min_length.inc(length);

    return false;
  } else if (length > config.max_genomic_length) {
    stats.filtered_max_length.inc(length);

    return false;
  }

  const auto max_n = config.max_ambiguous_bases;
  if (max_n < length && seq.count_ns() > max_n) {
    stats.filtered_ambiguous.inc(length);
    return false;
  }

  if (config.min_complexity > 0.0 && seq.complexity() < config.min_complexity) {
    stats.filtered_low_complexity.inc(length);
    return false;
  }

  return true;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// Implementations for `trimmed_reads`

trimmed_reads::trimmed_reads(const output_sample_files& map, const bool eof)
  : m_map(map)
  , m_chunks()
{
  const auto& filenames = map.filenames();
  const auto& pipeline_steps = map.pipeline_steps();
  AR_REQUIRE(filenames.size() == pipeline_steps.size());

  for (size_t i = 0; i < pipeline_steps.size(); ++i) {
    m_chunks.emplace_back(std::make_unique<fastq_output_chunk>(eof));
  }
}

void
trimmed_reads::add(const fastq& read, const read_type type)
{
  const size_t offset = m_map.offset(type);
  if (offset != output_sample_files::disabled) {
    m_chunks.at(offset)->add(read);
  }
}

chunk_vec
trimmed_reads::finalize()
{
  chunk_vec chunks;

  for (size_t i = 0; i < m_chunks.size(); ++i) {
    chunks.emplace_back(m_map.pipeline_steps().at(i),
                        std::move(m_chunks.at(i)));
  }

  return chunks;
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `reads_processor`

reads_processor::reads_processor(const userconfig& config,
                                 const output_sample_files& output,
                                 const size_t nth,
                                 trim_stats_ptr sink)
  : analytical_step(processing_order::unordered, "reads_processor")
  , m_config(config)
  , m_adapters(config.adapters.get_adapter_set(nth))
  , m_output(output)
  , m_nth(nth)
  , m_stats()
  , m_stats_sink(sink)
{
  AR_REQUIRE(m_stats_sink);

  for (size_t i = 0; i < m_config.max_threads; ++i) {
    m_stats.emplace_back(m_config.report_sample_rate);
  }
}

void
reads_processor::finalize()
{
  while (auto next = m_stats.try_acquire()) {
    *m_stats_sink += *next;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `se_reads_processor`

se_reads_processor::se_reads_processor(const userconfig& config,
                                       const output_sample_files& output,
                                       const size_t nth,
                                       trim_stats_ptr sink)
  : reads_processor(config, output, nth, sink)
{
}

chunk_vec
se_reads_processor::process(chunk_ptr chunk)
{
  auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  trimmed_reads chunks(m_output, read_chunk.eof);

  auto stats = m_stats.acquire();
  stats->adapter_trimmed_reads.resize_up_to(m_config.adapters.adapter_count());
  stats->adapter_trimmed_bases.resize_up_to(m_config.adapters.adapter_count());

  const auto aligner = sequence_aligner(m_adapters);

  for (auto& read : read_chunk.reads_1) {
    const size_t in_length = read.length();

    // Trim fixed number of bases from 5' and/or 3' termini
    pre_trim_read_termini(m_config, *stats, read, read_type::mate_1);
    // Trim poly-X tails for zero or more X
    pre_trim_poly_x_tail(m_config, *stats, read);

    const alignment_info alignment =
      aligner.align_single_end(read, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      const auto length = read.length();
      alignment.truncate_single_end(read);

      stats->adapter_trimmed_reads.inc(alignment.adapter_id);
      stats->adapter_trimmed_bases.inc(alignment.adapter_id,
                                       length - read.length());
    }

    // Add (optional) user specified prefixes to read names
    read.add_prefix_to_name(m_config.prefix_read_1);

    // Trim fixed number of bases from 5' and/or 3' termini
    post_trim_read_termini(m_config, *stats, read, read_type::mate_1);
    // Trim poly-X tails for zero or more X
    post_trim_poly_x_tail(m_config, *stats, read);
    // Sliding window trimming or single-base trimming of low quality bases
    post_trim_read_by_quality(m_config, *stats, read);

    stats->total_trimmed.inc_reads(in_length != read.length());
    stats->total_trimmed.inc_bases(in_length - read.length());

    if (is_acceptable_read(m_config, *stats, read)) {
      stats->read_1->process(read);
      chunks.add(read, read_type::mate_1);
    } else {
      stats->discarded->process(read);
      chunks.add(read, read_type::discarded);
    }
  }

  m_stats.release(stats);

  return chunks.finalize();
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `pe_reads_processor`

void
add_pe_statistics(const trimming_statistics& stats,
                  const fastq& read,
                  read_type type)
{
  switch (type) {
    case read_type::mate_1:
      stats.read_1->process(read);
      break;
    case read_type::mate_2:
      stats.read_2->process(read);
      break;

    case read_type::merged:
      stats.merged->process(read);
      break;

    case read_type::singleton:
      stats.singleton->process(read);
      break;

    case read_type::discarded:
      stats.discarded->process(read);
      break;

    case read_type::max:
    case read_type::unidentified_1:
    case read_type::unidentified_2:
      AR_FAIL("unhandled read type");

    default:
      AR_FAIL("invalid read type");
  }
}

pe_reads_processor::pe_reads_processor(const userconfig& config,
                                       const output_sample_files& output,
                                       const size_t nth,
                                       trim_stats_ptr sink)
  : reads_processor(config, output, nth, sink)
  , m_rngs()
{
  if (m_config.merge == merge_strategy::original) {
    std::mt19937 seed(config.merge_seed);
    for (size_t i = 0; i < m_config.max_threads; ++i) {
      m_rngs.emplace_back(seed());
    }
  }
}

chunk_vec
pe_reads_processor::process(chunk_ptr chunk)
{
  sequence_merger merger;
  merger.set_merge_strategy(m_config.merge);
  merger.set_max_recalculated_score(m_config.merge_quality_max);

  const auto aligner = sequence_aligner(m_adapters);

  auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  trimmed_reads chunks(m_output, read_chunk.eof);

  auto stats = m_stats.acquire();
  stats->adapter_trimmed_reads.resize_up_to(m_config.adapters.adapter_count());
  stats->adapter_trimmed_bases.resize_up_to(m_config.adapters.adapter_count());

  AR_REQUIRE(read_chunk.reads_1.size() == read_chunk.reads_2.size());
  std::unique_ptr<std::mt19937> rng;
  if (m_config.merge == merge_strategy::original) {
    rng = m_rngs.acquire();
    merger.set_rng(rng.get());
  }

  auto it_1 = read_chunk.reads_1.begin();
  auto it_2 = read_chunk.reads_2.begin();
  while (it_1 != read_chunk.reads_1.end()) {
    fastq& read_1 = *it_1++;
    fastq& read_2 = *it_2++;

    const size_t in_length_1 = read_1.length();
    const size_t in_length_2 = read_2.length();

    // Trim fixed number of bases from 5' and/or 3' termini
    pre_trim_read_termini(m_config, *stats, read_1, read_type::mate_1);
    pre_trim_read_termini(m_config, *stats, read_2, read_type::mate_2);

    // Trim poly-X tails for zero or more X
    pre_trim_poly_x_tail(m_config, *stats, read_1);
    pre_trim_poly_x_tail(m_config, *stats, read_2);

    // Reverse complement to match the orientation of read_1
    read_2.reverse_complement();

    const alignment_info alignment =
      aligner.align_paired_end(read_1, read_2, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      const size_t insert_size = alignment.insert_size(read_1, read_2);
      const bool can_merge_alignment = m_config.can_merge_alignment(alignment);
      if (can_merge_alignment) {
        // Insert size calculated from untrimmed reads
        stats->insert_sizes.resize_up_to(insert_size + 1);
        stats->insert_sizes.inc(insert_size);
        stats->overlapping_reads += 2;
      }

      const size_t pre_trimmed_bp = read_1.length() + read_2.length();
      const size_t n_adapters = alignment.truncate_paired_end(read_1, read_2);
      const size_t post_trimmed_bp = read_1.length() + read_2.length();

      stats->adapter_trimmed_reads.inc(alignment.adapter_id, n_adapters);
      stats->adapter_trimmed_bases.inc(alignment.adapter_id,
                                       pre_trimmed_bp - post_trimmed_bp);

      if (m_config.merge != merge_strategy::none && can_merge_alignment) {
        // Track if one or both source reads are trimmed post merging
        merged_reads mstats(read_1, read_2, insert_size);

        // Merge read_2 into read_1
        merger.merge(alignment, read_1, read_2);

        // Track amount of overlapping bases "lost" due to read merging
        stats->reads_merged.inc_reads(2);
        stats->reads_merged.inc_bases(post_trimmed_bp - read_1.length());

        // Add (optional) user specified prefix to read names
        read_1.add_prefix_to_name(m_config.prefix_merged);

        // Trim fixed number of bases from 5' and/or 3' termini
        post_trim_read_termini(
          m_config, *stats, read_1, read_type::merged, &mstats);

        if (!m_config.preserve5p) {
          // A merged read essentially consists of two 5p termini, so neither
          // end should be trimmed if --preserve5p is used.
          post_trim_read_by_quality(m_config, *stats, read_1, &mstats);
        }

        // Track total amount of bases trimmed/lost
        stats->total_trimmed.inc_reads(2);
        stats->total_trimmed.inc_bases(in_length_1 + in_length_2 -
                                       read_1.length());

        if (is_acceptable_read(m_config, *stats, read_1)) {
          stats->merged->process(read_1, 2);
          chunks.add(read_1, read_type::merged);
        } else {
          stats->discarded->process(read_1, 2);
          chunks.add(read_1, read_type::discarded);
        }

        continue;
      }
    }

    // Reads were not aligned or merging is not enabled
    // Undo reverse complementation (post truncation of adapters)
    read_2.reverse_complement();

    // Add (optional) user specified prefixes to read names
    read_1.add_prefix_to_name(m_config.prefix_read_1);
    read_2.add_prefix_to_name(m_config.prefix_read_2);

    // Trim fixed number of bases from 5' and/or 3' termini
    post_trim_read_termini(m_config, *stats, read_1, read_type::mate_1);
    post_trim_read_termini(m_config, *stats, read_2, read_type::mate_2);

    // Trim poly-X tails for zero or more X
    post_trim_poly_x_tail(m_config, *stats, read_1);
    post_trim_poly_x_tail(m_config, *stats, read_2);

    // Sliding window trimming or single-base trimming of low quality bases
    post_trim_read_by_quality(m_config, *stats, read_1);
    post_trim_read_by_quality(m_config, *stats, read_2);

    // Are the reads good enough? Not too many Ns?
    const bool is_ok_1 = is_acceptable_read(m_config, *stats, read_1);
    const bool is_ok_2 = is_acceptable_read(m_config, *stats, read_2);

    read_type type_1;
    read_type type_2;
    if (is_ok_1 && is_ok_2) {
      type_1 = read_type::mate_1;
      type_2 = read_type::mate_2;
    } else if (is_ok_1) {
      type_1 = read_type::singleton;
      type_2 = read_type::discarded;
    } else if (is_ok_2) {
      type_1 = read_type::discarded;
      type_2 = read_type::singleton;
    } else {
      type_1 = read_type::discarded;
      type_2 = read_type::discarded;
    }

    add_pe_statistics(*stats, read_1, type_1);
    add_pe_statistics(*stats, read_2, type_2);

    // Total number of reads/bases trimmed
    stats->total_trimmed.inc_reads((in_length_1 != read_1.length()) +
                                   (in_length_2 != read_2.length()));
    stats->total_trimmed.inc_bases((in_length_1 - read_1.length()) +
                                   (in_length_2 - read_2.length()));

    // Queue reads last, since this result in modifications to lengths
    chunks.add(read_1, type_1);
    chunks.add(read_2, type_2);
  }

  m_stats.release(stats);
  m_rngs.release(rng);

  return chunks.finalize();
}

} // namespace adapterremoval
