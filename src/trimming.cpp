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
#include <algorithm> // for min
#include <stdexcept> // for invalid_argument
#include <utility>   // for pair, move

#include "adapterset.hpp" // for adapter_set
#include "alignment.hpp"  // for alignment_info, sequence_merger, align_pai...
#include "counts.hpp"     // for counts
#include "debug.hpp"      // for AR_DEBUG_ASSERT
#include "fastq_io.hpp"   // for output_chunk_ptr, fastq_read_chunk, fastq_...
#include "statistics.hpp" // for trimming_statistics, fastq_statistics
#include "trimming.hpp"
#include "userconfig.hpp" // for userconfig, output_sample_files, output_sa...

////////////////////////////////////////////////////////////////////////////////
// Helper functions

/** Trims fixed numbers of bases from the 5' and/or 3' termini of reads. **/
void
trim_read_termini(const userconfig& config,
                  trimming_statistics& stats,
                  fastq& read,
                  read_type type)
{
  size_t trim_5p = 0;
  size_t trim_3p = 0;

  switch (type) {
    case read_type::mate_1:
      trim_5p = config.trim_fixed_5p.first;
      trim_3p = config.trim_fixed_3p.first;
      break;

    case read_type::mate_2:
      trim_5p = config.trim_fixed_5p.second;
      trim_3p = config.trim_fixed_3p.second;
      break;

    case read_type::merged:
      if (config.paired_ended_mode) {
        trim_5p = config.trim_fixed_5p.first;
        trim_3p = config.trim_fixed_5p.second;
      } else {
        trim_5p = config.trim_fixed_5p.first;
        trim_3p = config.trim_fixed_3p.first;
      }
      break;

    default:
      throw std::invalid_argument("Invalid read type in trim_read_termini");
  }

  const auto length = read.length();
  if (trim_5p || trim_3p) {
    if (trim_5p + trim_3p < read.length()) {
      read.truncate(trim_5p,
                    read.length() - std::min(read.length(), trim_5p + trim_3p));
    } else {
      read.truncate(0, 0);
    }

    stats.terminal_bases_trimmed += length - read.length();
  }
}

/** Trims a read if enabled, returning the total number of bases removed. */
bool
trim_sequence_by_quality(const userconfig& config,
                         trimming_statistics& stats,
                         fastq& read)
{
  fastq::ntrimmed trimmed;
  if (config.trim_window_length >= 0) {
    trimmed = read.trim_windowed_bases(config.trim_ambiguous_bases,
                                       config.low_quality_score,
                                       config.trim_window_length,
                                       config.preserve5p);
  } else if (config.trim_ambiguous_bases || config.trim_by_quality) {
    const char quality_score =
      config.trim_by_quality ? config.low_quality_score : -1;

    trimmed = read.trim_trailing_bases(
      config.trim_ambiguous_bases, quality_score, config.preserve5p);
  }

  if (trimmed.first || trimmed.second) {
    stats.low_quality_trimmed_reads++;
    stats.low_quality_trimmed_bases += trimmed.first + trimmed.second;

    return true;
  }

  return false;
}

bool
is_acceptable_read(const userconfig& config,
                   trimming_statistics& stats,
                   const fastq& seq)
{
  const auto length = seq.length();

  if (length < config.min_genomic_length) {
    stats.filtered_min_length_reads++;
    stats.filtered_min_length_bases += length;

    return false;
  } else if (length > config.max_genomic_length) {
    stats.filtered_max_length_reads++;
    stats.filtered_max_length_bases += length;

    return false;
  }

  const auto max_n = config.max_ambiguous_bases;
  if (max_n < length && seq.count_ns() > max_n) {
    stats.filtered_ambiguous_reads++;
    stats.filtered_ambiguous_bases += length;
    return false;
  }

  if (config.min_complexity > 0.0 && seq.complexity() < config.min_complexity) {
    stats.filtered_low_complexity_reads++;
    stats.filtered_low_complexity_bases += length;
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `trimmed_reads`

trimmed_reads::trimmed_reads(const output_sample_files& map, const bool eof)
  : m_map(map)
  , m_chunks()
{
  AR_DEBUG_ASSERT(map.filenames.size() == map.steps.size());

  for (size_t i = 0; i < map.steps.size(); ++i) {
    AR_DEBUG_ASSERT(map.filenames.at(i).size());
    AR_DEBUG_ASSERT(map.steps.at(i) != output_sample_files::disabled);

    m_chunks.emplace_back(new fastq_output_chunk(eof));
  }
}

void
trimmed_reads::add(fastq& read, const read_type type)
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
    chunks.emplace_back(m_map.steps.at(i), std::move(m_chunks.at(i)));
  }

  return chunks;
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `reads_processor`

reads_processor::reads_processor(const userconfig& config,
                                 const output_sample_files& output,
                                 size_t nth)
  : analytical_step(processing_order::unordered)
  , m_config(config)
  , m_adapters(config.adapters.get_adapter_set(nth))
  , m_stats()
  , m_output(output)
  , m_nth(nth)
{
  for (size_t i = 0; i < m_config.max_threads; ++i) {
    m_stats.emplace_back(m_config.report_sample_rate);
  }
}

trim_stats_ptr
reads_processor::get_final_statistics()
{
  auto stats = m_stats.acquire();
  while (!m_stats.empty()) {
    *stats += *m_stats.acquire();
  }

  return stats;
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `se_reads_processor`

se_reads_processor::se_reads_processor(const userconfig& config,
                                       const output_sample_files& output,
                                       size_t nth)
  : reads_processor(config, output, nth)
{}

chunk_vec
se_reads_processor::process(chunk_ptr chunk)
{
  auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  trimmed_reads chunks(m_output, read_chunk.eof);

  auto stats = m_stats.acquire();
  stats->adapter_trimmed_reads.resize_up_to(m_config.adapters.adapter_count());
  stats->adapter_trimmed_bases.resize_up_to(m_config.adapters.adapter_count());

  auto aligner = sequence_aligner(m_adapters);
  aligner.set_mismatch_threshold(m_config.mismatch_threshold);

  for (auto& read : read_chunk.reads_1) {
    const alignment_info alignment =
      aligner.align_single_end(read, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      const auto length = read.length();
      alignment.truncate_single_end(read);

      stats->adapter_trimmed_reads.inc(alignment.adapter_id);
      stats->adapter_trimmed_bases.inc(alignment.adapter_id,
                                       length - read.length());
    }

    trim_read_termini(m_config, *stats, read, read_type::mate_1);
    trim_sequence_by_quality(m_config, *stats, read);
    if (is_acceptable_read(m_config, *stats, read)) {
      stats->read_1->process(read);
      chunks.add(read, read_type::mate_1);
    } else {
      stats->discarded->process(read);
      chunks.add(read, read_type::discarded_1);
    }
  }

  m_stats.release(stats);

  return chunks.finalize();
}

////////////////////////////////////////////////////////////////////////////////
// Implementations for `pe_reads_processor`

pe_reads_processor::pe_reads_processor(const userconfig& config,
                                       const output_sample_files& output,
                                       size_t nth)
  : reads_processor(config, output, nth)
{}

chunk_vec
pe_reads_processor::process(chunk_ptr chunk)
{
  sequence_merger merger;
  merger.set_conservative(m_config.merge_conservatively);
  merger.set_max_recalculated_score(m_config.quality_max);

  auto aligner = sequence_aligner(m_adapters);
  aligner.set_mismatch_threshold(m_config.mismatch_threshold);

  auto& read_chunk = dynamic_cast<fastq_read_chunk&>(*chunk);
  trimmed_reads chunks(m_output, read_chunk.eof);

  auto stats = m_stats.acquire();
  stats->adapter_trimmed_reads.resize_up_to(m_config.adapters.adapter_count());
  stats->adapter_trimmed_bases.resize_up_to(m_config.adapters.adapter_count());

  AR_DEBUG_ASSERT(read_chunk.reads_1.size() == read_chunk.reads_2.size());

  auto it_1 = read_chunk.reads_1.begin();
  auto it_2 = read_chunk.reads_2.begin();
  while (it_1 != read_chunk.reads_1.end()) {
    fastq read_1 = *it_1++;
    fastq read_2 = *it_2++;

    // Reverse complement to match the orientation of read_1
    read_2.reverse_complement();

    const alignment_info alignment =
      aligner.align_paired_end(read_1, read_2, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      const size_t length = read_1.length() + read_2.length();
      const size_t n_adapters = alignment.truncate_paired_end(read_1, read_2);
      stats->adapter_trimmed_reads.inc(alignment.adapter_id, n_adapters);
      stats->adapter_trimmed_bases.inc(
        alignment.adapter_id, length - read_1.length() - read_2.length());

      if (m_config.can_merge_alignment(alignment)) {
        stats->overlapping_reads_merged += 2;
        // Merge read_2 into read_1
        merger.merge(alignment, read_1, read_2);

        trim_read_termini(m_config, *stats, read_1, read_type::merged);

        if (!m_config.preserve5p) {
          // A merged read essentially consists of two 5p termini, both
          // informative for PCR duplicate removal.
          trim_sequence_by_quality(m_config, *stats, read_1);
        }

        if (is_acceptable_read(m_config, *stats, read_1)) {
          stats->merged->process(read_1, 2);
          chunks.add(read_1, read_type::merged);
        } else {
          stats->discarded->process(read_1, 2);
          chunks.add(read_1, read_type::discarded_1);
        }

        continue;
      }
    }

    // Reads were not aligned or merging is not enabled
    // Undo reverse complementation (post truncation of adapters)
    read_2.reverse_complement();

    // Trim fixed number of bases from 5' and/or 3' termini
    trim_read_termini(m_config, *stats, read_1, read_type::mate_1);
    trim_read_termini(m_config, *stats, read_2, read_type::mate_2);

    // Sliding window trimming or single-base trimming
    trim_sequence_by_quality(m_config, *stats, read_1);
    trim_sequence_by_quality(m_config, *stats, read_2);

    // Are the reads good enough? Not too many Ns?
    const bool is_ok_1 = is_acceptable_read(m_config, *stats, read_1);
    const bool is_ok_2 = is_acceptable_read(m_config, *stats, read_2);

    read_type type_1;
    read_type type_2;
    if (is_ok_1 && is_ok_2) {
      type_1 = read_type::mate_1;
      type_2 = read_type::mate_2;
    } else if (is_ok_1) {
      type_1 = read_type::singleton_1;
      type_2 = read_type::discarded_2;
    } else if (is_ok_2) {
      type_1 = read_type::discarded_1;
      type_2 = read_type::singleton_2;
    } else {
      type_1 = read_type::discarded_1;
      type_2 = read_type::discarded_2;
    }

    if (is_ok_1) {
      stats->read_1->process(read_1);
    } else {
      stats->discarded->process(read_1);
    }

    if (is_ok_2) {
      stats->read_2->process(read_2);
    } else {
      stats->discarded->process(read_2);
    }

    // Queue reads last, since this result in modifications to lengths
    chunks.add(read_1, type_1);
    chunks.add(read_2, type_2);
  }

  m_stats.release(stats);

  return chunks.finalize();
}
