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
#include <random>

#include "debug.hpp"
#include "statistics.hpp"
#include "trimmed_reads.hpp"
#include "trimming.hpp"
#include "userconfig.hpp"

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

    case read_type::collapsed:
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

  return true;
}

reads_processor::reads_processor(const userconfig& config, size_t nth)
  : analytical_step(analytical_step::ordering::unordered)
  , m_config(config)
  , m_adapters(config.adapters.get_adapter_set(nth))
  , m_stats()
  , m_nth(nth)
{
  for (size_t i = 0; i < m_config.max_threads; ++i) {
    m_stats.emplace_back(m_config.report_sample_rate);
  }
}

statistics_ptr
reads_processor::get_final_statistics()
{
  auto stats = m_stats.acquire();
  while (!m_stats.empty()) {
    *stats += *m_stats.acquire();
  }

  return stats;
}

se_reads_processor::se_reads_processor(const userconfig& config, size_t nth)
  : reads_processor(config, nth)
{}

chunk_vec
se_reads_processor::process(analytical_chunk* chunk)
{
  const size_t offset = (m_nth + 1) * ai_analyses_offset;

  read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
  trimmed_reads chunks(m_config, offset, read_chunk->eof);
  auto stats = m_stats.acquire();

  for (auto& read : read_chunk->reads_1) {
    const alignment_info alignment =
      align_single_ended_sequence(read, m_adapters, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      const auto length = read.length();
      truncate_single_ended_sequence(alignment, read);

      stats->adapter_trimmed_reads.inc(alignment.adapter_id);
      stats->adapter_trimmed_bases.inc(alignment.adapter_id,
                                       length - read.length());
    }

    trim_read_termini(m_config, *stats, read, read_type::mate_1);
    trim_sequence_by_quality(m_config, *stats, read);
    if (is_acceptable_read(m_config, *stats, read)) {
      stats->read_1.process(read);
      chunks.add_mate_1_read(read, read_status::passed);
    } else {
      stats->discarded.process(read);
      chunks.add_mate_1_read(read, read_status::failed);
    }
  }

  m_stats.release(stats);

  return chunks.finalize();
}

pe_reads_processor::pe_reads_processor(const userconfig& config, size_t nth)
  : reads_processor(config, nth)
  , m_rngs()
{
  std::mt19937 seed(config.seed);
  for (size_t i = 0; i < m_config.max_threads; ++i) {
    m_rngs.emplace_back(seed());
  }
}

chunk_vec
pe_reads_processor::process(analytical_chunk* chunk)
{
  const size_t offset = (m_nth + 1) * ai_analyses_offset;
  const char mate_separator =
    m_config.combined_output ? '\0' : m_config.mate_separator;

  sequence_merger merger;
  merger.set_mate_separator(mate_separator);
  merger.set_conservative(m_config.collapse_conservatively);

  std::unique_ptr<std::mt19937> rng;
  if (!m_config.deterministic && !m_config.collapse_conservatively) {
    rng = m_rngs.acquire();
    merger.set_rng(rng.get());
  }

  read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
  trimmed_reads chunks(m_config, offset, read_chunk->eof);
  statistics_ptr stats = m_stats.acquire();

  AR_DEBUG_ASSERT(read_chunk->reads_1.size() == read_chunk->reads_2.size());

  auto it_1 = read_chunk->reads_1.begin();
  auto it_2 = read_chunk->reads_2.begin();
  while (it_1 != read_chunk->reads_1.end()) {
    fastq read_1 = *it_1++;
    fastq read_2 = *it_2++;

    // Throws if read-names or mate numbering does not match
    fastq::validate_paired_reads(read_1, read_2, m_config.mate_separator);

    // Reverse complement to match the orientation of read_1
    read_2.reverse_complement();

    const alignment_info alignment =
      align_paired_ended_sequences(read_1, read_2, m_adapters, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      const size_t length = read_1.length() + read_2.length();
      const size_t n_adapters =
        truncate_paired_ended_sequences(alignment, read_1, read_2);
      stats->adapter_trimmed_reads.inc(alignment.adapter_id, n_adapters);
      stats->adapter_trimmed_bases.inc(
        alignment.adapter_id, length - read_1.length() - read_2.length());

      if (m_config.is_alignment_collapsible(alignment)) {
        stats->overlapping_reads_merged += 2;
        fastq collapsed_read = merger.merge(alignment, read_1, read_2);

        trim_read_termini(
          m_config, *stats, collapsed_read, read_type::collapsed);

        bool trimmed = false;
        if (!m_config.preserve5p) {
          // A collapsed read essentially consists of two 5p termini, both
          // informative for PCR duplicate removal.
          trimmed = trim_sequence_by_quality(m_config, *stats, collapsed_read);
        }

        // If trimmed, the external coordinates are no longer reliable
        // for determining the size of the original template.
        collapsed_read.add_prefix_to_header(trimmed ? "MT_" : "M_");

        if (is_acceptable_read(m_config, *stats, collapsed_read)) {
          stats->merged.process(collapsed_read, 2);
          chunks.add_collapsed_read(collapsed_read, read_status::passed, 2);
        } else {
          stats->discarded.process(collapsed_read, 2);
          chunks.add_collapsed_read(collapsed_read, read_status::failed, 2);
        }

        if (m_config.combined_output) {
          // FIXME: Does this make sense?
          // Dummy read with read-count of zero; both mates have
          // already been accounted for in process_collapsed_read
          read_2.add_prefix_to_header(trimmed ? "MT_" : "M_");
          chunks.add_mate_2_read(read_2, read_status::failed, 0);
        }
        continue;
      }
    }

    // Reads were not aligned or collapsing is not enabled
    // Undo reverse complementation (post truncation of adapters)
    read_2.reverse_complement();

    // Trim fixed number of bases from 5' and/or 3' termini
    trim_read_termini(m_config, *stats, read_1, read_type::mate_1);
    trim_read_termini(m_config, *stats, read_2, read_type::mate_2);

    // Sliding window trimming or single-base trimming
    trim_sequence_by_quality(m_config, *stats, read_1);
    trim_sequence_by_quality(m_config, *stats, read_2);

    // Are the reads good enough? Not too many Ns?
    read_status state_1 = read_status::passed;
    if (is_acceptable_read(m_config, *stats, read_1)) {
      stats->read_1.process(read_1);
    } else {
      stats->discarded.process(read_1);
      state_1 = read_status::failed;
    }

    read_status state_2 = read_status::passed;
    if (is_acceptable_read(m_config, *stats, read_2)) {
      stats->read_2.process(read_2);
    } else {
      stats->discarded.process(read_2);
      state_2 = read_status::failed;
    }

    // Queue reads last, since this result in modifications to lengths
    chunks.add_pe_reads(read_1, state_1, read_2, state_2);
  }

  m_stats.release(stats);
  m_rngs.release(rng);

  return chunks.finalize();
}
