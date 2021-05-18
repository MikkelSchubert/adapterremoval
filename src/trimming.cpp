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
#include "trimmed_reads.hpp"
#include "trimming.hpp"
#include "userconfig.hpp"

namespace ar {

/** Trims fixed numbers of bases from the 5' and/or 3' termini of reads. **/
void
trim_read_termini_if_enabled(const userconfig& config,
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
      throw std::invalid_argument(
        "Invalid read type in trim_read_termini_if_enabled");
  }

  if (trim_5p || trim_3p) {
    if (trim_5p + trim_3p < read.length()) {
      read.truncate(trim_5p,
                    read.length() - std::min(read.length(), trim_5p + trim_3p));
    } else {
      read.truncate(0, 0);
    }
  }
}

/** Trims a read if enabled, returning the #bases removed from each end. */
fastq::ntrimmed
trim_sequence_by_quality_if_enabled(const userconfig& config, fastq& read)
{
  if (config.trim_window_length >= 0) {
    return read.trim_windowed_bases(config.trim_ambiguous_bases,
                                    config.low_quality_score,
                                    config.trim_window_length,
                                    config.preserve5p);
  } else if (config.trim_ambiguous_bases || config.trim_by_quality) {
    const char quality_score =
      config.trim_by_quality ? config.low_quality_score : -1;

    return read.trim_trailing_bases(
      config.trim_ambiguous_bases, quality_score, config.preserve5p);
  }

  return fastq::ntrimmed();
}

void
process_collapsed_read(const userconfig& config,
                       statistics& stats,
                       fastq& collapsed_read,
                       fastq* mate_read,
                       trimmed_reads& chunks)
{
  trim_read_termini_if_enabled(config, collapsed_read, read_type::collapsed);

  fastq::ntrimmed trimmed;
  if (!config.preserve5p) {
    // A collapsed read essentially consists of two 5p termini, both
    // informative for PCR duplicate removal.
    trimmed = trim_sequence_by_quality_if_enabled(config, collapsed_read);
  }

  // If trimmed, the external coordinates are no longer reliable
  // for determining the size of the original template.
  const bool was_trimmed = trimmed.first || trimmed.second;
  collapsed_read.add_prefix_to_header(was_trimmed ? "MT_" : "M_");
  if (mate_read) {
    mate_read->add_prefix_to_header(was_trimmed ? "MT_" : "M_");
  }

  if (config.is_acceptable_read(collapsed_read)) {
    stats.total_number_of_nucleotides += collapsed_read.length();
    stats.total_number_of_good_reads++;
    stats.inc_length_count(read_type::collapsed, collapsed_read.length());

    chunks.add_collapsed_read(collapsed_read, read_status::passed, 2);
    stats.number_of_collapsed++;
  } else {
    stats.discard1++;
    stats.discard2++;
    stats.inc_length_count(read_type::discarded, collapsed_read.length());

    chunks.add_collapsed_read(collapsed_read, read_status::failed, 2);
  }
}

reads_processor::reads_processor(const userconfig& config, size_t nth)
  : analytical_step(analytical_step::ordering::unordered)
  , m_config(config)
  , m_adapters(config.adapters.get_adapter_set(nth))
  , m_stats(config)
  , m_nth(nth)
{}

statistics_ptr
reads_processor::get_final_statistics()
{
  return m_stats.finalize();
}

reads_processor::stats_sink::stats_sink(const userconfig& config)
  : m_config(config)
{}

reads_processor::stats_sink::pointer
reads_processor::stats_sink::new_sink() const
{
  return m_config.create_stats();
}

void
reads_processor::stats_sink::reduce(pointer& dst, const pointer& src) const
{
  (*dst) += (*src);
}

se_reads_processor::se_reads_processor(const userconfig& config, size_t nth)
  : reads_processor(config, nth)
{}

chunk_vec
se_reads_processor::process(analytical_chunk* chunk)
{
  const size_t offset = m_nth * ai_analyses_offset;

  read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
  trimmed_reads chunks(m_config, offset, read_chunk->eof);
  stats_sink::pointer stats = m_stats.get_sink();

  for (auto& read : read_chunk->reads_1) {
    const alignment_info alignment =
      align_single_ended_sequence(read, m_adapters, m_config.shift);

    if (m_config.is_good_alignment(alignment)) {
      truncate_single_ended_sequence(alignment, read);
      stats->number_of_reads_with_adapter.at(alignment.adapter_id)++;
      stats->well_aligned_reads++;
    } else {
      stats->unaligned_reads++;
    }

    trim_read_termini_if_enabled(m_config, read, read_type::mate_1);
    trim_sequence_by_quality_if_enabled(m_config, read);
    if (m_config.is_acceptable_read(read)) {
      stats->keep1++;
      stats->total_number_of_good_reads++;
      stats->total_number_of_nucleotides += read.length();

      chunks.add_mate_1_read(read, read_status::passed);
      stats->inc_length_count(read_type::mate_1, read.length());
    } else {
      stats->discard1++;
      stats->inc_length_count(read_type::discarded, read.length());

      chunks.add_mate_1_read(read, read_status::failed);
    }
  }

  stats->records += read_chunk->reads_1.size();
  m_stats.return_sink(std::move(stats));

  return chunks.finalize();
}

rng_sink::rng_sink(unsigned seed)
  : m_seed(seed)
{}

rng_sink::pointer
rng_sink::new_sink() const
{
  return pointer(new std::mt19937(m_seed()));
}

void
rng_sink::reduce(pointer&, const pointer&) const
{
  // Intentionally left empty
}

pe_reads_processor::pe_reads_processor(const userconfig& config, size_t nth)
  : reads_processor(config, nth)
  , m_rngs(config.seed)
{}

chunk_vec
pe_reads_processor::process(analytical_chunk* chunk)
{
  const size_t offset = m_nth * ai_analyses_offset;
  const char mate_separator =
    m_config.combined_output ? '\0' : m_config.mate_separator;

  sequence_merger merger;
  merger.set_mate_separator(mate_separator);
  merger.set_conservative(m_config.collapse_conservatively);

  std::unique_ptr<std::mt19937> rng;
  if (!m_config.deterministic && !m_config.collapse_conservatively) {
    rng = m_rngs.get_sink();
    merger.set_rng(rng.get());
  }

  read_chunk_ptr read_chunk(dynamic_cast<fastq_read_chunk*>(chunk));
  trimmed_reads chunks(m_config, offset, read_chunk->eof);
  statistics_ptr stats = m_stats.get_sink();

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
      stats->well_aligned_reads++;
      const size_t n_adapters =
        truncate_paired_ended_sequences(alignment, read_1, read_2);
      stats->number_of_reads_with_adapter.at(alignment.adapter_id) +=
        n_adapters;

      if (m_config.is_alignment_collapsible(alignment)) {
        fastq collapsed_read = merger.merge(alignment, read_1, read_2);

        process_collapsed_read(m_config,
                               *stats,
                               collapsed_read,
                               // Make sure read_2 header is updated, if needed
                               m_config.combined_output ? &read_2 : nullptr,
                               chunks);

        if (m_config.combined_output) {
          // Dummy read with read-count of zero; both mates have
          // already been accounted for in process_collapsed_read
          chunks.add_mate_2_read(read_2, read_status::failed, 0);
        }
        continue;
      }
    } else {
      stats->unaligned_reads++;
    }

    // Reads were not aligned or collapsing is not enabled
    // Undo reverse complementation (post truncation of adapters)
    read_2.reverse_complement();

    // Trim fixed number of bases from 5' and/or 3' termini
    trim_read_termini_if_enabled(m_config, read_1, read_type::mate_1);
    trim_read_termini_if_enabled(m_config, read_2, read_type::mate_2);
    // Sliding window trimming or single-base trimming
    trim_sequence_by_quality_if_enabled(m_config, read_1);
    trim_sequence_by_quality_if_enabled(m_config, read_2);

    // Are the reads good enough? Not too many Ns?
    const bool read_1_acceptable = m_config.is_acceptable_read(read_1);
    const bool read_2_acceptable = m_config.is_acceptable_read(read_2);

    stats->total_number_of_nucleotides +=
      read_1_acceptable ? read_1.length() : 0u;
    stats->total_number_of_nucleotides +=
      read_2_acceptable ? read_2.length() : 0u;
    stats->total_number_of_good_reads += read_1_acceptable;
    stats->total_number_of_good_reads += read_2_acceptable;

    const read_status state_1 =
      read_1_acceptable ? read_status::passed : read_status::failed;
    const read_status state_2 =
      read_2_acceptable ? read_status::passed : read_status::failed;

    if (read_1_acceptable && read_2_acceptable) {
      stats->inc_length_count(read_type::mate_1, read_1.length());
      stats->inc_length_count(read_type::mate_2, read_2.length());
    } else {
      // Count singleton reads
      stats->keep1 += read_1_acceptable && !read_2_acceptable;
      stats->keep2 += read_2_acceptable && !read_1_acceptable;

      stats->discard1 += !read_1_acceptable;
      stats->discard2 += !read_2_acceptable;

      stats->inc_length_count(read_1_acceptable ? read_type::singleton
                                                : read_type::discarded,
                              read_1.length());
      stats->inc_length_count(read_2_acceptable ? read_type::singleton
                                                : read_type::discarded,
                              read_2.length());
    }

    // Queue reads last, since this result in modifications to lengths
    chunks.add_pe_reads(read_1, state_1, read_2, state_2);
  }

  stats->records += read_chunk->reads_1.size();
  m_stats.return_sink(std::move(stats));
  m_rngs.return_sink(std::move(rng));

  return chunks.finalize();
}

} // namespace ar