// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp"   // for string_vec
#include "sequence_sets.hpp" // for sample
#include <cstddef>           // for size_t
#include <functional>        // for function

namespace adapterremoval {

class buffer;
class fastq;
class read_group;
enum class output_format;

enum class fastq_flags
{
  //! SE read
  se,
  //! SE read that failed QC
  se_fail,
  //! Mate 1 read
  pe_1,
  //! Mate 1 read that failed QC
  pe_1_fail,
  //! Mate 1 read
  pe_2,
  //! Mate 1 read that failed QC
  pe_2_fail,

};

/** Struct containing metadata required to serialize a record */
struct serializer_settings
{
  const fastq_flags flags = fastq_flags::se;
  const char mate_separator = '\0';
  const bool demultiplexing_only = false;
};

class fastq_serializer
{
public:
  static void header(buffer& buf, const string_vec& args, const sample& s);
  static void record(buffer& buf,
                     const fastq& record,
                     const sample_sequences& sequences,
                     const serializer_settings& settings);
};

class sam_serializer
{
public:
  static void header(buffer& buf, const string_vec& args, const sample& s);
  static void record(buffer& buf,
                     const fastq& record,
                     const sample_sequences& sequences,
                     const serializer_settings& settings);
};

class bam_serializer
{
public:
  static void header(buffer& buf, const string_vec& args, const sample& s);
  static void record(buffer& buf,
                     const fastq& record,
                     const sample_sequences& sequences,
                     const serializer_settings& settings);
};

class serializer
{
public:
  explicit serializer(output_format format);

  /** Set the sample; used to write headers/records for SAM/BAM */
  void set_sample(const sample& s) { m_sample = s; }

  /** Set the mate separator; used to trim mate information for SAM/BAM */
  void set_mate_separator(char value) { m_mate_separator = value; }

  /** In demultiplexing only mode, barcodes must always be recorded */
  void set_demultiplexing_only(bool value) { m_demultiplexing_only = value; }

  void header(buffer& buf, const string_vec& args) const;
  void record(buffer& buf,
              const fastq& record,
              fastq_flags flags,
              size_t barcode) const;

private:
  //! Sample information to be added to SAM/BAM records
  sample m_sample{};
  //! Mate separator in processed reads
  char m_mate_separator = '\0';
  //! Indicates if output has only been demultiplexed but not trimmed
  bool m_demultiplexing_only = false;
  //! Function used to serialize file headers for processed reads
  std::function<decltype(fastq_serializer::header)> m_header{};
  //! Function used to serialize processed reads
  std::function<decltype(fastq_serializer::record)> m_record{};
};

} // namespace adapterremoval
