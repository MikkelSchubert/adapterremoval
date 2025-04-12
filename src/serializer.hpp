// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp"   // for string_vec, read_type, read_file
#include "sequence_sets.hpp" // for sample
#include <cstddef>           // for size_t
#include <iosfwd>            // for ostream

namespace adapterremoval {

class buffer;
class fastq;
class read_group;
enum class output_format;

/** Class for recording meta-information about reads for serialization */
class read_meta
{
public:
  /** Creates read meta-data for the specified read type */
  constexpr explicit read_meta(read_type type)
    : m_type(type)
  {
  }

  /** Returns the file type associated with the read type */
  [[nodiscard]] read_file get_file() const noexcept;

  /** Overwrites the current read type */
  constexpr auto& type(read_type v) noexcept
  {
    m_type = v;
    return *this;
  }

  /** Overwrites the current barcode */
  constexpr auto& barcode(size_t v) noexcept
  {
    m_barcode = v;
    return *this;
  }

  /** Creates debug representation of read meta data */
  friend std::ostream& operator<<(std::ostream& os, const read_meta& value);

private:
  friend class serializer;

  read_type m_type{};
  //! The barcode for this sample; 0 is also used if not demultiplexing
  size_t m_barcode{};
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
  void record(buffer& buf, const fastq& record, read_meta meta) const;

private:
  /** Write record to buffer in in FASTQ format */
  void fastq_record(buffer& buf,
                    const fastq& record,
                    const read_meta& meta,
                    const sample_sequences& sequences) const;

  /** Write header to buffer in in SAM format */
  static void sam_header(buffer& buf, const string_vec& args, const sample& s);
  void sam_record(buffer& buf,
                  const fastq& record,
                  const read_meta& meta,
                  const sample_sequences& sequences) const;

  /** Write header to buffer in in BAM format (uncompressed) */
  static void bam_header(buffer& buf, const string_vec& args, const sample& s);
  void bam_record(buffer& buf,
                  const fastq& record,
                  const read_meta& meta,
                  const sample_sequences& sequences) const;

  //! User specified output format
  output_format m_format = output_format::fastq;
  //! Sample information to be added to SAM/BAM records
  sample m_sample{};
  //! Mate separator in processed reads
  char m_mate_separator = '\0';
  //! Indicates if output has only been demultiplexed but not trimmed
  bool m_demultiplexing_only = false;
};

/** Creates debug representation of read meta data */
std::ostream&
operator<<(std::ostream& os, const read_file& value);

} // namespace adapterremoval
