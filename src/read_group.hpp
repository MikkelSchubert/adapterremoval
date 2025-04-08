// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "commontypes.hpp" // for barcode_orientation
#include <iosfwd>          // for ostream
#include <string>          // for string
#include <string_view>     // for string_view

namespace adapterremoval {

/** Contains SAM/BAM read-group information */
class read_group
{
public:
  /** Creates minimal read-group with an ID and nothing else */
  read_group();

  /**
   * Parses a read-group string in the form "ID:1\tSM:sample" (optionally
   * including a leading "@RG\t"). Throws std::invalid_argument if the value
   * is invalid.
   */
  explicit read_group(std::string_view value);

  /** Returns the read-group ID for use in per-read 'RG' tags */
  [[nodiscard]] std::string_view id() const { return m_id; }

  /** Returns the full @RG header, not including a trailing new-line */
  [[nodiscard]] std::string_view header() const { return m_header; }

  /** Adds/replaces the barcode (ID) tag */
  void set_id(std::string_view id);

  /** Adds/replaces the sample (SM) tag */
  void set_sample(std::string_view name) { update_tag("SM", name); }

  /** Adds/replaces the barcode (BC) tag */
  void set_barcodes(std::string_view value) { update_tag("BC", value); }

  /** Adds/replaces the description (DS) tag */
  void set_description(std::string_view value) { update_tag("DS", value); }

  /** Adds/replaces the custom orientation ('or') tag */
  void set_orientation(barcode_orientation value);

  /** Returns true if the header and ID are identical */
  [[nodiscard]] bool operator==(const read_group& other) const;

  friend std::ostream& operator<<(std::ostream& os, const read_group& value);

private:
  /** Updates or adds the specified tag; sets `m_id` if key is `ID` */
  void update_tag(std::string_view key, std::string_view value);

  //! The full read_group header, including leading `@RG\t`
  std::string m_header{};
  //! Value mapping reads (via `RG:Z:${ID}`) to the @RG header
  std::string m_id{};
};

} // namespace adapterremoval
