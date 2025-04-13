// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "errors.hpp"        // for parsing_error
#include "sequence_sets.hpp" // for read_group
#include "testing.hpp"       // for TEST_CASE, REQUIRE, ...
#include <initializer_list>  // for initializer_list
#include <string_view>       // for string_view

using Contains = Catch::Matchers::StdString::ContainsMatcher;

namespace adapterremoval {

namespace {

/** Convenience function for loading in PE mode */
sample_set
sample_set_pe(std::initializer_list<std::string_view> lines)
{
  return { lines, barcode_config().paired_end_mode(true) };
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Sample sets -- initializer list

TEST_CASE("Overlapping SE barcodes fail", "[sample_set]")
{
  const std::string message =
    "Duplicate mate 1 barcodes found in \'initializer_list\': \'ACGT\'. Even "
    "if these are associated with different mate 2 barcodes, it is not "
    "possible to distinguish between these in single-end mode!";

  const std::initializer_list<std::string_view> barcodes_1{
    "sample1 ACGT",
    "sample2 ACGT",
  };

  const std::initializer_list<std::string_view> barcodes_2{
    "sample1 ACGT CGTG",
    "sample2 ACGT TATA",
  };

  // SE mode
  CHECK_THROWS_MESSAGE(sample_set(barcodes_1), parsing_error, message);
  CHECK_THROWS_MESSAGE(sample_set(barcodes_2), parsing_error, message);

  // PE mode
  CHECK_THROWS_MESSAGE(
    sample_set_pe(barcodes_1),
    parsing_error,
    "Duplicate barcode pairs found in \'initializer_list\' with barcodes "
    "\'ACGT\' and \'\'. please verify correctness of the barcode table and "
    "remove any duplicate entries!");
  CHECK_NOTHROW(sample_set_pe(barcodes_2));
}

TEST_CASE("Overlapping PE barcodes fail", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(
    sample_set_pe({
      "sample1 ACGT CGTG",
      "sample2 ACGT CGTG",
    }),
    parsing_error,
    "Duplicate barcode pairs found in \'initializer_list\' with barcodes "
    "\'ACGT\' and \'CGTG\'. please verify correctness of the barcode table and "
    "remove any duplicate entries!");
}

TEST_CASE("Variable length SE barcodes fail #1", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(
    sample_set_pe({
      "sample1 ACGTA",
      "sample2 TGCT",
    }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 5 base-pairs "
    "long, but barcode \'TGCT\' is 4 base-pairs long. Variable length barcodes "
    "are not supported");
}

TEST_CASE("Variable length SE barcodes fail #2", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(
    sample_set_pe({
      "sample1 ACGT",
      "sample2 TGCTA",
    }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 4 base-pairs "
    "long, but barcode \'TGCTA\' is 5 base-pairs long. Variable length "
    "barcodes are not supported");
}

TEST_CASE("Variable length PE barcodes fail #1", "[sample_set]")
{

  REQUIRE_THROWS_MESSAGE(
    sample_set_pe({
      "sample1 ACGTT CGTG",
      "sample2 CGTG ACGT",
    }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 5 base-pairs "
    "long, but barcode \'CGTG\' is 4 base-pairs long. Variable length barcodes "
    "are not supported");
}

TEST_CASE("Variable length PE barcodes fail #2", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(
    sample_set_pe({
      "sample1 ACGT CGTGT",
      "sample2 CGTG ACGT",
    }),
    parsing_error,
    "Inconsistent mate 2 barcode lengths found: Last barcode was 5 base-pairs "
    "long, but barcode \'ACGT\' is 4 base-pairs long. Variable length barcodes "
    "are not supported");
}

TEST_CASE("Variable length PE barcodes fail #3", "[sample_set]")
{
  REQUIRE_THROWS_MESSAGE(
    sample_set_pe({
      "sample1 ACGT CGTG",
      "sample2 CGTGA ACGT",
    }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 4 base-pairs "
    "long, but barcode \'CGTGA\' is 5 base-pairs long. Variable length "
    "barcodes are not supported");
}

TEST_CASE("Variable length PE barcodes fail #4", "[sample_set]")
{

  REQUIRE_THROWS_MESSAGE(
    sample_set_pe({
      "sample1 ACGT CGTG",
      "sample2 CGTGA ACGTA",
    }),
    parsing_error,
    "Inconsistent mate 1 barcode lengths found: Last barcode was 4 base-pairs "
    "long, but barcode \'CGTGA\' is 5 base-pairs long. Variable length "
    "barcodes are not supported");
}

} // namespace adapterremoval
