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
#include "barcode_table.hpp" // for barcode_table
#include "errors.hpp"        // for parsing_error
#include "fastq.hpp"         // for fastq, sequence_pair_vec, fastq_pair
#include "sequence.hpp"      // for dna_sequence
#include "testing.hpp"       // for catch.hpp, StringMaker
#include <string>            // for basic_string, operator==, string

// Ignore nucleotide and quality strings
// spell-checker:ignoreRegExp /"[!-~]+"/g
// Ignore nucleotide comments
// spell-checker:ignoreRegExp /\W[acgtnACGTN]+\W/g

namespace adapterremoval {

TEST_CASE("what()", "[barcodes::errors]")
{
  parsing_error err("test error");

  REQUIRE(std::string(err.what()) == "test error");
}

TEST_CASE("copy constructor", "[barcodes::errors]")
{
  parsing_error err("test error");

  REQUIRE(std::string(parsing_error(err).what()) == "test error");
}

TEST_CASE("Empty barcode-table is OK", "[barcodes::constructor]")
{
  const sequence_pair_vec barcodes;

  barcode_table table(barcodes, 0, 0, 0);
}

TEST_CASE("Overlapping SE barcodes fail", "[barcodes::constructor]")
{
  const sequence_pair_vec barcodes = {
    { dna_sequence{ "ACGT" }, dna_sequence{} },
    { dna_sequence{ "ACGT" }, dna_sequence{} }
  };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

TEST_CASE("Overlapping PE barcodes fail", "[barcodes::constructor]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } },
    { dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } },
  };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

TEST_CASE("Partially overlapping PE barcodes are OK", "[barcodes::constructor]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } },
    { dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } },
  };

  barcode_table table(barcodes, 0, 0, 0);
}

TEST_CASE("Variable length SE barcodes fail #1", "[barcodes::construction]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACGT" }, dna_sequence{} },
                                 { dna_sequence{ "TGCTA" }, dna_sequence{} } };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

TEST_CASE("Variable length SE barcodes fail #2", "[barcodes::construction]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACGT" }, dna_sequence{} },
                                 { dna_sequence{ "TGCTA" }, dna_sequence{} } };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #1", "[barcodes::construction]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACGTT" }, dna_sequence{ "CGTG" } },
    { dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } },
  };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #2", "[barcodes::construction]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACGT" }, dna_sequence{ "CGTGT" } },
    { dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } },
  };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #3", "[barcodes::construction]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } },
    { dna_sequence{ "CGTGA" }, dna_sequence{ "ACGT" } },
  };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

TEST_CASE("Variable length PE barcodes fail #4", "[barcodes::construction]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } },
    { dna_sequence{ "CGTGA" }, dna_sequence{ "ACGTA" } },
  };

  REQUIRE_THROWS_AS(barcode_table(barcodes, 0, 0, 0), parsing_error);
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for SE reads, when using a table containing SE barcodes

TEST_CASE("Exact match among different SE barcodes for SE reads",
          "[barcodes::exact]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} },
                                 { dna_sequence{ "CACAC" }, dna_sequence{} },
                                 { dna_sequence{ "AATTC" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq()) == 0);
  REQUIRE(table.identify(fastq("B", "CACAC"), fastq()) == 1);
  REQUIRE(table.identify(fastq("C", "AATTC"), fastq()) == 2);
}

TEST_CASE("Exact match among similar SE barcodes for SE reads - differs at 5p",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} },
                                 { dna_sequence{ "TCCCA" }, dna_sequence{} },
                                 { dna_sequence{ "AGCCA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("B", "TCCCA")) == 1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("C", "AGCCA")) == 2);
}

TEST_CASE("Exact match among similar SE barcodes for SE reads - differs at 3p",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} },
                                 { dna_sequence{ "ACCCT" }, dna_sequence{} },
                                 { dna_sequence{ "ACCTA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("C", "ACCTA")) == 2);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("B", "ACCCT")) == 1);
}

TEST_CASE("Shorter and longer reads for SE barcodes and SE reads",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} } };
  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCAG")) == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for SE reads, when using a table containing PE barcodes

TEST_CASE("Exact match among different PE barcodes for SE reads",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "CACAC" }, dna_sequence{ "GATGC" } },
    { dna_sequence{ "AATTC" }, dna_sequence{ "TGCGG" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("C", "AATTC")) == 2);
  REQUIRE(table.identify(fastq("B", "CACAC")) == 1);
}

TEST_CASE("Exact match among similar PE barcodes for SE reads - differs at 5p",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "TCCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "AGCCA" }, dna_sequence{ "ATTTC" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("B", "TCCCA")) == 1);
  REQUIRE(table.identify(fastq("C", "AGCCA")) == 2);
}

TEST_CASE("Exact match among similar PE barcodes for SE reads - differs at 3p",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "ACCCT" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "AGCCA" }, dna_sequence{ "GTTTA" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("C", "AGCCA")) == 2);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("B", "ACCCT")) == 1);
}

TEST_CASE("Shorter and longer reads for PE barcodes and SE reads",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" },
                                   dna_sequence{ "TGATA" } } };
  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCAG")) == 0);
}

TEST_CASE(
  "Exact match among different PE barcodes for SE reads with ambiguous results",
  "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "AATTC" }, dna_sequence{ "GATGC" } },
    { dna_sequence{ "AATTC" }, dna_sequence{ "TGCGG" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("C", "AATTC")) == -2);
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for PE reads, when using a table containing SE barcodes

TEST_CASE("Exact match among different SE barcodes for PE reads",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} },
                                 { dna_sequence{ "CACAC" }, dna_sequence{} },
                                 { dna_sequence{ "AATTC" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("B", "CACAC"), fastq("E", "GATGC")) == 1);
  REQUIRE(table.identify(fastq("C", "AATTC"), fastq("F", "TGCGG")) == 2);
}

TEST_CASE("Exact match among similar SE barcodes for PE reads - differs at 5p",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} },
                                 { dna_sequence{ "TCCCA" }, dna_sequence{} },
                                 { dna_sequence{ "ATCCA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("B", "TCCCA"), fastq("E", "GTTTC")) == 1);
  REQUIRE(table.identify(fastq("C", "ATCCA"), fastq("F", "ATTTC")) == 2);
}

TEST_CASE("Exact match among similar SE barcodes for PE reads - differs at 3p",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} },
                                 { dna_sequence{ "ACCCT" }, dna_sequence{} },
                                 { dna_sequence{ "ACCGA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("B", "ACCCT"), fastq("E", "GTTTC")) == 1);
  REQUIRE(table.identify(fastq("C", "ACCGA"), fastq("F", "GTTTA")) == 2);
}

TEST_CASE("Shorter and longer reads for SE barcodes and PE reads",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} } };
  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATGA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATGA")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCAG"), fastq("B", "TGATGA")) == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for PE reads, when using a table containing PE barcodes

TEST_CASE("Exact match among different PE barcodes for PE reads",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "CACAC" }, dna_sequence{ "GATGC" } },
    { dna_sequence{ "AATTC" }, dna_sequence{ "TGCGG" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("B", "CACAC"), fastq("E", "GATGC")) == 1);
  REQUIRE(table.identify(fastq("C", "AATTC"), fastq("F", "TGCGG")) == 2);
}

TEST_CASE("Exact match among similar PE barcodes for PE reads - differs at 5p",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "TCCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "ACCCA" }, dna_sequence{ "ATTTC" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("B", "TCCCA"), fastq("E", "GTTTC")) == 1);
  REQUIRE(table.identify(fastq("C", "ACCCA"), fastq("F", "ATTTC")) == 2);
}

TEST_CASE("Exact match among similar PE barcodes for PE reads - differs at 3p",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "ACCCT" }, dna_sequence{ "GTTTC" } },
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTA" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("B", "ACCCT"), fastq("E", "GTTTC")) == 1);
  REQUIRE(table.identify(fastq("C", "ACCCA"), fastq("F", "GTTTA")) == 2);
}

TEST_CASE("Shorter and longer reads for PE barcodes and PE reads",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" },
                                   dna_sequence{ "TGATGA" } } };
  const barcode_table table(barcodes, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATG")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATGA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATGAT")) == -1);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATG")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATGA")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATGAG")) == 0);

  REQUIRE(table.identify(fastq("A", "ACCCAT"), fastq("B", "TGATG")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCAT"), fastq("B", "TGATGA")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCAT"), fastq("B", "TGATGAG")) == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for SE reads, when reads have mismatches

TEST_CASE("Exact matching reads with mismatches for SE barcodes",
          "[barcodes::exact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 0, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == -1);
  REQUIRE(table.identify(fastq("A", "AGCCA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACACA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCTA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCTT")) == -1);
  REQUIRE(table.identify(fastq("A", "NCCCA")) == -1);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "AGCCA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACACA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCTA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCTT"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) == -1);
}

TEST_CASE("Exact matching reads with mismatches for PE barcodes",
          "[barcodes::exact::pe]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(barcodes, 0, 0, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == -1);
  REQUIRE(table.identify(fastq("A", "AGCCA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACACA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCTA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCTT")) == -1);
  REQUIRE(table.identify(fastq("A", "NCCCA")) == -1);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTA")) == -1);
  REQUIRE(table.identify(fastq("A", "GCCCA"), fastq("B", "GATTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTN")) == -1);
}

///////////////////////////////////////////////////////////////////////////////
// Inexact matching with global limits

TEST_CASE("Global limits override local limits for SE barcodes",
          "[barcodes::inexact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 0, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == -1);
  REQUIRE(table.identify(fastq("A", "NCCCA")) == -1);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACACA"), fastq("B", "GTTTC")) == -1);
}

TEST_CASE("Global limits override local limits for PE barcodes",
          "[barcodes::inexact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(barcodes, 0, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCNA")) == -1);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACACA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTAG")) == -1);
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "NTTTC")) == -1);
}

TEST_CASE("Mismatches in R1 only with SE table", "[barcodes::inexact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 1, 1, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCT")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCN")) == 0);
  REQUIRE(table.identify(fastq("A", "ANCCN")) == -1);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "GCTCA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTG")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCNA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "NCCNA"), fastq("B", "GTTTC")) == -1);
}

TEST_CASE("Mismatches in R1 only with PE table", "[barcodes::inexact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(barcodes, 1, 1, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCT")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCN")) == 0);
  REQUIRE(table.identify(fastq("A", "ANCCN")) == -1);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "GCTCA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTG")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCN"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "NTTTC")) == -1);
}

TEST_CASE("Mismatches in R2 only with SE table", "[barcodes::inexact::pe]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(barcodes, 1, 0, 1);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTT")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTG")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ANTTG")) == 0);
}

TEST_CASE("Mismatches in R2 only with PE table", "[barcodes::inexact::pe]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(barcodes, 1, 0, 1);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTT")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTG")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTN")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCN"), fastq("B", "GTTTT")) == -1);
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTNTT")) == -1);
}

TEST_CASE("Mismatches in R1/R2 with PE table", "[barcodes::inexact::pe]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(barcodes, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTAC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTAC")) == -1);
  REQUIRE(table.identify(fastq("A", "TCCAA"), fastq("B", "GTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GGTGC")) == -1);
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTNC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTNTC")) == -1);
}

TEST_CASE("Multiple mismatches in R1/R2 with PE table",
          "[barcodes::inexact::pe]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(barcodes, 2, 2, 2);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTAC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTAC")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCAA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GGTGC")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCCT"), fastq("B", "CTTTC")) == -1);
  REQUIRE(table.identify(fastq("A", "TCCNA"), fastq("B", "GTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "NTTTC")) == 0);
  REQUIRE(table.identify(fastq("A", "NCCCN"), fastq("B", "GTTTA")) == -1);
}

TEST_CASE("Ambiguous matches in R1 for SE table", "[barcodes::inexact::se]")
{
  sequence_pair_vec barcodes = { { dna_sequence{ "ACCCA" }, dna_sequence{} },
                                 { dna_sequence{ "TCCCA" }, dna_sequence{} },
                                 { dna_sequence{ "ACCCG" }, dna_sequence{} } };

  const barcode_table table(barcodes, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == 1);
  REQUIRE(table.identify(fastq("A", "CCCCA")) == -2);
  REQUIRE(table.identify(fastq("A", "ACCCC")) == -2);
  REQUIRE(table.identify(fastq("A", "ACCCN")) == -2);
}

TEST_CASE("Ambiguous matches in SE R1 for PE table", "[barcodes::inexact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "TGCGT" } },
    { dna_sequence{ "TCCCA" }, dna_sequence{ "AAGTT" } },
    { dna_sequence{ "ACCCG" }, dna_sequence{ "CGGCA" } },
  };

  const barcode_table table(barcodes, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == 1);
  REQUIRE(table.identify(fastq("A", "CCCCA")) == -2);
  REQUIRE(table.identify(fastq("A", "ACCCC")) == -2);
  REQUIRE(table.identify(fastq("A", "ACCCN")) == -2);
}

TEST_CASE("Mismatch resulting in apparent match", "[barcodes::inexact::se]")
{
  sequence_pair_vec barcodes = {
    { dna_sequence{ "ACCCA" }, dna_sequence{ "TGCGT" } },
    { dna_sequence{ "TCCTT" }, dna_sequence{ "AAGTT" } },
  };

  const barcode_table table(barcodes, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "TCCTT")) == 1);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == 0);
  REQUIRE(table.identify(fastq("A", "ACCTT")) == 1);
}

} // namespace adapterremoval
