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
#include "sequence_sets.hpp" // for sample_set
#include "testing.hpp"       // for catch.hpp, StringMaker
#include <sstream>           // for ostringstream
#include <string>            // for basic_string, operator==, string

// Ignore nucleotide and quality strings
// spell-checker:ignoreRegExp /"[!-~]+"/g
// Ignore nucleotide comments
// spell-checker:ignoreRegExp /\W[acgtnACGTN]+\W/g

namespace Catch {

template<>
std::string
StringMaker<adapterremoval::barcode_key, void>::convert(
  adapterremoval::barcode_key const& value)
{
  std::ostringstream ss;

  ss << "barcode_key{ sample = " << value.sample
     << ", barcode = " << value.barcode << " }";

  return ss.str();
}

} // namespace Catch

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

TEST_CASE("Partially overlapping PE barcodes are OK", "[barcodes::constructor]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACGT" }, dna_sequence{ "CGTG" } },
    sample{ "sample2", dna_sequence{ "CGTG" }, dna_sequence{ "ACGT" } },
  };

  barcode_table table(samples, 0, 0, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for SE reads, when using a table containing SE barcodes

TEST_CASE("Exact match among different SE barcodes for SE reads",
          "[barcodes::exact]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{} },
    sample{ "sample2", dna_sequence{ "CACAC" }, dna_sequence{} },
    sample{ "sample3", dna_sequence{ "AATTC" }, dna_sequence{} }
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq()) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "CACAC"), fastq()) == barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "AATTC"), fastq()) == barcode_key{ 2, 0 });
}

TEST_CASE("Exact match among similar SE barcodes for SE reads - differs at 5p",
          "[barcodes::exact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{} },
    sample{ "sample2", dna_sequence{ "TCCCA" }, dna_sequence{} },
    sample{ "sample3", dna_sequence{ "AGCCA" }, dna_sequence{} }
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("B", "TCCCA")) == barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("C", "AGCCA")) == barcode_key{ 2, 0 });
}

TEST_CASE("Exact match among similar SE barcodes for SE reads - differs at 3p",
          "[barcodes::exact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{} },
    sample{ "sample2", dna_sequence{ "ACCCT" }, dna_sequence{} },
    sample{ "sample3", dna_sequence{ "ACCTA" }, dna_sequence{} }
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("C", "ACCTA")) == barcode_key{ 2, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "ACCCT")) == barcode_key{ 1, 0 });
}

TEST_CASE("Shorter and longer reads for SE barcodes and SE reads",
          "[barcodes::exact::se]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{} } };
  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC")) == barcode_key{});
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCAG")) == barcode_key{ 0, 0 });
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for SE reads, when using a table containing PE barcodes

TEST_CASE("Exact match among different PE barcodes for SE reads",
          "[barcodes::exact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample2", dna_sequence{ "CACAC" }, dna_sequence{ "GATGC" } },
    sample{ "sample3", dna_sequence{ "AATTC" }, dna_sequence{ "TGCGG" } },
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("C", "AATTC")) == barcode_key{ 2, 0 });
  REQUIRE(table.identify(fastq("B", "CACAC")) == barcode_key{ 1, 0 });
}

TEST_CASE("Exact match among similar PE barcodes for SE reads - differs at 5p",
          "[barcodes::exact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample2", dna_sequence{ "TCCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample3", dna_sequence{ "AGCCA" }, dna_sequence{ "ATTTC" } },
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "TCCCA")) == barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "AGCCA")) == barcode_key{ 2, 0 });
}

TEST_CASE("Exact match among similar PE barcodes for SE reads - differs at 3p",
          "[barcodes::exact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample2", dna_sequence{ "ACCCT" }, dna_sequence{ "GTTTC" } },
    sample{ "sample3", dna_sequence{ "AGCCA" }, dna_sequence{ "GTTTA" } },
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("C", "AGCCA")) == barcode_key{ 2, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "ACCCT")) == barcode_key{ 1, 0 });
}

TEST_CASE("Shorter and longer reads for PE barcodes and SE reads",
          "[barcodes::exact::se]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "TGATA" } } };
  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCAG")) == barcode_key{ 0, 0 });
}

TEST_CASE(
  "Exact match among different PE barcodes for SE reads with ambiguous results",
  "[barcodes::exact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample2", dna_sequence{ "AATTC" }, dna_sequence{ "GATGC" } },
    sample{ "sample3", dna_sequence{ "AATTC" }, dna_sequence{ "TGCGG" } },
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("C", "AATTC")) == barcode_key{ -2, -2 });
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for PE reads, when using a table containing SE barcodes

TEST_CASE("Exact match among different SE barcodes for PE reads",
          "[barcodes::exact::pe]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{} },
    sample{ "sample2", dna_sequence{ "CACAC" }, dna_sequence{} },
    sample{ "sample3", dna_sequence{ "AATTC" }, dna_sequence{} }
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "CACAC"), fastq("E", "GATGC")) ==
          barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "AATTC"), fastq("F", "TGCGG")) ==
          barcode_key{ 2, 0 });
}

TEST_CASE("Exact match among similar SE barcodes for PE reads - differs at 5p",
          "[barcodes::exact::pe]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{} },
    sample{ "sample2", dna_sequence{ "TCCCA" }, dna_sequence{} },
    sample{ "sample3", dna_sequence{ "ATCCA" }, dna_sequence{} }
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "TCCCA"), fastq("E", "GTTTC")) ==
          barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "ATCCA"), fastq("F", "ATTTC")) ==
          barcode_key{ 2, 0 });
}

TEST_CASE("Exact match among similar SE barcodes for PE reads - differs at 3p",
          "[barcodes::exact::pe]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{} },
    sample{ "sample2", dna_sequence{ "ACCCT" }, dna_sequence{} },
    sample{ "sample3", dna_sequence{ "ACCGA" }, dna_sequence{} }
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "ACCCT"), fastq("E", "GTTTC")) ==
          barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "ACCGA"), fastq("F", "GTTTA")) ==
          barcode_key{ 2, 0 });
}

TEST_CASE("Shorter and longer reads for SE barcodes and PE reads",
          "[barcodes::exact::se]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{} } };
  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATGA")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATGA")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCAG"), fastq("B", "TGATGA")) ==
          barcode_key{ 0, 0 });
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for PE reads, when using a table containing PE barcodes

TEST_CASE("Exact match among different PE barcodes for PE reads",
          "[barcodes::exact::pe]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample2", dna_sequence{ "CACAC" }, dna_sequence{ "GATGC" } },
    sample{ "sample3", dna_sequence{ "AATTC" }, dna_sequence{ "TGCGG" } },
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "CACAC"), fastq("E", "GATGC")) ==
          barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "AATTC"), fastq("F", "TGCGG")) ==
          barcode_key{ 2, 0 });
}

TEST_CASE("Exact match among similar PE barcodes for PE reads - differs at 5p",
          "[barcodes::exact::pe]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample2", dna_sequence{ "TCCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample3", dna_sequence{ "ACCCA" }, dna_sequence{ "ATTTC" } },
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "TCCCA"), fastq("E", "GTTTC")) ==
          barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "ACCCA"), fastq("F", "ATTTC")) ==
          barcode_key{ 2, 0 });
}

TEST_CASE("Exact match among similar PE barcodes for PE reads - differs at 3p",
          "[barcodes::exact::pe]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
    sample{ "sample2", dna_sequence{ "ACCCT" }, dna_sequence{ "GTTTC" } },
    sample{ "sample3", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTA" } },
  };

  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("D", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("B", "ACCCT"), fastq("E", "GTTTC")) ==
          barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("C", "ACCCA"), fastq("F", "GTTTA")) ==
          barcode_key{ 2, 0 });
}

TEST_CASE("Shorter and longer reads for PE barcodes and PE reads",
          "[barcodes::exact::pe]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "TGATGA" } } };
  const barcode_table table(samples, 0, 0, 0);

  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATG")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATGA")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCC"), fastq("B", "TGATGAT")) ==
          barcode_key{ -1, -1 });

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATG")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATGA")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "TGATGAG")) ==
          barcode_key{ 0, 0 });

  REQUIRE(table.identify(fastq("A", "ACCCAT"), fastq("B", "TGATG")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCAT"), fastq("B", "TGATGA")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCAT"), fastq("B", "TGATGAG")) ==
          barcode_key{ 0, 0 });
}

///////////////////////////////////////////////////////////////////////////////
// Exact matches for SE reads, when reads have mismatches

TEST_CASE("Exact matching reads with mismatches for SE barcodes",
          "[barcodes::exact::se]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(samples, 0, 0, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "AGCCA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACACA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCTA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCTT")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "NCCCA")) == barcode_key{ -1, -1 });

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "AGCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACACA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCTA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCTT"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Exact matching reads with mismatches for PE barcodes",
          "[barcodes::exact::pe]")
{
  const sample_set samples = {
    sample{ "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(samples, 0, 0, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "AGCCA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACACA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCTA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCTT")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "NCCCA")) == barcode_key{ -1, -1 });

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTA")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "GCCCA"), fastq("B", "GATTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTN")) ==
          barcode_key{ -1, -1 });
}

///////////////////////////////////////////////////////////////////////////////
// Inexact matching with global limits

TEST_CASE("Global limits override local limits for SE barcodes",
          "[barcodes::inexact::se]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(samples, 0, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "NCCCA")) == barcode_key{ -1, -1 });

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACACA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Global limits override local limits for PE barcodes",
          "[barcodes::inexact::se]")
{
  const sample_set samples = {
    sample{ "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(samples, 0, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCNA")) == barcode_key{ -1, -1 });

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACACA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTAG")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "NTTTC")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Mismatches in R1 only with SE table", "[barcodes::inexact::se]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(samples, 1, 1, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCT")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCN")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ANCCN")) == barcode_key{ -1, -1 });

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "GCTCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTG")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCNA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "NCCNA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Mismatches in R1 only with PE table", "[barcodes::inexact::se]")
{
  const sample_set samples = {
    sample{ "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(samples, 1, 1, 0);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCT")) == barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCN")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ANCCN")) == barcode_key{ -1, -1 });

  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "GCTCA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTG")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCN"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "NTTTC")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Mismatches in R2 only with SE table", "[barcodes::inexact::pe]")
{
  const sample_set samples = { sample{
    "sample", dna_sequence{ "ACCCA" }, dna_sequence{} } };

  const barcode_table table(samples, 1, 0, 1);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTT")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTG")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ANTTG")) ==
          barcode_key{ 0, 0 });
}

TEST_CASE("Mismatches in R2 only with PE table", "[barcodes::inexact::pe]")
{
  const sample_set samples = {
    sample{ "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(samples, 1, 0, 1);
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTT")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "ATTTG")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTN")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCN"), fastq("B", "GTTTT")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTNTT")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Mismatches in R1/R2 with PE table", "[barcodes::inexact::pe]")
{
  const sample_set samples = {
    sample{ "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(samples, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTAC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTAC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "TCCAA"), fastq("B", "GTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GGTGC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "ANCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTNC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACTCA"), fastq("B", "GTNTC")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Multiple mismatches in R1/R2 with PE table",
          "[barcodes::inexact::pe]")
{
  const sample_set samples = {
    sample{ "sample", dna_sequence{ "ACCCA" }, dna_sequence{ "GTTTC" } },
  };

  const barcode_table table(samples, 2, 2, 2);
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GTTAC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "GTTAC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCAA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCA"), fastq("B", "GGTGC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCT"), fastq("B", "CTTTC")) ==
          barcode_key{ -1, -1 });
  REQUIRE(table.identify(fastq("A", "TCCNA"), fastq("B", "GTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCCT"), fastq("B", "NTTTC")) ==
          barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "NCCCN"), fastq("B", "GTTTA")) ==
          barcode_key{ -1, -1 });
}

TEST_CASE("Ambiguous matches in R1 for SE table", "[barcodes::inexact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{} },
    sample{ "sample2", dna_sequence{ "TCCCA" }, dna_sequence{} },
    sample{ "sample3", dna_sequence{ "ACCCG" }, dna_sequence{} }
  };

  const barcode_table table(samples, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("A", "CCCCA")) == barcode_key{ -2, -2 });
  REQUIRE(table.identify(fastq("A", "ACCCC")) == barcode_key{ -2, -2 });
  REQUIRE(table.identify(fastq("A", "ACCCN")) == barcode_key{ -2, -2 });
}

TEST_CASE("Ambiguous matches in SE R1 for PE table", "[barcodes::inexact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "TGCGT" } },
    sample{ "sample2", dna_sequence{ "TCCCA" }, dna_sequence{ "AAGTT" } },
    sample{ "sample3", dna_sequence{ "ACCCG" }, dna_sequence{ "CGGCA" } },
  };

  const barcode_table table(samples, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("A", "CCCCA")) == barcode_key{ -2, -2 });
  REQUIRE(table.identify(fastq("A", "ACCCC")) == barcode_key{ -2, -2 });
  REQUIRE(table.identify(fastq("A", "ACCCN")) == barcode_key{ -2, -2 });
}

TEST_CASE("Mismatch resulting in apparent match", "[barcodes::inexact::se]")
{
  const sample_set samples = {
    sample{ "sample1", dna_sequence{ "ACCCA" }, dna_sequence{ "TGCGT" } },
    sample{ "sample2", dna_sequence{ "TCCTT" }, dna_sequence{ "AAGTT" } },
  };

  const barcode_table table(samples, 1, 1, 1);
  REQUIRE(table.identify(fastq("A", "ACCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "TCCTT")) == barcode_key{ 1, 0 });
  REQUIRE(table.identify(fastq("A", "TCCCA")) == barcode_key{ 0, 0 });
  REQUIRE(table.identify(fastq("A", "ACCTT")) == barcode_key{ 1, 0 });
}

} // namespace adapterremoval
