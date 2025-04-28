// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>

#include "errors.hpp" // for assertion_failed
#include "fastq.hpp"
#include "fastq_enc.hpp" // for fastq_encoding
#include "testing.hpp"   // for Catch

namespace adapterremoval {

TEST_CASE("fastq_encoding::phred_to_p")
{
  CHECK(fastq_encoding::phred_to_p(-0) == Approx(1.0));
  CHECK(fastq_encoding::phred_to_p(0) == Approx(1.0));
  CHECK(fastq_encoding::phred_to_p(1) == Approx(0.7943282));
  CHECK(fastq_encoding::phred_to_p(10) == Approx(0.1));
  CHECK(fastq_encoding::phred_to_p(20) == Approx(0.01));
  CHECK(fastq_encoding::phred_to_p(30) == Approx(0.001));
  CHECK(fastq_encoding::phred_to_p(90) == Approx(1e-9));
  // Supports higher phred scores than the 33 encoding
  CHECK(fastq_encoding::phred_to_p(100) == Approx(1e-10));

  CHECK_THROWS_AS(fastq_encoding::phred_to_p(-1), assert_failed);
}

TEST_CASE("fastq_encoding::p_to_phred")
{
  CHECK(fastq_encoding::p_to_phred(1.0) == Approx(0));
  CHECK(fastq_encoding::p_to_phred(0.7943282) == Approx(1));
  CHECK(fastq_encoding::p_to_phred(0.1) == Approx(10));
  CHECK(fastq_encoding::p_to_phred(0.01) == Approx(20));
  CHECK(fastq_encoding::p_to_phred(0.001) == Approx(30));
  CHECK(fastq_encoding::p_to_phred(1e-9) == Approx(90));
  // Supports higher phred scores than the 33 encoding
  CHECK(fastq_encoding::p_to_phred(1e-10) == Approx(100));

  CHECK_THROWS_AS(fastq_encoding::p_to_phred(-0.005), assert_failed);
}

TEST_CASE("fastq_encoding::p_to_phred_33")
{
  CHECK(fastq_encoding::p_to_phred_33(1.0) == '!');
  CHECK(fastq_encoding::p_to_phred_33(0.7943282) == '"');
  CHECK(fastq_encoding::p_to_phred_33(0.1) == '+');
  CHECK(fastq_encoding::p_to_phred_33(0.01) == '5');
  CHECK(fastq_encoding::p_to_phred_33(0.001) == '?');
  CHECK(fastq_encoding::p_to_phred_33(1e-9) == '{');
  // Supports higher phred scores than the 33 encoding
  CHECK(fastq_encoding::p_to_phred_33(1e-10) == '~');

  CHECK_THROWS_AS(fastq_encoding::p_to_phred_33(-0.005), assert_failed);
}

} // namespace adapterremoval