/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include "counts.hpp"
#include "testing.hpp"
#include <sstream>

namespace adapterremoval {

template<typename T>
std::ostream&
operator<<(std::ostream& os, counts_tmpl<T> const& value)
{
  os << "{";
  for (size_t i = 0; i < value.size(); ++i) {
    if (i) {
      os << ", ";
    }

    os << value.get(i);
  }

  return os << "}";
}

TEST_CASE("counts debug operator<<")
{
  std::ostringstream os;
  os << counts({ 7, 9, 13 });

  REQUIRE(os.str() == "{7, 9, 13}");
}

////////////////////////////////////////////////////////////////////////////////
// counts

TEST_CASE("counts constructors")
{
  SECTION("default")
  {
    counts c;
    REQUIRE(c.size() == 0);
    REQUIRE(c.product() == 0);
  }

  SECTION("with capacity")
  {
    const size_t N = 7;
    counts c = counts(N);
    REQUIRE(c.size() == 7);
    REQUIRE(c.product() == 0);
  }

  SECTION("from init list")
  {
    counts c = { 7, 9, 13 };
    REQUIRE(c.size() == 3);
    REQUIRE(c.get(0) == 7);
    REQUIRE(c.get(1) == 9);
    REQUIRE(c.get(2) == 13);
  }
}

TEST_CASE("counts resize_up_to")
{
  SECTION("empty")
  {
    counts c;
    REQUIRE(c.size() == 0);
    c.resize_up_to(7);
    REQUIRE(c.size() == 7);
    c.resize_up_to(1);
    REQUIRE(c.size() == 7);
  }

  SECTION("non-empty")
  {
    counts c = { 7, 9, 13 };
    REQUIRE(c.size() == 3);
    c.resize_up_to(7);
    REQUIRE(c == counts({ 7, 9, 13, 0, 0, 0, 0 }));
    c.resize_up_to(1);
    REQUIRE(c == counts({ 7, 9, 13, 0, 0, 0, 0 }));
  }
}

TEST_CASE("counts inc and get")
{
  counts c = counts(3);
  c.inc(1);
  REQUIRE(c.get(1) == 1);
  c.inc(1, 2);
  REQUIRE(c.get(1) == 3);
  REQUIRE(c.get(0) == 0);
  REQUIRE(c.get(2) == 0);
}

TEST_CASE("counts get beyond size")
{
  counts c = { 1 };
  REQUIRE(c.get(0) == 1);
  REQUIRE(c.get(1) == 0);
  REQUIRE(c.get(2) == 0);
}

TEST_CASE("counts sum")
{
  SECTION("empty")
  {
    const counts c;
    REQUIRE(c.sum() == 0);
  }

  SECTION("sum of values")
  {
    const counts c = { 1, 20, 300 };
    REQUIRE(c.sum() == 1 + 20 + 300);
  }

  SECTION("range")
  {
    const counts c = { 1, 20, 300 };
    REQUIRE(c.sum(1) == 20 + 300);
  }

  SECTION("range exceeding length")
  {
    const counts c = { 1, 20, 300 };
    REQUIRE(c.sum(1, 10) == 20 + 300);
  }

  SECTION("empty range")
  {
    const counts c = { 1, 20, 300 };
    REQUIRE(c.sum(1, 1) == 0);
  }
}

TEST_CASE("counts product")
{
  counts c = { 1, 20, 300 };
  REQUIRE(c.product() == 0 * 1 + 1 * 20 + 2 * 300);
}

TEST_CASE("counts trim empty")
{
  SECTION("empty")
  {
    counts c;
    REQUIRE(c.size() == 0);
    c = c.trim();
    REQUIRE(c.size() == 0);
  }

  SECTION("completely trimmed")
  {
    counts c = counts(10);
    REQUIRE(c.size() == 10);
    c = c.trim();
    REQUIRE(c.size() == 0);
  }

  SECTION("partially trimmed")
  {
    counts c = counts(10);
    REQUIRE(c.size() == 10);
    c.inc(5);
    c = c.trim();
    REQUIRE(c.size() == 6);
  }
}

TEST_CASE("counts normalize")
{
  using Catch::WithinAbs;

  SECTION("empty")
  {
    const counts c;
    const rates r = c.normalize();
    REQUIRE(r.size() == 0);
  }

  SECTION("non-empty")
  {
    const counts c = { 9, 7, 13 };
    const rates r = c.normalize();
    REQUIRE(r.size() == 3);
    REQUIRE_THAT(r.get(0), WithinAbs(9.0 / 29.0, 1e-6));
    REQUIRE_THAT(r.get(1), WithinAbs(7.0 / 29.0, 1e-6));
    REQUIRE_THAT(r.get(2), WithinAbs(13.0 / 29.0, 1e-6));
  }

  SECTION("zeros")
  {
    const counts c = { 0, 0, 0 };
    const rates r = c.normalize();
    REQUIRE(r.size() == 3);
    for (size_t i = 0; i < r.size(); ++i) {
      REQUIRE(std::isnan(r.get(i)));
    }
  }
}

TEST_CASE("counts operator+")
{
  REQUIRE(counts() + counts() == counts());
  REQUIRE(counts({ 7, 9, 13 }) + counts() == counts({ 7, 9, 13 }));
  REQUIRE(counts() + counts({ 7, 9, 13 }) == counts({ 7, 9, 13 }));
  REQUIRE(counts({ 100, 10 }) + counts({ 7, 9, 13 }) ==
          counts({ 107, 19, 13 }));
}

TEST_CASE("counts operator+=")
{
  SECTION("empty plus empty")
  {
    counts output;
    output += counts();
    REQUIRE(output == counts());
  }

  SECTION("non-empty plus empty")
  {
    counts output = counts({ 7, 9, 13 });
    output += counts();
    REQUIRE(output == counts({ 7, 9, 13 }));
  }

  SECTION("empty plus non-empty")
  {
    counts output;
    output += counts({ 7, 9, 13 });
    REQUIRE(output == counts({ 7, 9, 13 }));
  }

  SECTION("non-empty plus non-empty")
  {
    counts output = counts({ 100, 10 });
    output += counts({ 7, 9, 13 });
    REQUIRE(output == counts({ 107, 19, 13 }));
  }
}

TEST_CASE("counts operator/")
{
  using Catch::WithinAbs;

  SECTION("empty div empty")
  {
    REQUIRE(counts() / counts() == rates());
  }

  SECTION("non-empty div empty")
  {
    REQUIRE_THROWS_AS(counts({ 7, 9, 13 }) / counts(), assert_failed);
  }

  SECTION("empty div non-empty")
  {
    REQUIRE_THROWS_AS(counts() / counts({ 7, 9, 13 }), assert_failed);
  }

  SECTION("non-empty div non-empty")
  {
    const rates output = counts({ 7, 9, 13 }) / counts({ 1000, 100, 10 });
    REQUIRE(output.size() == 3);
    REQUIRE_THAT(output.get(0), WithinAbs(7.0 / 1000.0, 1e-6));
    REQUIRE_THAT(output.get(1), WithinAbs(9.0 / 100.0, 1e-6));
    REQUIRE_THAT(output.get(2), WithinAbs(13.0 / 10.0, 1e-6));
  }
}

TEST_CASE("counts operator==")
{
  const counts empty;
  const counts non_empty_1 = { 1, 2, 3 };
  const counts non_empty_2 = { 3, 2, 1 };

  REQUIRE(empty == empty);
  REQUIRE(non_empty_1 == non_empty_1);
  REQUIRE(!(empty == non_empty_1));
  REQUIRE(!(non_empty_1 == non_empty_2));
}

////////////////////////////////////////////////////////////////////////////////
// indexed_count

TEST_CASE("indexed_count<ACGT> constructors")
{
  indexed_count<ACGT> c;
  REQUIRE(c.size() == ACGT::size);
  REQUIRE(c.get('A') == 0);
  REQUIRE(c.get('C') == 0);
  REQUIRE(c.get('G') == 0);
  REQUIRE(c.get('T') == 0);
}

TEST_CASE("indexed_count<ACGT> inc and get")
{
  indexed_count<ACGT> c;
  c.inc('A', 2);
  c.inc('C', 30);
  c.inc('G', 500);
  c.inc('T', 9001);
  REQUIRE(c.get('A') == 2);
  REQUIRE(c.get('C') == 30);
  REQUIRE(c.get('G') == 500);
  REQUIRE(c.get('T') == 9001);
}

TEST_CASE("indexed_count<ACGTN> inc and get")
{
  indexed_count<ACGTN> c;
  c.inc('A', 2);
  c.inc('C', 30);
  c.inc('G', 500);
  c.inc('T', 9001);
  c.inc('N', 2022);
  REQUIRE(c.get('A') == 2);
  REQUIRE(c.get('C') == 30);
  REQUIRE(c.get('G') == 500);
  REQUIRE(c.get('T') == 9001);
  REQUIRE(c.get('N') == 2022);
}

TEST_CASE("indexed_count<ACGT> sum")
{
  SECTION("empty")
  {
    const indexed_count<ACGT> c;
    REQUIRE(c.sum() == 0);
  }

  SECTION("sum of values")
  {
    indexed_count<ACGT> c;
    c.inc('A', 2);
    c.inc('C', 30);
    c.inc('G', 500);
    c.inc('T', 9001);
    REQUIRE(c.sum() == 2 + 30 + 500 + 9001);
  }
}

TEST_CASE("indexed_count<ACGT> operator+=")
{
  indexed_count<ACGT> input_1;
  input_1.inc('A', 2);
  input_1.inc('C', 30);
  input_1.inc('G', 500);
  input_1.inc('T', 9001);

  indexed_count<ACGT> input_2;
  input_2.inc('A', 7);
  input_2.inc('C', 9);
  input_2.inc('G', 13);
  input_2.inc('T', 42);

  input_1 += input_2;
  REQUIRE(input_1.get('A') == 2 + 7);
  REQUIRE(input_1.get('C') == 30 + 9);
  REQUIRE(input_1.get('G') == 500 + 13);
  REQUIRE(input_1.get('T') == 9001 + 42);
}

TEST_CASE("indexed_count<ACGT> operator==")
{
  const indexed_count<ACGT> empty;
  indexed_count<ACGT> non_empty_1;
  indexed_count<ACGT> non_empty_2;

  non_empty_1.inc('A');
  non_empty_2.inc('T');

  REQUIRE(empty == empty);
  REQUIRE(non_empty_1 == non_empty_1);
  REQUIRE(!(empty == non_empty_1));
  REQUIRE(!(non_empty_1 == non_empty_2));
}

////////////////////////////////////////////////////////////////////////////////
// indexed_counts<ACGT>

TEST_CASE("indexed_counts<ACGT> constructors")
{
  SECTION("default")
  {
    indexed_counts<ACGT> c;
    REQUIRE(c.size() == 0);
  }

  SECTION("with capacity")
  {
    const size_t N = 7;
    indexed_counts<ACGT> c = indexed_counts<ACGT>(N);
    REQUIRE(c.size() == 7);
  }
}

TEST_CASE("indexed_counts<ACGT> resize_up_to")
{
  SECTION("empty")
  {
    indexed_counts<ACGT> c;
    REQUIRE(c.size() == 0);
    c.resize_up_to(7);
    REQUIRE(c.size() == 7);
    c.resize_up_to(1);
    REQUIRE(c.size() == 7);
  }

  SECTION("non-empty")
  {
    auto c = indexed_counts<ACGT>(2);
    REQUIRE(c.size() == 2);
    c.inc('A', 1);
    c.resize_up_to(5);
    REQUIRE(c.to_counts('A') == counts({ 0, 1, 0, 0, 0 }));
    REQUIRE(c.to_counts('C') == counts({ 0, 0, 0, 0, 0 }));
    REQUIRE(c.to_counts('G') == counts({ 0, 0, 0, 0, 0 }));
    REQUIRE(c.to_counts('T') == counts({ 0, 0, 0, 0, 0 }));
    c.resize_up_to(1);
    REQUIRE(c.to_counts('A') == counts({ 0, 1, 0, 0, 0 }));
    REQUIRE(c.to_counts('C') == counts({ 0, 0, 0, 0, 0 }));
    REQUIRE(c.to_counts('G') == counts({ 0, 0, 0, 0, 0 }));
    REQUIRE(c.to_counts('T') == counts({ 0, 0, 0, 0, 0 }));
  }
}

TEST_CASE("indexed_counts<ACGT> inc and get")
{
  indexed_counts<ACGT> c = indexed_counts<ACGT>(3);
  c.inc('A', 1);
  c.inc('C', 1, 7);
  c.inc('G', 2, 3);
  c.inc('T', 0, 6);
  c.inc('A', 1, 4);

  REQUIRE(c.get('A', 1) == 5);
  REQUIRE(c.get('C', 1) == 7);
  REQUIRE(c.get('G', 2) == 3);
  REQUIRE(c.get('T', 0) == 6);
}

TEST_CASE("indexed_counts<ACGT> get beyond size")
{
  indexed_counts<ACGT> c;
  REQUIRE(c.size() == 0);
  REQUIRE(c.get('A', 0) == 0);
  REQUIRE(c.get('A', 1) == 0);
  REQUIRE(c.get('A', 2) == 0);
}

TEST_CASE("indexed_counts<ACGT> to counts")
{
  indexed_counts<ACGT> c = indexed_counts<ACGT>(3);
  c.inc('A', 1);
  c.inc('C', 1, 7);
  c.inc('G', 2, 3);
  c.inc('T', 0, 6);
  REQUIRE(c.to_counts('A') == counts({ 0, 1, 0 }));
  REQUIRE(c.to_counts('C') == counts({ 0, 7, 0 }));
  REQUIRE(c.to_counts('G') == counts({ 0, 0, 3 }));
  REQUIRE(c.to_counts('T') == counts({ 6, 0, 0 }));
}

TEST_CASE("indexed_counts<ACGT> merge")
{
  indexed_counts<ACGT> c = indexed_counts<ACGT>(3);
  c.inc('A', 1);
  c.inc('C', 1, 7);
  c.inc('G', 2, 3);
  c.inc('T', 0, 6);
  REQUIRE(c.merge() == counts({ 6, 8, 3 }));
}

TEST_CASE("indexed_counts<ACGT> operator+=")
{
  auto non_empty = indexed_counts<ACGT>(3);
  non_empty.inc('A', 1, 5);
  non_empty.inc('C', 1, 7);
  non_empty.inc('G', 2, 3);
  non_empty.inc('T', 0, 6);

  SECTION("empty plus empty")
  {
    indexed_counts<ACGT> output;
    output += indexed_counts<ACGT>();
    REQUIRE(output == indexed_counts<ACGT>());
  }

  SECTION("non-empty plus empty")
  {
    indexed_counts<ACGT> output = non_empty;
    output += indexed_counts<ACGT>();
    REQUIRE(output == non_empty);
  }

  SECTION("empty plus non-empty")
  {
    indexed_counts<ACGT> output;
    output += non_empty;
    REQUIRE(output == non_empty);
  }

  SECTION("non-empty plus non-empty")
  {
    indexed_counts<ACGT> output = non_empty;
    output += non_empty;
    REQUIRE(output.to_counts('A') == counts({ 0, 10, 0 }));
    REQUIRE(output.to_counts('C') == counts({ 0, 14, 0 }));
    REQUIRE(output.to_counts('G') == counts({ 0, 0, 6 }));
    REQUIRE(output.to_counts('T') == counts({ 12, 0, 0 }));
  }
}

TEST_CASE("indexed_counts<ACGT> operator==")
{
  const indexed_counts<ACGT> empty;
  auto non_empty_1 = indexed_counts<ACGT>(1);
  auto non_empty_2 = indexed_counts<ACGT>(1);

  non_empty_1.inc('A', 0);
  non_empty_2.inc('T', 0);

  REQUIRE(empty == empty);
  REQUIRE(non_empty_1 == non_empty_1);
  REQUIRE(!(empty == non_empty_1));
  REQUIRE(!(non_empty_1 == non_empty_2));
}

} // namespace adapterremoval
