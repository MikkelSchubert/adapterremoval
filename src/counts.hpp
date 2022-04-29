/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2021 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#pragma once

#include <limits>
#include <vector>

#include "vecutils.hpp"

/** Helper class for working with simple count statistics. */
template<typename T>
class counts_tmpl
{
public:
  explicit counts_tmpl(size_t size = 0)
    : m_values(size, T())
  {
  }

  /** Increase the storage size to accomondate at least size items. */
  inline void resize_up_to(size_t size)
  {
    if (m_values.size() < size) {
      m_values.resize(size);
    }
  }

  /** Increment the specified value. */
  inline void inc(size_t n, T x = 1) { m_values.at(n) += x; }

  /** The sum of values in the specified range. */
  uint64_t sum(size_t from = 0,
               size_t to = std::numeric_limits<size_t>::max()) const
  {
    T total = 0;
    to = std::min(size(), to);
    for (size_t i = from; i < to; ++i) {
      total += m_values.at(i);
    }

    return total;
  }

  /** The sum of products from multiplying values with their 0-based indexes. */
  T product() const
  {
    T total = T();
    for (size_t i = 0; i < m_values.size(); ++i) {
      total += i * m_values.at(i);
    }

    return total;
  }

  /** Returns the count for n. */
  inline T get(size_t n) const
  {
    if (n >= m_values.size()) {
      return T();
    }

    return m_values.at(n);
  }

  /** Returns the number of values. */
  inline size_t size() const { return m_values.size(); }

  /** Return counts with trailing zero values trimmed. */
  counts_tmpl<T> trim() const
  {
    auto result = *this;

    while (result.m_values.size() && !result.m_values.back()) {
      result.m_values.pop_back();
    }

    return result;
  }

  /** + operator. */
  counts_tmpl<T> operator+(const counts_tmpl<T>& other) const
  {
    auto result = *this;
    result += other;

    return result;
  }

  /** += operator. */
  counts_tmpl<T>& operator+=(const counts_tmpl<T>& other)
  {
    merge_vectors(m_values, other.m_values);

    return *this;
  }

  /** Subtract a constants from the values. */
  counts_tmpl<T> operator-(const T value) const
  {
    auto result = *this;
    for (auto& it : result.m_values) {
      it -= value;
    }

    return result;
  }

  /** / operator. Always returns double values. */
  counts_tmpl<double> operator/(const counts_tmpl<T>& denom) const
  {
    const size_t length = std::max<size_t>(size(), denom.size());
    counts_tmpl<double> result(length);

    for (size_t i = 0; i < length; ++i) {
      result.inc(i, static_cast<double>(get(i)) / denom.get(i));
    }

    return result;
  }

  /** / operator for constants. Always returns double values. */
  counts_tmpl<double> operator/(T denom) const
  {
    counts_tmpl<double> result(size());

    for (size_t i = 0; i < size(); ++i) {
      result.inc(i, static_cast<double>(get(i)) / denom);
    }

    return result;
  }

private:
  std::vector<T> m_values;
};

/**
 * Counter for multiple classes of values, e.g. nucleotides. Allows for more
 * cache-efficient counting of several types along a sequence. However this
 * comes at the cost of making retrival of counter objects for each value type
 * more expensive, since a new instance has to be created.
 *
 * Template variable N defines the number of value classes.
 */
template<typename T, size_t N>
class indexed_counts
{
public:
  /* Creates a counter of the specified size for each value class. */
  explicit indexed_counts(size_t size = 0)
    : m_counts(size * N)
  {
  }

  /** Increase the storage size to accomondate at least size items. */
  inline void resize_up_to(size_t size) { m_counts.resize_up_to(size * N); }

  /** Increment the specified value for a given index. */
  inline void inc(size_t i, size_t n, T x = 1) { m_counts.inc(i + n * N, x); }

  /** Returns the count for n for a given index. */
  inline T get(size_t i, size_t n) const { return m_counts.get(i + n * N); }

  /** Create a standard counter object for a given value class */
  counts_tmpl<T> to_counts(size_t i) const
  {
    counts_tmpl<T> counts(size());
    for (size_t j = 0; j < counts.size(); ++j) {
      counts.inc(j, get(i, j));
    }

    return counts;
  }

  /** Merges all value classes into a single counter */
  counts_tmpl<T> merge() const
  {
    counts_tmpl<T> counts(size());
    for (size_t i = 0; i < m_counts.size(); ++i) {
      counts.inc(i / N, m_counts.get(i));
    }

    return counts;
  }

  /** Returns the number of values per index. */
  inline size_t size() const { return m_counts.size() / N; }

  /** += operator. */
  indexed_counts<T, N>& operator+=(const indexed_counts<T, N>& other)
  {
    m_counts += other.m_counts;

    return *this;
  }

private:
  counts_tmpl<T> m_counts;
};

//! Standard class for counts data
typedef counts_tmpl<int64_t> counts;
//! Standard class for rates, averages, etc.
typedef counts_tmpl<double> rates;

//! Counter indexed by ACGT nucleotides using ACGT_TO_IDX
using acgt_counts = indexed_counts<int64_t, 4>;
//! Counter indexed by ACGTN nucleotides using ACGTN_TO_IDX
using acgtn_counts = indexed_counts<int64_t, 5>;
