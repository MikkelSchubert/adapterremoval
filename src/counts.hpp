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
  {}

  /** Increase the storage size to accomondate at least size items. */
  void resize_up_to(size_t size)
  {
    if (m_values.size() < size) {
      m_values.resize(size);
    }
  }

  /** Increment the specified value. */
  void inc(size_t n, T x = 1) { m_values.at(n) += x; }

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
  T get(size_t n) const
  {
    if (n >= m_values.size()) {
      return T();
    }

    return m_values.at(n);
  }

  /** Returns the number of values. */
  size_t size() const { return m_values.size(); }

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

//! Standard class for counts data
typedef counts_tmpl<int64_t> counts;
//! Standard class for rates, averages, etc.
typedef counts_tmpl<double> rates;
