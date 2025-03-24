// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2021 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "debug.hpp"        // for AR_REQUIRE
#include "utilities.hpp"    // for merge
#include <algorithm>        // for max, min
#include <array>            // for operator==, array
#include <cstddef>          // for size_t
#include <cstdint>          // for int64_t
#include <initializer_list> // for initializer_list
#include <limits>           // for numeric_limits, numeric_limits<>::is_iec559
#include <ostream>          // for ostream
#include <type_traits>      // for is_same, is_floating_point, is_integral
#include <utility>          // for move
#include <vector>           // for operator==, vector

namespace adapterremoval {

/** Helper class for working with simple count statistics. */
template<typename T>
class counts_tmpl
{
  static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be number type");

public:
  counts_tmpl(std::initializer_list<T> values)
    : m_counts(std::move(values))
  {
  }

  explicit counts_tmpl(size_t size = 0)
    : m_counts(size, T())
  {
  }

  /** Increase the storage size to accommodate at least size items. */
  inline void resize_up_to(size_t size)
  {
    if (m_counts.size() < size) {
      m_counts.resize(size);
    }
  }

  /** Increment the specified value. */
  inline void inc(size_t n, T count = 1) { m_counts.at(n) += count; }

  /** The sum of values in the specified range. */
  T sum(size_t from = 0, size_t to = std::numeric_limits<size_t>::max()) const
  {
    AR_REQUIRE(from <= to);
    T total = T();
    to = std::min(size(), to);
    for (size_t i = from; i < to; ++i) {
      total += m_counts.at(i);
    }

    return total;
  }

  /** The sum of products from multiplying values with their 0-based indexes. */
  T product() const
  {
    T total = T();
    for (size_t i = 0; i < m_counts.size(); ++i) {
      total += i * m_counts.at(i);
    }

    return total;
  }

  /** Returns the count for n. */
  inline T get(size_t n) const
  {
    if (n >= m_counts.size()) {
      return T();
    }

    return m_counts.at(n);
  }

  /** Returns the number of values. */
  inline size_t size() const { return m_counts.size(); }

  /** Return counts with trailing zero values trimmed. */
  counts_tmpl<T> trim() const
  {
    auto result = *this;

    while (result.m_counts.size() && !result.m_counts.back()) {
      result.m_counts.pop_back();
    }

    return result;
  }

  /* Normalize counts to frequencies in the range [0; 1] */
  counts_tmpl<double> normalize() const { return *this / sum(); }

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
    merge(m_counts, other.m_counts);

    return *this;
  }

  /** / operator. Always returns double values. */
  counts_tmpl<double> operator/(const counts_tmpl<T>& denom) const
  {
    AR_REQUIRE(size() == denom.size());
    counts_tmpl<double> result(size());

    for (size_t i = 0; i < size(); ++i) {
      result.inc(i, static_cast<double>(get(i)) / denom.get(i));
    }

    return result;
  }

  /** / operator for constants. Always returns double values. */
  counts_tmpl<double> operator/(T denom) const
  {
    static_assert(std::numeric_limits<double>::is_iec559, "IEC 559 assumed");
    counts_tmpl<double> result(size());

    for (size_t i = 0; i < size(); ++i) {
      result.inc(i, static_cast<double>(get(i)) / denom);
    }

    return result;
  }

  bool operator==(const counts_tmpl<T>& other) const
  {
    return m_counts == other.m_counts;
  }

  [[nodiscard]] auto begin() const { return m_counts.begin(); }

  [[nodiscard]] auto end() const { return m_counts.end(); }

private:
  std::vector<T> m_counts;
};

template<typename I, typename T = int64_t>
class indexed_count
{
  static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                "T must be number type");

public:
  using value_type = typename I::value_type;

  /* Creates a counter of the specified size for each value class. */
  indexed_count() = default;

  /** Returns the number of values. */
  inline size_t size() const { return I::indices; }

  /** Increment the specified value. */
  template<typename V>
  inline void inc(V index, T count = 1)
  {
    // Assert to catch help mixups of arguments, pending stronger typed solution
    static_assert(std::is_same<V, value_type>::value, "probably a mistake");

    m_counts.at(I::to_index(index)) += count;
  }

  /** Returns the count for n for a given index. */
  template<typename V>
  inline T get(V index) const
  {
    // Assert to catch help mixups of arguments, pending stronger typed solution
    static_assert(std::is_same<V, value_type>::value, "probably a mistake");

    return m_counts.at(I::to_index(index));
  }

  T sum() const
  {
    T total = T();
    for (auto count : m_counts) {
      total += count;
    }

    return total;
  }

  indexed_count<I, T>& operator+=(const indexed_count<I, T>& other)
  {
    adapterremoval::merge(m_counts, other.m_counts);

    return *this;
  }

  bool operator==(const indexed_count<I, T>& other) const
  {
    return m_counts == other.m_counts;
  }

private:
  std::array<T, I::indices> m_counts = {};
};

/**
 * Counter for multiple classes of values, e.g. nucleotides. Allows for more
 * cache-efficient counting of several types along a sequence. However this
 * comes at the cost of making retrieval of counter objects for each value type
 * more expensive, since a new instance has to be created.
 */
template<typename I, typename T = int64_t>
class indexed_counts
{
public:
  //! Type of indexed values
  using value_type = typename I::value_type;

  /* Creates a counter of the specified size for each value class. */
  explicit indexed_counts(size_t size = 0)
    : m_counts(size)
  {
  }

  /** Increase the storage size to accommodate at least size items. */
  inline void resize_up_to(size_t size)
  {
    if (m_counts.size() < size) {
      m_counts.resize(size);
    }
  }

  /** Increment the specified value for a given index. */
  template<typename V>
  inline void inc(V index, size_t offset, T count = 1)
  {
    // Assert to catch help mixups of arguments, pending stronger typed solution
    static_assert(std::is_same<V, value_type>::value, "probably a mistake");

    m_counts.at(offset).inc(index, count);
  }

  /** Increment the specified value for a given index (unchecked). */
  template<typename V>
  inline void inc_unsafe(V index, size_t offset, T count = 1)
  {
    // Assert to catch help mixups of arguments, pending stronger typed solution
    static_assert(std::is_same<V, value_type>::value, "probably a mistake");

    m_counts[offset].inc(index, count);
  }

  /** Returns the count for n for a given index. */
  inline const indexed_count<I, T>& get(size_t offset) const
  {
    return m_counts.at(offset);
  }

  /** Returns the count for n for a given index. */
  template<typename V>
  inline T get(V index, size_t offset) const
  {
    // Assert to catch help mixups of arguments, pending stronger typed solution
    static_assert(std::is_same<V, value_type>::value, "probably a mistake");

    if (offset >= size()) {
      return T();
    }

    return m_counts.at(offset).get(index);
  }

  /** Create a standard counter object for a given value class */
  template<typename V>
  counts_tmpl<T> to_counts(V v) const
  {
    // Assert to catch help mixups of arguments, pending stronger typed solution
    static_assert(std::is_same<V, value_type>::value, "probably a mistake");

    counts_tmpl<T> counts(size());
    for (size_t j = 0; j < counts.size(); ++j) {
      counts.inc(j, get(v, j));
    }

    return counts;
  }

  /** Merges all value classes into a single counter */
  counts_tmpl<T> merge() const
  {
    counts_tmpl<T> counts(m_counts.size());
    for (size_t i = 0; i < m_counts.size(); ++i) {
      counts.inc(i, m_counts.at(i).sum());
    }

    return counts;
  }

  /** Returns the number of values per index. */
  inline size_t size() const { return m_counts.size(); }

  /** += operator. */
  indexed_counts<I, T>& operator+=(const indexed_counts<I, T>& other)
  {
    resize_up_to(other.size());

    for (size_t i = 0; i < other.size(); ++i) {
      m_counts.at(i) += other.m_counts.at(i);
    }

    return *this;
  }

  bool operator==(const indexed_counts<I, T>& other) const
  {
    return m_counts == other.m_counts;
  }

private:
  static_assert(sizeof(indexed_count<I, T>) == I::indices * sizeof(T),
                "assuming no padding");

  std::vector<indexed_count<I, T>> m_counts;
};

//! Standard class for counts data
using counts = counts_tmpl<int64_t>;
//! Standard class for rates, averages, etc.
using rates = counts_tmpl<double>;

/** Stream operator for debugging output */
template<typename T>
std::ostream&
operator<<(std::ostream& os, const counts_tmpl<T>& value)
{
  os << "counts{[";
  bool first = true;
  for (const auto v : value) {
    if (!first) {
      os << ", ";
    }

    os << v;
    first = false;
  }

  return os << "]}";
}

/** Stream operator for debugging output */
template<typename I, typename T>
std::ostream&
operator<<(std::ostream& os, const indexed_count<I, T>& value)
{
  os << "indexed_count{";
  bool first = true;
  for (const auto key : I::values) {
    if (!first) {
      os << ", ";
    }

    os << key << "=" << value.get(key);
    first = false;
  }

  return os << "}";
}

/** Stream operator for debugging output */
template<typename I, typename T>
std::ostream&
operator<<(std::ostream& os, const indexed_counts<I, T>& value)
{
  bool first = true;
  os << "indexed_counts{";
  for (const auto key : I::values) {
    if (!first) {
      os << ", ";
    }

    os << key << "=[";
    for (size_t i = 0; i < value.size(); ++i) {
      if (i) {
        os << ", ";
      }

      os << value.get(key, i);
    }

    os << "]";
    first = false;
  }

  return os << "}";
}

} // namespace adapterremoval
