// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "debug.hpp"    // for AR_REQUIRE, AR_FAIL
#include <algorithm>    // for copy, max, copy_backward
#include <memory>       // for shared_ptr
#include <mutex>        // for unique_lock
#include <shared_mutex> // for shared_lock, shared_mutex
#include <vector>       // for vector

namespace adapterremoval {

/**
 * Simple thread-safe storage backed by a vector: Allows threads to acquire,
 * modify, and then release pre-created, mutable values. This is used to reduce
 * the amount of locking that is required to update statistics and similar.
 */
template<typename T>
class threadstate
{
public:
  /** Creates empty container */
  threadstate() = default;

  //! The underlying data type stored in the container
  using value_type = T;
  //! Exclusive pointer to the underlying data type
  using pointer = std::unique_ptr<value_type>;

  /** Create a new state value with the given arguments **/
  template<class... Args>
  void emplace_back(Args&&... args)
  {
    emplace_back_n(1, std::forward<Args>(args)...);
  }

  /** Create N new state values using the given arguments **/
  template<class... Args>
  void emplace_back_n(size_t n, Args&&... args)
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    for (; n; n--) {
      m_values.emplace_back(std::make_unique<T>(std::forward<Args>(args)...));
    }
  }

  /** Attempt to acquire a value; returns null if no value are available **/
  pointer try_acquire()
  {
    {
      std::lock_guard<std::mutex> lock(m_mutex);
      if (!m_values.empty()) {
        pointer value = std::move(m_values.back());
        m_values.pop_back();
        return value;
      }
    }

    return pointer{};
  }

  /** Attempt to acquire a value; abort if no values are available */
  pointer acquire()
  {
    if (auto value = try_acquire()) {
      return value;
    }

    AR_FAIL("could not acquire thread state");
  }

  /** Release ownership of a value */
  void release(pointer value)
  {
    AR_REQUIRE(value);
    std::lock_guard<std::mutex> lock(m_mutex);

    m_values.push_back(std::move(value));
  }

  /** Merge also values into the arguments, using `+=`, then drop them */
  void merge_into(T& value)
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    while (!m_values.empty()) {
      value += *m_values.back();
      m_values.pop_back();
    }
  }

private:
  //! Mutex used to protect access to `m_values`
  std::mutex m_mutex{};
  //! Pre-created thread state
  std::vector<pointer> m_values{};
};

/**
 * This class encapsulates data that must be updated in a threaded context; the
 * `reader` class provides shared read-only access, while the `writer` class
 * provides exclusive read/write access.
 *
 * Any mutable (member) variables used by the wrapped data class must either be
 * atomic or guarded by their own mutexes.
 */
template<typename T>
class threadsafe_data
{
  struct inner
  {
    template<class... Args>
    explicit inner(Args&&... args)
      : mutex()
      , data(std::forward<Args>(args)...)
    {
    }

    //! Mutable mutex to allow const access to the data
    mutable std::shared_mutex mutex;
    T data;
  };

  using data_ptr = std::shared_ptr<inner>;

public:
  /** The `reader`class provides non-exclusive read access to the inner data */
  template<typename U>
  class reader
  {
  public:
    /** Default move constructor */
    reader(reader&& ptr) = default;
    /** Default move assignment operator */
    reader& operator=(reader&& ptr) = default;
    /** Default destructor */
    ~reader() = default;

    /** Returns (shared) reference to the inner data */
    const U& operator*() const
    {
      AR_REQUIRE(m_inner);
      return m_inner->data;
    }

    /** Returns (shared) pointer to the inner data */
    const U* operator->() const
    {
      AR_REQUIRE(m_inner);
      return &m_inner->data;
    }

    /** Copies are disallowed to prevent the creation of additional locks */
    reader(const reader& ptr) = delete;
    /** Copies are disallowed to prevent the creation of additional locks */
    reader& operator=(const reader& ptr) = delete;

  private:
    friend class threadsafe_data;

    /** Locks the mutex in ptr non-exclusively; ptr MUST be non-NULL  */
    explicit reader(data_ptr ptr)
      : m_inner(std::move(ptr))
      , m_lock(m_inner->mutex)
    {
    }

    //! Pointer to the inner data
    data_ptr m_inner;
    //! Non-exclusive lock of the mutex protecting the data
    std::shared_lock<std::shared_mutex> m_lock;
  };

  /** The `reader`class provides exclusive write access to the inner data */
  template<typename U>
  class writer
  {
  public:
    /** Default move constructor */
    writer(writer&& ptr) = default;
    /** Default move assignment operator */
    writer& operator=(writer&& ptr) = default;
    /** Default destructor */
    ~writer() = default;

    /** Returns exclusive reference to the inner data */
    U& operator*()
    {
      AR_REQUIRE(m_inner);
      return m_inner->data;
    }

    /** Returns exclusive pointer to the inner data */
    U* operator->()
    {
      AR_REQUIRE(m_inner);
      return &m_inner->data;
    }

    /** Copies are disallowed to prevent deadlocks */
    writer(const writer& ptr) = delete;
    /** Copies are disallowed to prevent deadlocks */
    writer& operator=(const writer& ptr) = delete;

  private:
    friend class threadsafe_data;

    /** Locks the mutex in ptr exclusively; ptr MUST be non-NULL  */
    explicit writer(data_ptr ptr)
      : m_inner(std::move(ptr))
      , m_lock(m_inner->mutex)
    {
    }

    //! Pointer to the inner data; shared to prevent dangling pointers
    data_ptr m_inner;
    //! Exclusive lock of the mutex protecting the data
    std::unique_lock<std::shared_mutex> m_lock;
  };

  /** Create default initialized thread-safe shared data */
  threadsafe_data()
    : m_inner(std::make_unique<inner>())
  {
  }

  /** Create thread-safe shared data */
  explicit threadsafe_data(T value)
    : m_inner(std::make_unique<inner>(std::move(value)))
  {
  }

  /** Move and copy construction behaves like a shared_ptr */
  threadsafe_data(const threadsafe_data&) = default;
  /** Move and copy construction behaves like a shared_ptr */
  threadsafe_data& operator=(threadsafe_data&&) = default;
  /** Move and copy construction behaves like a shared_ptr */
  threadsafe_data(threadsafe_data&&) = default;
  /** Move and copy construction behaves like a shared_ptr */
  threadsafe_data& operator=(const threadsafe_data&) = default;

  /** Default destructor */
  ~threadsafe_data() = default;

  /** Returns a smart pointer that provides shared (read) access */
  [[nodiscard]] reader<T> get_reader() const
  {
    AR_REQUIRE(m_inner);
    return reader<T>{ m_inner };
  }

  /** Returns a smart pointer that provides exclusive (write) access */
  [[nodiscard]] writer<T> get_writer() const
  {
    AR_REQUIRE(m_inner);
    return writer<T>{ m_inner };
  }

private:
  //! Shared pointer to inner data; shared_ptr is used to prevent use-after-free
  std::shared_ptr<inner> m_inner;
};

} // namespace adapterremoval
