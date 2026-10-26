// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "debug.hpp"    // for AR_REQUIRE, AR_FAIL
#include <algorithm>    // for copy, max, copy_backward
#include <memory>       // for shared_ptr
#include <mutex>        // for unique_lock
#include <shared_mutex> // for shared_lock, shared_mutex
#include <string>       // for string
#include <string_view>  // for string_view
#include <type_traits>  // for decay_t, enable_if_t
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

  //! The underlying datatype stored in the container
  using value_type = T;
  //! Exclusive pointer to the underlying datatype
  using pointer = std::unique_ptr<value_type>;

  /** Create a new state value with the given arguments **/
  template<class... Args>
  void emplace_back(Args&&... args)
  {
    emplace_back_n(1, std::forward<Args>(args)...);
  }

  /** Create N new state values using the given arguments **/
  template<class... Args>
  void emplace_back_n(size_t n, const Args&... args)
  {
    const std::unique_lock lock(m_mutex);
    for (; n; n--) {
      m_values.emplace_back(std::make_unique<T>(args...));
    }
  }

  /** Attempt to acquire a value; returns null if no value are available **/
  pointer try_acquire()
  {
    {
      const std::unique_lock lock(m_mutex);
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
    const std::unique_lock lock(m_mutex);

    m_values.push_back(std::move(value));
  }

  /** Merge also values into the arguments, using `+=`, then drop them */
  void merge_into(T& value)
  {
    const std::unique_lock lock(m_mutex);
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
  /** The `reader` class provides non-exclusive read access to the inner data */
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

  /** The `writer` class provides exclusive write access to the inner data */
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
    : m_inner(std::make_shared<inner>())
  {
  }

  /** Create thread-safe shared data; cannot be used with initializer lists */
  template<class Arg,
           class... Args,
           // exclude self to avoid overlap with copy constructor
           std::enable_if_t<!std::is_same_v<std::decay_t<Arg>, threadsafe_data>,
                            bool> = true>
  explicit threadsafe_data(Arg&& arg, Args&&... args)
    : m_inner(std::make_shared<inner>(std::forward<Arg>(arg),
                                      std::forward<Args>(args)...))
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
  [[nodiscard]] writer<T> get_writer()
  {
    AR_REQUIRE(m_inner);
    return writer<T>{ m_inner };
  }

private:
  //! Shared pointer to inner data; shared_ptr is used to prevent use-after-free
  std::shared_ptr<inner> m_inner;
};

/** Called to log and validate the final `reusable_objects` state */
void
finalize_reusable_objects(std::string_view description,
                          size_t n_allocated,
                          size_t n_released,
                          size_t n_acquired);

/**
 * Template for reusable objects. These must be default constructible, movable,
 * and have a `reserve`, `resize`, and `capacity` function, and users must take
 * care to clear any values that they might contain, before using them. This is
 * intended to reduce the number of allocations made in hot paths, which would
 * otherwise significantly impact overall throughput.
 *
 * CAPACITY specifies the default capacity for these objects; users should
 * avoid exceeding this capacity, to avoid increasing the memory overhead.
 *
 * All objects must be returned to the cache, before `finalize` is called.
 */
template<typename T>
class reusable_objects
{
public:
  /** Create object cache with the given description used for logging */
  reusable_objects(std::string_view description, size_t capacity)
    : m_description(description)
    , m_capacity(capacity)
  {
  }

  /**
   * Acquire object with the given size and (at least) the default capacity. May
   * return a new object or a previously released object. `size` may not exceed
   * the default capacity.
   */
  T acquire(const size_t size)
  {
    AR_REQUIRE(size <= m_capacity);
    T value;

    {
      auto writer = m_state.get_writer();

      writer->n_acquired++;
      if (writer->objects.empty()) {
        writer->n_allocated++;
      } else {
        value = std::move(writer->objects.back());
        writer->objects.pop_back();
      }
    }

    value.reserve(m_capacity);
    value.resize(size);

    return value;
  }

  /**
   * Return value to the cache; verifies that capacity is unchanged. Returning
   * an object with no capacity is a no-op.
   */
  void release(T&& value)
  {
    if (value.capacity()) {
      AR_REQUIRE(value.capacity() >= m_capacity);

      auto writer = m_state.get_writer();
      // catch objects created elsewhere being released
      AR_REQUIRE(writer->objects.size() < writer->n_allocated);
      writer->objects.emplace_back(std::move(value));
    }
  }

  /**
   * Log various performance metrics about the cache and verifies that all
   * allocated objects have been released.
   */
  void finalize() const
  {
    auto reader = m_state.get_reader();
    finalize_reusable_objects(m_description,
                              reader->n_allocated,
                              reader->objects.size(),
                              reader->n_acquired);
  }

private:
  /** Internal state of the reusable object cache */
  struct reusable_state
  {
    //! Reusable objects; order is undefined
    std::vector<T> objects{};

    //! The number of times `acquire()` has been called
    size_t n_acquired = 0;
    //! The number of times objects were allocated when `acquire()` was called
    size_t n_allocated = 0;
  };

  //! Description of this cache
  std::string m_description;
  //! Default capacity for the reusable objects
  size_t m_capacity;

  //! Thread-safe wrapper around internal state; only member allowed to change
  threadsafe_data<reusable_state> m_state{};
};

} // namespace adapterremoval
