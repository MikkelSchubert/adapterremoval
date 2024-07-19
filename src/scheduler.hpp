/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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

#include <algorithm>          // for copy, max, copy_backward
#include <atomic>             // for atomic_bool
#include <condition_variable> // for condition_variable
#include <cstddef>            // for size_t
#include <memory>             // for unique_ptr, shared_ptr, make_unique
#include <mutex>              // for mutex, lock_guard
#include <queue>              // for queue
#include <string>             // for string
#include <type_traits>        // for is_base_of
#include <utility>            // for pair
#include <vector>             // for vector

namespace adapterremoval {

class scheduler_step;
enum class threadtype;

/** Simple thread-safe storage backed by a vector. **/
template<typename T>
class threadstate
{
public:
  threadstate() = default;

  using pointer = std::unique_ptr<T>;

  /** Create new state value. **/
  template<class... Args>
  void emplace_back(Args&&... args)
  {
    std::lock_guard<std::mutex> lock(m_mutex);

    m_values.emplace_back(std::make_unique<T>(std::forward<Args>(args)...));
  }

  /** Acquire ownership of a value. **/
  pointer acquire()
  {
    auto ptr = try_acquire();
    AR_REQUIRE(ptr);

    return ptr;
  }

  /** Try to acquire ownership of a value, returning null if there are none. */
  pointer try_acquire()
  {
    std::lock_guard<std::mutex> lock(m_mutex);

    pointer value;
    if (!m_values.empty()) {
      value = std::move(m_values.back());
      m_values.pop_back();
    }

    return value;
  }

  /** Release ownership of a value. **/
  void release(pointer& value)
  {
    std::lock_guard<std::mutex> lock(m_mutex);

    m_values.push_back(std::move(value));
  }

private:
  mutable std::mutex m_mutex{};
  std::vector<pointer> m_values{};
};

/**
 * Base-class for data-chunks produced, processed and consumed by a pipeline.
 */
class analytical_chunk
{
public:
  /** Constructor; does nothing. */
  analytical_chunk() = default;

  /** Destructor; does nothing. */
  virtual ~analytical_chunk() = default;

  analytical_chunk(const analytical_chunk&) = delete;
  analytical_chunk(analytical_chunk&&) = delete;
  analytical_chunk& operator=(const analytical_chunk&) = delete;
  analytical_chunk& operator=(analytical_chunk&&) = delete;
};

using chunk_ptr = std::unique_ptr<analytical_chunk>;
using chunk_pair = std::pair<size_t, chunk_ptr>;
using chunk_vec = std::vector<chunk_pair>;

/** Ordering of input for analytical steps. */
enum class processing_order
{
  //! Data must be consumed in the input order. May involve IO
  ordered,
  //! Data must be consumed in the input order and only involves IO
  ordered_io,
  //! Data may be consumed in any order. May involve IO
  unordered
};

/**
 * Base class for analytical steps in a pipeline.
 *
 * Each step must implement the 'process' function as described below; note
 * that this function may be called simultaneously by multiple threads, and that
 * thread-safe storage (e.g. threadstate) must be used for writable
 * resources used by the step.
 */
class analytical_step
{
public:
  /**
   * @param step_order Indicates the expected ordering of chunks; processing
   *                   steps are typically unordered, while IO is typically
   *                   ordered in order to ensure that output order matches
   *                   input order.
   */
  analytical_step(processing_order step_order, std::string name);

  /** Destructor; does nothing in base class. **/
  virtual ~analytical_step() = default;

  /**
   * Function called by pipeline to generate / process / consume data chunks.
   *
   * Initially, the first step in the pipeline will receive nullptr; during
   * subsequent cycles, the pipeline will return the value output from the
   * last step to the initial step, which may re-use it to avoid allocations;
   * if this is not done, the chunk must be freed by the first step.
   *
   * Best performance is therefore obtained if the chunk contains buffers
   * for all steps, and these can be re-used across cycles, thereby reducing
   * the number of (de)allocations that must be performed.
   *
   * To terminate the pipeline, the first step must cease to return chunks;
   * however, any other step MUST return valid chunks, even if no input data
   * was provided. This is to ensure that tracking of chunk ordering can be
   * maintained across steps.
   *
   * The only exceptions to this rule are steps which ONLY has unordered
   * downstream steps (that is, no chunk generated by this step will be
   * processed by an ordered step later in the pipeline).
   */
  virtual chunk_vec process(chunk_ptr chunk) = 0;

  /**
   * Called once the pipeline has been run to completion; this function is
   * called on nodes in the same order as the pipeline.
   */
  virtual void finalize() {}

  /** Returns the expected ordering (ordered / unordered) for input data. **/
  processing_order ordering() const { return m_step_order; }

  /** Returns the name of the analytical step (type) */
  const std::string& name() const { return m_name; }

  analytical_step(const analytical_step&) = delete;
  analytical_step(analytical_step&&) = delete;
  analytical_step& operator=(const analytical_step&) = delete;
  analytical_step& operator=(analytical_step&&) = delete;

private:
  //! Stores the ordering of data chunks expected by the step
  const processing_order m_step_order;
  //! Human readable name for step; for debugging purposes
  const std::string m_name;
};

/**
 * Multithreaded scheduler.
 *
 * See 'analytical_step' for information on implementing analyses.
 */
class scheduler
{
public:
  scheduler() = default;
  ~scheduler() = default;

  /**
   * Adds a step to the pipeline.
   *
   * @param args Arguments passed to the analytical step constructor.
   * @return The unique ID of the newly added step.
   *
   * The ID specified here is specified as the first value of 'chunk_pair's
   * in order to determine to which analytical step a chunk is assigned.
   **/
  template<typename T, typename... Args>
  size_t add(Args&&... args)
  {
    static_assert(std::is_base_of<analytical_step, T>(),
                  "requires analytical_step sub-class");
    return add_step(std::make_unique<T>(std::forward<Args>(args)...));
  }

  /** Runs the pipeline with n threads; return false on error. */
  bool run(int nthreads);

  scheduler(const scheduler&) = delete;
  scheduler(scheduler&&) = delete;
  scheduler& operator=(const scheduler&) = delete;
  scheduler& operator=(scheduler&&) = delete;

private:
  using step_ptr = std::shared_ptr<scheduler_step>;
  using step_queue = std::queue<step_ptr>;
  using pipeline = std::vector<step_ptr>;

  size_t add_step(std::unique_ptr<analytical_step> step);

  /** Wrapper function which calls do_run on the provided thread. */
  static void run_wrapper(scheduler*, threadtype);
  /** Work function for calculation heavy threads; invoked by each thread. */
  void run_calc_loop();
  /** Work function for IO threads; invoked by each thread. */
  void run_io_loop();

  //! Analytical steps
  pipeline m_steps{};

  //! Condition used to signal the (potential) availability of IO work
  std::condition_variable m_condition_io{};
  //! Condition used to signal the (potential) availability of calculation work
  std::condition_variable m_condition_calc{};

  //! Counter used for sequential processing of data
  size_t m_chunk_counter = 0;
  //! The current number of running/queued tasks
  size_t m_tasks = 0;
  //! The maximum number of tasks to process simultaneously
  size_t m_tasks_max = 0;

  //! Lock used to control access to member variables
  std::mutex m_queue_lock{};
  //! Queue used for currently runnable steps involving calculations
  step_queue m_queue_calc{};
  //! Queue used for currently runnable steps involving only IO
  step_queue m_queue_io{};

  //! Set to indicate if errors have occurred
  std::atomic_bool m_errors{ false };
};

} // namespace adapterremoval
