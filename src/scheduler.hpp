/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <memory>             // for unique_ptr, shared_ptr
#include <mutex>              // for mutex, lock_guard
#include <queue>              // for queue
#include <stddef.h>           // for size_t
#include <string>             // for string
#include <utility>            // for pair
#include <vector>             // for vector

#include "debug.hpp" // for AR_DEBUG_ASSERT

struct scheduler_step;
struct data_chunk;

/** Simple thread-safe storage backed by a vector. **/
template<typename T>
class threadstate
{
public:
  threadstate()
    : m_mutex()
    , m_values()
  {}

  typedef std::unique_ptr<T> pointer;

  /** Create new state value. **/
  template<class... Args>
  void emplace_back(Args&&... args)
  {
    std::lock_guard<std::mutex> lock(m_mutex);

    m_values.emplace_back(new T(args...));
  }

  /** Acquire ownership of a value. **/
  pointer acquire()
  {
    auto ptr = try_acquire();
    AR_DEBUG_ASSERT(ptr);

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

  /** Returns true if there are no values. **/
  bool empty() const
  {
    std::lock_guard<std::mutex> lock(m_mutex);

    return m_values.empty();
  }

private:
  mutable std::mutex m_mutex;
  std::vector<pointer> m_values;
};

/**
 * Base-class for data-chunks produced, processed and consumed by a pipeline.
 */
class analytical_chunk
{
public:
  /** Constructor; does nothing. */
  analytical_chunk();

  /** Destructor; does nothing. */
  virtual ~analytical_chunk();
};

typedef std::unique_ptr<analytical_chunk> chunk_ptr;
typedef std::pair<size_t, chunk_ptr> chunk_pair;
typedef std::vector<chunk_pair> chunk_vec;

/** Ordering of input for analytical steps. */
enum class processing_order
{
  //! Data must be consumed in the input order
  ordered,
  //! Data must be consumed in the input order and involves disk IO
  ordered_io,
  //! Data may be consumed in any order
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
  analytical_step(processing_order step_order);

  /** Destructor; does nothing in base class. **/
  virtual ~analytical_step();

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
  virtual chunk_vec process(analytical_chunk* chunk) = 0;

  /**
   * Called once the pipeline has been run to completion; this function is
   * called on nodes in the same order as the pipeline.
   */
  virtual void finalize();

  /** Returns the expected ordering (ordered / unordered) for input data. **/
  processing_order ordering() const;

  //! Copy construction not supported
  analytical_step(const analytical_step&) = delete;
  //! Assignment not supported
  analytical_step& operator=(const analytical_step&) = delete;

private:
  //! Stores the ordering of data chunks expected by the step
  const processing_order m_step_order;
};

/**
 * Multithreaded scheduler.
 *
 * See 'analytical_step' for information on implementing analyses.
 */
class scheduler
{
public:
  /** Constructor. */
  scheduler();

  /** Frees any object passed via 'add_step'. **/
  ~scheduler();

  /**
   * Adds a step to the pipeline.
   *
   * @param name Textual name for the (type) of step being added.
   * @param step A analytical step; is deleted when scheduler is destroyed.
   * @return The unique ID of the newly added step.
   *
   * The ID specified here is specified as the first value of 'chunk_pair's
   * in order to determine to which analytical step a chunk is assigned.
   **/
  size_t add_step(const std::string& name, analytical_step* step);

  /** Runs the pipeline with n threads; return false on error. */
  bool run(int nthreads);

  //! Copy construction not supported
  scheduler(const scheduler&) = delete;
  //! Assignment not supported
  scheduler& operator=(const scheduler&) = delete;

private:
  typedef std::shared_ptr<scheduler_step> step_ptr;
  typedef std::queue<step_ptr> runables;
  typedef std::vector<step_ptr> pipeline;

  /** Wrapper function which calls do_run on the provided thread. */
  static void run_wrapper(scheduler*);
  /** Work function; invoked by each thread. */
  void do_run();

  /** Returns true if an error has occurred, and the run should terminate. */
  bool errors_occured();
  /** Mark that an error has occurred, and that the run should terminate. */
  void set_errors_occured();

  //! Analytical steps
  pipeline m_steps;

  //! Condition used to signal the (potential) availability of work
  std::condition_variable m_condition;

  //! Counter used for sequential processing of data
  size_t m_chunk_counter;
  //! Count of currently live tasks
  size_t m_live_tasks;
  //! The maximum number of tasks to process simultanously
  size_t m_live_tasks_max;

  //! Lock used to control access to chunks
  std::mutex m_queue_lock;
  //! Queue used for currently runnable steps involving only calculations
  runables m_queue_calc;
  //! Queue used for currently runnable steps involving IO
  runables m_queue_io;

  //! Indicates if a thread is doing IO; access control through 'm_queue_lock'
  bool m_io_active;
  //! Set to indicate if errors have occurred
  std::atomic_bool m_errors;
};

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'analytical_step'

inline void
analytical_step::finalize()
{}

inline processing_order
analytical_step::ordering() const
{
  return m_step_order;
}

///////////////////////////////////////////////////////////////////////////////
// Implementations for 'scheduler'

inline bool
scheduler::errors_occured()
{
  return m_errors.load();
}

inline void
scheduler::set_errors_occured()
{
  m_errors.store(true);
}
