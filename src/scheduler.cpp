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
#include <algorithm>    // for max
#include <exception>    // for exception
#include <iostream>     // for operator<<, basic_ostream, endl, cerr
#include <system_error> // for system_error
#include <thread>       // for thread

#include "debug.hpp" // for AR_DEBUG_ASSERT
#include "scheduler.hpp"
#include "strutils.hpp" // for cli_formatter
#include "threads.hpp"  // for print_locker, thread_abort

///////////////////////////////////////////////////////////////////////////////
// analytical_chunk

analytical_chunk::analytical_chunk() {}

analytical_chunk::~analytical_chunk() {}

///////////////////////////////////////////////////////////////////////////////
// analytical_step

analytical_step::analytical_step(processing_order step_order)
  : m_step_order(step_order)
{}

analytical_step::~analytical_step() {}

///////////////////////////////////////////////////////////////////////////////
// scheduler_step

typedef std::pair<size_t, chunk_ptr> data_chunk;

/** Class wrapping a analytical_step. */
class scheduler_step
{
public:
  /** Constructor; name is used in error messages related to bugs */
  scheduler_step(analytical_step* value, const std::string& name)
    : m_ptr(value)
    , m_next_chunk(0)
    , m_last_chunk(0)
    , m_queue()
    , m_name(name)
  {}

  /** Returns true if the chunk can be executed now; requires locking */
  bool can_run(size_t chunk_id) const
  {
    return ordering() == processing_order::unordered ||
           m_next_chunk == chunk_id;
  }

  /**
   * Returns true if the next chunk is queued; only applicable to ordered tasks,
   * since unordered tasks are always added to a scheduler queue and looping on
   * the same task will result in queues going out of sync. Requires locking.
   */
  bool can_run_next() const
  {
    return ordering() != processing_order::unordered && !m_queue.empty() &&
           m_queue.front().first == m_next_chunk;
  }

  /** FIFO queues a task; requires locking */
  void push_chunk(size_t id, chunk_ptr ptr)
  {
    m_queue.emplace_back(id, std::move(ptr));
    std::push_heap(m_queue.begin(), m_queue.end(), std::greater<chunk_pair>{});
  }

  /** Returns the oldest task; requires locking */
  data_chunk pop_chunk()
  {
    std::pop_heap(m_queue.begin(), m_queue.end(), std::greater<chunk_pair>{});
    auto task = std::move(m_queue.back());
    m_queue.pop_back();

    return task;
  }

  /** Name of the analytical task; for error messages when bugs are detected */
  const std::string name() const { return m_name; }
  /** Name of the analytical task; for error messages when bugs are detected */
  processing_order ordering() const { return m_ptr->ordering(); }
  /** Number of chunks waiting to be processed by this task; requires locking */
  size_t chunks_queued() const { return m_queue.size(); }
  /** Returns new ID for (sparse) output from ordered tasks; requires locking */
  size_t new_chunk_id() { return m_last_chunk++; }

  /** Processes a data chunk and returns the output chunks; assumes `can_run` */
  chunk_vec process(chunk_ptr chunk)
  {
    return m_ptr->process(std::move(chunk));
  }

  /** Increment the next chunk to be processed; requires locking */
  void increment_next_chunk() { m_next_chunk++; }
  /** Perform any final cleanup assosited with the task; requires locking */
  void finalize() { return m_ptr->finalize(); }

private:
  //! Copy construction not supported
  scheduler_step(const scheduler_step&) = delete;
  //! Assignment not supported
  scheduler_step& operator=(const scheduler_step&) = delete;

  //! Analytical step implementation
  std::unique_ptr<analytical_step> m_ptr;
  //! The next chunk to be processed
  size_t m_next_chunk;
  //! Running counter of output; for (possibly) sparse output of ordered steps
  size_t m_last_chunk;
  //! (Ordered) vector of chunks to be processed
  std::vector<chunk_pair> m_queue;
  //! Short name for step used for error reporting
  std::string m_name;
};

///////////////////////////////////////////////////////////////////////////////
// scheduler

scheduler::scheduler()
  : m_steps()
  , m_condition()
  , m_chunk_counter(0)
  , m_tasks(0)
  , m_tasks_max(0)
  , m_queue_lock()
  , m_queue_calc()
  , m_queue_io()
  , m_io_active(false)
  , m_errors(false)
{}

scheduler::~scheduler() {}

size_t
scheduler::add_step(const std::string& name, analytical_step* step)
{
  AR_DEBUG_ASSERT(step);

  const size_t step_id = m_steps.size();
  m_steps.emplace_back(new scheduler_step(step, name));

  return step_id;
}

bool
scheduler::run(int nthreads)
{
  AR_DEBUG_ASSERT(!m_steps.empty());
  AR_DEBUG_ASSERT(nthreads >= 1);
  AR_DEBUG_ASSERT(!m_chunk_counter);

  m_tasks_max = static_cast<size_t>(nthreads) * 3;

  std::vector<std::thread> threads;

  try {
    for (int i = 0; i < nthreads - 1; ++i) {
      threads.emplace_back(run_wrapper, this);
    }
  } catch (const std::system_error& error) {
    print_locker lock;
    std::cerr << "ERROR: Failed to create threads:\n"
              << cli_formatter::fmt(error.what()) << std::endl;

    set_errors_occured();
  }

  // Run the main thread (the only thread in case of non-threaded mode)
  run_wrapper(this);

  for (auto& thread : threads) {
    try {
      thread.join();
    } catch (const std::system_error& error) {
      std::cerr << "ERROR: Failed to join thread: " << error.what()
                << std::endl;
      set_errors_occured();
    }
  }

  if (errors_occured()) {
    return false;
  }

  for (auto& step : m_steps) {
    if (step->chunks_queued()) {
      print_locker lock;
      std::cerr << "ERROR: Not all parts run for step " << step->name() << "; "
                << step->chunks_queued() << " parts left" << std::endl;

      set_errors_occured();
    }
  }

  if (errors_occured()) {
    return false;
  }

  // Steps are added in reverse order
  for (auto it = m_steps.rbegin(); it != m_steps.rend(); ++it) {
    try {
      (*it)->finalize();
    } catch (const std::exception&) {
      std::cerr << "ERROR: Failed to finalize " << (*it)->name() << ":\n";
      throw;
    }
  }

  return true;
}

void
scheduler::run_wrapper(scheduler* sch)
{
  try {
    return sch->do_run();
  } catch (const thread_abort&) {
    print_locker lock;
    std::cerr << "Aborting thread due to error." << std::endl;
  } catch (const std::exception& error) {
    print_locker lock;
    std::cerr << "ERROR: Unhandled exception in thread:\n"
              << cli_formatter::fmt(error.what()) << std::endl;
  } catch (...) {
    print_locker lock;
    std::cerr << "ERROR: Unhandled, non-standard exception in thread"
              << std::endl;
  }

  sch->set_errors_occured();
  sch->m_condition.notify_all();
}

void
scheduler::do_run()
{
  std::unique_lock<std::mutex> lock(m_queue_lock);

  while (!errors_occured()) {
    step_ptr step;
    if (!m_io_active && !m_queue_io.empty()) {
      // Try to keep the disk busy by preferring IO tasks
      step = m_queue_io.front();
      m_queue_io.pop();
      m_io_active = true;
    } else if (!m_queue_calc.empty()) {
      // Otherwise try do do some non-IO work
      step = m_queue_calc.front();
      m_queue_calc.pop();
    } else if (!m_io_active && m_tasks < m_tasks_max) {
      // If all (or no) tasks are running and IO is idle then read another chunk
      m_tasks++;
      m_io_active = true;
      step = m_steps.back();
      step->push_chunk(m_chunk_counter++, chunk_ptr());
    } else if (m_tasks) {
      // There are either tasks running (which may produce new tasks) or tasks
      // that cannot yet be run due to IO already being active.
      m_condition.wait(lock);
      continue;
    } else {
      break;
    }

    do {
      const bool wake_thread =
        // CPU bound tasks can always be run (if there are any idle threads)
        m_queue_calc.size() ||
        // IO bound tasks can be run if IO is idle or we are low on tasks
        (!m_io_active && (m_queue_io.size() || m_tasks < m_tasks_max));

      data_chunk chunk = step->pop_chunk();

      lock.unlock();
      if (wake_thread) {
        m_condition.notify_one();
      }

      chunk_vec chunks = step->process(std::move(chunk.second));
      lock.lock();

      if (chunks.empty() && step == m_steps.back()) {
        // The source has stopped producing chunks; nothing more to do
        m_tasks_max = 0;
      }

      // Schedule each of the resulting blocks
      for (auto& result : chunks) {
        AR_DEBUG_ASSERT(result.first < m_steps.size());
        const step_ptr& recipient = m_steps.at(result.first);

        // Inherit reference count from source chunk
        auto next_id = chunk.first;
        if (step->ordering() != processing_order::unordered) {
          // Ordered steps are allowed to not return results, so the chunk
          // numbering is remembered for down-stream steps
          next_id = recipient->new_chunk_id();
        }

        recipient->push_chunk(next_id, std::move(result.second));
        if (recipient->can_run(next_id)) {
          if (recipient->ordering() == processing_order::ordered_io) {
            m_queue_io.push(recipient);
          } else {
            m_queue_calc.push(recipient);
          }
        }

        m_tasks++;
      }

      if (step->ordering() != processing_order::unordered) {
        // Indicate that the next chunk can be processed
        step->increment_next_chunk();
      }

      // One less task in memory
      m_tasks--;

      // If possible continue processing this task using the same thread
    } while (step->can_run_next() && !errors_occured());

    // Unlock use of IO steps after finishing processing
    if (step->ordering() == processing_order::ordered_io) {
      m_io_active = false;
    }
  }

  // Signal any waiting threads
  m_condition.notify_all();
}
