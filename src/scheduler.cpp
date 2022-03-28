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
// scheduler

struct data_chunk
{
  explicit data_chunk(size_t chunk_id_ = 0)
    : chunk_id(chunk_id_)
    , data()
  {}

  explicit data_chunk(size_t chunk_id_, chunk_ptr& data_)
    : chunk_id(chunk_id_)
    , data(std::move(data_))
  {}

  /** Sorts by id from high to low. **/
  bool operator<(const data_chunk& other) const
  {
    return chunk_id > other.chunk_id;
  }

  //! Strictly increasing counter; used to sort chunks for 'ordered' tasks
  size_t chunk_id;
  //! Use generated data; is normally not freed by this struct
  chunk_ptr data;
};

/** Simple priority queue; needed to allow moving values out of queue. */
class chunk_queue
{
public:
  chunk_queue()
    : m_chunks()
  {}

  template<class... Args>
  void emplace_back(Args&&... args)
  {
    m_chunks.emplace_back(args...);
    std::push_heap(m_chunks.begin(), m_chunks.end());
  }

  data_chunk pop()
  {
    std::pop_heap(m_chunks.begin(), m_chunks.end());
    data_chunk value = std::move(m_chunks.back());
    m_chunks.pop_back();

    return value;
  }

  const data_chunk& top() const { return m_chunks.front(); }

  bool empty() const { return m_chunks.empty(); }

  size_t size() const { return m_chunks.size(); }

private:
  typedef std::vector<data_chunk> chunk_vec;

  chunk_vec m_chunks;
};

struct scheduler_step
{
  scheduler_step(analytical_step* value, const std::string& name_)
    : ptr(value)
    , next_chunk(0)
    , last_chunk(0)
    , queue()
    , name(name_)
  {}

  bool can_run(size_t chunk_id) const
  {
    if (ptr->ordering() != processing_order::unordered) {
      return (next_chunk == chunk_id);
    }

    return true;
  }

  bool has_next() const
  {
    return queue.size() && queue.top().chunk_id == next_chunk;
  }

  //! Analytical step implementation
  std::unique_ptr<analytical_step> ptr;
  //! The next chunk to be processed
  size_t next_chunk;
  //! The last chunk queued to the step;
  //! Used to correct numbering for sparse output from sequential steps
  size_t last_chunk;
  //! (Ordered) vector of chunks to be processed
  chunk_queue queue;
  //! Short name for step used for error reporting
  std::string name;

  //! Copy construction not supported
  scheduler_step(const scheduler_step&) = delete;
  //! Assignment not supported
  scheduler_step& operator=(const scheduler_step&) = delete;
};

scheduler::scheduler()
  : m_steps()
  , m_condition()
  , m_chunk_counter(0)
  , m_tasks(0)
  , m_tasks_max(0)
  , m_live_tasks(0)
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
    if (!step->queue.empty()) {
      print_locker lock;
      std::cerr << "ERROR: Not all parts run for step " << step->name << "; "
                << step->queue.size() << " parts left" << std::endl;

      set_errors_occured();
    }
  }

  if (errors_occured()) {
    return false;
  }

  // Steps are added in reverse order
  for (auto it = m_steps.rbegin(); it != m_steps.rend(); ++it) {
    try {
      (*it)->ptr->finalize();
    } catch (const std::exception&) {
      std::cerr << "ERROR: Failed to finalizing task " << (*it)->name << ":\n";
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
      m_live_tasks++;
      m_io_active = true;
      step = m_steps.back();
      step->queue.emplace_back(m_chunk_counter++);
    } else if (m_live_tasks) {
      // There are either tasks running (which may produce new tasks) or tasks
      // that cannot yet be run due to IO already being active.
      m_condition.wait(lock);
      continue;
    } else {
      // We break if there is no live tasks to prevent bugs in the scheduling
      // from causing infinite loops. Leftover tasks are caught in `run`.
      break;
    }

    do {
      data_chunk chunk = step->queue.pop();

      lock.unlock();
      chunk_vec chunks = step->ptr->process(std::move(chunk.data));
      lock.lock();

      if (chunks.empty() && step == m_steps.back()) {
        // The source has stopped producing chunks; nothing more to do
        m_tasks_max = 0;
      }

      // Schedule each of the resulting blocks
      for (auto& result : chunks) {
        AR_DEBUG_ASSERT(result.first < m_steps.size());
        step_ptr& recipient = m_steps.at(result.first);

        // Inherit reference count from source chunk
        auto next_id = chunk.chunk_id;
        if (step->ptr->ordering() != processing_order::unordered) {
          // Ordered steps are allowed to not return results, so the chunk
          // numbering is remembered for down-stream steps
          next_id = recipient->last_chunk++;
        }

        recipient->queue.emplace_back(next_id, std::move(result.second));

        if (recipient->can_run(next_id)) {
          if (recipient->ptr->ordering() == processing_order::ordered_io) {
            m_queue_io.push(recipient);
          } else {
            m_queue_calc.push(recipient);
          }

          m_live_tasks++;
        }

        m_tasks++;
      }

      if (chunks.size()) {
        // This worked better than trying to be smart about it
        m_condition.notify_all();
      }

      if (step->ptr->ordering() != processing_order::unordered) {
        // Indicate that the next chunk can be processed
        step->next_chunk++;
      }

      // One less task in memory
      m_tasks--;

      // If possible continue processing this task using the same thread
    } while (step->ptr->ordering() != processing_order::unordered &&
             step->has_next() && !errors_occured());

    // Decrement number of running/runnable tasks
    m_live_tasks--;

    // Unlock use of IO steps after finishing processing
    if (step->ptr->ordering() == processing_order::ordered_io) {
      m_io_active = false;
    }
  }

  // Signal any waiting threads
  m_condition.notify_all();
}
