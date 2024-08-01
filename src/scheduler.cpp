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
#include "buffer.hpp"   // for buffer
#include "debug.hpp"    // for AR_REQUIRE, AR_FAIL
#include "fastq.hpp"    // for fastq
#include "logging.hpp"  // for error, log_stream, debug
#include "strutils.hpp" // for indent_lines
#include <algorithm>    // for max
#include <exception>    // for exception
#include <functional>   // for greater
#include <memory>       // for __shared_ptr_access, operator<, unique...
#include <system_error> // for system_error
#include <thread>       // for thread

#include "scheduler.hpp"

namespace adapterremoval {

/** Indicates the role of a worker thread */
enum class threadtype
{
  //! Mainly compute threads; may optionally perform IO
  cpu,
  //! Pure IO threads; compute work must be minimized
  io
};

///////////////////////////////////////////////////////////////////////////////
// analytical_chunk

void
analytical_chunk::add(const fastq& read)
{
  nucleotides += read.length();
  read.into_string(reads);
}

///////////////////////////////////////////////////////////////////////////////
// analytical_step

analytical_step::analytical_step(processing_order step_order, std::string name)
  : m_step_order(step_order)
  , m_name(std::move(name))
{
}

///////////////////////////////////////////////////////////////////////////////
// scheduler_step

using data_chunk = std::pair<size_t, chunk_ptr>;

/** Class wrapping a analytical_step. */
class scheduler_step
{
public:
  /** Constructor; name is used in error messages related to bugs */
  explicit scheduler_step(std::unique_ptr<analytical_step> value)
    : m_ptr(std::move(value))
  {
    AR_REQUIRE(m_ptr);
  }

  ~scheduler_step() = default;

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
  bool push_chunk(size_t chunk_id, chunk_ptr ptr)
  {
    m_queue.emplace_back(chunk_id, std::move(ptr));
    std::push_heap(m_queue.begin(), m_queue.end(), std::greater<>{});

    return m_next_chunk == chunk_id ||
           ordering() == processing_order::unordered;
  }

  /** Returns the oldest task; requires locking */
  data_chunk pop_chunk()
  {
    std::pop_heap(m_queue.begin(), m_queue.end(), std::greater<>{});
    auto task = std::move(m_queue.back());
    m_queue.pop_back();

    return task;
  }

  /** Name of the analytical task; for error messages when bugs are detected */
  const std::string& name() const { return m_ptr->name(); }

  /** Name of the analytical task; for error messages when bugs are detected */
  processing_order ordering() const { return m_ptr->ordering(); }

  /** Number of chunks waiting to be processed by this task; requires locking */
  size_t chunks_queued() const { return m_queue.size(); }

  /** Returns new ID for (sparse) output from ordered tasks; requires locking */
  size_t new_chunk_id() { return m_last_chunk++; }

  /** Processes a data chunk and returns the output chunks; assumes `can_run` */
  chunk_vec process(chunk_ptr chunk) const
  {
    return m_ptr->process(std::move(chunk));
  }

  /** Increment the next chunk to be processed; requires locking */
  bool increment_next_chunk()
  {
    m_next_chunk++;
    return can_run_next();
  }

  /** Perform any final cleanup associated with the task; requires locking */
  void finalize() const { return m_ptr->finalize(); }

  scheduler_step(const scheduler_step&) = delete;
  scheduler_step(scheduler_step&&) = delete;
  scheduler_step& operator=(const scheduler_step&) = delete;
  scheduler_step& operator=(scheduler_step&&) = delete;

private:
  //! Analytical step implementation
  std::unique_ptr<analytical_step> m_ptr = nullptr;
  //! The next chunk to be processed
  size_t m_next_chunk = 0;
  //! Running counter of output; for (possibly) sparse output of ordered steps
  size_t m_last_chunk = 0;
  //! (Ordered) vector of chunks to be processed
  std::vector<chunk_pair> m_queue{};
};

///////////////////////////////////////////////////////////////////////////////
// scheduler

size_t
scheduler::add_step(std::unique_ptr<analytical_step> step)
{
  AR_REQUIRE(step);

  const size_t step_id = m_steps.size();
  m_steps.push_back(std::make_shared<scheduler_step>(std::move(step)));

  return step_id;
}

bool
scheduler::run(int nthreads)
{
  AR_REQUIRE(!m_steps.empty());
  AR_REQUIRE(nthreads >= 1);
  AR_REQUIRE(!m_chunk_counter);
  // The last step added is assumed to be the initial/producing step
  AR_REQUIRE(m_steps.back()->ordering() == processing_order::ordered ||
             m_steps.back()->ordering() == processing_order::ordered_io);

  m_tasks_max = static_cast<size_t>(nthreads) * 3;

  std::vector<std::thread> threads;

  try {
    // CPU bound threads
    for (int i = 0; i < nthreads - 1; ++i) {
      threads.emplace_back(run_wrapper, this, threadtype::cpu);
    }

    // IO only threads
    for (int i = 0; i < 2; ++i) {
      threads.emplace_back(run_wrapper, this, threadtype::io);
    }
  } catch (const std::system_error& error) {
    log::error() << "Failed to create threads:\n" << indent_lines(error.what());
    m_errors = true;
  }

  // Run the main thread; the only calculation thread when "single-threaded"
  run_wrapper(this, threadtype::cpu);

  for (auto& thread : threads) {
    try {
      thread.join();
    } catch (const std::system_error& error) {
      log::error() << "Failed to join thread: " << error.what();
      m_errors = true;
    }
  }

  if (m_errors) {
    return false;
  }

  for (const auto& step : m_steps) {
    if (step->chunks_queued()) {
      log::error() << "Not all parts run for step " << step->name() << "; "
                   << step->chunks_queued() << " parts left";

      m_errors = true;
    }
  }

  if (m_errors) {
    return false;
  }

  // Steps are added in reverse order
  for (auto it = m_steps.rbegin(); it != m_steps.rend(); ++it) {
    try {
      (*it)->finalize();
    } catch (const std::exception&) {
      log::error() << "Failed to finalize task: " << (*it)->name();
      throw;
    }
  }

  return true;
}

void
scheduler::run_wrapper(scheduler* sch, threadtype thread_type)
{
  try {
    switch (thread_type) {
      case threadtype::cpu:
        sch->run_calc_loop();
        break;
      case threadtype::io:
        sch->run_io_loop();
        break;
      default:
        AR_FAIL("unexpected threadtype value");
    }
  } catch (const std::exception& error) {
    sch->m_errors = true;
    log::error() << error.what();
  } catch (...) {
    AR_FAIL("Unhandled, non-standard exception in thread");
  }

  // Signal any waiting threads
  sch->m_condition_calc.notify_all();
  sch->m_condition_io.notify_all();
}

void
scheduler::run_io_loop()
{
  std::unique_lock<std::mutex> lock(m_queue_lock);

  while (!m_errors) {
    step_ptr step;
    if (!m_queue_io.empty()) {
      step = m_queue_io.front();
      m_queue_io.pop();
    } else if (m_tasks || m_tasks_max) {
      m_condition_io.wait(lock);
      continue;
    } else {
      break;
    }

    const bool wake_thread = !m_queue_io.empty();
    data_chunk chunk = step->pop_chunk();

    lock.unlock();
    if (wake_thread) {
      m_condition_io.notify_one();
    }

    chunk_vec chunks = step->process(std::move(chunk.second));
    // Currently only (final) output steps
    AR_REQUIRE(chunks.empty());
    lock.lock();

    // Indicate that the next chunk can be processed
    if (step->increment_next_chunk()) {
      m_queue_io.push(step);
    }

    // One less task in memory
    m_tasks--;

    // Queue additional read tasks
    if (!m_queue_calc.empty() || m_tasks < m_tasks_max) {
      m_condition_calc.notify_one();
    }
  }

  m_condition_io.notify_all();
  m_condition_calc.notify_all();
}

void
scheduler::run_calc_loop()
{
  std::unique_lock<std::mutex> lock(m_queue_lock);

  while (!m_errors) {
    step_ptr step;
    if (!m_queue_calc.empty()) {
      // Otherwise try do do some non-IO work
      step = m_queue_calc.front();
      m_queue_calc.pop();
    } else if (m_tasks < m_tasks_max) {
      // If all (or no) tasks are running and IO is idle then read another chunk
      m_tasks++;
      step = m_steps.back();
      if (!step->push_chunk(m_chunk_counter++, chunk_ptr())) {
        continue;
      }
    } else if (m_tasks || m_tasks_max) {
      // There are either tasks running (which may produce new tasks) or tasks
      // that cannot yet be run due to IO already being active.
      m_condition_calc.wait(lock);
      continue;
    } else {
      break;
    }

    const bool wake_thread = !m_queue_calc.empty() || (m_tasks < m_tasks_max);
    data_chunk chunk = step->pop_chunk();

    lock.unlock();
    if (wake_thread) {
      m_condition_calc.notify_one();
    }

    chunk_vec chunks = step->process(std::move(chunk.second));
    lock.lock();

    if (chunks.empty() && step == m_steps.back()) {
      // The source has stopped producing chunks; nothing more to do
      m_tasks_max = 0;
    }

    // Schedule each of the resulting blocks
    for (auto& result : chunks) {
      AR_REQUIRE(result.first < m_steps.size());
      const step_ptr& recipient = m_steps.at(result.first);

      // Inherit reference count from source chunk
      auto next_id = chunk.first;
      if (step->ordering() != processing_order::unordered) {
        // Ordered steps are allowed to not return results, so the chunk
        // numbering is remembered for down-stream steps
        next_id = recipient->new_chunk_id();
      }

      if (recipient->push_chunk(next_id, std::move(result.second))) {
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
      if (step->increment_next_chunk()) {
        m_queue_calc.push(step);
      }
    }

    // One less task in memory
    m_tasks--;

    if (!m_queue_io.empty()) {
      m_condition_io.notify_one();
    }
  }

  m_condition_io.notify_all();
  m_condition_calc.notify_all();
}

} // namespace adapterremoval
