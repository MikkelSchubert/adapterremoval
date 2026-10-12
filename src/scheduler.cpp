// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#include "scheduler.hpp" // declarations
#include "debug.hpp"     // for AR_REQUIRE, AR_FAIL
#include "logging.hpp"   // for error, log_stream, debug
#include <algorithm>     // for pop_heap, push_heap
#include <cstddef>       // for size_t
#include <exception>     // for exception
#include <functional>    // for greater
#include <memory>        // for make_shared, unique_ptr
#include <mutex>         // for mutex, unique_lock
#include <ostream>       // for ostream
#include <string>        // for string
#include <system_error>  // for system_error
#include <thread>        // for thread
#include <utility>       // for move
#include <variant>       // for variant
#include <vector>        // for vector

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

// Out of line implementation to avoid adding dependency on fastq to every file
analytical_chunk::~analytical_chunk() = default;

///////////////////////////////////////////////////////////////////////////////
// analytical_step

analytical_step::analytical_step(processing_order step_order,
                                 std::string name,
                                 const step_id next_step)
  : m_step_order(step_order)
  , m_name(std::move(name))
  , m_next_step(next_step)
{
}

chunk_vec
analytical_step::signal()
{
  AR_FAIL("signal() called on " + m_name);
}

chunk_vec
analytical_step::process(chunk_ptr /* data */)
{
  AR_FAIL("process() called on " + m_name);
}

///////////////////////////////////////////////////////////////////////////////
// scheduler_step

using data_chunk = std::pair<size_t, chunk_ptr>;

/** Class wrapping an analytical_step */
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

  /** Adds chunk, returning true if the caller can queue it. Requires locking */
  [[nodiscard]] bool push_chunk(size_t chunk_id, chunk_ptr ptr)
  {
    AR_REQUIRE(ptr);
    m_queue.emplace_back(chunk_id, std::move(ptr));
    std::push_heap(m_queue.begin(), m_queue.end(), std::greater<>{});

    return chunk_id == m_next_chunk || !is_ordered();
  }

  /** Adds pending signal. Requires locking */
  void push_signal()
  {
    AR_REQUIRE(is_ordered() && !m_signals_suspended);
    m_pending_signals++;
  }

  /** Returns the oldest chunk or a signal; requires locking */
  [[nodiscard]] data_chunk pop_chunk()
  {
    if (!is_ordered() ||
        (!m_queue.empty() && m_queue.front().first == m_next_chunk)) {
      std::pop_heap(m_queue.begin(), m_queue.end(), std::greater<>{});
      auto task = std::move(m_queue.back());
      m_queue.pop_back();
      m_next_chunk++;

      return task;
    } else {
      AR_REQUIRE(m_pending_signals);
      m_pending_signals--;
      return { 0, nullptr };
    }
  }

  /** Name of the analytical task; for error messages when bugs are detected */
  [[nodiscard]] const std::string& name() const { return m_ptr->name(); }

  /** Returns the ordering of the analytical step */
  [[nodiscard]] processing_order ordering() const { return m_ptr->ordering(); }

  /** Returns true if the analytical step is ordered (IO or otherwise) */
  [[nodiscard]] bool is_ordered() const
  {
    return m_ptr->ordering() != processing_order::unordered;
  }

  /** Number of chunks waiting to be processed by this task; requires locking */
  [[nodiscard]] size_t chunks_queued() const { return m_queue.size(); }

  /** Returns new ID for (sparse) output from ordered tasks; requires locking */
  [[nodiscard]] size_t new_chunk_id() { return m_last_chunk++; }

  /** Processes a data chunk and returns the output chunks; assumes `can_run` */
  [[nodiscard]] chunk_vec process(chunk_ptr chunk) const
  {
    chunk_vec chunks;
    if (chunk) {
      chunks = m_ptr->process(std::move(chunk));
    } else {
      AR_REQUIRE(is_ordered());
      chunks = m_ptr->signal();
    }

    return chunks;
  }

  /** Suspend signalling and return N signals dropped; requires locking */
  [[nodiscard]] size_t suspend_signals(bool value)
  {
    AR_REQUIRE(m_is_queued);
    m_signals_suspended = value;

    if (m_signals_suspended && m_pending_signals) {
      const auto was_pending = m_pending_signals;
      m_pending_signals = 0;
      return was_pending;
    } else {
      return 0;
    }
  }

  /** Returns true if an ordered step has been queued; requires locking */
  [[nodiscard]] bool is_queued() const { return m_is_queued; }

  /** Sets flag indicating if the step has been queued. This has no effect on
   * unordered steps. Requires locking */
  void set_is_queued(bool value)
  {
    if (is_ordered()) {
      AR_REQUIRE(m_is_queued != value);
      m_is_queued = value;
    }
  }

  /** True if signalling of this step has been suspended; requires locking */
  [[nodiscard]] bool is_signalling_suspended() const
  {
    return m_signals_suspended;
  }

  /** Returns true if the step can be requeued after a task. Requires locking */
  [[nodiscard]] bool can_requeue_step() const
  {
    AR_REQUIRE(!m_is_queued);
    return is_ordered() &&
           (m_pending_signals ||
            (!m_queue.empty() && m_queue.front().first == m_next_chunk));
  }

  /** Perform any final cleanup associated with the task; requires locking */
  void finalize() const { m_ptr->finalize(); }

  /** Returns the ID of the step to receive the output from this step */
  [[nodiscard]] step_id next_step() const { return m_ptr->next_step(); }

  /** Set the signalling target for this step */
  void set_signal_target(step_id step, read_mate read)
  {
    if (read == read_mate::_1) {
      m_ptr->m_signal_targets.first = step;
    } else {
      m_ptr->m_signal_targets.second = step;
    }
  }

  [[nodiscard]] step_id get_signal_target(read_mate read) const
  {
    return m_ptr->signal_target(read);
  }

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
  //! Indicates if the step has been queued already (ordered only)
  bool m_is_queued = false;
  //! The number of pending signals; these are always processed last
  size_t m_pending_signals = 0;
  //! Is signalling of this step suspended?
  bool m_signals_suspended = false;
  //! (Ordered) vector of chunks to be processed
  std::vector<data_chunk> m_queue{};
};

///////////////////////////////////////////////////////////////////////////////
// scheduler

step_id
scheduler::add_step(std::unique_ptr<analytical_step> step)
{
  AR_REQUIRE(step);
  AR_REQUIRE(!m_tasks_max, "pipeline has already been started");

  m_steps.push_back(std::make_shared<scheduler_step>(std::move(step)));
  return step_id(m_steps.size() - 1);
}

void
scheduler::set_signal_target(const step_id src, step_id dst, read_mate read)
{
  AR_REQUIRE(src != step_id::UNDEFINED && m_steps.at(src.get())->is_ordered());
  AR_REQUIRE(dst == step_id::UNDEFINED ||
             (dst != src && m_steps.at(dst.get())->is_ordered()));
  AR_REQUIRE(!m_tasks_max, "pipeline has already been started");

  m_steps.at(src.get())->set_signal_target(dst, read);
}

void
scheduler::set_source_step(step_id step)
{
  AR_REQUIRE(!m_tasks_max, "pipeline has already been started");
  if (step == step_id::UNDEFINED) {
    m_source_step.reset();
  } else {
    m_source_step = m_steps.at(step.get());
    AR_REQUIRE(m_source_step->is_ordered(), "source step must be ordered");
  }
}

bool
scheduler::run(unsigned worker_threads, unsigned io_threads)
{
  AR_REQUIRE(!m_steps.empty());
  AR_REQUIRE(worker_threads >= 1);
  AR_REQUIRE(!m_tasks_max, "pipeline has already been started");

  if (!m_source_step) {
    // The last step added is defaults to being the initial/producing step
    m_source_step = m_steps.back();
  }

  AR_REQUIRE(m_source_step->is_ordered(), "source step must be ordered");

  m_tasks_max = static_cast<size_t>(worker_threads) * 3;

  std::vector<std::thread> threads;
  const auto run_wrapper = [](scheduler* sch, threadtype type) {
    sch->main_loop(type);
  };

  try {
    // CPU bound threads. The current thread executes one of them
    for (unsigned i = 0; i < worker_threads - 1; ++i) {
      threads.emplace_back(run_wrapper, this, threadtype::cpu);
    }

    // IO only threads (optional)
    for (unsigned i = 0; i < io_threads; ++i) {
      threads.emplace_back(run_wrapper, this, threadtype::io);
    }
  } catch (const std::system_error& error) {
    log::error() << "Failed to create worker threads";
    throw;
  }

  // Run the main thread; the only calculation thread when "single-threaded"
  main_loop(threadtype::cpu);

  for (auto& thread : threads) {
    try {
      thread.join();
    } catch (const std::system_error& error) {
      log::error() << "Failed to join worker threads";
      throw;
    }
  }

  for (const auto& step : m_steps) {
    if (step->chunks_queued()) {
      log::error() << "Not all parts run for step " << step->name() << "; "
                   << step->chunks_queued() << " parts left";
      return false;
    }
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
scheduler::main_loop(threadtype work_type)
{
  std::unique_lock lock{ m_queue_lock };
  const bool may_signal_source =
    work_type == threadtype::cpu ||
    m_source_step->ordering() == processing_order::ordered_io;

  while (true) {
    step_ptr step;
    if (work_type == threadtype::cpu && !m_queue_calc.empty()) {
      // Try to do some non-IO/mixed work
      step = m_queue_calc.front();
      m_queue_calc.pop();
    } else if (!m_queue_io.empty() &&
               (work_type == threadtype::io || m_sleeping_io_threads == 0)) {
      // Otherwise try to do some IO work; worker threads are allowed to steal
      // tasks from this queue if there are no available IO threads
      step = m_queue_io.front();
      m_queue_io.pop();
    } else if (may_signal_source && can_signal_source()) {
      // If all (or no) tasks are running and IO is idle then read another chunk
      m_tasks++;
      step = m_source_step;
      step->set_is_queued(true);
      step->push_signal();
    } else if (m_tasks || !m_source_step->is_signalling_suspended()) {
      // There are either tasks running (which may produce new tasks) and/or the
      // pipeline may signal for the source to produce more work
      if (work_type == threadtype::cpu) {
        if (can_wake_io_thread()) {
          m_condition_io.notify_one();
        }

        m_condition_calc.wait(lock);
      } else {
        if (can_wake_calc_thread()) {
          m_condition_calc.notify_one();
        }

        m_sleeping_io_threads++;
        m_condition_io.wait(lock);
        m_sleeping_io_threads--;
      }
      continue;
    } else {
      break;
    }

    const bool wake_io = can_wake_io_thread();
    const bool wake_calc = can_wake_calc_thread();
    data_chunk chunk = step->pop_chunk();
    AR_REQUIRE(chunk.second || !step->is_signalling_suspended());

    lock.unlock();
    if (wake_io) {
      m_condition_io.notify_one();
    }
    if (wake_calc) {
      m_condition_calc.notify_one();
    }

    auto chunks = step->process(std::move(chunk.second));
    lock.lock();

    // Schedule each of the resulting blocks
    for (auto& output : chunks) {
      step_id recipient_id = step->next_step();
      chunk_ptr output_chunk;
      if (auto* output_ptr = std::get_if<chunk_ptr>(&output)) {
        AR_REQUIRE(recipient_id != step_id::UNDEFINED);
        AR_REQUIRE(*output_ptr);

        output_chunk = std::move(*output_ptr);
      } else if (auto* recipient_id_ptr = std::get_if<step_id>(&output)) {
        recipient_id = *recipient_id_ptr;
        AR_REQUIRE(recipient_id == step->get_signal_target(read_mate::_1) ||
                   recipient_id == step->get_signal_target(read_mate::_2));
      } else if (auto* signal_ptr = std::get_if<scheduler_signal>(&output)) {
        AR_REQUIRE(step->is_ordered());
        m_tasks -=
          step->suspend_signals(*signal_ptr == scheduler_signal::suspend);
        continue;
      } else if (auto* output_pair = std::get_if<chunk_pair>(&output)) {
        AR_REQUIRE(recipient_id == step_id::UNDEFINED);
        AR_REQUIRE(output_pair->second);

        recipient_id = output_pair->first;
        output_chunk = std::move(output_pair->second);
      } else {
        AR_FAIL("unexpected variant");
      }

      AR_REQUIRE(recipient_id.get() < m_steps.size());
      const step_ptr& recipient = m_steps.at(recipient_id.get());

      bool queue_step = true;
      if (output_chunk) {
        // Ordered steps are allowed to not return results, so chunks are
        // re-numbered per recipient. Unordered tasks must return 1:1 results,
        // or never return results, so numbering is kept from upstream
        auto next_id =
          step->is_ordered() ? recipient->new_chunk_id() : chunk.first;

        queue_step = recipient->push_chunk(next_id, std::move(output_chunk));
      } else if (!recipient->is_signalling_suspended()) {
        recipient->push_signal();
      } else {
        continue;
      }

      m_tasks++;
      if (queue_step && !recipient->is_queued()) {
        queue_pending_step(recipient);
      }
    }

    // One less task in memory
    step->set_is_queued(false);
    m_tasks--;

    if (step->can_requeue_step()) {
      queue_pending_step(step);
    } else if (step == m_source_step && can_signal_source()) {
      queue_pending_step(step);
      step->push_signal();
      m_tasks++;
    }
  }

  m_condition_io.notify_all();
  m_condition_calc.notify_all();
}

bool
scheduler::can_signal_source() const
{
  return m_tasks < m_tasks_max && !m_source_step->is_queued() &&
         !m_source_step->is_signalling_suspended();
}

bool
scheduler::can_wake_calc_thread() const
{
  return !m_queue_calc.empty() || can_signal_source() ||
         (m_sleeping_io_threads == 0 && !m_queue_io.empty());
}

bool
scheduler::can_wake_io_thread() const
{
  return !m_queue_io.empty() ||
         (m_source_step->ordering() == processing_order::ordered_io &&
          can_signal_source());
}

void
scheduler::queue_pending_step(step_ptr step)
{
  step->set_is_queued(true);
  if (step->ordering() == processing_order::ordered_io) {
    m_queue_io.push(std::move(step));
  } else {
    m_queue_calc.push(std::move(step));
  }
}

std::ostream&
operator<<(std::ostream& os, step_id step)
{
  if (step == step_id::UNDEFINED) {
    return os << "step_id::UNDEFINED";
  } else {
    return os << "step_id{" << step.get() << "}";
  }
}

} // namespace adapterremoval
