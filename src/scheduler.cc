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

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <unistd.h>

#include "debug.h"
#include "scheduler.h"
#include "strutils.h"

namespace ar
{

///////////////////////////////////////////////////////////////////////////////
// analytical_chunk

analytical_chunk::analytical_chunk()
{
}

analytical_chunk::~analytical_chunk()
{
}


///////////////////////////////////////////////////////////////////////////////
// analytical_step

analytical_step::analytical_step(ordering step_order, bool file_io)
    : m_step_order(step_order)
    , m_file_io(file_io)
{
}

analytical_step::~analytical_step()
{
}


///////////////////////////////////////////////////////////////////////////////
// scheduler


struct data_chunk
{
    explicit data_chunk(size_t chunk_id_ = 0)
      : chunk_id(chunk_id_)
      , data()
      , counter(new bool())
    {
    }

    explicit data_chunk(const data_chunk& parent, chunk_ptr data_)
      : chunk_id(parent.chunk_id)
      , data(std::move(data_))
      , counter(parent.counter)
    {
    }

    /** Sorts by counter, data, type, in that order. **/
    bool operator<(const data_chunk& other) const
    {
        if (chunk_id != other.chunk_id) {
            return chunk_id > other.chunk_id;
        } else if (data != other.data) {
            return data > other.data;
        } else if (counter != other.counter) {
            return counter > other.counter;
        }

        return false;
    }

    bool unique() const
    {
        return counter.unique();
    }

    //! Strictly increasing counter; used to sort chunks for 'ordered' tasks
    size_t chunk_id;
    //! Use generated data; is normally not freed by this struct
    chunk_ptr data;

private:
    //! Reference counts
    std::shared_ptr<bool> counter;
};


/** Simple priority queue; needed to allow moving values out of queue. */
class chunk_queue
{
public:
    chunk_queue()
      : m_chunks()
    {
    }

    void push(data_chunk value)
    {
        m_chunks.push_back(std::move(value));
        std::push_heap(m_chunks.begin(), m_chunks.end());
    }

    data_chunk pop()
    {
        std::pop_heap(m_chunks.begin(), m_chunks.end());
        data_chunk value = std::move(m_chunks.back());
        m_chunks.pop_back();

        return value;
    }

    const data_chunk& top() const
    {
        return m_chunks.front();
    }

    bool empty() const
    {
        return m_chunks.empty();
    }

    size_t size() const
    {
        return m_chunks.size();
    }

private:
    typedef std::vector<data_chunk> chunk_vec;

    chunk_vec m_chunks;
};


struct scheduler_step
{
    scheduler_step(analytical_step* value, const std::string& name_)
      : lock()
      , ptr(value)
      , current_chunk(0)
      , last_chunk(0)
      , queue()
      , name(name_)
    {
    }

    bool can_run(size_t next_chunk)
    {
        if (ptr->get_ordering() == analytical_step::ordered) {
            return (current_chunk == next_chunk);
        }

        return true;
    }

    //! Mutex used to control access to step
    std::mutex lock;
    //! Analytical step implementation
    std::unique_ptr<analytical_step> ptr;
    //! The current chunk to be processed
    size_t current_chunk;
    //! The last chunk queued to the step;
    //! Used to correct numbering for sparse output from sequential steps
    size_t last_chunk;
    //! (Ordered) vector of chunks to be processed
    chunk_queue queue;
    //! Short name for step used for error reporting
    std::string name;

private:
    //! Not implemented
    scheduler_step(const scheduler_step&);
    //! Not implemented
    scheduler_step& operator=(const scheduler_step&);
};


scheduler::scheduler()
  : m_steps()
  , m_condition()
  , m_chunk_counter(0)
  , m_live_chunks(0)
  , m_queue_lock()
  , m_queue_calc()
  , m_queue_io()
  , m_io_active(false)
  , m_errors(false)
{
}


scheduler::~scheduler()
{
}


void scheduler::add_step(size_t step_id, const std::string& name,
                         analytical_step* step)
{
    if (m_steps.size() <= step_id) {
        m_steps.resize(step_id + 1);
    }

    AR_DEBUG_ASSERT(step);
    AR_DEBUG_ASSERT(!m_steps.at(step_id));

    m_steps.at(step_id) = step_ptr(new scheduler_step(step, name));
}


bool scheduler::run(int nthreads)
{
    AR_DEBUG_ASSERT(!m_steps.empty());
    AR_DEBUG_ASSERT(m_steps.front());
    AR_DEBUG_ASSERT(nthreads >= 1);
    AR_DEBUG_ASSERT(!m_chunk_counter);

    for (size_t task = 3 * static_cast<size_t>(nthreads); task; --task) {
        m_steps.front()->queue.push(data_chunk(m_chunk_counter++));
    }

    queue_analytical_step(m_steps.front(), 0);

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

    for (auto& thread: threads) {
        try {
            thread.join();
        } catch (const std::system_error& error) {
            std::cerr << "ERROR: Failed to join thread: " << error.what() << std::endl;
            set_errors_occured();
        }
    }

    if (errors_occured()) {
        return false;
    }

    for (auto& step: m_steps) {
        if (step && !step->queue.empty()) {
            print_locker lock;
            std::cerr << "ERROR: Not all parts run for step " << step->name
                          << "; " << step->queue.size() << " left ..." << std::endl;

            set_errors_occured();
        }
    }

    if (errors_occured()) {
        return false;
    }

    for (auto step: m_steps) {
        if (step) {
            try {
                step->ptr->finalize();
            } catch (const std::exception&) {
                std::cerr << "ERROR: Failed to finalizing task " << step->name << ":\n";
                throw;
            }
        }
    }

    return true;
}


void scheduler::run_wrapper(scheduler* sch)
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
        std::cerr << "ERROR: Unhandled, non-standard exception in thread" << std::endl;
    }

    sch->set_errors_occured();
    sch->m_condition.notify_all();
}


void scheduler::do_run()
{
    std::unique_lock<std::mutex> lock(m_queue_lock);

    while (!errors_occured()) {
        // Try to keep the disk busy by preferring IO chunks
        step_ptr current_step;
        if (m_io_active || m_queue_io.empty()) {
            if (!m_queue_calc.empty()) {
                current_step = m_queue_calc.front();
                m_queue_calc.pop();
            } else if (!m_live_chunks) {
                // Nothing left to do at all
                break;
            }
        } else {
            current_step = m_queue_io.front();
            m_queue_io.pop();
            m_io_active = true;
        }

        if (current_step) {
            lock.unlock();
            execute_analytical_step(current_step);
            lock.lock();
        } else {
            m_condition.wait(lock);
        }
    }

    // Signal any waiting threads
    m_condition.notify_all();
}


void scheduler::execute_analytical_step(const step_ptr& step)
{
    data_chunk chunk;

    {
        std::lock_guard<std::mutex> lock(step->lock);
        chunk = step->queue.pop();
    }

    chunk_vec chunks = step->ptr->process(chunk.data.release());

    std::lock_guard<std::mutex> lock(m_queue_lock);

    // Schedule each of the resulting blocks
    for (auto& result: chunks) {
        step_ptr& other_step = m_steps.at(result.first);
        AR_DEBUG_ASSERT(other_step != NULL);

        std::lock_guard<std::mutex> step_lock(other_step->lock);
        // Inherit reference count from source chunk
        data_chunk next_chunk(chunk, std::move(result.second));
        if (step->ptr->get_ordering() == analytical_step::ordered) {
            // Ordered steps are allowed to not return results, so the chunk
            // numbering is remembered for down-stream steps
            next_chunk.chunk_id = other_step->last_chunk++;
        }

        other_step->queue.push(std::move(next_chunk));
        queue_analytical_step(other_step, next_chunk.chunk_id);
    }

    // Unlock use of IO steps after finishing processing
    if (step->ptr->file_io()) {
        m_io_active = false;
        if (!m_queue_io.empty()) {
            m_condition.notify_all();
        }
    }

    // Reschedule current step if ordered and next chunk is available
    if (step->ptr->get_ordering() == analytical_step::ordered) {
        std::lock_guard<std::mutex> step_lock(step->lock);

        step->current_chunk++;
        if (!step->queue.empty()) {
            queue_analytical_step(step, step->queue.top().chunk_id);
        }
    }

    // End of the line for this chunk; re-schedule first step
    if (chunks.empty() && chunk.unique() && step != m_steps.front()) {
        step_ptr other_step = m_steps.front();

        std::lock_guard<std::mutex> step_lock(other_step->lock);
        other_step->queue.push(data_chunk(m_chunk_counter));

        queue_analytical_step(other_step, m_chunk_counter);

        m_chunk_counter++;
    }

    // Decrement counters before releasing lock
    chunk = data_chunk(0);
    m_live_chunks--;
}


void scheduler::queue_analytical_step(const step_ptr& step, size_t current)
{
    if (step->can_run(current)) {
        if (step->ptr->file_io()) {
            m_queue_io.push(step);
        } else {
            m_queue_calc.push(step);
        }

        m_live_chunks++;
        m_condition.notify_one();
    }
}

} // namespace ar
