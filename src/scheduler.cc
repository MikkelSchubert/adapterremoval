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

#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <unistd.h>
#include <cstdlib>

#include "scheduler.h"
#include "strutils.h"


///////////////////////////////////////////////////////////////////////////////
// exceptions

thread_abort::thread_abort()
  : thread_error("abort thread")
{
}


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
    data_chunk(unsigned chunk_id_ = 0,
               analytical_chunk* data_ = NULL)
      : chunk_id(chunk_id_)
      , data(data_)
      , nrefs(NULL)
    {
        increment_refs();
    }

    data_chunk(const data_chunk& parent, analytical_chunk* data = NULL)
      : chunk_id(parent.chunk_id)
      , data(data ? data : parent.data)
      , nrefs(parent.nrefs)
    {
        increment_refs();
    }

    /** Destructor; simply decrements reference count **/
    ~data_chunk()
    {
        decrement_refs();
    }

    data_chunk& operator=(const data_chunk& other)
    {
        other.increment_refs();
        decrement_refs();

        nrefs = other.nrefs;
        chunk_id = other.chunk_id;
        data = other.data;

        return *this;
    }

    /** Sorts by counter, data, type, in that order. **/
    bool operator<(const data_chunk& other) const
    {
        if (chunk_id > other.chunk_id) {
            return true;
        } else if (chunk_id == other.chunk_id) {
            if (data > other.data) {
                return true;
            }
        }

        return false;
    }

    bool unique()
    {
        return nrefs->current() == 1;
    }

    //! Strictly increasing counter; used to sort chunks for 'ordered' tasks
    unsigned chunk_id;
    //! Use generated data; is not freed by this struct
    analytical_chunk* data;

private:
    void increment_refs() const
    {
        if (nrefs) {
            nrefs->increment();
        } else {
            nrefs = new atomic_counter(1);
        }
    }

    void decrement_refs() const
    {
        if (nrefs->decrement() == 0) {
            delete nrefs;
            nrefs = NULL;
        }
    }

    //! Reference counts
    mutable atomic_counter* nrefs;
};


struct scheduler_step
{
    scheduler_step(analytical_step* value)
      : lock()
      , ptr(value)
      , current_chunk(0)
      , queue()
    {
    }

    /** Deletes any remaining chunks. */
    ~scheduler_step() {
        reset();
    }

    /** Cleans up after previous runs; deleting any remaining chunks. */
    void reset() {
        while (!queue.empty()) {
            delete queue.top().data;
            queue.pop();
        }
    }

    bool can_run(size_t next_chunk)
    {
        if (ptr->get_ordering() == analytical_step::ordered) {
            return (current_chunk == next_chunk);
        }

        return true;
    }


    //! Mutex used to control access to step
    mutex lock;
    //! Analytical step implementation
    std::auto_ptr<analytical_step> ptr;
    //! The current chunk to be processed
    unsigned current_chunk;
    //! (Ordered) vector of chunks to be processed
    std::priority_queue<data_chunk> queue;

private:
    //! Not implemented
    scheduler_step(const scheduler_step&);
    //! Not implemented
    scheduler_step& operator=(const scheduler_step&);
};


/** Simple structure used to pass parameters to threads. */
struct thread_info
{
    thread_info(unsigned seed_, scheduler* sch_)
      : seed(seed_)
      , sch(sch_)
    {
    }

    //! Per thread seed
    unsigned seed;
    //! Pointer to current scheduler
    scheduler* sch;
};


scheduler::scheduler()
  : m_steps()
  , m_running()
  , m_errors(false)
  , m_condition()
  , m_chunk_counter(0)
#ifdef AR_PTHREAD_SUPPORT
  , m_threads()
#endif
  , m_queue_lock()
  , m_queue_calc()
  , m_queue_io()
  , m_io_active(false)
  , m_live_chunks(0)
{
}


scheduler::~scheduler()
{
    mutex_locker lock(m_running);

    for (pipeline::iterator it = m_steps.begin(); it != m_steps.end(); ++it) {
        delete *it;
    }

    m_steps.clear();
}


void scheduler::add_step(size_t step_id, analytical_step* step)
{
    mutex_locker lock(m_running);
    if (!step) {
        throw std::invalid_argument("scheduler::scheduler: NULL not allowed");
    } else if (m_steps.size() <= step_id) {
        m_steps.resize(step_id + 1);
    } else if (m_steps.at(step_id)) {
        throw std::invalid_argument("scheduler::scheduler: id used multiple times");
    }

    m_steps.at(step_id) = new scheduler_step(step);
}



bool scheduler::run(int nthreads, unsigned seed)
{
    if (m_steps.empty()) {
        throw std::invalid_argument("pipeline must contain at least one step");
    } else if (!m_steps.at(0)) {
        throw std::invalid_argument("first step has not been specified");
    } else if (nthreads <= 0) {
        throw std::invalid_argument("scheduler::run: nthreads <= 0");
    }

    mutex_locker lock(m_running);

    m_chunk_counter = 0;
    for (pipeline::iterator it = m_steps.begin(); it != m_steps.end(); ++it) {
        if (*it) {
            (*it)->reset();
        }
    }

    for (unsigned task = 3 * nthreads; task; --task) {
        m_steps.front()->queue.push(data_chunk(m_chunk_counter++));
    }

    queue_analytical_step(m_steps.front(), 0);

    m_io_active = false;
    m_errors = !initialize_threads(nthreads - 1, seed + 1);

    // Signal for threads to start, or terminate in case of errors
    signal_threads();

    thread_info* info = new thread_info(seed, this);
    m_errors = !run_wrapper(info) || m_errors;
    m_errors = !join_threads() || m_errors;

    if (!m_errors) {
        for (pipeline::iterator it = m_steps.begin(); it != m_steps.end(); ++it) {
            if (*it) {
                (*it)->ptr->finalize();
            }
        }

        for (pipeline::iterator it = m_steps.begin(); it != m_steps.end(); ++it) {
            if (*it && !(*it)->queue.empty()) {
                print_locker lock;
                std::cerr << "ERROR: Not all parts run for step " << it - m_steps.begin()
                          << "; " << (*it)->queue.size() << " left ..." << std::endl;
                m_errors = true;
            }
        }
    }

    return !m_errors;
}


void* scheduler::run_wrapper(void* ptr)
{
    std::auto_ptr<thread_info> info(reinterpret_cast<thread_info*>(ptr));
    scheduler* sch = info->sch;

    // Set seed for RNG; rand is used in collapse_paired_ended_sequences()
    srandom(info->seed);

    try {
        return sch->do_run();
    } catch (const thread_abort&) {
        // Error messaging is assumed to have been done by thrower
    } catch (const std::exception& error) {
        print_locker lock;
        std::cerr << "Error in thread:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;
    } catch (...) {
        print_locker lock;
        std::cerr << "Unhandled exception in thread" << std::endl;
    }

    sch->m_errors = true;
    sch->signal_threads();

    return reinterpret_cast<void*>(false);
}


void* scheduler::do_run()
{
    // Wait to allow early termination in case of errors during setup
    m_condition.wait();

    while (!m_errors) {
        scheduler_step* current_step = NULL;

        {
            mutex_locker lock(m_queue_lock);

            // Try to keep the disk busy by preferring IO chunks
            if (m_io_active || m_queue_io.empty()) {
                if (!m_queue_calc.empty()) {
                    current_step = m_queue_calc.front();
                    m_queue_calc.pop_front();
                } else if (!m_live_chunks) {
                    // Nothing left to do at all
                    break;
                }
            } else {
                current_step = m_queue_io.front();
                m_queue_io.pop_front();
                m_io_active = true;
            }
        }

        if (current_step) {
            execute_analytical_step(current_step);
        } else {
            // Nothing to do yet ...
            m_condition.wait();
        }
    }

    // Signal any waiting threads
    m_condition.signal();

    return reinterpret_cast<void*>(true);
}


void scheduler::execute_analytical_step(scheduler_step* step)
{
    data_chunk chunk;

    {
        mutex_locker lock(step->lock);
        chunk = step->queue.top();
        step->queue.pop();
    }

    chunk_list chunks = step->ptr->process(chunk.data);

    // Unlock use of IO steps immediately after finishing processing
    if (step->ptr->file_io()) {
        mutex_locker lock(m_queue_lock);
        m_io_active = false;
        if (!m_queue_io.empty()) {
            m_condition.signal();
        }
    }

    // Schedule each of the resulting blocks
    for (chunk_list::iterator it = chunks.begin(); it != chunks.end(); ++it) {
        scheduler_step* other_step = m_steps.at(it->first);

        mutex_locker lock(other_step->lock);
        // Inherit reference count from source chunk
        other_step->queue.push(data_chunk(chunk, it->second));

        queue_analytical_step(other_step, chunk.chunk_id);
    }

    // Reschedule current step if ordered and next chunk is available
    {
        mutex_locker lock(step->lock);
        if (step->ptr->get_ordering() == analytical_step::ordered) {
            step->current_chunk++;

            if (!step->queue.empty()) {
                queue_analytical_step(step, step->queue.top().chunk_id);
            }
        }
    }

    // End of the line for this chunk; re-schedule first step
    if (chunks.empty() && chunk.unique() && step != m_steps.front()) {
        scheduler_step* other_step = m_steps.front();

        mutex_locker lock(other_step->lock);
        other_step->queue.push(data_chunk(m_chunk_counter));

        queue_analytical_step(other_step, m_chunk_counter);

        m_chunk_counter++;
    }

    // Counter is decremented last, so that threads do not exit while new
    // parts are being scheduled.
    {
        mutex_locker lock(m_queue_lock);
        m_live_chunks--;
    }
}


void scheduler::queue_analytical_step(scheduler_step* step, size_t current)
{
    if (step->can_run(current)) {
        mutex_locker lock(m_queue_lock);
        if (step->ptr->file_io()) {
            m_queue_io.push_back(step);
        } else {
            m_queue_calc.push_back(step);
        }

        m_live_chunks++;
        m_condition.signal();
    }
}


bool scheduler::initialize_threads(int nthreads, unsigned seed)
{
#ifdef AR_PTHREAD_SUPPORT
    if (!m_threads.empty()) {
        throw std::invalid_argument("scheduler::initialize_threads: threads must be empty");
    } else if (nthreads < 0) {
        throw std::invalid_argument("scheduler::initialize_threads: negative thread count");
    }

    try {
        for (int i = 0; i < nthreads; ++i) {
            m_threads.push_back(pthread_t());
            // Each thread is assigned a unique seed, based on the (user) seed
            thread_info* info = new thread_info(seed + i, this);
            switch (pthread_create(&m_threads.back(), NULL, &run_wrapper, info)) {
                case 0:
                    break;

                case EAGAIN:
                    throw thread_error("pthread_create: insufficient resources to create thread");

                case EINVAL:
                    throw thread_error("pthread_create: invalid attributes");

                case EPERM:
                    throw thread_error("pthread_create: insufficient permissions");

                default:
                    throw thread_error("pthread_create: unknown error");
            }
        }
    } catch (const thread_error& error) {
        print_locker lock;
        std::cerr << "Error creating threads:\n"
                  << cli_formatter::fmt(error.what()) << std::endl;

        m_threads.pop_back();
        return false;
    }
#else
    (void)nthreads;
    (void)seed;
#endif
    return true;
}


void scheduler::signal_threads()
{
#ifdef AR_PTHREAD_SUPPORT
    // Signal the main and all other threads
    for (size_t i = 0; i < m_threads.size() + 1; ++i) {
        m_condition.signal();
    }
#endif
}


bool scheduler::join_threads()
{
    bool join_result = true;

#ifdef AR_PTHREAD_SUPPORT
    for (thread_vector::iterator it = m_threads.begin(); it != m_threads.end(); ++it) {
        void* run_result = NULL;
        const int join_error = pthread_join(*it, &run_result);

        if (join_error) {
            join_result = false;

            print_locker lock;
            switch (join_error) {
                case EINVAL:
                    std::cerr << "Error in pthread_join: invalid thread" << std::endl;
                    break;

                case ESRCH:
                    std::cerr << "Error in pthread_join: thread not joinable" << std::endl;
                    break;

                case EDEADLK:
                    std::cerr << "Error in pthread_join: deadlock detected" << std::endl;
                    break;

                default:
                    std::cerr << "Error in pthread_join: unknown error: " << join_error << std::endl;
            }
        } else {
            join_result &= static_cast<bool>(run_result);
        }
    }

    m_threads.clear();
#endif

    return join_result;
}


