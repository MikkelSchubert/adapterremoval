// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2015 Mikkel Schubert <mikkelsch@gmail.com>
#pragma once

#include "buffer.hpp"         // for buffer
#include "commontypes.hpp"    // for read_mate
#include "fastq.hpp"          // for fastq
#include <condition_variable> // for condition_variable
#include <cstddef>            // for size_t
#include <cstdint>            // for uint32_t
#include <iosfwd>             // for ostream
#include <memory>             // for unique_ptr, shared_ptr, make_unique
#include <mutex>              // for mutex
#include <queue>              // for queue
#include <string>             // for string
#include <type_traits>        // for is_base_of
#include <utility>            // for pair
#include <variant>            // for variant
#include <vector>             // for vector

namespace adapterremoval {

class buffer;
class fastq;
class scheduler_step;
enum class threadtype;

/** Identifier returned by scheduler when adding a step */
class step_id
{
public:
  /** Special value representing unset/undefined steps */
  static const step_id UNDEFINED;

  /** Returns the underlying value representing the step */
  [[nodiscard]] constexpr size_t get() const { return m_id; }

  /** Equality operator */
  [[nodiscard]] constexpr bool operator==(const step_id& other) const
  {
    return this->m_id == other.m_id;
  }

  /** Inequality operator */
  [[nodiscard]] constexpr bool operator!=(const step_id& other) const
  {
    return !(*this == other);
  }

private:
  /** Constructs step from a numerical ID; only done by `scheduler` */
  explicit constexpr step_id(size_t id)
    : m_id(id)
  {
  }

  //! Numerical ID of existing step, or -1 if undefined
  size_t m_id;

  friend class scheduler;
  friend step_id step_id_test_helper(size_t);
};

inline constexpr step_id step_id::UNDEFINED{ static_cast<size_t>(-1) };

/** Converts step_id to string for debug logging */
std::ostream&
operator<<(std::ostream&, step_id);

/** Signal used to control step and scheduler behavior */
enum class scheduler_signal
{
  //! The step emitting this signal will no longer accept signals going forward
  suspend,
  //! The step emitting this signal will accept signals going forward
  resume,
};

/**
 * Base-class for data-chunks produced, processed and consumed by a pipeline.
 */
class analytical_chunk
{
public:
  /** Constructor; does nothing. */
  analytical_chunk() = default;

  /** Destructor; typically does nothing */
  virtual ~analytical_chunk() = 0;

  analytical_chunk(const analytical_chunk&) = delete;
  analytical_chunk(analytical_chunk&&) = delete;
  analytical_chunk& operator=(const analytical_chunk&) = delete;
  analytical_chunk& operator=(analytical_chunk&&) = delete;
};

/** Chunk containing unprocessed or partially processed reads */
class fastq_chunk : public analytical_chunk
{
public:
  fastq_chunk() = default;
  ~fastq_chunk() override = default;

  //! Indicates that this is the first chunk
  bool first = false;
  //! Indicates that this is the last chunk/that EOF has been reached
  bool eof = false;

  //! Original mate separator for (paired) reads
  char mate_separator = '\0';

  //! FASTQ records read from the mate 1 files
  std::vector<fastq> reads_1{};
  //! FASTQ records read from the mate 2 files
  std::vector<fastq> reads_2{};

  //! Barcode identified for the corresponding read (pair)
  std::vector<uint32_t> barcodes{};

  /** Total number of reads in the chunk */
  [[nodiscard]] size_t reads() const noexcept
  {
    return this->reads_1.size() + this->reads_2.size();
  }

  fastq_chunk(const fastq_chunk&) = delete;
  fastq_chunk(fastq_chunk&&) = delete;
  fastq_chunk& operator=(const fastq_chunk&) = delete;
  fastq_chunk& operator=(fastq_chunk&&) = delete;
};

/** Chunk containing processed reads to (maybe) be compressed and streamed */
class output_chunk : public analytical_chunk
{
public:
  output_chunk() = default;

  ~output_chunk() override = default;

  //! Indicates that this is the last chunk/that EOF has been reached
  bool eof = false;

  //! Buffers of (compressed) FASTQ reads
  std::vector<buffer> buffers{};

  //! CRC32 of (uncompressed) data; only set if EOF is true
  uint32_t crc32 = 0;

  //! Size of (uncompressed) data in buffers;
  size_t uncompressed_size = 0;

  output_chunk(const output_chunk&) = delete;
  output_chunk(output_chunk&&) = delete;
  output_chunk& operator=(const output_chunk&) = delete;
  output_chunk& operator=(output_chunk&&) = delete;
};

using chunk_ptr = std::unique_ptr<analytical_chunk>;
using chunk_pair = std::pair<step_id, chunk_ptr>;
using chunk_variant =
  std::variant<chunk_ptr, chunk_pair, step_id, scheduler_signal>;
using chunk_vec = std::vector<chunk_variant>;

using fastq_chunk_ptr = std::unique_ptr<fastq_chunk>;
using fastq_chunk_ptr_vec = std::vector<fastq_chunk_ptr>;

using output_chunk_ptr = std::unique_ptr<output_chunk>;
using output_chunk_ptr_vec = std::vector<output_chunk_ptr>;

/** Ordering of input for analytical steps. */
enum class processing_order
{
  //! Data must be consumed in the input order. May involve IO
  ordered,
  //! Data must be consumed in the input order and only involves IO
  ordered_io,
  //! Data may be consumed in any order. May not involve IO
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
   * @param name Human readable name of the step, used for debug messages.
   * @param next_step If specified, this step is used as the recipient for
   *                  output from `process`. Otherwise, a recipient must be
   *                  specified for each output from `process()`.
   */
  analytical_step(processing_order step_order,
                  std::string name,
                  step_id next_step);

  /** Destructor; does nothing in base class. **/
  virtual ~analytical_step() = default;

  /**
   * Function used to signal to a step that it should generate more data; for a
   * reader this means reading data from disk, while a processing step might in
   * turn signal a prior step to feed it more data. Only ordered steps are
   * allowed to send and/or receive signals.
   *
   * Ordered steps may produce any number of outputs per input, while unordered
   * steps must produce one output per input, or no output whatsoever. A ordered
   * step may receive (direct) input from multiple ordered steps, in which case
   * results are received interleaved from each source, in the order in which it
   * is produced by those direct sources. No other combination involving inputs
   * from multiple sources are allowed.
   *
   * The scheduler will signal the first step in the pipeline, i.e. the last
   * step added. The pipeline terminates once this step has returned `suspend`
   * and no pending tasks remain.
   *
   * This function may return one of several variants:
   *  - A chunk_ptr, which will be sent to the step returned by `next_step()`
   *  - A step_id recipient and the chunk_ptr it should receive
   *  - A step_id which should receive a `signal()` call
   *  - A scheduler_signal turning on/off signals for this step
   *
   * As a rule of thumb, a signal should only be sent if another step has
   * signalled this step to emit data, but this step is unable to produce more
   * output. A simple implementation processes data in `process` and then call
   * `signal()` to return it (if enough was processed) or to request more data
   * from upstream (if not).
   */
  virtual chunk_vec signal();

  /**
   * Function called by pipeline to process data chunks. This function differs
   * only from `signal()` in that it receives a chunk produced from a prior step
   */
  virtual chunk_vec process(chunk_ptr data);

  /**
   * Called once the pipeline has been run to completion; this function is
   * called on analytical steps in reverse order of addition to the pipeline.
   */
  virtual void finalize() {}

  /** Returns the ID of the step to receive the output from `process()` */
  [[nodiscard]] step_id next_step() const { return m_next_step; }

  [[nodiscard]] step_id signal_target(read_mate read = read_mate::_1) const
  {
    return read == read_mate::_1 ? m_signal_targets.first
                                 : m_signal_targets.second;
  }

  /** Returns the expected ordering (ordered / unordered) for input data. **/
  [[nodiscard]] processing_order ordering() const { return m_step_order; }

  /** Returns the name of the analytical step (type) */
  [[nodiscard]] const std::string& name() const { return m_name; }

  analytical_step(const analytical_step&) = delete;
  analytical_step(analytical_step&&) = delete;
  analytical_step& operator=(const analytical_step&) = delete;
  analytical_step& operator=(analytical_step&&) = delete;

private:
  //! Stores the ordering of data chunks expected by the step
  const processing_order m_step_order;
  //! Human readable name for step; for debugging purposes
  const std::string m_name;
  //! The step that should receive output from this step
  const step_id m_next_step;
  //! Two possible targets (mate 1 and 2) for signals sent by this step
  std::pair<step_id, step_id> m_signal_targets = {
    step_id::UNDEFINED,
    step_id::UNDEFINED,
  };

  friend class scheduler_step;
};

/**
 * Multithreaded scheduler.
 *
 * See 'analytical_step' for information on implementing analyzes.
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
   * This ID may be used to specify the recipient of output from other steps.
   **/
  template<typename T, typename... Args>
  step_id add(Args&&... args)
  {
    static_assert(std::is_base_of<analytical_step, T>(),
                  "requires analytical_step sub-class");
    return add_step(std::make_unique<T>(std::forward<Args>(args)...));
  }

  /** Set the signalling target for a step */
  void set_signal_target(step_id src,
                         step_id dst,
                         read_mate read = read_mate::_1);

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

  step_id add_step(std::unique_ptr<analytical_step> step);

  /** Wrapper function which calls do_run on the provided thread. */
  static void run_wrapper(scheduler* /* sch */, threadtype /* thread_type */);
  /** Work function for calculation heavy threads; invoked by each thread. */
  void run_calc_loop();
  /** Work function for IO threads; invoked by each thread. */
  void run_io_loop();

  /** Returns true if a signal can be queued for the source step */
  [[nodiscard]] bool can_signal_source() const;

  //! Analytical steps
  pipeline m_steps{};

  //! Condition used to signal the (potential) availability of IO work
  std::condition_variable m_condition_io{};
  //! Condition used to signal the (potential) availability of calculation work
  std::condition_variable m_condition_calc{};

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
};

} // namespace adapterremoval
