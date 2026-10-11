// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "commontypes.hpp" // for read_mate
#include "errors.hpp"      // for assert_failed
#include "scheduler.hpp"   // for step_id
#include "testing.hpp"     // for CHECK, CHECK_FALSE, TEST_CASE
#include <functional>      // for function
#include <memory>          // for make_unique
#include <sstream>         // for ostringstream
#include <utility>         // for move

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// helper functions and classes

step_id
step_id_test_helper(size_t i)
{
  return step_id{ i };
}

namespace {

using process_function =
  std::function<void(analytical_step&, chunk_vec&, chunk_ptr ptr)>;
using signal_function = std::function<void(analytical_step&, chunk_vec&)>;

class test_step : public analytical_step
{
public:
  test_step(processing_order order,
            std::string name,
            step_id next_step,
            process_function process = {},
            signal_function signal = {})
    : analytical_step(order, std::move(name), next_step)
    , m_process(std::move(process))
    , m_signal(std::move(signal))
  {
  }

  test_step(processing_order order,
            std::string name,
            step_id next_step,
            signal_function signal)
    : test_step(order, std::move(name), next_step, nullptr, std::move(signal))
  {
  }

  chunk_vec process(chunk_ptr data) override
  {
    AR_REQUIRE(!!m_process);
    chunk_vec vec;
    m_process(*this, vec, std::move(data));
    return vec;
  }

  chunk_vec signal() override
  {
    AR_REQUIRE(!!m_signal);
    chunk_vec vec;
    m_signal(*this, vec);
    return vec;
  }

private:
  process_function m_process{};
  signal_function m_signal{};
};

} // namespace

////////////////////////////////////////////////////////////////////////////////
// step_id

TEST_CASE("step_id equality", "[step_id]")
{
  CHECK(step_id::UNDEFINED == step_id::UNDEFINED);
  CHECK_FALSE(step_id::UNDEFINED != step_id::UNDEFINED);

  CHECK_FALSE(step_id_test_helper(1) == step_id::UNDEFINED);
  CHECK(step_id_test_helper(1) != step_id::UNDEFINED);

  CHECK_FALSE(step_id::UNDEFINED == step_id_test_helper(1));
  CHECK(step_id::UNDEFINED != step_id_test_helper(1));

  CHECK_FALSE(step_id_test_helper(2) == step_id_test_helper(1));
  CHECK(step_id_test_helper(2) != step_id_test_helper(1));

  const auto s = step_id_test_helper(1);
  CHECK(s == s);
  CHECK_FALSE(s != s);
}

TEST_CASE("step_id get", "[step_id]")
{
  CHECK(step_id_test_helper(0).get() == 0);
  CHECK(step_id_test_helper(12345).get() == 12345);
  CHECK(step_id_test_helper(9223372036854775808LLU).get() ==
        9223372036854775808LLU);
  CHECK(step_id::UNDEFINED.get() == static_cast<size_t>(-1));
}

TEST_CASE("step_id id to string", "[step_id]")
{
  std::ostringstream os;
  os << step_id_test_helper(12345);

  CHECK(os.str() == "step_id{12345}");
}

TEST_CASE("step_id undefined to string", "[step_id]")
{
  std::ostringstream os;
  os << step_id::UNDEFINED;

  CHECK(os.str() == "step_id::UNDEFINED");
}

////////////////////////////////////////////////////////////////////////////////
// scheduler: signalling

TEST_CASE("initial step must not be unordered", "[scheduler]")
{
  scheduler sch;
  sch.add<test_step>(processing_order::unordered,
                     "test_step",
                     step_id::UNDEFINED);
  REQUIRE_THROWS_CONTAINS(sch.run(1),
                          assert_failed,
                          "source step must be ordered");
}

TEST_CASE("scheduler signals source until suspended", "[scheduler]")
{
  scheduler sch;

  size_t signals_received = 0;
  sch.add<test_step>(processing_order::ordered,
                     "test",
                     step_id::UNDEFINED,
                     [&](analytical_step&, chunk_vec& vec) {
                       if (++signals_received >= 3) {
                         vec.emplace_back(scheduler_signal::suspend);
                       }
                     });

  REQUIRE(sch.run(1));
  REQUIRE(signals_received == 3);
}

TEST_CASE("signals from downstream to source are filtered if suspended",
          "[scheduler]")
{
  scheduler sch;
  size_t chunks_received = 0;
  const auto step_2 =
    sch.add<test_step>(processing_order::ordered,
                       "dst",
                       step_id::UNDEFINED,
                       [&](analytical_step& step, chunk_vec& vec, chunk_ptr) {
                         chunks_received++;
                         vec.emplace_back(step.signal_target());
                       });

  size_t signals_received = 0;
  const auto step_1 =
    sch.add<test_step>(processing_order::ordered,
                       "src",
                       step_2,
                       [&](analytical_step&, chunk_vec& vec) {
                         signals_received++;
                         vec.emplace_back(std::make_unique<fastq_chunk>());
                         vec.emplace_back(scheduler_signal::suspend);
                       });

  sch.set_signal_target(step_2, step_1);

  REQUIRE(sch.run(1));
  REQUIRE(signals_received == 1);
  REQUIRE(chunks_received == 1);
}

TEST_CASE("signals from downstream to non-source", "[scheduler]")
{
  scheduler sch;
  size_t chunks_received_3 = 0;
  const auto step_3 =
    sch.add<test_step>(processing_order::ordered,
                       "dst",
                       step_id::UNDEFINED,
                       [&](analytical_step& step, chunk_vec& vec, chunk_ptr) {
                         if (++chunks_received_3 <= 5) {
                           vec.emplace_back(step.signal_target());
                         }
                       });

  size_t chunks_received_2 = 0;
  size_t signals_received_2 = 0;
  const auto step_2 = sch.add<test_step>(
    processing_order::ordered,
    "dst",
    step_3,
    [&](analytical_step&, chunk_vec& vec, chunk_ptr chunk) {
      chunks_received_2++;
      vec.emplace_back(std::move(chunk));
    },
    [&](analytical_step&, chunk_vec&) { ++signals_received_2; });

  size_t signals_received_1 = 0;
  sch.add<test_step>(processing_order::ordered,
                     "src",
                     step_2,
                     [&](analytical_step&, chunk_vec& vec) {
                       if (++signals_received_1 == 10) {
                         vec.emplace_back(scheduler_signal::suspend);
                       }
                       vec.emplace_back(std::make_unique<fastq_chunk>());
                     });

  sch.set_signal_target(step_3, step_2);

  REQUIRE(sch.run(1));
  REQUIRE(signals_received_1 == 10);
  REQUIRE(signals_received_2 == 5);
  REQUIRE(chunks_received_2 == 10);
  REQUIRE(chunks_received_3 == 10);
}

TEST_CASE("signals from downstream to non-source are filtered if suspended",
          "[scheduler]")
{
  scheduler sch;
  size_t chunks_received_3 = 0;
  const auto step_3 =
    sch.add<test_step>(processing_order::ordered,
                       "dst",
                       step_id::UNDEFINED,
                       [&](analytical_step& step, chunk_vec& vec, chunk_ptr) {
                         chunks_received_3++;
                         vec.emplace_back(step.signal_target());
                       });

  size_t chunks_received_2 = 0;
  size_t signals_received_2 = 0;
  const auto step_2 = sch.add<test_step>(
    processing_order::ordered,
    "dst",
    step_3,
    [&](analytical_step&, chunk_vec& vec, chunk_ptr chunk) {
      chunks_received_2++;
      vec.emplace_back(std::move(chunk));
    },
    [&](analytical_step&, chunk_vec& vec) {
      if (++signals_received_2 == 5) {
        vec.emplace_back(scheduler_signal::suspend);
      }
    });

  size_t signals_received_1 = 0;
  sch.add<test_step>(processing_order::ordered,
                     "src",
                     step_2,
                     [&](analytical_step&, chunk_vec& vec) {
                       if (++signals_received_1 == 10) {
                         vec.emplace_back(scheduler_signal::suspend);
                       }
                       vec.emplace_back(std::make_unique<fastq_chunk>());
                     });

  sch.set_signal_target(step_3, step_2);

  REQUIRE(sch.run(1));
  REQUIRE(signals_received_1 == 10);
  REQUIRE(signals_received_2 == 5);
  REQUIRE(chunks_received_2 == 10);
  REQUIRE(chunks_received_3 == 10);
}

TEST_CASE("signal source and recipient must be ordered", "[scheduler]")
{
  scheduler sch;
  auto u2 = sch.add<analytical_step>(processing_order::unordered,
                                     "unordered step 2",
                                     step_id::UNDEFINED);
  auto o2 = sch.add<analytical_step>(processing_order::ordered,
                                     "ordered step 2",
                                     step_id::UNDEFINED);
  auto u1 = sch.add<analytical_step>(processing_order::unordered,
                                     "unordered step 1",
                                     step_id::UNDEFINED);
  auto o1 = sch.add<analytical_step>(processing_order::ordered,
                                     "ordered step 1",
                                     step_id::UNDEFINED);

  auto mate = GENERATE(read_mate::_1, read_mate::_2);
  CHECK_NOTHROW(sch.set_signal_target(o1, o2, mate));
  CHECK_NOTHROW(sch.set_signal_target(o2, o1, mate));
  CHECK_THROWS_AS(sch.set_signal_target(u1, u2, mate), assert_failed);
  CHECK_THROWS_AS(sch.set_signal_target(u1, o2, mate), assert_failed);
  CHECK_THROWS_AS(sch.set_signal_target(o1, u1, mate), assert_failed);
}

} // namespace adapterremoval
