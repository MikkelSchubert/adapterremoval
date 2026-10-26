// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "buffer.hpp"    // for buffer
#include "errors.hpp"    // for assertion_failed
#include "logging.hpp"   // for log_capture
#include "strutils.hpp"  // for log_escape
#include "testing.hpp"   // for TEST_CASE, REQUIRE, ...
#include "threading.hpp" // for threadsafe_data
#include <atomic>        // for atomic_int
#include <chrono>        // for milliseconds
#include <cstddef>       // for size_t
#include <cstdint>       // for uint8_t
#include <memory>        // for unique_ptr
#include <string>        // for string
#include <thread>        // for std::thread
#include <utility>       // for move, ignore

namespace Catch {

template<>
std::string
fallbackStringifier(const std::unique_ptr<int>& value)
{
  ReusableStringStream ss;

  if (value) {
    ss << "unique_ptr({" << *value << "})";
  } else {
    ss << "unique_ptr(null)";
  }

  return ss.str();
}

template<>
std::string
fallbackStringifier(const std::unique_ptr<std::string>& value)
{
  ReusableStringStream ss;

  if (value) {
    ss << "unique_ptr({" << adapterremoval::log_escape(*value) << "})";
  } else {
    ss << "unique_ptr(null)";
  }

  return ss.str();
}

} // namespace Catch

namespace adapterremoval {

////////////////////////////////////////////////////////////////////////////////
// fallbackStringifier

TEST_CASE("unique_ptr fallback stringifier")
{
  REQUIRE(Catch::fallbackStringifier(std::make_unique<int>(123)) ==
          "unique_ptr({123})");
  REQUIRE(Catch::fallbackStringifier(std::unique_ptr<int>{}) ==
          "unique_ptr(null)");

  REQUIRE(Catch::fallbackStringifier(std::make_unique<std::string>("123")) ==
          "unique_ptr({'123'})");
  REQUIRE(Catch::fallbackStringifier(std::unique_ptr<std::string>{}) ==
          "unique_ptr(null)");
}

////////////////////////////////////////////////////////////////////////////////
// threadstate

TEST_CASE("empty threadstate", "[scheduler::threadstate]")
{
  threadstate<int> state;

  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
  REQUIRE(state.try_acquire() == threadstate<int>::pointer());
}

TEST_CASE("threadstate emplace back", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back(1234);

  SECTION("acquire")
  {
    auto value = state.acquire();
    REQUIRE(value);
    REQUIRE(*value == 1234);
    REQUIRE_THROWS_AS(state.acquire(), assert_failed);
  }

  SECTION("try acquire")
  {
    auto value = state.try_acquire();
    REQUIRE(value);
    REQUIRE(*value == 1234);
    REQUIRE(state.try_acquire() == threadstate<int>::pointer());
  }
}

TEST_CASE("threadstate emplace back n", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back_n(3, 1234);

  for (size_t i = 0; i < 3; ++i) {
    auto value = state.acquire();
    REQUIRE(value);
    REQUIRE(*value == 1234);
  }

  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate emplace back n (moved)", "[scheduler::threadstate]")
{
  threadstate<std::string> state;
  state.emplace_back_n(3, std::string{ __PRETTY_FUNCTION__ });

  for (size_t i = 0; i < 3; ++i) {
    auto value = state.acquire();
    REQUIRE(value);
    REQUIRE(*value == __PRETTY_FUNCTION__);
  }

  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate release", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back(1234);

  auto value = state.acquire();
  REQUIRE(value);
  REQUIRE(*value == 1234);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);

  state.release(std::move(value));

  value = state.acquire();
  REQUIRE(value);
  REQUIRE(*value == 1234);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate release rvalue", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.release(std::make_unique<int>(1234));

  auto value = state.acquire();
  REQUIRE(value);
  REQUIRE(*value == 1234);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate merge into", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.emplace_back(1);
  state.emplace_back(3);
  state.emplace_back(5);

  int destination = 10;
  state.merge_into(destination);
  REQUIRE(destination == 19);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

////////////////////////////////////////////////////////////////////////////////
// threadsafe_data

TEST_CASE("threadsafe_data default construction")
{
  threadsafe_data<std::string> data;

  {
    auto reader = data.get_reader();
    REQUIRE(reader->empty());
  }

  {
    auto writer = data.get_writer();
    REQUIRE(writer->empty());
  }
}

TEST_CASE("threadsafe_data with single argument")
{
  threadsafe_data<std::string> data{ "test data" };

  {
    auto reader = data.get_reader();
    REQUIRE(*reader == "test data");
  }

  {
    auto writer = data.get_writer();
    REQUIRE(*writer == "test data");
  }
}

TEST_CASE("threadsafe_data with multiple arguments")
{
  threadsafe_data<std::string> data{ 13, 'O' };

  {
    auto reader = data.get_reader();
    REQUIRE(*reader == "OOOOOOOOOOOOO");
  }

  {
    auto writer = data.get_writer();
    REQUIRE(*writer == "OOOOOOOOOOOOO");
  }
}

TEST_CASE("threadsafe_data are smart pointers")
{
  threadsafe_data<std::string> data_1{ "example" };
  threadsafe_data<std::string> data_2{ data_1 };
  threadsafe_data<std::string> data_3{ "other" };
  data_3 = data_2;

  {
    auto reader_1 = data_1.get_reader();
    REQUIRE(*reader_1 == "example");
    REQUIRE(reader_1->size() == 7);
    auto reader_2 = data_2.get_reader();
    REQUIRE(*reader_2 == "example");
    REQUIRE(reader_2->size() == 7);
    auto reader_3 = data_3.get_reader();
    REQUIRE(*reader_3 == "example");
    REQUIRE(reader_3->size() == 7);
  }

  *data_2.get_writer() = "modified";

  {
    auto reader_1 = data_1.get_reader();
    REQUIRE(*reader_1 == "modified");
    REQUIRE(reader_1->size() == 8);
    auto reader_2 = data_2.get_reader();
    REQUIRE(*reader_2 == "modified");
    REQUIRE(reader_2->size() == 8);
    auto reader_3 = data_3.get_reader();
    REQUIRE(*reader_3 == "modified");
    REQUIRE(reader_3->size() == 8);
  }
}

TEST_CASE("threadsafe_data are smart pointers from const")
{
  const threadsafe_data<std::string> data_1{ "example" };
  threadsafe_data<std::string> data_2{ data_1 };
  threadsafe_data<std::string> data_3{ "other" };
  data_3 = data_2;

  {
    auto reader_1 = data_1.get_reader();
    REQUIRE(*reader_1 == "example");
    REQUIRE(reader_1->size() == 7);
    auto reader_2 = data_2.get_reader();
    REQUIRE(*reader_2 == "example");
    REQUIRE(reader_2->size() == 7);
    auto reader_3 = data_3.get_reader();
    REQUIRE(*reader_3 == "example");
    REQUIRE(reader_3->size() == 7);
  }

  *data_2.get_writer() = "modified";

  {
    auto reader_1 = data_1.get_reader();
    REQUIRE(*reader_1 == "modified");
    REQUIRE(reader_1->size() == 8);
    auto reader_2 = data_2.get_reader();
    REQUIRE(*reader_2 == "modified");
    REQUIRE(reader_2->size() == 8);
    auto reader_3 = data_3.get_reader();
    REQUIRE(*reader_3 == "modified");
    REQUIRE(reader_3->size() == 8);
  }
}

TEST_CASE("threadsafe_data moves are checked")
{
  threadsafe_data<std::string> data_1{ "example" };
  const threadsafe_data<std::string> data_2{ std::move(data_1) };
  // NOLINTBEGIN(bugprone-use-after-move,hicpp-invalid-access-moved)
  threadsafe_data<std::string> data_3{ data_1 };

  REQUIRE_THROWS_AS(data_1.get_reader(), assert_failed);
  REQUIRE_THROWS_AS(data_1.get_writer(), assert_failed);

  REQUIRE_THROWS_AS(data_3.get_reader(), assert_failed);
  REQUIRE_THROWS_AS(data_3.get_writer(), assert_failed);
  // NOLINTEND(bugprone-use-after-move,hicpp-invalid-access-moved)

  REQUIRE(*data_2.get_reader() == "example");
}

TEST_CASE("threadsafe_data readers are shared")
{
  const threadsafe_data<std::string> data{ "example" };

  auto reader_1 = data.get_reader();
  auto reader_2 = data.get_reader();
  auto reader_3 = data.get_reader();

  REQUIRE(*reader_1 == "example");
  REQUIRE(*reader_2 == "example");
  REQUIRE(*reader_3 == "example");
}

TEST_CASE("threadsafe_data writers are exclusive #1")
{
  threadsafe_data<std::string> data{ "example" };
  std::atomic_int waiting = 2;

  const auto write_func = [&data, &waiting]() {
    auto writer = data.get_writer();

    for (waiting -= 1; waiting;) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    REQUIRE(*writer == "example");
    *writer = "modified";
  };

  const auto read_func = [&data, &waiting]() {
    for (waiting -= 1; waiting;) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    const auto reader = data.get_reader();
    REQUIRE(*reader == "modified");
  };

  std::thread thread_1{ write_func };
  std::thread thread_2{ read_func };
  thread_1.join();
  thread_2.join();
}

TEST_CASE("threadsafe_data writers are exclusive #2")
{
  threadsafe_data<std::string> data{ "example" };
  std::atomic_int waiting = 2;

  const auto write_func_1 = [&data, &waiting]() {
    auto writer = data.get_writer();

    for (waiting -= 1; waiting;) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    writer->push_back('!');
  };

  const auto write_func_2 = [&data, &waiting]() {
    for (waiting -= 1; waiting;) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    auto writer = data.get_writer();
    writer->push_back('?');
  };

  std::thread thread_1{ write_func_1 };
  std::thread thread_2{ write_func_2 };
  thread_1.join();
  thread_2.join();

  REQUIRE(*data.get_reader() == "example!?");
}

TEST_CASE("threadsafe_data writers are exclusive #3")
{
  threadsafe_data<std::string> data_1{ "example" };
  auto data_2 = data_1;
  std::atomic_int waiting = 2;

  const auto write_func_1 = [&data_1, &waiting]() {
    auto writer = data_1.get_writer();

    for (waiting -= 1; waiting;) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    writer->push_back('!');
  };

  const auto write_func_2 = [&data_2, &waiting]() {
    for (waiting -= 1; waiting;) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    auto writer = data_2.get_writer();
    writer->push_back('?');
  };

  std::thread thread_1{ write_func_1 };
  std::thread thread_2{ write_func_2 };
  thread_1.join();
  thread_2.join();

  REQUIRE(*data_1.get_reader() == "example!?");
  REQUIRE(*data_2.get_reader() == "example!?");
}

TEST_CASE("readers can be moved")
{
  const threadsafe_data<std::string> data{ "example" };
  auto reader_1 = data.get_reader();

  // NOLINTBEGIN(bugprone-use-after-move,hicpp-invalid-access-moved)
  const decltype(reader_1) reader_2{ std::move(reader_1) };
  REQUIRE_THROWS_AS(*reader_1, assert_failed);
  REQUIRE_THROWS_AS(reader_1->size(), assert_failed);
  // NOLINTEND(bugprone-use-after-move,hicpp-invalid-access-moved)

  REQUIRE(*reader_2 == "example");
  REQUIRE(reader_2->size() == 7);
}

TEST_CASE("writers can be moved")
{
  threadsafe_data<std::string> data{ "example" };

  auto writer_1 = data.get_writer();

  // NOLINTBEGIN(bugprone-use-after-move,hicpp-invalid-access-moved)
  decltype(writer_1) writer_2{ std::move(writer_1) };
  REQUIRE_THROWS_AS(*writer_1, assert_failed);
  REQUIRE_THROWS_AS(writer_1->size(), assert_failed);
  // NOLINTEND(bugprone-use-after-move,hicpp-invalid-access-moved)

  REQUIRE(*writer_2 == "example");
  REQUIRE(writer_2->size() == 7);
}

////////////////////////////////////////////////////////////////////////////////
// reusable_objects

TEST_CASE("reusable_objects::acquire returns previously released object",
          "[reusable_objects]")
{
  reusable_objects<buffer> rs{ "test", 64 };
  buffer buf = rs.acquire(64);
  const auto* data = buf.data();
  rs.release(std::move(buf));
  REQUIRE(buf.data() == nullptr);
  buf = rs.acquire(GENERATE(0, 1, 63, 64));
  REQUIRE(buf.data() == data);
}

TEST_CASE("reusable_objects::acquire returns object of requested size/capacity",
          "[reusable_objects]")
{
  reusable_objects<buffer> rs{ "test", 64 };

  const size_t size = GENERATE(0, 1, 63, 64);
  auto buf = rs.acquire(size);
  REQUIRE(buf.size() == size);
  REQUIRE(buf.capacity() == 64);
}

TEST_CASE("reusable_objects::acquire asserts if requested size is too big",
          "[reusable_objects]")
{
  reusable_objects<buffer> rs{ "test", 64 };

  const size_t size = GENERATE(65, 66, 1024);
  REQUIRE_THROWS_AS(rs.acquire(size), assert_failed);
}

TEST_CASE("reusable_objects::release ignores release of 0-capacity objects",
          "[reusable_objects]")
{
  reusable_objects<buffer> rs{ "test", 64 };

  // empty or non-empty storage
  if (GENERATE(true, false)) {
    auto buf = rs.acquire(64);
    rs.release(std::move(buf));
  }

  buffer buf;
  REQUIRE(buf.capacity() == 0);
  REQUIRE_NOTHROW(rs.release(std::move(buf)));
}

TEST_CASE("reusable_objects::release asserts if capacity is too low",
          "[reusable_objects]")
{
  reusable_objects<buffer> rs{ "test", 64 };

  std::ignore = rs.acquire(64);
  buffer buf{ GENERATE(1, 2, 62, 63) };
  REQUIRE_THROWS_AS(rs.release(std::move(buf)), assert_failed);
}

TEST_CASE("reusable_objects::release permits equal or higher capacity",
          "[reusable_objects]")
{
  reusable_objects<buffer> rs{ "test", 64 };

  std::ignore = rs.acquire(64);
  buffer buf{ GENERATE(64, 65, 66, 1024) };
  REQUIRE_NOTHROW(rs.release(std::move(buf)));
}

TEST_CASE("reusable_objects::release asserts if too many objects are returned",
          "[reusable_objects]")
{
  reusable_objects<buffer> rs{ "test", 64 };

  // empty or non-empty storage
  const int storage_size = GENERATE(0, 1, 2);
  for (int i = 0; i < storage_size; ++i) {
    rs.release(rs.acquire(64));
  }

  // empty objects may be returned
  buffer buf;
  REQUIRE_NOTHROW(rs.release(std::move(buf)));

  // non-empty objects may not; capacity matters, not size
  if (GENERATE(true, false)) {
    buf.reserve(64);
  } else {
    buf.resize(64);
  }

  REQUIRE_THROWS_AS(rs.release(std::move(buf)), assert_failed);
}

TEST_CASE("reusable_objects::finalize asserts if too few objects are returned",
          "[reusable_objects]")
{
  log::log_capture cap;
  reusable_objects<buffer> rs{ "test", 64 };

  const int storage_size = GENERATE(1, 2, 3);
  std::vector<buffer> storage;
  for (int i = 0; i < storage_size; ++i) {
    storage.emplace_back(rs.acquire(64));
  }

  // this results in some redundant tests, but is easier than writing them out
  const int n_returned = std::min<int>(GENERATE(0, 1, 2), storage_size - 1);
  for (int i = 0; i < n_returned; ++i) {
    rs.release(std::move(storage.at(i)));
  }

  // empty objects are not counted
  if (GENERATE(true, false)) {
    rs.release(buffer{});
  }

  REQUIRE_THROWS_AS(rs.finalize(), assert_failed);
}

TEST_CASE("reusable_objects::finalize logs utilization (defaults)",
          "[reusable_objects]")
{
  log::log_capture cap;
  reusable_objects<buffer> rs{ "test", 64 };

  REQUIRE_NOTHROW(rs.finalize());
  REQUIRE(cap.str() == "[DEBUG] Cache statistics for test: objects created = "
                       "0, objects released = 0, allocations saved = 0\n");
}

TEST_CASE("reusable_objects::finalize logs utilization (non-defaults)",
          "[reusable_objects]")
{
  log::log_capture cap;
  reusable_objects<buffer> rs{ "test", 64 };

  std::vector<buffer> storage;
  for (int i = 0; i < 7; ++i) {
    storage.emplace_back(rs.acquire(64));
  }

  rs.release(std::move(storage.back()));
  storage.pop_back();

  for (int i = 0; i < 20; ++i) {
    rs.release(rs.acquire(64));
  }

  for (auto& it : storage) {
    rs.release(std::move(it));
  }

  REQUIRE_NOTHROW(rs.finalize());
  REQUIRE(cap.str() == "[DEBUG] Cache statistics for test: objects created = "
                       "7, objects released = 7, allocations saved = 20\n");
}

TEST_CASE("reusable_objects::finalize logs error if any objects are missing",
          "[reusable_objects]")
{
  log::log_capture cap;
  reusable_objects<buffer> rs{ "test", 64 };

  std::vector<buffer> storage;
  for (int i = 0; i < 7; ++i) {
    storage.emplace_back(rs.acquire(64));
  }

  rs.release(std::move(storage.back()));
  storage.pop_back();
  storage.pop_back(); // drop one object

  for (int i = 0; i < 20; ++i) {
    rs.release(rs.acquire(64));
  }

  for (auto& it : storage) {
    rs.release(std::move(it));
  }

  REQUIRE_THROWS_AS(rs.finalize(), assert_failed);
  REQUIRE(cap.str() == "[ERROR] Cache statistics for test: objects created = "
                       "7, objects released = 6, allocations saved = 20\n");
}

TEST_CASE("reusable_objects share state", "[reusable_objects]")
{
  reusable_objects<buffer> rs_1{ "test_1", 64 };
  reusable_objects<buffer> rs_2{ rs_1 };

  uint8_t* data = nullptr;
  {
    auto buf = rs_1.acquire(64);
    data = buf.data();
    rs_1.release(std::move(buf));
  }

  auto buf = rs_2.acquire(64);
  REQUIRE(buf.data() == data);
}

} // namespace adapterremoval
