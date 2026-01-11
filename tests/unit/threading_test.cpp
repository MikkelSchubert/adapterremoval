// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
#include "catch.hpp"
#include "errors.hpp"    // for assertion_failed
#include "testing.hpp"   // for TEST_CASE, REQUIRE, ...
#include "threading.hpp" // for threadsafe_data
#include <memory>        // for unique_ptr
#include <mutex>         // for std::mutex
#include <thread>        // for std::thread

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

TEST_CASE("threadstate release lvalue", "[scheduler::threadstate]")
{
  threadstate<int> state;
  state.release(std::make_unique<int>(1234));

  auto value = state.acquire();
  REQUIRE(value);
  REQUIRE(*value == 1234);
  REQUIRE_THROWS_AS(state.acquire(), assert_failed);
}

TEST_CASE("threadstate merg into", "[scheduler::threadstate]")
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
    REQUIRE(*reader == "");
  }

  {
    auto writer = data.get_writer();
    REQUIRE(*writer == "");
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

TEST_CASE("threadsafe_data moves are checked")
{
  threadsafe_data<std::string> data_1{ "example" };
  threadsafe_data<std::string> data_2{ std::move(data_1) };
  threadsafe_data<std::string> data_3{ data_1 };

  REQUIRE_THROWS_AS(data_1.get_reader(), assert_failed);
  REQUIRE_THROWS_AS(data_1.get_writer(), assert_failed);

  REQUIRE_THROWS_AS(data_3.get_reader(), assert_failed);
  REQUIRE_THROWS_AS(data_3.get_writer(), assert_failed);

  REQUIRE(*data_2.get_reader() == "example");
}

TEST_CASE("threadsafe_data readers are shared")
{
  threadsafe_data<std::string> data{ "example" };

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
  std::mutex lock;
  lock.lock();

  const auto write_func = [&data, &lock]() {
    auto writer = data.get_writer();
    lock.unlock();
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    REQUIRE(*writer == "example");
    *writer = "modified";
  };

  const auto read_func = [&data, &lock]() {
    lock.lock();
    const auto reader = data.get_reader();
    REQUIRE(*reader == "modified");
  };

  std::vector<std::thread> threads;
  threads.emplace_back(write_func);
  std::this_thread::sleep_for(std::chrono::milliseconds(1));
  threads.emplace_back(read_func);

  for (auto& it : threads) {
    it.join();
  }
}

TEST_CASE("threadsafe_data writers are exclusive #2")
{
  threadsafe_data<std::string> data{ "example" };

  const auto write_func_1 = [&data]() {
    auto writer = data.get_writer();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    writer->push_back('!');
  };

  const auto write_func_2 = [&data]() {
    auto writer = data.get_writer();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    writer->push_back('?');
  };

  std::vector<std::thread> threads;
  threads.emplace_back(write_func_1);
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  threads.emplace_back(write_func_2);

  for (auto& it : threads) {
    it.join();
  }

  REQUIRE(*data.get_reader() == "example!?");
}

TEST_CASE("threadsafe_data writers are exclusive #3")
{
  threadsafe_data<std::string> data_1{ "example" };
  auto data_2 = data_1;

  const auto write_func_1 = [&data_1]() {
    auto writer = data_1.get_writer();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    writer->push_back('!');
  };

  const auto write_func_2 = [&data_2]() {
    auto writer = data_2.get_writer();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    writer->push_back('?');
  };

  std::vector<std::thread> threads;
  threads.emplace_back(write_func_1);
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  threads.emplace_back(write_func_2);

  for (auto& it : threads) {
    it.join();
  }

  REQUIRE(*data_1.get_reader() == "example!?");
  REQUIRE(*data_2.get_reader() == "example!?");
}

TEST_CASE("readers can be moved")
{
  threadsafe_data<std::string> data{ "example" };

  auto reader_1 = data.get_reader();
  decltype(reader_1) reader_2{ std::move(reader_1) };

  REQUIRE_THROWS_AS(*reader_1, assert_failed);
  REQUIRE_THROWS_AS(reader_1->size(), assert_failed);
  REQUIRE(*reader_2 == "example");
  REQUIRE(reader_2->size() == 7);
}

TEST_CASE("writers can be moved")
{
  threadsafe_data<std::string> data{ "example" };

  auto writer_1 = data.get_writer();
  decltype(writer_1) writer_2{ std::move(writer_1) };

  REQUIRE_THROWS_AS(*writer_1, assert_failed);
  REQUIRE_THROWS_AS(writer_1->size(), assert_failed);
  REQUIRE(*writer_2 == "example");
  REQUIRE(writer_2->size() == 7);
}

} // namespace adapterremoval
