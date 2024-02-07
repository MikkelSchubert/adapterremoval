#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "debug.hpp" // for assert_failed
#include <string>    // for string

namespace adapterremoval {

[[noreturn]] void
terminate(const std::string& message)
{
  throw assert_failed(message);
}

} // namespace adapterremoval
