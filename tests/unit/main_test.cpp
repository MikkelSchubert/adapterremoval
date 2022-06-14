#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "debug.hpp"

namespace adapterremoval {

[[noreturn]] void
terminate(const std::string& message)
{
  throw assert_failed(message);
}

} // namespace adapterremoval
