// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "main.hpp"       // declarations
#include "argparse.hpp"   // for parse_result, parse_result::error, parse_...
#include "debug.hpp"      // for AR_FAIL
#include "logging.hpp"    // for log_stream, error
#include "reports.hpp"    // for print_terminal_postamble, print_terminal_...
#include "userconfig.hpp" // for ar_command, userconfig, ar_command::demul...
#include <cstdlib>        // for abort, size_t
#include <exception>      // for set_terminate
#include <ios>            // for ios_base
#include <vector>         // for vector

namespace adapterremoval {

// See main_adapter_rm.cpp
int
remove_adapter_sequences(const userconfig& config);
// See main_adapter_id.cpp
int
identify_adapter_sequences(const userconfig& config);
// See main_demultiplex.cpp
int
demultiplex_sequences(const userconfig& config);
// See main_fastq_ro.cpp
int
fastq_report_only(const userconfig& config);
// See main_benchmark.cpp
int
benchmark(const userconfig& config);

[[noreturn]] void
terminate(const std::string& message)
{
  log::error()
    << message << "\n"
    << "This should not happen! Please file a bug-report at\n"
    << "    https://github.com/MikkelSchubert/adapterremoval/issues/new\n\n"
    << "Please specify the AdapterRemoval version used (" << VERSION << ")\n"
    << "and if possible describe the steps taken to reproduce this problem.";

  std::abort();
}

[[noreturn]] void
terminate_on_exception()
{
  try {
    std::exception_ptr eptr{ std::current_exception() };

    if (eptr) {
      std::rethrow_exception(eptr);
    } else {
      terminate("Aborted due to unknown error");
    }
  } catch (const std::exception& ex) {
    log::error()
      << ex.what()
      << "\nAdapterRemoval did not run to completion; do NOT make use of the "
         "resulting reads!";
  } catch (...) {
    terminate("Caught unknown exception type");
  }

  std::exit(1);
}

} // namespace adapterremoval

int
main(int argc, char* argv[])
{
  using namespace adapterremoval;

  std::set_terminate(&terminate_on_exception);
  std::ios_base::sync_with_stdio(false);

  userconfig config;

  const std::vector<std::string> argvec(argv, argv + argc);
  switch (config.parse_args(argvec)) {
    case argparse::parse_result::error:
      return 1;

    case argparse::parse_result::exit:
      // --version, --help, or similar used.
      return 0;

    case argparse::parse_result::ok:
      // Ok
      break;

    default:
      AR_FAIL("invalid argparse::parse_result");
  }

  print_terminal_preamble(config);

  auto returncode = 0;
  switch (config.run_type) {
    case ar_command::demultiplex_only:
    case ar_command::trim_adapters: {
      returncode = remove_adapter_sequences(config);
      break;
    }

    case ar_command::benchmark: {
      returncode = benchmark(config);
      break;
    }

    case ar_command::report_only: {
      returncode = identify_adapter_sequences(config);
      break;
    }

    default: {
      log::error() << "Unknown run-type: "
                   << static_cast<size_t>(config.run_type);
      return 1;
    }
  }

  print_terminal_postamble(config, returncode);

  return returncode;
}
