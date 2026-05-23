// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2011 Stinus Lindgreen <stinus@binf.ku.dk>
// SPDX-FileCopyrightText: 2014 Mikkel Schubert <mikkelsch@gmail.com>
#include "main.hpp"       // declarations
#include "argparse.hpp"   // for parse_result, parse_result::error, parse_...
#include "debug.hpp"      // for AR_FAIL, terminate_on_exception
#include "errors.hpp"     // for assert_failed, program_failure
#include "logging.hpp"    // for log_stream, error
#include "reports.hpp"    // for print_terminal_postamble, print_terminal_...
#include "userconfig.hpp" // for ar_command, userconfig, ar_command::demul...
#include "version.hpp"    // for version
#include <cstdlib>        // for abort, size_t, quick_exit
#include <exception>      // for set_terminate
#include <filesystem>     // for filesystem_error
#include <ios>            // for ios_base
#include <system_error>   // for system_error
#include <vector>         // for vector

namespace adapterremoval {

[[noreturn]] void
terminate_on_bug(std::string_view message)
{
  try {
    log::error()
      << message << "\n"
      << "This should not happen! Please file a bug-report at\n"
      << "    https://github.com/MikkelSchubert/adapterremoval/issues/new\n\n"
      << "Please specify the AdapterRemoval version ("
      << program::long_version()
      << ")\nand if possible describe the steps taken to reproduce this "
         "problem.";
  } catch (...) { // NOLINT(bugprone-empty-catch)
    // function is called by std::terminate, so cannot rethrow exceptions
  }

  std::abort();
}

[[noreturn]] void
terminate_on_failure(std::string_view message)
{
  try {
    log::error()
      << message
      << "\nAdapterRemoval did not run to completion; do NOT make use of the "
         "resulting reads!";
  } catch (...) { // NOLINT(bugprone-empty-catch)
    // function is called by std::terminate, so cannot rethrow exceptions
  }

  std::quick_exit(1);
}

[[noreturn]] void
terminate_on_exception()
{
  try {
    std::exception_ptr eptr{ std::current_exception() };

    if (eptr) {
      std::rethrow_exception(eptr);
    } else {
      terminate_on_bug("Aborted due to unknown error");
    }
  } catch (const program_failure& ex) {
    terminate_on_failure(ex.what()); // probably not a bug
  } catch (const std::filesystem::filesystem_error& ex) {
    terminate_on_failure(ex.what()); // probably not a bug
  } catch (const std::system_error& ex) {
    terminate_on_failure(ex.what()); // may or may not be a bug
  } catch (const std::exception& ex) {
    terminate_on_bug(ex.what()); // probably a bug
  } catch (...) {
    terminate_on_bug("Caught unknown exception type");
  }

  std::abort();
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
      break;

    default:
      AR_FAIL("invalid argparse::parse_result");
  }

  print_terminal_preamble(config);

  auto returncode = 0;
  switch (config.run_type) {
    case ar_command::demultiplex_only:
      [[fallthrough]];
    case ar_command::trim_adapters: {
      returncode = remove_adapter_sequences(config);
      break;
    }

    case ar_command::benchmark: {
      returncode = benchmark(config);
      break;
    }

    case ar_command::report_only: {
      returncode = generate_reports(config);
      break;
    }

    default: {
      log::error() << "Unknown run-type: "
                   << static_cast<size_t>(config.run_type);
      return 1;
    }
  }

  // Print final summary, or warnings if the run failed
  print_terminal_postamble(config, returncode);

  return returncode;
}
