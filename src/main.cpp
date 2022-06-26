/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * Schubert et al. (2016). AdapterRemoval v2: rapid adapter trimming,    *
 * identification, and read merging. BMC Research Notes, 12;9(1):88      *
 * https://doi.org/10.1186/s13104-016-1900-2                             *
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
#include <stddef.h> // for size_t

#include "argparse.hpp"   // for parse_result, parse_result::error
#include "debug.hpp"      // for AR_FAIL
#include "logging.hpp"    // for log
#include "main.hpp"       // header
#include "reports.hpp"    // for print_terminal_preamble
#include "userconfig.hpp" // for userconfig, ar_command, ar_command::demult...

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

[[noreturn]] void
terminate(const std::string& message)
{
  log::error() << message << "\n"
               << "This should not happen! Please file a bug-report at\n    "
               << "https://github.com/MikkelSchubert/adapterremoval/issues/new";

  std::abort();
}

} // namespace adapterremoval

int
main(int argc, char* argv[])
{
  using namespace adapterremoval;

  std::ios_base::sync_with_stdio(false);

  userconfig config(NAME, VERSION, HELPTEXT);

  switch (config.parse_args(argc, argv)) {
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
    case ar_command::trim_adapters: {
      returncode = remove_adapter_sequences(config);
      break;
    }

    case ar_command::demultiplex_sequences: {
      returncode = demultiplex_sequences(config);
      break;
    }

    case ar_command::identify_adapters: {
      returncode = identify_adapter_sequences(config);
      break;
    }

    case ar_command::report_only: {
      returncode = fastq_report_only(config);
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
