/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include "debug.hpp"      // for AR_FAIL
#include "logging.hpp"    // for log
#include "main.hpp"       // for NAME, VERSION
#include "userconfig.hpp" // for userconfig

namespace adapterremoval {

namespace {

__attribute__((target("sse2"))) bool
supports_sse2()
{
  return true;
}

__attribute__((target("default"))) bool
supports_sse2()
{
  return false;
}

__attribute__((target("avx2"))) bool
supports_avx2()
{
  return true;
}

__attribute__((target("default"))) bool
supports_avx2()
{
  return false;
}

} // namespace

void
print_terminal_preamble(const userconfig& config)
{
  if (supports_sse2() || supports_avx2()) {
    log::info() << NAME << " " << VERSION << " ("
                << (supports_sse2() ? "SSE2" : "")
                << (supports_avx2() ? " AVX2" : "") << ")";
  } else {
    log::info() << NAME << " " << VERSION;
    log::warn() << "Hardware acceleration (SSE2/AVX2) disabled!";
  }

  switch (config.run_type) {
    case ar_command::trim_adapters:
      if (config.paired_ended_mode) {
        log::info() << "Trimming adapters from PE reads";
      } else {
        log::info() << "Trimming adapters from SE reads";
      }
      break;
    case ar_command::demultiplex_sequences:
      log::info() << "Demultiplexing reads";
      break;
    case ar_command::identify_adapters:
      log::info() << "Attempting to identify adapter sequences";
      break;
    case ar_command::report_only:
      log::info() << "Generating FASTQ quality report";
      break;
    default:
      AR_FAIL("invalid run type");
  }
}

void
print_terminal_postamble(const userconfig& config, bool any_errors)
{
  if (any_errors) {
    switch (config.run_type) {
      case ar_command::trim_adapters:
      case ar_command::demultiplex_sequences:
        log::error() << "AdapterRemoval did not run to completion;\n"
                     << "    do NOT make use of the resulting reads!";
      case ar_command::identify_adapters:
      case ar_command::report_only:
        break;
      default:
        AR_FAIL("invalid run type");
    }

    return;
  }

  switch (config.run_type) {
    case ar_command::trim_adapters:
      log::info() << "Adapter trimming complete";
      break;
    case ar_command::demultiplex_sequences:
      log::info() << "Demultiplexing complete";
      break;
    case ar_command::identify_adapters:
      log::info() << "Adapter identification complete";
      break;
    case ar_command::report_only:
      log::info() << "FASTQ quality report generation complete";
      break;
    default:
      AR_FAIL("invalid run type");
  }
}

} // namespace adapterremoval
