// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "commontypes.hpp" // for trimming_strategy, trimming_strategy::mott
#include "debug.hpp"       // for AR_FAIL
#include "fastq_enc.hpp"   // for fastq_encoding
#include "logging.hpp"     // for info, log_stream, error, warn
#include "reports.hpp"     // for print_terminal_postamble, print_terminal_...
#include "simd.hpp"        // for name, instruction_set, instruction_set::none
#include "userconfig.hpp"  // for userconfig, ar_command, ar_command::demul...
#include <iomanip>         // for setprecision
#include <ios>             // for fixed

namespace adapterremoval {

void
print_trimming_parameters(const userconfig& config)
{
  switch (config.trim) {
    case trimming_strategy::mott:
      log::info() << "  - Mott based trimming with Phred encoded quality score "
                  << std::fixed << std::setprecision(1)
                  << fastq_encoding::p_to_phred(config.trim_mott_rate);
      break;
    case trimming_strategy::window:
      log::info() << "  - Window based quality based trimming with window size "
                  << config.trim_window_length << "and minimum quality score "
                  << config.trim_quality_score
                  << (config.trim_ambiguous_bases ? " (including Ns)" : "");
      break;
    case trimming_strategy::per_base:
      if (config.trim_low_quality_bases) {
        log::info() << "  - Per-base based quality based trimming with minimum "
                    << "quality score " << config.trim_quality_score
                    << (config.trim_ambiguous_bases ? " (including Ns)" : "");
      } else if (config.trim_ambiguous_bases) {
        log::info() << "  - Per-base based trimming of Ns";
      } else {
        AR_FAIL("this should not be possible");
      }
      break;
    case trimming_strategy::none:
      log::info() << "  - Quality based trimming disabled";
      break;

    default:
      AR_FAIL("not implemented");
  }
}

void
print_terminal_preamble(const userconfig& config)
{
  log::log_preamble();

  if (config.simd == simd::instruction_set::none) {
    log::warn() << "Hardware accelerated alignments disabled!";
  } else {
    log::info() << "Using " << simd::name(config.simd)
                << " accelerated alignments";
  }

  switch (config.run_type) {
    case ar_command::trim_adapters:
      if (config.paired_ended_mode) {
        log::info() << "Trimming adapters from PE reads:";
      } else {
        log::info() << "Trimming adapters from SE reads:";
      }

      print_trimming_parameters(config);
      break;
    case ar_command::benchmark:
      log::info() << "Benchmarking sub-systems";
      break;
    case ar_command::demultiplex_only:
      log::info() << "Demultiplexing reads";
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
      case ar_command::demultiplex_only:
        log::error() << "AdapterRemoval did not run to completion;\n"
                     << "    do NOT make use of the resulting reads!";
      case ar_command::benchmark:
      case ar_command::report_only:
        break;
      default:
        AR_FAIL("invalid run type");
    }

    return;
  }

  switch (config.run_type) {
    case ar_command::benchmark:
      log::info() << "Benchmarking complete";
      break;
    case ar_command::trim_adapters:
      log::info() << "Adapter trimming complete";
      break;
    case ar_command::demultiplex_only:
      log::info() << "Demultiplexing complete";
      break;
    case ar_command::report_only:
      log::info() << "FASTQ quality report generation complete";
      break;
    default:
      AR_FAIL("invalid run type");
  }
}

} // namespace adapterremoval
