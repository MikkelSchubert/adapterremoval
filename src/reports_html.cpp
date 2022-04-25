/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
#include <cstring>  // for size_t, strerror
#include <errno.h>  // for errno
#include <fstream>  // for ofstream
#include <iostream> // for ofstream, operator<<, basic_ostream, endl
#include <string>   // for operator+, string, operator<<
#include <vector>   // for vector

#include "adapterset.hpp"            // for adapter_set
#include "counts.hpp"                // for counts, counts_tmpl
#include "debug.hpp"                 // for AR_FAIL
#include "fastq.hpp"                 // for fastq_pair_vec, IDX_TO_ACGT, fastq
#include "json.hpp"                  // for json_writer, json_section
#include "main.hpp"                  // for NAME, VERSION
#include "reports_template_html.hpp" // for template strings
#include "statistics.hpp"            // for fastq_statistics, ...
#include "strutils.hpp"              // for cli_formatter
#include "userconfig.hpp"            // for userconfig, ...
#include "utilities.hpp"             // for make_shared

bool
write_html_report(const userconfig& config,
                  const statistics& stats,
                  const std::string& filename)
{
  try {
    std::ofstream output(filename, std::ofstream::out);
    if (!output.is_open()) {
      throw std::ofstream::failure(std::strerror(errno));
    }

    output.exceptions(std::ofstream::failbit | std::ofstream::badbit);

    HTMLTmplHead head;
    head.set_name(NAME);
    head.set_version(VERSION);
    head.write(output);

    HTMLTmplBody body;
    body.write(output);

  } catch (const std::ios_base::failure& error) {
    std::cerr << "Error writing JSON report to '" << filename << "':\n"
              << cli_formatter::fmt(error.what()) << std::endl;
    return false;
  }

  return true;
}
