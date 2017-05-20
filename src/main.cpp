/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2011 by Stinus Lindgreen - stinus@binf.ku.dk            *
 * Copyright (C) 2014 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#include <iostream>

#include "debug.hpp"
#include "main.hpp"
#include "userconfig.hpp"

namespace ar
{

// See main_adapter_rm.cpp
int remove_adapter_sequences(const userconfig& config);
// See main_adapter_id.cpp
int identify_adapter_sequences(const userconfig& config);
// See main_demultiplex.cpp
int demultiplex_sequences(const userconfig& config);

} // namespace ar


int main(int argc, char *argv[])
{
    using namespace ar;
    std::ios_base::sync_with_stdio(false);

    userconfig config(NAME, VERSION, HELPTEXT);
    switch (config.parse_args(argc, argv)) {
        case argparse::parse_result::error: {
            return 1;
        }

        case argparse::parse_result::exit: {
            // --version, --help, or similar used.
            return 0;
        }

        default: {
            // Ok
        }
    }

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
            return identify_adapter_sequences(config);
        }

        default: {
            std::cerr << "ERROR: Unknown run-type: "
                      << static_cast<size_t>(config.run_type)
                      << std::endl;
            return 1;
        }
    }

    if (returncode) {
        std::cerr << "ERROR: AdapterRemoval did not run to completion;\n"
                  << "       do NOT make use of resulting trimmed reads!"
                  << std::endl;
    }

    return returncode;
}
