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
#include <cstdlib>

#include "debug.h"
#include "main.h"
#include "userconfig.h"

namespace ar
{

// See main_adapter_rm.cc
int remove_adapter_sequences(const userconfig& config);
// See main_adapter_id.cc
int identify_adapter_sequences(const userconfig& config);

} // namespace ar


int main(int argc, char *argv[])
{
    using namespace ar;
    std::ios_base::sync_with_stdio(false);

    userconfig config(NAME, VERSION, HELPTEXT);
    const argparse::parse_result result = config.parse_args(argc, argv);
    if (result == argparse::pr_error) {
        return 1;
    } else if (result == argparse::pr_exit) {
        // --version, --help, or similar used.
        return 0;
    }

    // User-specified or timestamp-based seed;
    // note, however, that threads are reseeded upon initialization
    srandom(config.seed);

    if (config.identify_adapters) {
        return identify_adapter_sequences(config);
    } else {
        return remove_adapter_sequences(config);
    }

    return 1;
}
