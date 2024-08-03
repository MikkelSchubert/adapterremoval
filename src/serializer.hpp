/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2024 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#pragma once

#include <cstddef> // for size_t

namespace adapterremoval {

class fastq;
class buffer;
enum class output_format;

enum class fastq_flags
{
  //! SE read
  se,
  //! SE read that failed QC
  se_fail,
  //! Mate 1 read
  pe_1,
  //! Mate 1 read that failed QC
  pe_1_fail,
  //! Mate 1 read
  pe_2,
  //! Mate 1 read that failed QC
  pe_2_fail,

};

class fastq_serializer
{
public:
  static void header(buffer& buf, output_format format);
  static void record(buffer& buf,
                     const fastq& record,
                     output_format format,
                     fastq_flags flags);
};

} // namespace adapterremoval
