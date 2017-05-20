/*************************************************************************\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2015 by Mikkel Schubert - mikkelsch@gmail.com           *
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
#ifndef UTILS_H
#define UTILS_H

#include <vector>

namespace ar
{

/**
 * Merge two vectors by adding each value in src to each value in 'dst'.
 *
 * If 'dst' is shorter than 'src', it is resized to the same length as 'src'.
 */
template <typename T>
void merge_vectors(std::vector<T>& dst, const std::vector<T>& src)
{
    if (dst.size() < src.size()) {
        dst.resize(src.size());
    }

    typename std::vector<T>::iterator dst_it = dst.begin();
    typename std::vector<T>::const_iterator src_it = src.begin();
    while (src_it != src.end()) {
        *dst_it++ += *src_it++;
    }
}


/**
 * Merge each pair of value in vectors of vectors using merge_vector.
 *
 * If 'dst' is shorter than 'src', it is resized to the same length as 'src'.
 */
template <typename T>
void merge_sub_vectors(std::vector<std::vector<T> >& dst, const std::vector<std::vector<T> >& src)
{
    if (dst.size() < src.size()) {
        dst.resize(src.size());
    }

    typename std::vector<std::vector<T> >::iterator dst_it = dst.begin();
    typename std::vector<std::vector<T> >::const_iterator src_it = src.begin();
    while (src_it != src.end()) {
        merge_vectors(*dst_it++, *src_it++);
    }
}

} // namespace ar

#endif
