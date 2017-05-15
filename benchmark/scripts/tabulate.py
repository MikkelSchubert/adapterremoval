#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import sys
import argparse
import collections

import numpy


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('filter', choices=("basic", "throughput"),
                        type=str.lower,
                        help="TODO")
    parser.add_argument('table', help="Benchmarking results", nargs="+")

    return parser.parse_args(argv)


def read_table(fpath, filter_func):
    results = []
    with open(fpath) as handle:
        header = handle.readline().rstrip().split("\t")
        for line in handle:
            fields = line.rstrip().split("\t")
            assert len(fields) == len(header), (header, fields)
            row = dict(zip(header, fields))

            name = row["Name"]
            if name.endswith("_gz") or name.endswith("_bz2"):
                row["Name"], row["Compression"] = name.rsplit("_", 1)
            else:
                row["Compression"] = "NA"

            if filter_func(row):
                results.append(row)

    return results


def _filter_basic(row):
    return int(row["ReadLen"]) == 100 \
        and int(row["Threads"]) == 1


def _filter_throughput(row):
    return int(row["ReadLen"]) >= 100 \
        and int(row["Threads"]) >= 1


_FILTERS = {
    "basic": _filter_basic,
    "throughput": _filter_throughput
}


def _collapse_basic(rows):
    max_nth = 0

    keys = ("Name", "ReadLen", "Reads", "Metric", "Compression")
    results = collections.defaultdict(dict)
    for row in rows:
        key_row = tuple(row[key] for key in keys)
        key_nth = int(row["Nth"])

        max_nth = max(max_nth, key_nth)

        subdd = results[key_row]
        assert key_nth >= 1, (row, subdd)
        assert key_nth not in subdd, (key_nth, subdd)

        subdd[key_nth] = row

    for key, values in results.iteritems():
        assert len(values) == max_nth and max(values) == max_nth, values

    collapsed_rows = []
    value_keys = ("SEN", "SPC", "PPV", "NPV", "MCC", "Rate")
    for rows in results.itervalues():
        result = {}

        for key in value_keys:
            if any((key in row) for row in rows.itervalues()):
                values = tuple(float(row[key]) for row in rows.itervalues())
                result[key] = "%.3f" % (numpy.mean(values),)

        for key in keys:
            result[key] = rows[1][key]

        collapsed_rows.append(result)

    return ("ReadLen", "Reads", "Name", "Metric", "Compression") + value_keys, \
        collapsed_rows


def _collapse_throughput(rows):
    max_nth = 0

    keys = ("Name", "Reads", "Metric", "Compression")
    rate_keys = set()

    results = collections.defaultdict(dict)
    for row in rows:
        key_row = tuple(row[key] for key in keys)
        key_nth = int(row["Nth"])
        assert key_nth >= 1, (row, key_row)

        max_nth = max(max_nth, key_nth)

        subdd = results[key_row]
        subdd.setdefault(key_nth, {})
        subdd = subdd[key_nth]

        if int(row["Threads"]) == 1 and int(row["ReadLen"]) == 100:
            subdd.update(row)

        key = "%03i/%i" % (int(row["ReadLen"]), int(row["Threads"]))
        rate_keys.add(key)

        subdd[key] = row["Rate"]

    for key, values in results.iteritems():
        assert len(values) == max_nth and max(values) == max_nth, values

    collapsed_rows = []
    value_keys = tuple(sorted(rate_keys))
    for rows in results.itervalues():
        result = {}

        for key in value_keys:
            if any((key in row) for row in rows.itervalues()):
                values = tuple(float(row[key]) for row in rows.itervalues())
                result[key] = "%.1f" % (numpy.mean(values),)

        for key in keys:
            result[key] = rows[1][key]

        collapsed_rows.append(result)

    return ("Name", "Metric", "Compression") + value_keys, \
        collapsed_rows


_COLLAPSERS = {
    "basic": _collapse_basic,
    "throughput": _collapse_throughput,
}


def get_output_keys(rows):
    keys = ["Name", "ReadLen", "Reads", "Metric", "Compression",
            "SEN", "SPC", "PPV", "NPV", "MCC"]
    max_threads = 1

    for row in rows:
        rate_key = "Rate(%i)" % (max_threads,)
        if rate_key in row:
            keys.append(rate_key)
            max_threads += 1

    return tuple(keys)


def main(argv):
    args = parse_args(argv)

    rows = []
    for fpath in args.table:
        rows.extend(read_table(fpath, _FILTERS[args.filter]))

    keys, rows = _COLLAPSERS[args.filter](rows)

    output_rows = []
    for row in sorted(rows):
        output_rows.append(tuple(str(row.get(key, "NA")) for key in keys))

    print "\t".join(keys)
    for row in sorted(output_rows):
        print "\t".join(row)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
