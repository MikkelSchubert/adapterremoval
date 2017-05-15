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
import bz2
import gzip
import argparse


def open_ro(filename):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(filename)
    try:
        header = handle.read(2)

        if header == "\x1f\x8b":
            handle.close()
            # TODO: Re-use handle (fileobj)
            handle = gzip.open(filename)
        elif header == "BZ":
            handle.close()
            handle = bz2.BZ2File(filename)
        else:
            handle.seek(0)

        return handle
    except:
        handle.close()
        raise


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('file')

    return parser.parse_args(argv)


def try_cast(value):
    try:
        return int(value)
    except ValueError:
        return value


def main(argv):
    args = parse_args(argv)

    rows = []
    in_table = False
    in_table_header = None
    with open_ro(args.file) as file_in:
        for line in file_in:
            line_s = line.rstrip()

            if in_table:
                if not in_table_header:
                    in_table_header = line_s.split('\t')
                    if line_s.startswith("#"):
                        in_table_header[0] = in_table_header[0][1:]

                    if "Cycle" not in in_table_header:
                        in_table = False
                        in_table_header = None
                elif not line_s or line_s.startswith("#") or line_s == "<<END":
                    if rows:
                        for row in sorted(rows):
                            print "\t".join(map(str, row))
                    rows = []
                    in_table = False
                    in_table_header = None
                else:
                    row = dict(zip(in_table_header,
                                   map(try_cast, line_s.split('\t'))))
                    row["Cycle"] = row["Cycle"] * 2
                    rows.append(tuple(row[key] for key in in_table_header))
                    row["Cycle"] = row["Cycle"] - 1
                    rows.append(tuple(row[key] for key in in_table_header))
                    continue
            elif line.startswith("["):
                in_table = True

            print line,

        for row in sorted(rows):
            print "\t".join(map(str, row))

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
