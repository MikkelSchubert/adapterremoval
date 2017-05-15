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
import os
import sys
import argparse


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('root')

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    keys = []
    rows = []

    for run_info in os.listdir(args.root):
        readlen, replicate = run_info.split("_")

        for filename in os.listdir(os.path.join(args.root, run_info)):
            if not filename.endswith(".table"):
                continue

            fpath = os.path.join(args.root, run_info, filename)
            with open(fpath) as handle:
                line = handle.readline()
                header = line.rstrip().split('\t')
                if keys is None:
                    keys = header
                elif header != keys:
                    if set(header).issubset(keys):
                        pass  # Fewer keys
                    elif set(keys).issubset(header):
                        keys = header  # More keys
                    else:
                        assert False, fpath

                for line in handle:
                    rows.append(dict(zip(header, line.rstrip().split('\t'))))
                    rows[-1]["Nth"] = str(int(replicate))  # Strip leading zeros
                    rows[-1]["ReadLen"] = readlen

    keys.insert(0, "Nth")
    keys.insert(0, "ReadLen")
    print "\t".join(keys)
    rows = ["\t".join(row.get(key, "NA") for key in keys) for row in rows]
    rows.sort()

    for row in rows:
        print row

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
