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
import argparse
import bz2
import gzip
import os
import sys


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


def read_barcodes(filename):
    header = None
    results = {}
    for line in open_ro(filename):
        if line.startswith("@"):
            record = dict(zip(header, line.rstrip().split('\t')))

            name = record["readId"].split("/")[0]
            results[name] = "_".join((record["barcode_1"], record["barcode_2"]))
        elif line.startswith("#readId"):
            header = line[1:].rstrip().split('\t')
    return dict(results)


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('info')
    parser.add_argument('files', nargs="+")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    expectations = read_barcodes(args.info)

    total = 0
    correct = 0

    for filename in args.files:
        basename = os.path.basename(filename)
        # Implicitly handle unidentified, will be "unidentified_1", etc.
        barcode = basename.split(".")[1]

        with open_ro(filename) as handle:
            header = handle.readline()
            while header:
                name = header.split(None, 1)[0].split("/")[0]
                if expectations[name] == barcode:
                    correct += 1
                else:
                    print expectations[name], barcode
                total += 1

                handle.readline()
                handle.readline()
                handle.readline()
                header = handle.readline()

    print total, correct

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
