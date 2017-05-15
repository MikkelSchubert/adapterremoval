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
import random
import argparse


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument('--seed', default=None, type=int)
    parser.add_argument('--replicates', default=5, type=int)

    return parser.parse_args(argv)


def read_fasta(fpath):
    with open(fpath) as handle:
        header = handle.readline().rstrip()
        sequence = handle.readline().rstrip()
        assert not handle.readline()

    assert header.startswith(">")

    header = header.split("/", 1)[0].split(None, 1)[0]

    return header, sequence


def shuffle(rng, seq):
    seq = list(seq)
    rng.shuffle(seq)
    return "".join(seq)


def main(argv):
    args = parse_args(argv)
    if args.seed is None:
        args.seed = random.randint(0, sys.maxint)
    rng = random.Random(args.seed)
    name, sequence = read_fasta(args.file)

    print "Seed =", args.seed
    print "Seq  =", sequence

    for _ in xrange(args.replicates - 1):
        print "New  =", shuffle(rng, sequence)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
