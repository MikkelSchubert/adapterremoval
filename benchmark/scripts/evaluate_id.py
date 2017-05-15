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
import collections

import numpy


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    parser.add_argument('--adapter1', default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTG")
    parser.add_argument('--adapter2', default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")

    return parser.parse_args(argv)


def cmp_pcr(obs, ref):
    ncorrect = 0
    for ref_nt, obs_nt in zip(ref, obs):
        if ref_nt == obs_nt:
            ncorrect += 1
        else:
            break

    return ncorrect


class evalulate_minion(object):
    def __init__(self, pcr1, pcr2):
        self._pcr1 = pcr1
        self._pcr2 = pcr2
        self._results = collections.defaultdict(dict)

    def __call__(self, readlen, insert_mean, replicate, fpath):
        prefix = os.path.join(fpath, "mate")

        for method, result in self._parse_report_minion(prefix).items():
            self._results[(readlen, insert_mean, method)][replicate] = result

    def finalize(self):
        return dict(self._results)

    def _parse_report_minion(self, prefix):
        results = collections.defaultdict(dict)
        for result in self._parse_report_minion_part(prefix + "1.txt"):
            key = "minion:%s" % (result["criterion"],)
            results[key]["PCR1"] = max(results[key].get("PCR1", 0),
                                       cmp_pcr(result["sequence"], self._pcr1))

        for result in self._parse_report_minion_part(prefix + "2.txt"):
            key = "minion:%s" % (result["criterion"],)
            results[key]["PCR2"] = max(results[key].get("PCR2", 0),
                                       cmp_pcr(result["sequence"], self._pcr2))

        return dict(results)

    @classmethod
    def _parse_report_minion_part(cls, filename):
        results = []
        with open(filename) as handle:
            result = {}
            for line in handle:
                line = line.strip()
                if not line:
                    if result:
                        results.append(result)
                        result = {}
                else:
                    key, value = map(str.strip, line.split("=", 1))
                    result[key] = value

            if result:
                results.append(result)

        return results


class evalulate_adapterrm(object):
    def __init__(self, pcr1, pcr2):
        self._pcr = [pcr1, pcr2]
        self._results = collections.defaultdict(dict)

    def __call__(self, readlen, insert_mean, replicate, fpath):
        prefix = os.path.join(fpath, "mates.txt")

        self._results[(readlen, insert_mean, "adapterremovalv2")][replicate] \
            = self._parse_report(prefix)

    def finalize(self):
        return dict(self._results)

    def _parse_report(self, fpath):
        result = {}
        with open(fpath) as handle:
            key = "Consensus:"
            counter = 1
            for line in handle:
                line = line.strip()
                if line.startswith(key):
                    adapter = line[line.index(key) + len(key):].strip()
                    result["PCR%i" % (counter,)] \
                        = cmp_pcr(adapter, self._pcr[counter - 1])

                    counter += 1
                    if counter > 2:
                        break
        return result


def print_rows(program, results):
    rows = []
    for (readlen, insert_mean, method), replicates in sorted(results.items()):
        mean_pcr1 = numpy.mean([repl["PCR1"] for repl in replicates.values()])
        mean_pcr2 = numpy.mean([repl["PCR2"] for repl in replicates.values()])

        rows.append((method, readlen, insert_mean, mean_pcr1, mean_pcr2))

    for row in sorted(rows):
        print "\t".join(map(str, row))


def main(argv):
    args = parse_args(argv)
    eval_funcs = {
        "minion": evalulate_minion(args.adapter1, args.adapter2),
        "adapterremovalv2": evalulate_adapterrm(args.adapter1, args.adapter2),
    }

    for run in os.listdir(args.root):
        readlen, insert_mean, replicate = map(int, run.split("_"))

        for program in os.listdir(os.path.join(args.root, run)):
            if program in eval_funcs:
                fpath = os.path.join(args.root, run, program)
                eval_funcs[program](readlen, insert_mean, replicate, fpath)

    print "\t".join(("Method", "ReadLen", "InsertMean", "PCR1", "PCR2"))
    for program, func in sorted(eval_funcs.items()):
        print_rows(program, func.finalize())

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
