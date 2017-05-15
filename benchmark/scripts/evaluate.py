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
import bz2
import math
import gzip
import argparse
import itertools
import collections


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


def parse_timelog(filename):
    if not os.path.exists(filename):
        return float("nan"), float("nan")

    total_time, total_ram = 0.0, 0.0
    with open(filename) as handle:
        for line in handle:
            key, value = line.lower().strip().split(": ")
            if "(wall clock)" in key:
                mult = 1.0
                for field in reversed(value.split(":")):
                    total_time += float(field) * mult
                    mult *= 60.0
            elif key == "maximum resident set size (kbytes)":
                total_ram += float(value)
    return total_time, total_ram


class Summary(object):
    def __init__(self):
        self.true_pos = 0
        self.false_pos = 0
        self.true_neg = 0
        self.false_neg = 0
        self.cpu_time = 0.0
        self.max_ram = 0.0

    @property
    def total(self):
        return self.true_pos + self.false_pos \
            + self.true_neg + self.false_neg

    def print_summary(self, args, stattype, header=True):
        stat = self
        _summary_stats = [
            ("Reads", args.read_mode, None),
            ("Metric", stattype, None),
            ("Name", args.name, None),
            ("Threads", args.threads, None),
            ("Total", stat.total, None),
            ("TP", stat.true_pos, None),
            ("FP", stat.false_pos, None),
            ("TN", stat.true_neg, None),
            ("FN", stat.false_neg, None),

            ("SEN", stat.true_pos, stat.true_pos + stat.false_neg),
            ("SPC", stat.true_neg, stat.true_neg + stat.false_pos),
            ("PPV", stat.true_pos, stat.true_pos + stat.false_pos),
            ("NPV", stat.true_neg, stat.true_neg + stat.false_neg),
            ("MCC", (stat.true_pos * stat.true_neg -
                     stat.false_pos * stat.false_neg),
                    math.sqrt((stat.true_pos + stat.false_pos) *
                              (stat.true_pos + stat.false_neg) *
                              (stat.true_neg + stat.false_pos) *
                              (stat.true_neg + stat.false_neg))),
            ("WallTime", "%.3f" % (stat.cpu_time,), None),
            ("MaxRAM", "%.3f" % (stat.max_ram / (2.0 ** 10),), None),
            ("Rate", ("%.1f"  % ((stat.total / (stat.cpu_time or float("nan"))) / 1e3)).rjust(6), None)
        ]

        row = []
        if header:
            args.out.write("\t".join(name for (name, _, _) in _summary_stats)
                           + "\n")
        for (_, numerator, denominator) in _summary_stats:
            if denominator is None:
                value = str(numerator)
            elif not denominator:
                value = "NA"
            else:
                value = "%.4f" % (float(numerator) / float(denominator),)
            row.append(value)

        args.out.write("\t".join(row) + "\n")


def missing_reads(expectations, se_mode):
    if se_mode:
        reads = [(key, [pair[0], None])
                 for key, pair in expectations.iteritems()
                 if pair[0] is not None]
    else:
        reads = [(key, pair) for key, pair in expectations.iteritems()
                 if pair != [None, None]]

    max_errors = 0
    for name, mates in reads:
        for mate, record in enumerate(mates, start=1):
            if record is not None:
                if max_errors > 4:
                    sys.stderr.write("    ...\n")
                    return True

                sys.stderr.write("    Read not found: %s/%i\n" % (name, mate))
                max_errors += 1

    return bool(reads)


def pop_expected_read(expectations, record):
    is_collapsed = False
    name = record["name"]
    if name.startswith("M_") or name.startswith("MT_"):
        name = name.split("_", 1)[1]
        is_collapsed = True

    expected = expectations[name]
    if is_collapsed:
        expectations[name] = [None, None]
        return expected[0] or expected[1], True

    mate = expected[record["mate"] - 1]
    expected[record["mate"] - 1] = None
    return mate, False


def _perf_trimming(args, expectations, stat, record):
    ins_length, collapsed = pop_expected_read(expectations, record)
    obs_length = record["obs_len"]

    if args.barcode_len and obs_length == args.read_len:
        # Barcode not correctly identified / trimmed
        stat.false_neg += 1
        return

    read_len = args.read_len - args.barcode_len
    if not collapsed:
        if ins_length >= read_len:
            if obs_length != read_len:
                stat.false_pos += 1
            else:
                stat.true_neg += 1
        elif obs_length == ins_length:
            stat.true_pos += 1
        elif obs_length < ins_length:
            stat.false_pos += 1
        elif obs_length > ins_length:
            stat.false_neg += 1
        else:
            assert False, record
    else:
        # If collapsed, the read accounts for 2 sequences
        if obs_length == ins_length:
            # Merged reads containing adapters are considered positives
            if ins_length >= read_len:
                stat.true_neg += 2
            else:
                stat.true_pos += 2
        else:
            stat.false_pos += 2


def _perf_collapsed(args, expectations, stat, record):
    ins_length, collapsed = pop_expected_read(expectations, record)
    obs_length = record["obs_len"]

    if ins_length and ins_length < args.read_len * 2:
        if collapsed:
            if obs_length == ins_length:
                stat.true_pos += 2
            else:
                stat.false_pos += 2
        else:
            stat.false_neg += 1
    elif collapsed:
        stat.false_pos += 2
    else:
        stat.true_neg += 1

    return collapsed


def evaluate(args, expectations, trimmed_reads):
    stat = Summary()

    if args.collapsed:
        eval_func = _perf_collapsed
        eval_title = "Collapsed"
    else:
        eval_func = _perf_trimming
        eval_title = "Trimmed"

    for record in trimmed_reads:
        eval_func(args, expectations, stat, record)

    if args.dimers_are_discarded:
        for key, pair in expectations.items():
            if args.read_mode in ("SE", "MIXED_SE"):
                pair[1] = None

            for mate, read in enumerate(pair, start=1):
                if read is not None:
                    record = {"name": key, "mate": mate, "obs_len": 0}
                    eval_func(args, expectations, stat, record)

    cpu_time, max_ram = parse_timelog(args.time)
    stat.cpu_time += cpu_time
    stat.max_ram += max_ram

    stat.print_summary(args, eval_title)

    return not missing_reads(expectations, args.read_mode in ("SE", "MIXED_SE"))


def find_fastq_files(root):
    results = []
    for filename in os.listdir(root):
        fpath = os.path.join(root, filename)
        if os.path.isfile(fpath):
            with open_ro(fpath) as handle:
                header = handle.readline()
                handle.readline()
                seperator = handle.readline()

                if header.startswith("@") and seperator.startswith("+"):
                    results.append(fpath)

    return results


def read_fastq(filename):
    if not os.path.exists(filename):
        sys.stderr.write("    Skipping %r ...\n" % (filename,))
        return

    counter = 0
    sys.stderr.write("    Reading %r ...\n" % (filename,))
    with open_ro(filename) as handle:
        while True:
            header = handle.readline()
            sequence = handle.readline()
            _ = handle.readline()  # Seperator
            _ = handle.readline()  # Qualities
            if not header:
                break

            fields = header.split(None, 1)
            result = {}

            # Split is needed for ARv1x discarded
            name = fields[0].split("__")[0]
            if name[-2] == "/" and not name.startswith("M_"):
                result["name"] = name[1:-2]
                result["mate"] = int(name[-1])
                assert result["mate"] in (1, 2), result["mate"]
            else:
                result["name"] = name[1:]
                result["mate"] = None
            result["obs_len"] = len(sequence.rstrip())

            yield result
            counter += 1

    sys.stderr.write("        Read %i records\n" % (counter,))


def read_expected_lengths(args):
    header = None
    results = collections.defaultdict(lambda: [None, None])
    for line in open_ro(args.read_info):
        if line.startswith("@"):
            record = dict(zip(header, line.rstrip().split('\t')))

            name = record["readId"][1:]
            name, mate = name.split("/")

            results[name][int(mate) - 1] = int(record["insertSize"])
        elif line.startswith("#read length: "):
            args.read_len = int(line[14:])
        elif line.startswith("#readId"):
            header = line[1:].rstrip().split('\t')
    return dict(results)


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    parser.add_argument('read_info')
    parser.add_argument('--read-mode', default="PE",
                        choices=("SE", "PE", "MIXED_SE", "MIXED_PE"))
    parser.add_argument('--time', default=None)
    parser.add_argument('--name', default="NA")
    parser.add_argument('--collapsed', default=False, action="store_true")
    parser.add_argument('--threads', default=1, type=int)
    parser.add_argument("--dimers-are-discarded", default=False,
                        action="store_true")
    parser.add_argument("--keep", default=False, action="store_true")
    parser.add_argument("--barcode-len", default=0, type=int)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    args.name = args.root.rstrip("/").rsplit("/", 1)[-1]

    name = args.name.rsplit("_", 1)
    if name[-1].startswith("t") and name[-1][1:].isdigit():
        args.name = name[0]
        args.threads = int(name[1][1:])

    if args.time is None:
        args.time = os.path.join(args.root, "time")

    args.out = open(os.path.join(args.root, "table"), "w")

    expected = read_expected_lengths(args)
    fastq_files = find_fastq_files(args.root)

    trimmed_reads = []
    for filename in fastq_files:
        trimmed_reads.append(read_fastq(filename))

    if not evaluate(args, expected, itertools.chain(*trimmed_reads)):
        return 1

    args.out.close()
    os.rename(os.path.join(args.root, "table"), args.root + ".table")

    if not args.keep:
        for fpath in fastq_files:
            os.remove(fpath)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
