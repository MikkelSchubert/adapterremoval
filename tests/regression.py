#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2022 Mikkel Schubert <MikkelSch@gmail.com>
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
from collections import namedtuple
import copy
import difflib
import errno
import gzip
import io
import itertools
import json
import math
import multiprocessing
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import traceback
import typing
from itertools import groupby, islice


#############################################################################
_COLORS_ENABLED = not os.getenv("NO_COLOR")


def _do_print_color(*vargs_, **kwargs):
    """Utility function: Prints using shell colors."""
    colorcode = kwargs.pop("colorcode")
    destination = kwargs.pop("file", sys.stdout)

    # No colors if output is redirected (e.g. less, file, etc.)
    if _COLORS_ENABLED and destination.isatty():
        vargs = list(vargs_)
        for index, varg in enumerate(vargs):
            varg_lines = []
            # Newlines terminate the color-code for e.g. 'less', so ensure that
            # each line is color-coded, while preserving the list of arguments
            for line in str(varg).split("\n"):
                varg_lines.append("\033[00;%im%s\033[00m" % (colorcode, line))
            vargs[index] = "\n".join(varg_lines)
    else:
        vargs = vargs_

    print(*vargs, file=destination, **kwargs)

    if "\n" in kwargs.get("end", "\n"):
        destination.flush()


def print_ok(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (green)."""
    _do_print_color(*vargs, colorcode=32, **kwargs)


def print_warn(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (green)."""
    _do_print_color(*vargs, colorcode=33, **kwargs)


def print_err(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (red)."""
    _do_print_color(*vargs, colorcode=31, **kwargs)


#############################################################################


def read_file(filename, mode="rt"):
    try:
        with open(filename, mode) as handle:
            return handle.read()
    except OSError as error:
        raise TestError("ERROR while reading data:\n    {}".format(error))


def read_json(filename):
    data = read_file(filename)

    try:
        return data, json.loads(data)
    except json.JSONDecodeError as error:
        raise TestError("ERROR while reading {!r}:\n    {}".format(filename, error))


def read_lines(filename):
    return read_file(filename).splitlines(keepends=True)


def read_and_decompress_lines(filename):
    value = read_file(filename, "rb")

    header = value[:2]
    if header == b"\x1f\x8b":
        if not filename.endswith(".gz"):
            raise TestError("{} is gzipped, but lacks .gz extension".format(filename))

        value = gzip.decompress(value)
    elif filename.endswith(".gz"):
        raise TestError("{} has .gz extension, but is not compressed".format(filename))

    return io.StringIO(value.decode("utf-8")).readlines()


def escape_whitespace(line):
    out = []
    whitespace = {
        "\n": "\\n",
        "\r": "\\r",
        "\t": "\\t",
        "\b": "\\b",
        "\f": "\\f",
        " ": " ",
        "\\": "\\\\",
    }

    for char in line:
        escaped = whitespace.get(char)
        if escaped is None:
            if char.isspace():
                escaped = "\\{:o}".format(ord(char))
            else:
                escaped = char

        out.append(escaped)

    return "".join(out)


def truncate_lines(lines, max_lines=float("inf")):
    lines = list(islice(lines, max_lines + 1))
    if len(lines) > max_lines:
        lines[-1] = "..."

    return lines


def pretty_output(lines, max_lines, padding=0):
    result = []
    prefix = " " * padding
    lines = truncate_lines(lines, max_lines)
    for line in lines:
        result.append("%s>  %s" % (prefix, escape_whitespace(line)))

    return "\n".join(result)


def filter_errors(
    text,
    _re=re.compile(r"(\[(?:ERROR|WARNING)\] .*)"),
):
    text = text.split("\n")
    errors = []
    for line in text:
        match = _re.search(line)
        if match:
            errors.append(match.group())

    return errors if errors else text


def write_data(path, data):
    open_ = gzip.open if path.endswith(".gz") else open
    with open_(path, "wt") as handle:
        handle.write(data)


def cmd_to_s(cmd):
    return " ".join(shlex.quote(field) for field in cmd)


def path_to_s(path):
    return ".".join(str(value) for value in path)


def classname(v):
    if v is None:
        return "None"

    return v.__class__.__name__


################################################################################

# Wildcards with (limited) type checking
JSON_WILDCARDS = {
    "...": lambda _: True,
    "...str": lambda it: isinstance(it, str),
    "...float": lambda it: isinstance(it, float) and it == it,
    "...[str]": lambda it: isinstance(it, list) and all(isinstance(v, str) for v in it),
}

# JSON path elements that do not require quoting
JSON_PATH = re.compile("[a-z0-9_]", re.I)


def diff_json(reference, observed, path=("",)):
    def _err(what, ref, obs):
        return "{} at {}: {} is not {}".format(what, path_to_s(path), obs, ref)

    if reference == observed:
        return

    if isinstance(reference, typing.Hashable) and reference in JSON_WILDCARDS:
        if not JSON_WILDCARDS[reference](observed):
            yield _err("wildcard mismatch", repr(reference), repr(observed))
    elif type(reference) != type(observed):
        yield _err("type mismatch", classname(reference), classname(observed))
    elif isinstance(reference, list) and isinstance(observed, list):
        if len(reference) != len(observed):
            yield _err("length mismatch", len(reference), len(observed))

        for idx, (ref, obs) in enumerate(zip(reference, observed)):
            yield from diff_json(ref, obs, path + ("[{}]".format(idx),))
    elif isinstance(reference, dict):
        for key in reference.keys() | observed.keys():
            if key not in reference:
                yield _err("unexpected key", "in reference", key)
            elif key not in observed:
                yield _err("unexpected key", "in observation", key)
            else:
                if not JSON_PATH.match(key):
                    key = repr(key)

                yield from diff_json(reference[key], observed[key], path + (key,))
    elif isinstance(reference, float):
        # There should be no NANs in this data, so this test also catches those
        if not (abs(reference - observed) < 1e-6):
            yield _err("mismatch", repr(observed), repr(reference))
    else:
        yield _err("mismatch", repr(reference), repr(observed))


################################################################################


class TestError(RuntimeError):
    pass


_INPUT_FILES_FQ = {
    "input_1",
    "input_2",
}

# Files supplied by the end-user
_INPUT_FILES = {
    "barcodes",
    "adapters",
} | _INPUT_FILES_FQ

_OUTPUT_FILES_FQ = {
    "output",
}

# Files generated by AdapterRemoval
_OUTPUT_FILES = {
    "json",
    "html",
    "ignore",
} | _OUTPUT_FILES_FQ

_TEST_FILES = _INPUT_FILES | _OUTPUT_FILES


# Option value in regression test specification
Default = namedtuple("Default", ("spec", "value"))


# Layout of regression tests in JSON format
_TEST_SPECIFICATION = {
    "arguments": [str],
    "skip": Default(bool, False),
    "exhaustive": Default((bool, type(None)), None),
    "return_code": Default(int, 0),
    "stdout": Default([str], []),
    "stderr": Default([str], []),
    "files": Default(
        {
            "input_1": Default([str], []),
            "input_2": Default([str], []),
            "adapters": Default([str], []),
            "barcodes": Default([str], []),
            "output": Default([str], []),
            "json": Default([str], []),
            "html": Default([str], []),
            "ignore": Default([str], []),
        },
        {},
    ),
}


def validate(spec, value, path=("{}",)):
    def _check_type(exp, obs, path=path):
        if not isinstance(obs, exp):
            raise TestError(
                "expected {} at {}, found {!r}".format(exp, path_to_s(path), obs)
            )

        return obs

    if isinstance(spec, Default):
        if value is None:
            value = copy.deepcopy(spec.value)

        return validate(spec.spec, value, path)
    elif isinstance(spec, dict):
        _check_type(dict, value)

        out = {}
        for key, subspec in spec.items():
            out[key] = validate(subspec, value.pop(key, None), path + (key,))

        if value:
            raise TestError(
                "unexpected values at {}: {}".format(path_to_s(path), value)
            )

        return out
    elif isinstance(spec, list):
        _check_type(list, value)

        (spec,) = spec
        return [validate(spec, value, path + ("[]",)) for value in value]
    elif isinstance(spec, (set, frozenset)):
        if value not in spec:
            raise TestError("unexpected value at {}: {}".format(path_to_s(path), value))

        return value
    else:
        return _check_type(spec, value)


class TestFile(namedtuple("TestFile", ("name", "kind", "data", "json"))):
    @classmethod
    def parse(cls, root, name, kind):
        data = None
        json_data = None
        # HTML files are not compared in detail
        if kind not in ("html", "ignore"):
            filename = os.path.join(root, name)
            if kind == "json":
                data, json_data = read_json(filename)
            else:
                # Reference data is not expected to be compresses, regardless of extension
                data = read_lines(filename)

        return TestFile(name=name, kind=kind, data=data, json=json_data)

    def compare_with_file(self, expected, observed: str):
        if not os.path.isfile(observed):
            raise TestError("file {} not created".format(observed))

        # Some files just need to exist; HTML reports, etc.
        if self.kind == "json":
            return self._compare_json(expected, observed)
        elif self.data is None:
            return

        return self._compare_text(expected, observed)

    def _compare_json(self, expected: str, observed: str):
        _, data = read_json(observed)
        differences = diff_json(reference=self.json, observed=data)
        differences = truncate_lines(differences, 4)

        if differences:
            differences = "\n".join(
                "    {}. {}".format(idx, line)
                for idx, line in enumerate(differences, start=1)
            )

            self._raise_test_error(
                label="JSON",
                expected=expected,
                observed=observed,
                differences=differences,
            )

    def _compare_text(self, expected: str, observed: str):
        data = read_and_decompress_lines(observed)
        if data == self.data:
            return

        lines = difflib.unified_diff(self.data, data, "expected", "observed")

        self._raise_test_error(
            label="output",
            expected=expected,
            observed=observed,
            differences=pretty_output(lines, 10),
        )

    def _raise_test_error(self, label, expected, observed, differences):
        raise TestError(
            (
                "Mismatches in {} file:\n"
                "  Expected   = {}\n"
                "  Observed   = {}\n"
                "  Mismatches =\n{}"
            ).format(
                label,
                shlex.quote(expected),
                shlex.quote(observed),
                differences,
            )
        )


class TestConfig(
    namedtuple(
        "TestConfig",
        (
            "path",
            "name",
            "variant",
            "skip",
            "exhaustive",
            "arguments",
            "return_code",
            "stdout",
            "stderr",
            "files",
        ),
    )
):
    @classmethod
    def load(cls, name, filepath):
        _, data = read_json(filepath)
        data = validate(_TEST_SPECIFICATION, data)
        assert isinstance(data, dict)

        root = os.path.dirname(filepath)
        files = []
        for key, values in data.pop("files").items():
            for filename in values:
                files.append(TestFile.parse(root, filename, key))

        self = TestConfig(
            path=filepath,
            name=name,
            variant=(),
            arguments=tuple(data.pop("arguments")),
            skip=data.pop("skip"),
            exhaustive=data.pop("exhaustive"),
            return_code=data.pop("return_code"),
            stdout=tuple(data.pop("stdout")),
            stderr=tuple(data.pop("stderr")),
            files=tuple(files),
        )

        assert not data, data

        return self

    def build_command(self, executable):
        command = [executable] + list(self.arguments)
        input_1 = []
        input_2 = []

        for it in self.files:
            if it.kind == "input_1":
                input_1.append(it.name)
            elif it.kind == "input_2":
                input_2.append(it.name)
            elif it.kind == "barcodes":
                command += ["--barcode-list", it.name]
            elif it.kind == "adapters":
                command += ["--adapter-list", it.name]
            elif it.kind == "json":
                command += ["--report-sample-rate", "1"]
            elif it.kind not in ("output", "json", "html", "ignore"):
                raise NotImplementedError(it.kind)

        if input_1:
            command += ["--file1"] + input_1

        if input_2:
            command += ["--file2"] + input_2

        return command

    def get_files(self, keys):
        if set(keys) - _TEST_FILES:
            raise ValueError(keys)

        for item in self.files:
            if item.kind in keys:
                yield item


class TestMutator:
    @classmethod
    def create_exhaustive_tests(cls, it, exhaustive):
        if it.exhaustive is not None:
            exhaustive = it.exhaustive

        if not (exhaustive and it.files):
            yield it
            return

        # For simplicity's sake, filenames are not checked in exhaustive tests
        for file in it.files:
            if file.kind == "json":
                masked_stats = file._replace(json=cls._mask_filenames(file.json))
                break
        else:
            masked_stats = None

        # Test with split or interleaved input (if possible)
        for it in cls._interleave_input(it, masked_stats):
            # Test with gzip compression of input and output files
            for compress_in, compress_out in itertools.product((False, True), repeat=2):
                variant = []
                arguments = []
                compressed_files = set()

                if compress_in:
                    variant.append("igz")
                    compressed_files |= _INPUT_FILES_FQ

                if compress_out:
                    variant.append("ogz")
                    arguments.append("--gzip")
                    compressed_files |= _OUTPUT_FILES_FQ

                if variant:
                    yield it._replace(
                        variant=it.variant + tuple(variant),
                        arguments=it.arguments + tuple(arguments),
                        files=cls._compress_files(
                            it.files,
                            compressed_files,
                            masked_stats,
                        ),
                    )
                else:
                    yield it

    @classmethod
    def _interleave_input(cls, test, masked_stats):
        yield test

        input_files = {}
        other_files = []
        for it in test.files:
            if it.kind in _INPUT_FILES_FQ:
                # Ensure exactly one trailing newline
                lines = list(it.data)
                while lines and not lines[-1].rstrip():
                    lines.pop()
                lines.append("\n")

                input_files.setdefault(it.kind, []).append(lines)
            elif it.kind == "json":
                if masked_stats is not None:
                    other_files.append(masked_stats)
            else:
                other_files.append(it)

        if len(input_files) == 2:
            files_1, files_2 = input_files.values()

            n_lines_1 = sum(map(len, files_1))
            n_lines_2 = sum(map(len, files_2))
            if n_lines_1 == n_lines_2 and n_lines_1 % 4 == 0:
                for idx, (lines_1, lines_2) in enumerate(zip(files_1, files_2)):
                    records_1 = cls._lines_to_records(lines_1)
                    records_2 = cls._lines_to_records(lines_2)

                    data = []
                    for record_1, record_2 in zip(records_1, records_2):
                        data.extend(record_1)
                        data.extend(record_2)

                    other_files.append(
                        TestFile(
                            name="input_{}.fastq".format(idx),
                            kind="input_1",
                            data=data,
                            json=None,
                        )
                    )

                yield test._replace(
                    variant=test.variant + ("intl",),
                    files=other_files,
                    arguments=test.arguments + ("--interleaved-input",),
                )

    @staticmethod
    def _lines_to_records(lines):
        assert len(lines) % 4 == 0, lines
        for idx in range(0, len(lines), 4):
            yield lines[idx : idx + 4]

    @classmethod
    def _compress_files(cls, files, to_compress, masked_stats):
        if not to_compress:
            return files

        updated = []
        for it in files:
            if it.kind in to_compress:
                it = it._replace(name=it.name + ".gz")
            elif it.kind == "json":
                if masked_stats is not None:
                    it = masked_stats

            updated.append(it)

        return tuple(updated)

    @classmethod
    def _mask_filenames(cls, data):
        if isinstance(data, dict):
            updates = []
            for key, value in data.items():
                if key == "filenames":
                    updates.append((key, "...[str]"))
                else:
                    new_value = cls._mask_filenames(value)
                    if new_value is not value:
                        updates.append((key, new_value))
        elif isinstance(data, list):
            updates = []
            for idx, value in enumerate(data):
                new_value = cls._mask_filenames(value)
                if new_value is not value:
                    updates.append((idx, new_value))
        else:
            return data

        if updates:
            data = copy.copy(data)
            for idx, value in updates:
                data[idx] = value

        return data


class TestRunner:
    def __init__(
        self,
        test: TestConfig,
        root,
        executable,
        keep_all: bool = False,
    ):
        self._test = test

        self.keep_all = keep_all
        self.executable = executable
        self.name = " / ".join(test.name)

        variant = test.variant if test.variant else ("basic",)
        self.path = os.path.join(root, "_".join(test.name + variant))

        self._test_path = os.path.join(self.path, "test")
        self._exp_path = os.path.join(self.path, "expected")

    @property
    def spec_path(self):
        return self._test.path

    @property
    def skip(self):
        return self._test.skip

    @property
    def command(self):
        return self._test.build_command(self.executable)

    def run(self):
        if self.skip:
            return

        # Create folder containing test input data
        self._setup(self._test_path, _INPUT_FILES)
        # Attempt to execute the command; errors here are probably logic errors
        returncode, stdout, stderr = self._execute()

        try:
            self._evaluate_return_code(returncode, stderr)
            self._evaluate_terminal_output(stdout, stderr)
            self._evaluate_output_files()
        except:
            # Create folder containing reference output data for comparison
            self._setup(self._exp_path, _OUTPUT_FILES)

            # Write observed terminal output
            write_data(os.path.join(self._test_path, "_stdout.txt"), stdout)
            write_data(os.path.join(self._test_path, "_stderr.txt"), stderr)
            write_data(
                os.path.join(self._test_path, "_run.sh"),
                "#!/bin/bash\n\n{}\n".format(cmd_to_s(self.command)),
            )

            raise

        # Cleanup after successful test
        if not self.keep_all:
            shutil.rmtree(self.path)

    def _setup(self, root, keys, verbose=False):
        os.makedirs(root)

        for it in self._test.get_files(keys):
            # Ignored files are not written in order to make diffing simpler
            if it.data is not None:
                write_data(os.path.join(root, it.name), "".join(it.data))

    def _execute(self):
        proc = subprocess.Popen(
            self.command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            preexec_fn=os.setsid,
            cwd=self._test_path,
        )

        try:
            stdout, stderr = proc.communicate(timeout=3)
        except subprocess.TimeoutExpired:
            raise TestError("Test took too long to run")
        finally:
            proc.terminate()

        stdout = stdout.decode("utf-8")
        stderr = stderr.decode("utf-8")

        return proc.returncode, stdout, stderr

    def _evaluate_return_code(self, returncode, stderr):
        if returncode != self._test.return_code:
            raise TestError(
                "Expected return-code {}, but AdapterRemoval returned {}:\n{}".format(
                    self._test.return_code,
                    returncode,
                    pretty_output(filter_errors(stderr), 4),
                )
            )

    def _evaluate_terminal_output(self, stdout, stderr):
        # Output to STDOUT is not normally expected
        if stdout and not self._test.stdout:
            raise TestError("Unexpected output to STDOUT: %r" % (stdout,))

        self._compare_terminal_output("stdout", self._test.stdout, stdout)
        self._compare_terminal_output("stderr", self._test.stderr, stderr)

    @staticmethod
    def _compare_terminal_output(
        pipe,
        expected_lines,
        actual_text,
    ):
        expected_lines = list(expected_lines[::-1])
        actual_lines = actual_text.split("\n")[::-1]

        # Lines are expected to be found in the specified order, but additional
        # lines of output are permitted between required lines
        while expected_lines and actual_lines:
            actual = actual_lines[-1].strip()
            expected = expected_lines[-1].strip()

            if expected in actual:
                expected_lines.pop()
            actual_lines.pop()

        if expected_lines:
            raise TestError("expected {} not found: {}".format(pipe, expected_lines[0]))

    def _evaluate_output_files(self):
        expected_files = set()
        for it in self._test.get_files(_OUTPUT_FILES | _INPUT_FILES):
            expected_files.add(it.name)

            it.compare_with_file(
                expected=os.path.join(self._exp_path, it.name),
                observed=os.path.join(self._test_path, it.name),
            )

        observed_files = set(os.listdir(self._test_path))
        unexpected_files = observed_files - expected_files
        if unexpected_files:
            raise TestError("unexpected files created {!r}".format(unexpected_files))


class TestUpdater:
    def __init__(
        self,
        test: TestConfig,
        executable: str,
    ):
        self.executable = executable
        self.name = " / ".join(test.name)
        self.path = os.path.dirname(test.path)
        self.spec_path = test.path
        self.skip = test.skip or test.return_code or not test.files
        self._test = test

    @property
    def command(self):
        return self._test.build_command(self.executable)

    def run(self):
        if not self.skip:
            self._run()
            self._update_files()

    def _run(self):
        proc = subprocess.Popen(
            self.command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            close_fds=True,
            preexec_fn=os.setsid,
            cwd=self.path,
        )

        try:
            proc.wait(timeout=3)
        except subprocess.TimeoutExpired:
            raise TestError("Test took too long to run")
        finally:
            proc.terminate()

    def _update_files(self):
        for it in self._test.files:
            filename = os.path.join(self.path, it.name)

            if it.kind in ("ignore", "html"):
                try:
                    os.remove(filename)
                except OSError as error:
                    if error.errno != errno.ENOENT:
                        raise
            elif it.kind == "json":
                metadata = [
                    (b'    "version":', b' "...str",\n'),
                    (b'    "command":', b' "...[str]",\n'),
                    (b'    "runtime":', b' "...float"\n'),
                ]

                lines = []
                with open(filename, "rb") as handle:
                    for line in handle:
                        for key, value in list(metadata):
                            if line.startswith(key):
                                metadata.remove((key, value))
                                line = key + value
                                break

                        lines.append(line)

                with open(filename, "wb") as handle:
                    handle.writelines(lines)

            else:
                with open(filename, "rb") as handle:
                    data = handle.read()

                # Some files may intentionally be written compressed
                if data.startswith(b"\x1f\x8b"):
                    with open(filename, "wb") as handle:
                        handle.write(gzip.decompress(data))


################################################################################


def test_name(root, current, filename):
    postfix = current[len(root) :]
    name = [v.strip() for v in postfix.split("/") if v.strip()]
    if filename != "test.json":
        name.append(filename[4:-5].strip("_"))
    elif not name:
        name = [os.path.basename(current)]

    return tuple(name)


def collect_tests(root):
    tests = []
    for dirname, _, filenames in os.walk(root):
        for filename in filenames:
            normalized = filename.lower()

            if normalized.startswith("test") and normalized.endswith(".json"):
                name = test_name(root, dirname, normalized)
                filepath = os.path.join(dirname, filename)

                try:
                    test = TestConfig.load(name=name, filepath=filepath)
                except TestError as error:
                    print_err(
                        "Error while loading test {!r}: {}".format(filepath, error)
                    )
                    sys.exit(1)

                tests.append(test)

    return tests


def build_test_runners(args, tests):
    result = []
    for template in tests:
        for test in TestMutator.create_exhaustive_tests(template, args.exhaustive):
            test = TestRunner(
                test=test,
                root=args.work_dir,
                executable=args.executable,
                keep_all=args.keep_all,
            )

            result.append(test)

    return result


def build_test_updaters(args, tests):
    return [TestUpdater(test, args.executable) for test in tests if not test.skip]


def collect_unused_files(root, tests):
    observed_files = set()
    for dirname, _, filenames in os.walk(root):
        for filename in filenames:
            observed_files.add(os.path.join(dirname, filename))

    recorded_files = set()
    for test in tests:
        recorded_files.add(test.path)
        dirname = os.path.dirname(test.path)
        for it in test.files:
            recorded_files.add(os.path.join(dirname, it.name))

    return observed_files - recorded_files


def run_test(runner: TestRunner):
    try:
        runner.run()
        return runner, None
    except TestError as error:
        return runner, error
    except Exception:
        return runner, TestError(traceback.format_exc())


################################################################################


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("work_dir", help="Directory in which to run test-cases.")
    parser.add_argument("source_dir", help="Directory containing test-cases.")
    parser.add_argument(
        "--executable",
        type=os.path.abspath,
        default=os.path.abspath("./build/adapterremoval3"),
        help="Path to AdapterRemoval executable",
    )
    parser.add_argument(
        "--max-failures",
        type=int,
        default=5,
        help="Maximum number of failures to report. Set to zero or less for no limit",
    )
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="Keep both completed and failed test folders, "
        "not just the folders of failed tests.",
    )
    parser.add_argument(
        "--exhaustive",
        action="store_true",
        help="Perform exhaustive testing by varying input/output formats",
    )
    parser.add_argument(
        "--create-updated-reference",
        action="store_true",
        help="Creates updated reference files by running test commands in the "
        "specified work folder, mirroring the structure of the input folder. "
        "Commands that are expected to fail/not generate output are not run.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=min(4, multiprocessing.cpu_count()),
        help="Size of thread-pool for concurrent execution of tests",
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    print("Testing executable %r" % (args.executable,))
    if not os.path.exists(args.executable):
        print_err("ERROR: Executable does not exist")
        return 1

    if args.max_failures <= 0:
        args.max_failures = float("inf")

    if args.create_updated_reference:
        # test-sets may generate (partially) overlapping files
        args.threads = 1

    if args.create_updated_reference:
        os.makedirs(os.path.dirname(args.work_dir), exist_ok=True)
        shutil.copytree(
            src=args.source_dir,
            dst=args.work_dir,
            symlinks=True,
            dirs_exist_ok=True,
        )

        args.source_dir = args.work_dir
        args.keep_all = True

    print("Reading test-cases from %r" % (args.source_dir,))
    tests = collect_tests(args.source_dir)
    tests.sort(key=lambda it: it.path)
    print("  {:,} tests found".format(len(tests)))

    unused_files = collect_unused_files(args.source_dir, tests)
    if unused_files:
        print_err("Found {} unused files:".format(len(unused_files)))
        for idx, filepath in enumerate(unused_files, start=1):
            print("  {}. {}".format(idx, shlex.quote(filepath)))

        return 1

    if args.create_updated_reference:
        exhaustive_tests = build_test_updaters(args, tests)
    else:
        os.makedirs(args.work_dir, exist_ok=True)
        args.work_dir = tempfile.mkdtemp(dir=args.work_dir)
        exhaustive_tests = build_test_runners(args, tests)

    print("  {:,} test variants generated".format(len(exhaustive_tests)))
    print("Writing test-cases results to %r" % (args.work_dir,))

    n_failures = 0
    n_successes = 0
    n_skipped = 0
    print("\nRunning tests:")

    with multiprocessing.Pool(args.threads) as pool:
        results = pool.imap(run_test, exhaustive_tests)
        grouped_results = groupby(results, lambda it: it[0].name)

        for idx, (name, test_results) in enumerate(grouped_results, start=1):
            print("  %i of %i: %s " % (idx, len(tests), name), end="")

            last_test = None
            last_error = None
            test_skipped = False
            for test, error in test_results:
                if error is not None:
                    print_err("X", end="")
                    if last_error is None:
                        last_error = error
                        last_test = test
                elif test.skip:
                    print_warn(".", end="")
                    test_skipped = True
                else:
                    print_ok(".", end="")

                sys.stdout.flush()

            if test_skipped:
                print_warn(" [SKIPPED]")
                n_skipped += 1
            elif last_error is None:
                n_successes += 1
                print_ok(" [OK]")
            else:
                n_failures += 1
                assert last_test is not None
                print_err("\nTest {} failed:".format(last_test.name))
                print_err("  Specification = {}".format(last_test.spec_path))
                print_err("  Directory     = {}".format(last_test.path))
                print_err("  Command       = {}".format(cmd_to_s(last_test.command)))

                error = "\n                  ".join(str(last_error).split("\n"))
                print_err("  Error         = {}".format(error))

                if n_failures >= args.max_failures:
                    pool.terminate()
                    break

    if not (n_failures or args.create_updated_reference):
        os.rmdir(args.work_dir)

    if n_failures >= args.max_failures:
        print_err(
            "\nAborted after %i of %i tests failed."
            % (n_failures, n_failures + n_successes + n_skipped)
        )
    elif n_failures:
        print_err("\n%i of %i tests failed." % (n_failures, len(exhaustive_tests)))
    else:
        print_ok("\nAll %i tests succeeded." % (len(exhaustive_tests),))

    if n_skipped:
        print_warn("Skipped {} tests!".format(n_skipped))

    return 1 if n_failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
