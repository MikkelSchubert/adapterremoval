#!/usr/bin/env python3
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
from __future__ import annotations

import argparse
import difflib
import errno
import gzip
import itertools
import json
import math
import multiprocessing
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import traceback
from dataclasses import dataclass
from itertools import groupby, islice
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    NoReturn,
    Sequence,
    Union,
    cast,
    overload,
)

try:
    import jsonschema
except ImportError:
    JSON_SCHEMA_VALIDATION = False  # pyright: ignore[reportConstantRedefinition]
else:
    JSON_SCHEMA_VALIDATION = True


if TYPE_CHECKING:
    from typing_extensions import Literal, TypeAlias

    JSON: TypeAlias = Union[
        Dict[str, "JSON"], List["JSON"], str, int, float, bool, None
    ]

#############################################################################
_COLORS_ENABLED = not os.getenv("NO_COLOR")


def _do_print_color(
    *vargs_: object,
    colorcode: int,
    end: str,
) -> None:
    """Utility function: Prints using shell colors."""
    destination = sys.stdout

    # No colors if output is redirected (e.g. less, file, etc.)
    if _COLORS_ENABLED and destination.isatty():
        vargs = list(vargs_)
        for index, varg in enumerate(vargs):
            varg_lines: list[str] = []
            # Newlines terminate the color-code for e.g. 'less', so ensure that
            # each line is color-coded, while preserving the list of arguments
            for line in str(varg).split("\n"):
                varg_lines.append("\033[00;%im%s\033[00m" % (colorcode, line))
            vargs[index] = "\n".join(varg_lines)
    else:
        vargs = vargs_

    print(*vargs, file=destination, end=end)
    if "\n" in end:
        destination.flush()


def print_ok(*vargs: object, end: str = "\n") -> None:
    """Equivalent to print, but prints using shell colorcodes (green)."""
    _do_print_color(*vargs, colorcode=32, end=end)


def print_warn(*vargs: object, end: str = "\n") -> None:
    """Equivalent to print, but prints using shell colorcodes (green)."""
    _do_print_color(*vargs, colorcode=33, end=end)


def print_err(*vargs: object, end: str = "\n") -> None:
    """Equivalent to print, but prints using shell colorcodes (red)."""
    _do_print_color(*vargs, colorcode=31, end=end)


#############################################################################


@overload
def read_file(filename: Path, mode: Literal["rt"] = "rt") -> str: ...


@overload
def read_file(filename: Path, mode: Literal["rb"] = "rb") -> bytes: ...


def read_file(filename: Path, mode: Literal["rb", "rt"] = "rt") -> str | bytes:
    try:
        with filename.open(mode) as handle:
            return handle.read()
    except OSError as error:
        raise TestError(f"ERROR while reading data:\n    {error}") from error


def read_json(filename: Path) -> tuple[str, JSON]:
    text = read_file(filename)

    try:
        return text, json.loads(text)
    except json.JSONDecodeError as error:
        raise TestError(f"ERROR while reading {quote(filename)}:") from error


def read_and_decompress_file(filename: Path) -> str:
    value = read_file(filename, "rb")

    header = value[:2]
    if header == b"\x1f\x8b":
        if not filename.suffix == ".gz":
            raise TestError(f"{filename} is gzipped, but lacks .gz extension")

        value = gzip.decompress(value)
    elif filename.suffix == ".gz":
        raise TestError(f"{filename} has .gz extension, but is not compressed")

    return value.decode("utf-8")


def escape_whitespace(line: str) -> str:
    out: list[str] = []
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
            escaped = f"\\{ord(char):o}" if char.isspace() else char

        out.append(escaped)

    return "".join(out)


def truncate_lines(lines: Iterable[str], max_lines: int) -> list[str]:
    lines = list(islice(lines, max_lines + 1))
    if len(lines) > max_lines:
        lines[-1] = "..."

    return lines


def pretty_output(lines: Iterable[str], max_lines: int, padding: int = 0) -> str:
    result: list[str] = []
    prefix = " " * padding
    lines = truncate_lines(lines, max_lines)
    for line in lines:
        result.append(f"{prefix}>  {escape_whitespace(line)}")

    return "\n".join(result)


def pretty_returncode(returncode: int) -> str:
    if returncode >= 0:
        return str(returncode)

    try:
        return f"{returncode} ({signal.Signals(-returncode).name})"
    except ValueError:
        return str(returncode)


def filter_errors(
    text: str,
    _re: re.Pattern[str] = re.compile(r"(\[(?:ERROR|WARNING)\] .*)"),
) -> list[str]:
    lines = text.split("\n")
    errors: list[str] = []
    for line in lines:
        match = _re.search(line)
        if match:
            errors.append(match.group())

    return errors if errors else lines


def write_data(path: Path, data: str) -> None:
    open_ = gzip.open if path.suffix == ".gz" else open
    with open_(path, "wt") as handle:
        handle.write(data)


def quote(path: Path | str) -> str:
    if isinstance(path, Path):
        path = str(path)

    return shlex.quote(path)


def cmd_to_s(cmd: list[str | Path]) -> str:
    return " ".join(quote(field) for field in cmd)


def path_to_s(path: Iterable[str]) -> str:
    return ".".join(str(value) for value in path)


def classname(v: object) -> str:
    if v is None:
        return "None"

    return v.__class__.__name__


################################################################################

# Wildcards with (limited) type checking
JSON_WILDCARDS: dict[str, Callable[[object], bool]] = {
    "...": lambda _: True,
    "...str": lambda it: isinstance(it, str),
    "...float": lambda it: isinstance(it, float) and not math.isnan(it),
    "...[str]": lambda it: isinstance(it, list)
    and all(isinstance(v, str) for v in cast(List[object], it)),
}

# JSON path elements that do not require quoting
JSON_PATH = re.compile("[a-z0-9_]", re.I)


def diff_json(
    reference: JSON,
    observed: JSON,
    path: tuple[str, ...] = ("",),
) -> Iterator[str]:
    def _err(what: str, ref: object, obs: object) -> str:
        return f"{what} at {path_to_s(path)}: {obs} is not {ref}"

    if reference == observed:
        return

    if not isinstance(reference, (dict, list)) and reference in JSON_WILDCARDS:
        if not JSON_WILDCARDS[reference](observed):
            yield _err("wildcard mismatch", repr(reference), repr(observed))
    elif type(reference) != type(observed):
        yield _err("type mismatch", classname(reference), classname(observed))
    elif isinstance(reference, list) and isinstance(observed, list):
        if len(reference) != len(observed):
            yield _err("length mismatch", len(reference), len(observed))

        for idx, (ref, obs) in enumerate(zip(reference, observed)):
            yield from diff_json(ref, obs, (*path, f"[{idx}]"))
    elif isinstance(reference, dict):
        assert isinstance(observed, dict)

        for key in reference.keys() | observed.keys():
            if key not in reference:
                yield _err("unexpected key", "in reference", key)
            elif key not in observed:
                yield _err("unexpected key", "in observation", key)
            else:
                if not JSON_PATH.match(key):
                    key = repr(key)

                yield from diff_json(reference[key], observed[key], (*path, key))
    elif isinstance(reference, float):
        assert isinstance(observed, float)

        # There should be no NANs in this data, so this test also catches those
        if not (abs(reference - observed) < 1e-6):
            yield _err("mismatch", repr(observed), repr(reference))
    else:
        yield _err("mismatch", repr(reference), repr(observed))


################################################################################


class TestError(RuntimeError):
    pass


class JSONSchemaError(TestError):
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


class JSONValidator:
    _schema: JSON

    __slots__ = [
        "_schema",
    ]

    def __init__(self, filepath: Path | None) -> None:
        if filepath is None:
            self._schema = {}
        elif not JSON_SCHEMA_VALIDATION:
            print_warn("JSON schema provided, but `jsonschema` could not be loaded")

            self._schema = {}
        else:
            print(f"Reading JSON schema from {quote(filepath)}")
            _, self._schema = read_json(filepath)

    def __call__(self, data: JSON) -> None:
        if self._schema:
            try:
                jsonschema.validate(instance=data, schema=self._schema)
            except jsonschema.ValidationError as error:
                raise JSONSchemaError("Invalid JSON file") from error

    def __bool__(self) -> bool:
        return bool(self._schema)


def json_pop_value(data: JSON, path: tuple[str, ...], default: JSON) -> JSON:
    if not isinstance(data, dict):
        raise TestError("test specification did not contain a dict")

    return data.pop(path[-1], default)


def json_pop_dict(data: JSON, path: tuple[str, ...]) -> dict[str, JSON]:
    value = json_pop_value(data, path, {})
    if isinstance(value, dict):
        return value

    raise TestError(f"expected dict at {path_to_s(path)}, found {value!r}")


def json_pop_tuple_of_str(
    data: JSON,
    path: tuple[str, ...],
    default: JSON,
) -> tuple[str, ...]:
    value = json_pop_value(data, path, default)
    if isinstance(value, list) and all(isinstance(it, str) for it in value):
        return tuple(cast(List[str], value))

    raise TestError(f"expected string list at {path_to_s(path)}, found {value!r}")


def json_pop_int(data: JSON, path: tuple[str, ...], *, default: JSON) -> int:
    value = json_pop_value(data, path, default)
    if isinstance(value, int):
        return value

    raise TestError(f"expected int at {path_to_s(path)}, found {value!r}")


def json_pop_bool(data: JSON, path: tuple[str, ...], *, default: JSON) -> bool:
    value = json_pop_value(data, path, default)
    if isinstance(value, bool):
        return value

    raise TestError(f"expected bool at {path_to_s(path)}, found {value!r}")


@dataclass
class TestFile:
    name: str
    kind: str

    def compare_with_file(self, expected: Path, observed: Path) -> None:
        raise NotImplementedError

    @classmethod
    def parse(
        cls,
        *,
        root: Path,
        name: str,
        kind: str,
        json_validator: JSONValidator | None,
    ) -> TestFile:
        # HTML files are not compared in detail
        if kind in ("html", "ignore"):
            return TestMiscFile(name=name, kind=kind)

        filename = root / name
        if kind == "json":
            text, data = read_json(filename)
            return TestJsonFile(
                name=name,
                kind=kind,
                text=text,
                data=data,
                json_validator=json_validator,
            )

        # Reference data is not compresses, regardless of extension
        return TestTextFile(name=name, kind=kind, text=read_file(filename))

    def _raise_test_error(
        self,
        label: str,
        expected: Path,
        observed: Path,
        differences: str,
    ) -> NoReturn:
        raise TestError(
            f"Mismatches in {label} file:\n"
            f"  Expected   = {quote(expected)}\n"
            f"  Observed   = {quote(observed)}\n"
            f"  Mismatches =\n{differences}"
        )


@dataclass
class TestTextFile(TestFile):
    text: str

    def compare_with_file(self, expected: Path, observed: Path) -> None:
        if not observed.is_file():
            raise TestError(f"file {observed} not created")

        text = read_and_decompress_file(observed)
        if text == self.text:
            return

        lines = difflib.unified_diff(
            self.text.splitlines(keepends=True),
            text.splitlines(keepends=True),
            "expected",
            "observed",
        )

        self._raise_test_error(
            label="output",
            expected=expected,
            observed=observed,
            differences=pretty_output(lines, 10),
        )


@dataclass
class TestJsonFile(TestFile):
    text: str
    data: JSON
    json_validator: JSONValidator | None

    def compare_with_file(self, expected: Path, observed: Path) -> None:
        _, data = read_json(observed)
        differences = diff_json(reference=self.data, observed=data)
        differences = truncate_lines(differences, 4)

        if differences:
            differences = "\n".join(
                f"    {idx}. {line}" for idx, line in enumerate(differences, start=1)
            )

            self._raise_test_error(
                label="JSON",
                expected=expected,
                observed=observed,
                differences=differences,
            )

        if self.json_validator is not None:
            self.json_validator(data)

    def mask_filenames(self) -> TestJsonFile:
        def _mask_filenames(data: JSON) -> JSON:
            if isinstance(data, dict):
                data = dict(data)
                for key, value in data.items():
                    if key == "filenames":
                        data[key] = "...[str]"
                    else:
                        new_value = _mask_filenames(value)
                        if new_value is not value:
                            data[key] = new_value
            elif isinstance(data, list):
                data = [_mask_filenames(value) for value in data]

            return data

        return TestJsonFile(
            name=self.name,
            kind=self.kind,
            text=self.text,
            data=_mask_filenames(self.data),
            json_validator=self.json_validator,
        )


class TestMiscFile(TestFile):
    def compare_with_file(self, expected: Path, observed: Path) -> None: ...


class TestConfig(NamedTuple):
    path: Path
    name: tuple[str, ...]
    variant: tuple[str, ...]
    skip: bool
    exhaustive: bool | None
    arguments: tuple[str, ...]
    return_code: int
    stdout: tuple[str, ...]
    stderr: tuple[str, ...]
    files: tuple[TestFile, ...]

    @classmethod
    def load(
        cls,
        name: tuple[str, ...],
        filepath: Path,
        json_validator: JSONValidator | None,
    ) -> TestConfig:
        root = filepath.parent
        _, data = read_json(filepath)

        files = json_pop_dict(data, ("files",))
        test_files: list[TestFile] = []
        for key in _TEST_FILES:
            for filename in json_pop_tuple_of_str(files, ("files", key), []):
                test_files.append(
                    TestFile.parse(
                        root=root,
                        name=filename,
                        kind=key,
                        json_validator=json_validator,
                    )
                )

        if files:
            raise TestError(f"Unexpected files: {files}")

        self = TestConfig(
            path=filepath,
            name=name,
            variant=(),
            arguments=json_pop_tuple_of_str(data, ("arguments",), default=[]),
            skip=json_pop_bool(data, ("skip",), default=False),
            exhaustive=json_pop_bool(data, ("exhaustive",), default=False),
            return_code=json_pop_int(data, ("return_code",), default=0),
            stdout=json_pop_tuple_of_str(data, ("stdout",), default=[]),
            stderr=json_pop_tuple_of_str(data, ("stderr",), default=[]),
            files=tuple(test_files),
        )

        if data:
            raise TestError(f"Unexpected JSON values: {data}")

        return self

    def build_command(self, executable: Path) -> list[str | Path]:
        command: list[str | Path] = [executable, *self.arguments]
        input_1: list[str] = []
        input_2: list[str] = []

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
            command += ["--file1", *input_1]

        if input_2:
            command += ["--file2", *input_2]

        return command

    def get_files(self, keys: Iterable[str]) -> Iterator[TestFile]:
        if set(keys) - _TEST_FILES:
            raise ValueError(keys)

        for item in self.files:
            if item.kind in keys:
                yield item


class TestMutator:
    @classmethod
    def create_exhaustive_tests(
        cls,
        it: TestConfig,
        *,
        exhaustive: bool,
    ) -> Iterator[TestConfig]:
        if it.exhaustive is not None:
            exhaustive = it.exhaustive

        if not (exhaustive and it.files):
            yield it
            return

        # For simplicity's sake, filenames are not checked in exhaustive tests
        for file in it.files:
            if isinstance(file, TestJsonFile):
                masked_stats = file.mask_filenames()
                break
        else:
            masked_stats = None

        # Test with split or interleaved input (if possible)
        for it in cls._interleave_input(it, masked_stats):
            # Test with gzip compression of input and output files
            for compress_in, compress_out in itertools.product((False, True), repeat=2):
                variant: list[str] = []
                arguments: list[str] = []
                compressed_files: set[str] = set()

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
    def _interleave_input(
        cls,
        test: TestConfig,
        masked_stats: TestFile | None,
    ) -> Iterator[TestConfig]:
        yield test

        input_files: dict[str, list[list[str]]] = {}
        other_files: list[TestFile] = []
        for it in test.files:
            if it.kind in _INPUT_FILES_FQ:
                assert isinstance(it, TestTextFile)

                # Ensure exactly one trailing newline
                lines = it.text.splitlines(keepends=True)
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

                    data: list[str] = []
                    for record_1, record_2 in zip(records_1, records_2):
                        data.extend(record_1)
                        data.extend(record_2)

                    other_files.append(
                        TestTextFile(
                            name=f"input_{idx}.fastq",
                            kind="input_1",
                            text="".join(data),
                        )
                    )

                yield test._replace(
                    variant=(*test.variant, "intl"),
                    files=other_files,
                    arguments=(*test.arguments, "--interleaved-input"),
                )

    @staticmethod
    def _lines_to_records(lines: list[str]) -> Iterator[list[str]]:
        assert len(lines) % 4 == 0, lines
        for idx in range(0, len(lines), 4):
            yield lines[idx : idx + 4]

    @classmethod
    def _compress_files(
        cls,
        files: tuple[TestFile, ...],
        to_compress: set[str],
        masked_stats: TestFile | None,
    ) -> tuple[TestFile, ...]:
        if not to_compress:
            return files

        updated: list[TestFile] = []
        for it in files:
            if it.kind in to_compress:
                assert isinstance(it, TestTextFile)

                # Tests may manually specify gzip output
                if not it.name.endswith(".gz"):
                    it = TestTextFile(name=f"{it.name}.gz", kind=it.kind, text=it.text)
            elif it.kind == "json" and masked_stats is not None:
                it = masked_stats

            updated.append(it)

        return tuple(updated)


class TestRunner:
    def __init__(
        self,
        *,
        test: TestConfig,
        root: Path,
        executable: Path,
        keep_all: bool = False,
    ) -> None:
        self._test = test

        self.keep_all = keep_all
        self.executable = executable
        self.name = " / ".join(test.name)

        variant = test.variant if test.variant else ("basic",)
        self.path = root / "_".join(test.name + variant)

        self._test_path = self.path / "test"
        self._exp_path = self.path / "expected"

    @property
    def spec_path(self) -> Path:
        return self._test.path

    @property
    def skip(self) -> bool:
        return self._test.skip

    @property
    def command(self) -> list[str | Path]:
        return self._test.build_command(self.executable)

    def run(self) -> None:
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
            write_data(self._test_path / "_stdout.txt", stdout)
            write_data(self._test_path / "_stderr.txt", stderr)
            write_data(
                self._test_path / "_run.sh",
                f"#!/bin/bash\n\n{cmd_to_s(self.command)}\n",
            )

            raise

        # Cleanup after successful test
        if not self.keep_all:
            shutil.rmtree(self.path)

    def _setup(self, root: Path, keys: Iterable[str]) -> None:
        root.mkdir(parents=True)

        for it in self._test.get_files(keys):
            # Ignored files are not written in order to make diffing simpler
            if isinstance(it, (TestTextFile, TestJsonFile)):
                write_data(root / it.name, it.text)

    def _execute(self) -> tuple[int, str, str]:
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
        except subprocess.TimeoutExpired as error:
            raise TestError("Test took too long to run") from error
        finally:
            proc.terminate()

        stdout = stdout.decode("utf-8")
        stderr = stderr.decode("utf-8")

        return proc.returncode, stdout, stderr

    def _evaluate_return_code(self, returncode: int, stderr: str) -> None:
        if returncode != self._test.return_code:
            raise TestError(
                f"Expected return-code {self._test.return_code}, but AdapterRemoval "
                f"returned {pretty_returncode(returncode)}:\n"
                f"{pretty_output(filter_errors(stderr), 4)}"
            )

    def _evaluate_terminal_output(self, stdout: str, stderr: str) -> None:
        # Output to STDOUT is not normally expected
        if stdout and not self._test.stdout:
            raise TestError(f"Unexpected output to STDOUT: {stdout!r}")

        self._compare_terminal_output("stdout", self._test.stdout, stdout)
        self._compare_terminal_output("stderr", self._test.stderr, stderr)

    @staticmethod
    def _compare_terminal_output(
        pipe: str,
        expected_lines: Sequence[str],
        actual_text: str,
    ) -> None:
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
            raise TestError(f"expected {pipe} not found: {expected_lines[0]}")

    def _evaluate_output_files(self) -> None:
        expected_files: set[str] = set()
        for it in self._test.get_files(_OUTPUT_FILES | _INPUT_FILES):
            expected_files.add(it.name)

            it.compare_with_file(
                expected=self._exp_path / it.name,
                observed=self._test_path / it.name,
            )

        observed_files = set(os.listdir(self._test_path))
        unexpected_files = observed_files - expected_files
        if unexpected_files:
            raise TestError(f"unexpected files created {unexpected_files!r}")


class TestUpdater:
    def __init__(
        self,
        test: TestConfig,
        executable: Path,
    ) -> None:
        self.executable = executable
        self.name = " / ".join(test.name)
        self.path = test.path.parent
        self.spec_path = test.path
        self.skip = test.skip or test.return_code or not test.files
        self._test = test

    @property
    def command(self) -> list[str | Path]:
        return self._test.build_command(self.executable)

    def run(self) -> None:
        if not self.skip:
            self._run()
            self._update_files()

    def _run(self) -> None:
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
        except subprocess.TimeoutExpired as error:
            raise TestError("Test took too long to run") from error
        finally:
            proc.terminate()

    def _update_files(self) -> None:
        for it in self._test.files:
            filename = self.path / it.name

            if it.kind in ("ignore", "html"):
                try:
                    filename.unlink()
                except OSError as error:
                    if error.errno != errno.ENOENT:
                        raise
            elif it.kind == "json":
                metadata: list[tuple[bytes, bytes]] = [
                    (b'    "version":', b' "...str",\n'),
                    (b'    "command":', b' "...[str]",\n'),
                    (b'    "runtime":', b' "...float",\n'),
                    (b'    "timestamp":', b' "...str"\n'),
                ]

                lines: list[bytes] = []
                with filename.open("rb") as handle:
                    for line in handle:
                        for key, value in list(metadata):
                            if line.startswith(key):
                                metadata.remove((key, value))
                                line = key + value
                                break

                        lines.append(line)

                with filename.open("wb") as handle:
                    handle.writelines(lines)

            else:
                # Some files may intentionally be written compressed
                data = read_file(filename, "rb")
                if data.startswith(b"\x1f\x8b"):
                    with filename.open("wb") as handle:
                        handle.write(gzip.decompress(data))


################################################################################


def test_name(root: Path, current: Path, filename: str) -> tuple[str, ...]:
    root_p = root.parts
    current_p = current.parts
    if len(root_p) > len(current_p) or current_p[: len(root_p)] != root_p:
        raise ValueError((root, current))

    name = list(current_p[len(root_p) :])

    if filename != "test.json":
        name.append(filename[4:-5].strip("_"))
    elif not name:
        name = [current.name]

    return tuple(name)


def collect_tests(
    *,
    root: Path,
    json_validator: JSONValidator | None,
) -> list[TestConfig]:
    tests: list[TestConfig] = []
    for dirname, _, filenames in os.walk(root):
        dirname = Path(dirname)
        for filename in filenames:
            normalized = filename.lower()

            if normalized.startswith("test") and normalized.endswith(".json"):
                name = test_name(root, dirname, normalized)
                filepath = dirname / filename

                try:
                    test = TestConfig.load(
                        name=name,
                        filepath=filepath,
                        json_validator=json_validator,
                    )
                except TestError as error:
                    print_err(f"Error while loading test {quote(filepath)}: {error}")
                    sys.exit(1)

                tests.append(test)

    return tests


def build_test_runners(args: Args, tests: list[TestConfig]) -> list[TestRunner]:
    result: list[TestRunner] = []
    for template in tests:
        for test in TestMutator.create_exhaustive_tests(
            template, exhaustive=args.exhaustive
        ):
            test = TestRunner(
                test=test,
                root=args.work_dir,
                executable=args.executable,
                keep_all=args.keep_all,
            )

            result.append(test)

    return result


def build_test_updaters(args: Args, tests: list[TestConfig]) -> list[TestUpdater]:
    return [TestUpdater(test, args.executable) for test in tests if not test.skip]


def collect_unused_files(root: Path, tests: list[TestConfig]) -> set[Path]:
    observed_files: set[Path] = set()
    for dirname, _, filenames in os.walk(root):
        for filename in filenames:
            observed_files.add(Path(dirname) / filename)

    recorded_files: set[Path] = set()
    for test in tests:
        recorded_files.add(test.path)
        dirname = test.path.parent
        for it in test.files:
            recorded_files.add(dirname / it.name)

    return observed_files - recorded_files


def run_test(
    runner: TestRunner | TestUpdater,
) -> tuple[TestRunner | TestUpdater, Exception | None]:
    try:
        runner.run()
    except TestError as error:
        return runner, error
    except Exception:  # noqa: BLE001
        return runner, TestError(traceback.format_exc())
    else:
        return runner, None


################################################################################


class Args(argparse.Namespace):
    work_dir: Path
    source_dir: Path
    executable: Path
    max_failures: float
    keep_all: bool
    exhaustive: bool
    json_schema: Path | None
    schema_validation_required: bool
    create_updated_reference: bool
    threads: int


def parse_args(argv: list[str]) -> Args:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "work_dir",
        type=Path,
        help="Directory in which to run test-cases",
    )
    parser.add_argument(
        "source_dir",
        type=Path,
        help="Directory containing test-cases",
    )
    parser.add_argument(
        "--executable",
        type=lambda value: Path(value).absolute(),
        default=Path("./build/adapterremoval3").absolute(),
        help="Path to AdapterRemoval executable",
    )
    parser.add_argument(
        "--max-failures",
        type=int,
        default=0,
        help="Maximum number of failures to report. Set to zero or less for no limit",
    )
    parser.add_argument(
        "--keep-all",
        default=False,
        action="store_true",
        help="Keep both completed and failed test folders, "
        "not just the folders of failed tests.",
    )
    parser.add_argument(
        "--exhaustive",
        default=False,
        action="store_true",
        help="Perform exhaustive testing by varying input/output formats",
    )
    parser.add_argument(
        "--json-schema",
        type=Path,
        help="Validate output JSON files using the supplied schema; requires Python "
        "module `jsonschema`; validating the schema is only performed if the "
        "required is installed, otherwise it is skipped",
    )
    parser.add_argument(
        "--schema-validation-required",
        default=False,
        action="store_true",
        help="Terminate if JSON schema validation cannot be performed",
    )
    parser.add_argument(
        "--create-updated-reference",
        default=False,
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

    return parser.parse_args(argv, namespace=Args())


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    print(f"Testing executable {quote(args.executable)}")
    if not args.executable.is_file():
        print_err("ERROR: Executable does not exist")
        return 1
    elif args.schema_validation_required and not args.json_schema:
        print_err("ERROR: No JSON schema was provided")
        return 1
    elif args.schema_validation_required and not JSON_SCHEMA_VALIDATION:
        print_err("ERROR: Required module `jsonschema` could not be loaded")
        return 1

    if args.max_failures <= 0:
        args.max_failures = float("inf")

    if args.create_updated_reference:
        # test-sets may generate (partially) overlapping files
        args.threads = 1

    if args.create_updated_reference:
        # `dirs_exist_ok` for `copytree` requires Python >= 3.8
        if args.work_dir.exists():
            try:
                args.work_dir.rmdir()
            except OSError as error:
                print_err(f"ERROR: Could not remove existing work-dir: {error}")
                return 1

        args.work_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(
            src=args.source_dir,
            dst=args.work_dir,
            symlinks=True,
        )

        args.source_dir = args.work_dir
        args.keep_all = True

    json_validator = JSONValidator(args.json_schema)

    print(f"Reading test-cases from {quote(args.source_dir)}")
    tests = collect_tests(root=args.source_dir, json_validator=json_validator)
    tests.sort(key=lambda it: it.path)
    print(f"  {len(tests):,} tests found")

    unused_files = collect_unused_files(args.source_dir, tests)
    if unused_files:
        print_err(f"Found {len(unused_files)} unused files:")
        for idx, filepath in enumerate(unused_files, start=1):
            print(f"  {idx}. {quote(filepath)}")

        return 1

    if args.create_updated_reference:
        exhaustive_tests = build_test_updaters(args, tests)
    else:
        args.work_dir.mkdir(parents=True, exist_ok=True)
        args.work_dir = Path(tempfile.mkdtemp(dir=args.work_dir))
        exhaustive_tests = build_test_runners(args, tests)

    print(f"  {len(exhaustive_tests):,} test variants generated")
    print(f"Writing test-cases results to {quote(args.work_dir)}")

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
                print_err(f"\nTest {last_test.name} failed:")
                print_err(f"  Specification = {last_test.spec_path}")
                print_err(f"  Directory     = {last_test.path}")
                print_err(f"  Command       = {cmd_to_s(last_test.command)}")

                error = "\n                  ".join(str(last_error).split("\n"))
                print_err(f"  Error         = {error}")

                if n_failures >= args.max_failures:
                    pool.terminate()
                    break

    if not (n_failures or args.create_updated_reference):
        args.work_dir.rmdir()

    if n_failures >= args.max_failures:
        print_err(
            "\nAborted after %i of %i tests failed."
            % (n_failures, n_failures + n_successes + n_skipped)
        )
    elif n_failures:
        print_err("\n%i of %i tests failed." % (n_failures, len(exhaustive_tests)))
    else:
        print_ok("\nAll %i tests succeeded." % (len(exhaustive_tests),))

    if json_validator:
        print_ok("JSON file validation performed.")

    if n_skipped:
        print_warn(f"Skipped {n_skipped} tests!")

    if not json_validator:
        print_warn("JSON files were not validated!")

    return 1 if n_failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
