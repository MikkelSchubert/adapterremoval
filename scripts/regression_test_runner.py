#!/usr/bin/env python3
# /// script
# requires-python = ">=3.7"
# dependencies = [
#     "fastjsonschema==2.21.1",
# ]
# [tool.uv]
# exclude-newer = "2025-02-26T00:00:00Z"
# ///
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import argparse
import difflib
import errno
import gzip
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
    Any,
    Callable,
    ClassVar,
    Dict,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    NoReturn,
    Sequence,
    Tuple,
    Union,
    cast,
    overload,
)

try:
    import fastjsonschema
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
_COLORS_ENABLED = False


def _do_print_color(
    *vargs_: object,
    color_code: int,
    end: str,
) -> None:
    """Utility function: Prints using shell colors."""
    destination = sys.stdout

    # No colors if output is redirected (e.g. less, file, etc.)
    if _COLORS_ENABLED:
        vargs = list(vargs_)
        for index, varg in enumerate(vargs):
            varg_lines: list[str] = []
            # Newlines terminate the color-code for e.g. 'less', so ensure that
            # each line is color-coded, while preserving the list of arguments
            for line in str(varg).split("\n"):
                varg_lines.append("\033[00;%im%s\033[00m" % (color_code, line))
            vargs[index] = "\n".join(varg_lines)
    else:
        vargs = vargs_

    print(*vargs, file=destination, end=end)
    if "\n" in end:
        destination.flush()


def print_ok(*vargs: object, end: str = "\n") -> None:
    """Equivalent to print, but prints using shell color codes (green)."""
    _do_print_color(*vargs, color_code=32, end=end)


def print_warn(*vargs: object, end: str = "\n") -> None:
    """Equivalent to print, but prints using shell color codes (green)."""
    _do_print_color(*vargs, color_code=33, end=end)


def print_err(*vargs: object, end: str = "\n") -> None:
    """Equivalent to print, but prints using shell color codes (red)."""
    _do_print_color(*vargs, color_code=31, end=end)


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
        raise TestError(f"ERROR while reading data: {error}") from error


def read_json(filename: Path) -> tuple[str, JSON]:
    text = read_file(filename)

    try:
        return text, json.loads(text)
    except json.JSONDecodeError as error:
        raise TestError(f"ERROR while reading {quote(filename)}: {error}") from error


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


def truncate_lines(
    lines: Iterable[str],
    *,
    max_lines: int,
    tail: bool = False,
) -> list[str]:
    if tail:
        lines = list(lines)[-(max_lines + 1) :]
        if len(lines) > max_lines:
            lines[0] = "..."
    else:
        lines = list(islice(lines, max_lines + 1))
        if len(lines) > max_lines:
            lines[-1] = "..."

    return lines


def pretty_output(
    lines: Iterable[str],
    *,
    max_lines: int,
    padding: int = 0,
    tail: bool = False,
) -> str:
    result: list[str] = []
    prefix = " " * padding
    lines = truncate_lines(lines, max_lines=max_lines, tail=tail)
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
    lines = text.rstrip().split("\n")
    errors: list[str] = []
    for line in lines:
        match = _re.search(line)
        if match:
            errors.append(match.group())

    return errors if errors else lines


def print_long_err(
    prefix: str,
    error: str,
    *,
    max_lines: int | None = None,
    tail: bool = False,
):
    lines = filter_errors(error)
    if max_lines is not None:
        lines = truncate_lines(lines, max_lines=max_lines, tail=tail)

    for line in lines:
        print_err(prefix, line)
        prefix = " " * len(prefix)


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


def simplify_cmd(cmd: list[str | Path]) -> str:
    cmd = list(cmd)
    if cmd:
        cmd[0] = Path(cmd[0]).name

    return " ".join(quote(field) for field in cmd)


def path_to_s(path: Iterable[str | int]) -> str:
    return ".".join(str(value) for value in path)


def relative_to(root: Path, path: Path) -> Path:
    if path.is_absolute():
        try:
            return path.relative_to(root)
        except ValueError:
            pass

    return path


def classname(v: object) -> str:
    if v is None:
        return "None"

    return v.__class__.__name__


################################################################################


def _is_str_list(value: object) -> bool:
    return isinstance(value, list) and all(
        isinstance(v, str) for v in cast(List[object], value)
    )


# Wildcards with (limited) type checking
JSON_WILDCARDS: dict[str, Callable[[object], bool]] = {
    "...": lambda _: True,
    "...str": lambda it: isinstance(it, str),
    "...float": lambda it: isinstance(it, float) and not math.isnan(it),
    "...[str]": _is_str_list,
    "...[str] | None": lambda it: it is None or _is_str_list(it),
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
    elif type(reference) is not type(observed):
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
                key_s = key if JSON_PATH.match(key) else repr(key)

                yield from diff_json(reference[key], observed[key], (*path, key_s))
    elif isinstance(reference, float):
        assert isinstance(observed, float)

        # There should be no NANs in this data, so this test also catches those
        if not (abs(reference - observed) < 1e-6):
            yield _err("mismatch", repr(observed), repr(reference))
    else:
        yield _err("mismatch", repr(reference), repr(observed))


################################################################################


class TestError(RuntimeError):
    def __init__(self, message: str, *, output: str | None = None) -> None:
        super().__init__(message)
        self.output = output


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

_IGNORED_FILES = {
    "meson.build",
}


class JSONValidator:
    # Validator is stored a a class-var, to support multiprocessing (pickling)
    _validator: ClassVar[Any] = None

    @classmethod
    def load_schema(cls, filepath: Path) -> None:
        _, schema = read_json(filepath)

        cls._validator = fastjsonschema.compile(schema, use_default=False)

    @classmethod
    def validate(cls, data: JSON) -> None:
        if cls._validator is not None:
            try:
                cls._validator(data)
            except fastjsonschema.JsonSchemaValueException as error:
                error_path = path_to_s(error.path)
                raise JSONSchemaError(
                    f"Schema error at {error_path}: {error.message}"
                ) from error
            except fastjsonschema.JsonSchemaException as error:
                raise JSONSchemaError(str(error)) from error

    @classmethod
    def can_validate(cls) -> bool:
        return cls._validator is not None


def json_pop_value(data: JSON, path: tuple[str, ...], default: JSON) -> JSON:
    if not isinstance(data, dict):
        raise TestError("test specification did not contain a dict")

    return data.pop(path[-1], default)


def json_pop_dict(data: JSON, path: tuple[str, ...]) -> dict[str, JSON]:
    value = json_pop_value(data, path, {})
    if isinstance(value, dict):
        return value

    raise TestError(f"expected dict at {path_to_s(path)}, found {value!r}")


def json_pop_optional_str(
    data: JSON,
    path: tuple[str, ...],
    *,
    default: JSON = None,
) -> str | None:
    value = json_pop_value(data, path, default)
    if value is None or isinstance(value, str):
        return value

    raise TestError(f"expected string or null at {path_to_s(path)}, found {value!r}")


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
    def parse(cls, *, root: Path, name: str, kind: str) -> TestFile:
        # HTML files are not compared in detail
        if kind in ("html", "ignore"):
            return TestMiscFile(name=name, kind=kind)

        filename = root / name
        if kind == "json":
            text, data = read_json(filename)
            return TestJsonFile(name=name, kind=kind, text=text, data=data)

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
        expected_text = self.text.splitlines(keepends=True)
        observed_text = read_and_decompress_file(observed).splitlines(keepends=True)

        if self._is_sam():
            self._preprocess_sam(expected_text)
            self._preprocess_sam(observed_text)

        if observed_text == expected_text:
            return

        lines = difflib.unified_diff(
            expected_text,
            observed_text,
            "expected",
            "observed",
        )

        self._raise_test_error(
            label="output",
            expected=expected,
            observed=observed,
            differences=pretty_output(lines, max_lines=10),
        )

    def _is_sam(self) -> bool:
        # Simple check for SAM header, which should always be present in reference files
        return self.text.startswith("@")

    @classmethod
    def _preprocess_sam(cls, lines: list[str]) -> None:
        for line_idx, line in enumerate(lines):
            if line.startswith("@PG") and "ID:adapterremoval" in line:
                fields = line.split("\t")
                for idx, value in enumerate(fields):
                    if value.startswith("CL:"):
                        # The location of the executable must be normalized
                        command = value[3:].split(" ")
                        command[0] = os.path.basename(command[0])
                        fields[idx] = "CL:" + " ".join(command)
                    elif value.startswith("VN:"):
                        # Program version must be masked
                        fields[idx] = "VN:3.x.x"

                lines[line_idx] = "\t".join(fields)


@dataclass
class TestJsonFile(TestFile):
    text: str
    data: JSON

    def compare_with_file(self, expected: Path, observed: Path) -> None:
        _, data = read_json(observed)
        differences = diff_json(reference=self.data, observed=data)
        differences = truncate_lines(differences, max_lines=4)

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

        JSONValidator.validate(data)

    def mask_filenames(self) -> TestJsonFile:
        def _mask_filenames(data: JSON) -> JSON:
            if isinstance(data, dict):
                data = dict(data)
                for key, value in data.items():
                    if key == "filenames":
                        data[key] = "...[str] | None"
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
        )


class TestMiscFile(TestFile):
    def compare_with_file(self, expected: Path, observed: Path) -> None:
        pass


class TestConfig(NamedTuple):
    description: str | None
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
    def load(cls, name: tuple[str, ...], filepath: Path) -> TestConfig:
        root = filepath.parent
        _, data = read_json(filepath)
        assert isinstance(data, dict)

        while True:
            inherit_from = json_pop_optional_str(data, ("inherit_from",))
            if inherit_from is not None:
                _, template = read_json(filepath.parent / inherit_from)
                assert isinstance(template, dict)
                for key, value in data.items():
                    if key in template and template[key] == value:
                        print_warn(
                            f"{{Key: value}} pair {{{key!r}: {value!r}}} in "
                            f"{quote(filepath)} is redundant"
                        )

                template.update(data)
                data = template
            else:
                break

        files = json_pop_dict(data, ("files",))
        test_files: list[TestFile] = []
        for key in _TEST_FILES:
            for filename in json_pop_tuple_of_str(files, ("files", key), []):
                test_files.append(TestFile.parse(root=root, name=filename, kind=key))

        if files:
            raise TestError(f"Unexpected files: {files}")

        self = TestConfig(
            description=json_pop_optional_str(data, ("description",)),
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

        def has_filetype(filetype: str) -> bool:
            return any(it.kind == filetype for it in self.files)

        if has_filetype("barcodes") and "--barcode-list" not in self.arguments:
            raise TestError(f"Missing --barcode-list in {quote(filepath)}")
        if has_filetype("adapters") and "--adapter-list" not in self.arguments:
            raise TestError(f"Missing --adapter-list in {quote(filepath)}")

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
            elif it.kind == "json":
                # Ensure that report generation is deterministic
                command += ["--report-sample-rate", "1"]
            elif it.kind not in (
                "adapters",
                "barcodes",
                "output",
                "json",
                "html",
                "ignore",
            ):
                raise NotImplementedError(it.kind)

        if input_1:
            command += ["--in-file1", *input_1]

        if input_2:
            command += ["--in-file2", *input_2]

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
            # Test without/with gzip compression of input files
            yield it
            yield it._replace(
                variant=it.variant + ("gz",),
                arguments=it.arguments,
                files=cls._compress_files(
                    files=it.files,
                    to_compress=_INPUT_FILES_FQ,
                    masked_stats=masked_stats,
                ),
            )

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

                # Remove any trailing newlines
                lines = it.text.splitlines(keepends=True)
                while lines and not lines[-1].rstrip():
                    lines.pop()

                input_files.setdefault(it.kind, []).append(lines)
            elif it.kind == "json":
                if masked_stats is not None:
                    other_files.append(masked_stats)
            else:
                other_files.append(it)

        if len(input_files) == 2:
            files_1 = input_files["input_1"]
            files_2 = input_files["input_2"]

            n_lines_1 = sum(map(len, files_1))
            n_lines_2 = sum(map(len, files_1))

            if n_lines_1 == n_lines_2 and n_lines_1 % 4 == 0:
                for idx, (file_1, file_2) in enumerate(zip(files_1, files_2)):
                    records_1 = cls._lines_to_records(file_1)
                    records_2 = cls._lines_to_records(file_2)

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
            else:
                print_warn("Could not interleave input files in", test.path)

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
    def description(self) -> str | None:
        return self._test.description

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

            command = self.command
            command = [os.path.abspath(command[0]), *command[1:]]
            command_s = cmd_to_s(command)

            write_data(
                self._test_path / "_run.sh",
                f"#!/bin/bash\n\n{command_s}\n",
            )

            write_data(
                self._test_path / "_gdb.sh",
                "#!/bin/bash\n"
                "gdb \\\n"
                "  --eval-command='break std::terminate' \\\n"
                "  --eval-command='run' \\\n"
                f"  --cd={quote(self._test_path)} \\\n"
                f"  --args {command_s}\n",
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
                f"returned {pretty_returncode(returncode)}",
                output=stderr,
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
            raise TestError(
                f"expected {pipe} not found: {expected_lines[0]}",
                output=actual_text,
            )

    def _evaluate_output_files(self) -> None:
        observed_files = set(os.listdir(self._test_path))
        expected_files: set[str] = {
            it.name for it in self._test.get_files(_OUTPUT_FILES | _INPUT_FILES)
        }

        if observed_files != expected_files:
            changes: list[str] = []
            for name in sorted(observed_files.symmetric_difference(expected_files)):
                changes.append(f"-{name}" if name in expected_files else f"+{name}")

            raise TestError(f"files do not match expectations: {', '.join(changes)}")

        for it in self._test.get_files(_OUTPUT_FILES | _INPUT_FILES):
            expected_files.add(it.name)

            it.compare_with_file(
                expected=self._exp_path / it.name,
                observed=self._test_path / it.name,
            )


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
    def description(self) -> str | None:
        return self._test.description

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
                assert isinstance(it, TestJsonFile)
                proc = subprocess.Popen(
                    [
                        sys.executable,
                        Path(__file__).parent / "update_regression_test.py",
                        "/dev/stdin",
                        filename,
                        "--save",
                    ],
                    stdin=subprocess.PIPE,
                    close_fds=True,
                    preexec_fn=os.setsid,
                )

                proc.communicate(it.text.encode())
                if proc.returncode != 0:
                    raise TestError(f"error dating test {self.path}")
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


def collect_tests(*, root: Path) -> list[TestConfig]:
    tests: list[TestConfig] = []
    for dirname, _, filenames in os.walk(root):
        dirname = Path(dirname)
        for filename in filenames:
            normalized = filename.lower()

            if normalized.startswith("test") and normalized.endswith(".json"):
                name = test_name(root, dirname, normalized)
                filepath = dirname / filename

                try:
                    test = TestConfig.load(name=name, filepath=filepath)
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
            if filename not in _IGNORED_FILES:
                observed_files.add(Path(dirname) / filename)

    recorded_files: set[Path] = set()
    for test in tests:
        recorded_files.add(test.path)
        dirname = test.path.parent
        for it in test.files:
            if it.kind not in ("html", "ignore"):
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
    verbose: bool
    use_colors: str


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
        "module `fastjsonschema`; validating the schema is only performed if the "
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
        default=min(8, multiprocessing.cpu_count()),
        help="Size of thread-pool for concurrent execution of tests",
    )
    parser.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="Print the name of each test that is executed",
    )
    parser.add_argument(
        "--use-colors",
        type=str.lower,
        choices=("auto", "always", "never"),
        default="auto",
        help="Enable color output",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=os.getcwd(),
        help="Source root used to construct relative paths",
    )

    return parser.parse_args(argv, namespace=Args())


def main(argv: list[str]) -> int:
    global _COLORS_ENABLED

    args = parse_args(argv)
    if args.use_colors == "auto":
        _COLORS_ENABLED = sys.stdout.isatty() and not os.getenv("NO_COLOR")
    elif args.use_colors == "always":
        _COLORS_ENABLED = True

    print(f"Testing executable {quote(args.executable)}")
    if not args.executable.is_file():
        print_err("ERROR: Executable does not exist")
        return 1
    elif args.schema_validation_required and not args.json_schema:
        print_err("ERROR: No JSON schema was provided")
        return 1
    elif args.schema_validation_required and not JSON_SCHEMA_VALIDATION:
        print_err("ERROR: Required module `fastjsonschema` could not be loaded")
        return 1
    elif not args.source_dir.is_dir():
        print_err(f"ERROR: {quote(args.source_dir)} is not a directory")
        return 1

    if args.max_failures <= 0:
        args.max_failures = float("inf")

    if args.create_updated_reference:
        # test-sets may generate (partially) overlapping files
        args.threads = 1

    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir = Path(tempfile.mkdtemp(dir=args.work_dir))

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

    if args.json_schema is not None and JSON_SCHEMA_VALIDATION:
        print(f"Reading JSON schema from {quote(args.json_schema)}")
        JSONValidator.load_schema(args.json_schema)

    print(f"Reading test-cases from {quote(args.source_dir)}")
    tests = collect_tests(root=args.source_dir)
    tests.sort(key=lambda it: it.path)

    if not tests:
        print_err(f"No tests found in {quote(args.source_dir)}")
        return 1

    unused_files = collect_unused_files(args.source_dir, tests)
    if unused_files:
        print_err(f"Found {len(unused_files)} unused files:")
        for idx, filepath in enumerate(unused_files, start=1):
            print(f"  {idx}. {quote(filepath)}")

        return 1

    if args.create_updated_reference:
        exhaustive_tests = build_test_updaters(args, tests)
    else:
        exhaustive_tests = build_test_runners(args, tests)

    print(f"Found {len(tests)} specifications; generated {len(exhaustive_tests)} tests")
    print(f"Using work-dir {quote(args.work_dir)}")

    n_failures = 0
    n_successes = 0
    n_skipped = 0
    failures: List[Tuple[Union[TestRunner, TestUpdater], Exception]] = []
    with multiprocessing.Pool(args.threads) as pool:
        results = pool.imap(run_test, exhaustive_tests)
        grouped_results = groupby(results, lambda it: it[0].name)

        for idx, (name, test_results) in enumerate(grouped_results, start=1):
            if args.verbose:
                print("  %i of %i: %s " % (idx, len(tests), name), end="")

            last_error = None
            test_skipped = False
            for test, error in test_results:
                if error is not None:
                    print_err("x", end="")
                    failures.append((test, error))
                    last_error = error
                elif test.skip:
                    print_warn("s", end="")
                    test_skipped = True
                else:
                    print_ok(".", end="")

                sys.stdout.flush()

            if test_skipped:
                if args.verbose:
                    print_warn(" [SKIPPED]")
                n_skipped += 1
            elif last_error is None:
                if args.verbose:
                    print_ok(" [OK]")
                n_successes += 1
            else:
                if args.verbose:
                    print_err(" [FAILED]")
                n_failures += 1
                if n_failures >= args.max_failures:
                    break

    for test, error in failures:
        print_err(f"\nTest {test.name} failed:")
        if test.description is not None:
            print_err("  Description   =", test.description)
        print_err("  Specification =", relative_to(args.source_root, test.spec_path))
        print_err("  Directory     =", relative_to(args.source_root, test.path))
        print_err("  Command       =", simplify_cmd(test.command))

        print_long_err("  Error         =", str(error))
        if isinstance(error, TestError) and error.output is not None:
            print_long_err("  Output        =", error.output, max_lines=4, tail=True)

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
        print("\nAll %i tests succeeded." % (len(exhaustive_tests),))

    if JSONValidator.can_validate():
        print("JSON file validation performed.")
    else:
        fastjsonschema_was = "available" if JSON_SCHEMA_VALIDATION else "unavailable"
        schema_was = "not provided" if args.json_schema is None else "provided"

        print_warn(
            f"JSON files not validated: `fastjsonschema` was {fastjsonschema_was}, "
            f"and JSON schema was {schema_was}"
        )

    if n_skipped:
        print_warn(f"Skipped {n_skipped} tests!")

    if args.create_updated_reference:
        print_ok(f"Updated regression tests written to {quote(args.work_dir)}")

    return 1 if n_failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
