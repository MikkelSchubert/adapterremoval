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
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import argparse
import functools
import json
import sys
from pathlib import Path
from typing import Any

import fastjsonschema


def load_json(filepath: Path) -> object:
    with filepath.open("rb") as handle:
        return json.load(handle)


class Args(argparse.Namespace):
    schema: Path
    files: list[Path]


def parse_args(argv: list[str]) -> Args:
    parser = argparse.ArgumentParser(
        formatter_class=functools.partial(
            argparse.ArgumentDefaultsHelpFormatter,
            width=79,
        ),
        allow_abbrev=False,
    )

    parser.add_argument("schema", type=Path, help="JSON schema")
    parser.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="One or more JSON files to be validated using the supplied schema",
    )

    return parser.parse_args(argv, namespace=Args())


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    schema = load_json(args.schema)
    validator: Any = fastjsonschema.compile(schema, use_default=False)
    exit_code = 0

    for filepath in args.files:
        try:
            validator(load_json(filepath))
        except fastjsonschema.JsonSchemaValueException as error:
            path = ".".join(map(str, error.path))
            print(f"ERROR in {filepath} at {path}: {error.message}", file=sys.stderr)
            exit_code = 1
        else:
            print(f"Successfully validated {filepath}", filepath, file=sys.stderr)

    return exit_code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
