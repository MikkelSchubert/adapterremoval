#!/usr/bin/env python3
# pyright: strict
from __future__ import annotations

import argparse
import functools
import json
import sys
from pathlib import Path

import jsonschema
import jsonschema.exceptions


def load_json(filepath: Path) -> object:
    with filepath.open("rb") as handle:
        return json.load(handle)


def error_path(error: jsonschema.exceptions.ValidationError) -> str:
    path = ".".join(map(str, error.path))
    if error.parent is not None:
        assert isinstance(error.parent, jsonschema.exceptions.ValidationError)
        return f"{error_path(error.parent)}.{path}"

    return path


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
    exit_code = 0

    for filepath in args.files:
        try:
            jsonschema.validate(load_json(filepath), schema=schema)
        except jsonschema.exceptions.ValidationError as error:
            print(
                f"ERROR in {filepath} at {error_path(error)}: {error.message}",
                file=sys.stderr,
            )
            exit_code = 1
        else:
            print(f"Successfully validated {filepath}", filepath, file=sys.stderr)

    return exit_code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
