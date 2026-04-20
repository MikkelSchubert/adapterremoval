#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import argparse
import functools
import re
import sys
from pathlib import Path

SKIP_LIST: set[tuple[Path, str]] = {
    (Path("README.md"), "2.3.4"),  # Link to current stable version
}

# Captures subset of https://semver.org/, based on AdapterRemoval usage
SEMVER_RE = re.compile(
    r"(0|[1-9]\d*)"  # major version
    r"\.(0|[1-9]\d*)"  # minor version
    r"\.(0|[1-9]\d*)"  # patch version
    r"(?:-?((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)*))?"  # pre-release
)


def collect_version_strings(filepath: Path, n: int | None = None) -> list[str]:
    text = filepath.read_text()
    matches: list[str] = []

    skip_line = False
    for line in text.splitlines():
        if "NO_VERSION_CHECK" in line:
            skip_line = True

        if not skip_line:
            for match in SEMVER_RE.finditer(line):
                matches.append(match.group(0))

        skip_line = False
        if "NO_VERSION_CHECK_NEXT_LINE" in line:
            skip_line = True

    if n is not None and n < len(matches):
        matches = matches[:n]

    return sorted(set(matches))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=functools.partial(
            argparse.ArgumentDefaultsHelpFormatter,
            width=79,
        ),
        allow_abbrev=False,
    )

    _ = parser.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="Text files in which to check version strings",
    )
    _ = parser.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="Print detected version strings",
    )

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    # the version specified in meson.build is considered the authoritative version; this
    # is a bit arbitrary, but the check just needs to trigger if any version has changed
    reference = Path("meson.build")
    versions = collect_version_strings(reference, n=1)
    if len(versions) != 1:
        print(f"ERROR getting version from {reference}: {versions}", file=sys.stderr)
        return 1

    (expected,) = versions
    any_errors = False

    for filepath in args.files:
        for match in collect_version_strings(filepath):
            if args.verbose:
                print(match, filepath, file=sys.stderr)

            if match != expected and (filepath, match) not in SKIP_LIST:
                any_errors = True

                print(
                    f"Found version {match!r} in {filepath}, expected {expected!r}",
                    file=sys.stderr,
                )

    return 1 if any_errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
