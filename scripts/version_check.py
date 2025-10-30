#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2025 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import argparse
import functools
import re
import sys
from pathlib import Path

# Based on https://semver.org/
SEMVER_RE = re.compile(
    r"(0|[1-9]\d*)"  # major version
    r"\.(0|[1-9]\d*)"  # minor version
    r"\.(0|[1-9]\d*)"  # patch version
    r"(?:-?((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)*))?"  # pre-release
)

IGNORED_FILES = (
    "CHANGES.md",
    "Containerfile",
    "uv.lock",
)


def collect_version_strings(filepath: Path) -> list[str]:
    text = filepath.read_text()
    matches: set[str] = set()

    for match in SEMVER_RE.finditer(text):
        value = match.group(0)

        # Trim what are typically file extensions
        while True:
            parts = value.rsplit(".", 1)
            if len(parts) > 1 and parts[-1].isalpha():
                value = parts[0]
            else:
                break

        if value:
            matches.add(value)

    return sorted(matches)


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

    # the version specified in main.hpp is considered the authoritative version; this is
    # a bit arbitrary, but the check just needs to trigger if any version has changed
    reference = Path("src") / "main.hpp"
    versions = collect_version_strings(reference)
    if len(versions) != 1:
        print(f"ERROR getting version from {reference}: {versions}", file=sys.stderr)
        return 1

    (expected,) = versions
    any_errors = False

    for filepath in args.files:
        if (
            # Some specific files contain non-AR version strings
            filepath.name not in IGNORED_FILES
            # Reference SAM files are expected to contain out-dated version strings
            and filepath.suffix != ".sam"
            # Reference JSON reports are expected to contain out-dated version strings
            and (filepath.suffix != ".json" or not filepath.is_relative_to("tests"))
            # workflows contains numerous other software versions, but no AR versions
            and not filepath.is_relative_to(".github")
        ):
            for match in collect_version_strings(filepath):
                if args.verbose:
                    print(match, filepath, file=sys.stderr)

                if "3." <= match <= "4.":
                    if match != expected:
                        any_errors = True

                        print(
                            f"Found version {match!r} in {filepath}, expected",
                            f"{expected!r}",
                            file=sys.stderr,
                        )
        elif args.verbose:
            print("skipping", filepath, file=sys.stderr)

    return 1 if any_errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
