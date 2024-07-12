#!/usr/bin/env python3
"""
/*************************************************************************\\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Mikkel Schubert - mikkelsch@gmail.com           *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
\\*************************************************************************/
"""

from __future__ import annotations

import argparse
import functools
import io
import re
import sys
from enum import Enum
from pathlib import Path
from typing import NamedTuple, NoReturn

_RE_SECTION = re.compile(r"<!--\s+template:\s+([a-z0-9_]+)\s+-->", re.I)
_RE_FIELD = re.compile(
    r"""
    (
        {{\w+}}                           # Single template value
        | \[\[\w+\]\]                     # One or more template values (repeated line)
        | JS_TEMPLATE_[a-z0-9_]+          # Single value embedded in JS
        | JS_DEFAULT_\w+\s*=\s*\w+        # Value with default token embedded in JS
        | JS_DEFAULT_\w+\s*=\s*".*[^\\]"  # Value with quoted default embedded in JS
    )
    """,
    re.IGNORECASE | re.VERBOSE | re.ASCII,
)
_CPP_ENCODED = {
    "\\": "\\\\",
    '"': '\\"',
    "\a": "\\a",
    "\b": "\\b",
    "\f": "\\f",
    "\r": "\\r",
    "\t": "\\t",
    "\v": "\\v",
    "\n": "\\n",
}

_BUILTIN_VARS = {
    "id": "id",
}


_BASE_CLASS_HEADER = """
class html_template
{
public:
  html_template() = default;
  virtual ~html_template() = default;
  virtual void write(std::ostream& out) = 0;
};"""


def abort(fmt: str, *args: object) -> NoReturn:
    print("ERROR:", fmt.format(*args), file=sys.stderr)
    sys.exit(1)


class FieldType(Enum):
    REQUIRED = "required"
    REPEATED = "repeated"
    DEFAULT = "default"


class Field(NamedTuple):
    name: str
    kind: FieldType
    default: str | None


class Section(NamedTuple):
    lines: list[str]
    variables: list[Field]


def read_template(filepath: Path) -> dict[str, Section]:
    current: list[str] | None = None
    sections: dict[str, list[str]] = {}
    with filepath.open() as handle:
        for line in handle:
            match = _RE_SECTION.search(line)
            if match is not None:
                (name,) = match.groups()
                name = name.upper()
                if name in sections:
                    abort("Duplicate section name found: {!r}", name)
                current = sections[name.upper()] = []
            elif current is not None:
                current.append(line)

    result: dict[str, Section] = {}
    for key, lines in sections.items():
        text = "".join(lines)
        variables: dict[str, Field] = {}
        for value in _RE_FIELD.findall(text):
            default = None

            if value.lower().startswith("js_template_"):
                name = value[12:]
                kind = FieldType.REQUIRED
            elif value.lower().startswith("js_default_"):
                name, default = value[11:].split("=", 1)
                kind = FieldType.DEFAULT

                default = default.strip()
                if default.startswith('"'):
                    default = default[1:-1]
            elif value.startswith("{"):
                name = value[2:-2]
                kind = FieldType.REQUIRED
            else:
                name = value[2:-2]
                kind = FieldType.REPEATED

            name = name.strip().lower()
            if name not in _BUILTIN_VARS:
                variable = variables.get(name)
                if variable is not None and variable.kind != kind:
                    abort("{!r} is both {!r} and {!r}", name, kind, variable.kind)

                variables[name] = Field(
                    name=name,
                    kind=kind,
                    default=default,
                )

        result[key] = Section(
            lines=lines,
            variables=sorted(variables.values(), key=lambda it: it.name),
        )

    return result


def quote(value: str) -> str:
    result = ['"']
    for char in value:
        encoded = _CPP_ENCODED.get(char)
        if encoded is not None:
            result.append(encoded)
        else:
            result.append(char)
    result.append('"')

    return "".join(result)


def inject_variables(value: str) -> str:
    result = [" ", "out"]

    repeater = None
    for field in _RE_FIELD.split(value):
        lc_field = field.lower()

        result.append("<<")
        if lc_field.startswith("js_template_"):
            name = field[12:].lower()
            result.append(_BUILTIN_VARS.get(name, f"m_{name}"))
        elif lc_field.startswith("js_default_"):
            name, _ = field[11:].lower().split("=", 1)
            result.append(_BUILTIN_VARS.get(name, f"m_{name}"))
        elif field.startswith("{{") and field.endswith("}}"):
            name = field[2:-2].lower()
            result.append(_BUILTIN_VARS.get(name, f"m_{name}"))
        elif field.startswith("[[") and field.endswith("]]"):
            if repeater is not None:
                abort("multiple repeater values: {!r} and {!r}", repeater, field[2:-2])
            repeater = field[2:-2].lower()
            result.append("value")
        else:
            result.append(quote(field))

    result = " ".join(result) + ";"

    if repeater is None:
        return result

    return f"  for (const auto& value : m_{repeater}) {{\n  {result}\n  }}"


def to_classname(name: str) -> str:
    return f"html_{name.lower()}"


def write_header(sections: dict[str, Section]) -> str:
    handle = io.StringIO()

    def tprint(line: str, *args: object) -> None:
        print(line.format(*args), file=handle)

    tprint(__doc__.strip())
    tprint("#pragma once\n")
    tprint("#include <ostream>")
    tprint("#include <string>")
    tprint("#include <vector>")
    tprint("")
    tprint("namespace adapterremoval {{")

    tprint("{}", _BASE_CLASS_HEADER)

    for key, props in sections.items():
        classname = to_classname(key)

        tprint("\nclass {} : public html_template", classname)
        tprint("{{")
        tprint("public:")
        tprint("  {}();", classname)
        tprint("  ~{}() override;", classname)

        tprint("")
        tprint("  {}(const {}&) = delete;", classname, classname)
        tprint("  {}& operator=(const {}&) = delete;", classname, classname)

        if props.variables:
            tprint("")

        for field in props.variables:
            if field.kind == FieldType.REPEATED:
                tprint("  {}& add_{}(const std::string& value);", classname, field.name)
            else:
                tprint("  {}& set_{}(const std::string& value);", classname, field.name)

        tprint("\n  void write(std::ostream& out) override;")

        tprint("\nprivate:")
        tprint("  bool m_written;")

        for field in props.variables:
            if field.kind == FieldType.REPEATED:
                tprint("  std::vector<std::string> m_{};", field.name)
            else:
                tprint("  std::string m_{};", field.name)

            if field.kind != FieldType.DEFAULT:
                tprint("  bool m_{}_is_set;", field.name)

        tprint("}};")

    tprint("")
    tprint("}} // namespace adapterremoval")

    return handle.getvalue()


def write_implementations(sections: dict[str, Section], header_name: str) -> str:
    handle = io.StringIO()

    def tprint(line: str, *args: object) -> None:
        print(line.format(*args), file=handle)

    tprint(__doc__.strip())
    tprint('#include "{}"', header_name)
    tprint('#include "debug.hpp" // for AR_REQUIRE')
    tprint("#include <cstddef>   // for size_t")
    tprint("")
    tprint("namespace adapterremoval {{")
    tprint("")
    tprint("size_t g_html_id = 1;")

    for key, props in sections.items():
        classname = to_classname(key)

        tprint("\n{}::{}()", classname, classname)
        tprint("  : m_written()")
        for field in props.variables:
            if field.kind == FieldType.DEFAULT:
                tprint('  , m_{}("{}")', field.name, field.default)
            else:
                tprint("  , m_{}()", field.name)
                tprint("  , m_{}_is_set()", field.name)
        # Dummy comment to prevent re-formatting depending on num. of variables
        tprint("{{\n  //\n}}")

        tprint("\n{}::~{}()", classname, classname)
        tprint("{{")
        tprint('  AR_REQUIRE(m_written, "template {} was not written");', classname)
        tprint("}}\n")

        for field in props.variables:
            tprint("{}&", classname)
            if field.kind == FieldType.REPEATED:
                tprint("{}::add_{}(const std::string& value)", classname, field.name)
                tprint("{{")
                tprint("  m_{}.push_back(value);", field.name)
            else:
                tprint("{}::set_{}(const std::string& value)", classname, field.name)
                tprint("{{")
                tprint("  m_{} = value;", field.name)

            if field.kind != FieldType.DEFAULT:
                tprint("  m_{}_is_set = true;", field.name)

            tprint("  return *this;")
            tprint("}}\n")

        tprint("void")
        tprint("{}::write(std::ostream& out)", classname)
        tprint("{{")
        tprint('  AR_REQUIRE(!m_written, "template {} already written");', classname)

        for field in props.variables:
            if field.kind != FieldType.DEFAULT:
                tprint(
                    '  AR_REQUIRE(m_{0}_is_set, "{1}::{0} not set");',
                    field.name,
                    classname,
                )

        # prevent clang-format from adding linebreaks
        tprint("  // clang-format off")
        # cast to void to silence unused-variable warnings when ID isn't used
        tprint("  auto id = g_html_id; ++g_html_id; (void)id;")

        for line in props.lines:
            tprint("{}", inject_variables(line))

        tprint("  // clang-format on")
        tprint("  m_written = true;")
        tprint("}}")

    tprint("")
    tprint("}} // namespace adapterremoval")

    return handle.getvalue()


class Args(argparse.Namespace):
    template: Path
    output_prefix: Path


def parse_args(argv: list[str]) -> Args:
    parser = argparse.ArgumentParser(
        formatter_class=functools.partial(
            argparse.ArgumentDefaultsHelpFormatter,
            width=79,
        ),
        allow_abbrev=False,
    )

    parser.add_argument("template", type=Path, help="Path to HTML template")
    parser.add_argument("output_prefix", type=Path, help="Path prefix for output files")

    return parser.parse_args(argv, namespace=Args())


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    header_path = args.output_prefix.with_suffix(".hpp")
    impl_path = args.output_prefix.with_suffix(".cpp")

    sections = read_template(args.template)

    header = write_header(sections)
    implementations = write_implementations(sections, header_path.name)

    header_path.write_text(header)
    impl_path.write_text(implementations)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
