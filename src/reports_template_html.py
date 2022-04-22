#!/usr/bin/env python3
# -*- coding: utf8 -*-
"""
/*************************************************************************\\
 * AdapterRemoval - cleaning next-generation sequencing reads            *
 *                                                                       *
 * Copyright (C) 2022 by Stinus Lindgreen - stinus@binf.ku.dk            *
 *                                                                       *
 * If you use the program, please cite the paper:                        *
 * S. Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation *
 * Sequencing Reads, BMC Research Notes, 5:337                           *
 * http://www.biomedcentral.com/1756-0500/5/337/                         *
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
import argparse
import io
import re
import sys
from pathlib import Path

_RE_SECTION = re.compile(r"<!--\s+template:\s+([a-z_]+)\s+-->", re.I)
_RE_FIELD = re.compile(r"{{([a-z0-9_]+)}}", re.I)
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


def read_template(filepath):
    current = None
    sections = {}
    with filepath.open() as handle:
        for line in handle:
            match = _RE_SECTION.search(line)
            if match is not None:
                (name,) = match.groups()
                assert name not in sections, name
                current = sections[name.upper()] = []
            elif current is not None:
                current.append(line)

    result = {}
    for key, lines in sections.items():
        text = "".join(lines)
        fields = _RE_FIELD.findall(text)
        result[key] = {
            "lines": lines,
            "fields": sorted(frozenset(map(str.lower, fields))),
        }

    return result


def quote(value):
    result = ['"']
    for char in value:
        encoded = _CPP_ENCODED.get(char)
        if encoded is not None:
            result.append(encoded)
        else:
            result.append(char)
    result.append('"')

    return "".join(result)


def inject_variables(value):
    def _to_variable(match):
        return '" << m_{} << "'.format(match.group(1).lower())

    return _RE_FIELD.sub(_to_variable, value)


def to_classname(name):
    return "HTMLTmpl{}".format(name.title().replace("_", ""))


def write_header(sections):
    handle = io.StringIO()

    def tprint(line, *args):
        print(line.format(*args), file=handle)

    tprint(__doc__.strip())
    tprint("#pragma once\n")
    tprint("#include <fstream>")
    tprint("#include <string>")
    for key, props in sections.items():
        classname = to_classname(key)

        tprint("\nclass {}", classname)
        tprint("{{")
        tprint("public:")
        tprint("  {}();", classname)
        tprint("  ~{}();", classname)

        if props["fields"]:
            tprint("")
            for key in props["fields"]:
                tprint("  void set_{}(const std::string& value);", key)

        tprint("\n  void write(std::ofstream& out);")

        tprint("\nprivate:")
        tprint("  bool m_written;")
        for key in props["fields"]:
            tprint("  std::string m_{};", key)
            tprint("  bool m_{}_is_set;", key)

        tprint("}};")

    return handle.getvalue()


def write_implementations(sections, header_name):
    handle = io.StringIO()

    def tprint(line, *args):
        print(line.format(*args), file=handle)

    tprint(__doc__.strip())
    tprint('#include "{}"', header_name)
    tprint('#include "debug.hpp"', header_name)
    for key, props in sections.items():
        classname = to_classname(key)

        tprint("\n{}::{}()", classname, classname)
        tprint("  : m_written()")
        for field in props["fields"]:
            tprint("  , m_{}()", field)
            tprint("  , m_{}()", field + "_is_set")

        # Dummy comment to prevent re-formatting depending on num. of variables
        tprint("{{\n  //\n}}")

        tprint("\n{}::~{}()", classname, classname)
        tprint("{{")
        tprint("  AR_DEBUG_ASSERT(m_written);")
        tprint("}}\n")

        for field in props["fields"]:
            tprint("void")
            tprint("{}::set_{}(const std::string& value)", classname, field)
            tprint("{{")
            tprint("  m_{} = value;", field)
            tprint("  m_{}_is_set = true;", field)
            tprint("}}\n")

        tprint("void")
        tprint("{}::write(std::ofstream& out)", classname)
        tprint("{{")
        tprint("  AR_DEBUG_ASSERT(!m_written);")

        for field in props["fields"]:
            tprint("  AR_DEBUG_ASSERT(m_{}_is_set);", field)

        # prevent clang-format from adding linebreaks
        tprint("  // clang-format off")

        for line in props["lines"]:
            tprint("  out << {};", inject_variables(quote(line)))

        tprint("  // clang-format on")
        tprint("  m_written = true;")
        tprint("}}")

    return handle.getvalue()


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("template", type=Path, help="Path to HTML template")
    parser.add_argument("output_prefix", type=Path, help="Path prefix for output files")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    args.header = args.output_prefix.with_suffix(".hpp")
    args.impl = args.output_prefix.with_suffix(".cpp")

    sections = read_template(args.template)

    header = write_header(sections)
    implementations = write_implementations(sections, args.header.name)

    args.header.write_text(header)
    args.impl.write_text(implementations)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
