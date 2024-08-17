#!/usr/bin/env python3
# pyright: strict
from __future__ import annotations

import argparse
import functools
import sys
from dataclasses import dataclass, field
from enum import Enum, auto
from itertools import takewhile
from pathlib import Path
from typing import Callable, Iterable, Iterator, TypeAlias, Union


class JSONError(Exception):
    pass


class TextIter:
    def __init__(self, text: str) -> None:
        self._text = text
        self._pos = 0

    def get(self) -> str | None:
        if self._pos < len(self._text):
            return self._text[self._pos]

        return None

    def next(self) -> None:
        self._pos += 1

    def rewind(self) -> None:
        if self._pos <= 0:
            raise AssertionError("Rewind at starting position")
        self._pos -= 1

    def __iter__(self) -> Iterator[str]:
        return self

    def __next__(self) -> str:
        if self._pos < len(self._text):
            self._pos += 1
            return self._text[self._pos - 1]

        raise StopIteration


JSON: TypeAlias = Union["JSONValue", "JSONList", "JSONDict"]


@dataclass(frozen=True)
class JSONValue:
    token: list[str]
    value: str | int | float | bool | None

    def update(self, other: JSONValue) -> JSONValue:
        token = "".join(other.token).strip()
        # Preserve current whitespace
        leading_whitespace = list(takewhile(str.isspace, other.token))
        trailing_whitespace = list(takewhile(str.isspace, reversed(other.token)))

        return JSONValue(
            token=leading_whitespace + [token] + trailing_whitespace[::-1],
            value=other.value,
        )

    def __str__(self) -> str:
        return f"{''.join(self.token)}"

    def __repr__(self) -> str:
        return f"<{self}>"

    def __hash__(self) -> int:
        return hash(self.value)


@dataclass
class JSONList:
    leading_token: list[str] = field(default_factory=list)
    trailing_token: list[str] = field(default_factory=list)
    values: list[JSON] = field(default_factory=list)

    def update(self, other: JSONValue) -> JSONValue:
        token = "".join(other.token).strip()
        return JSONValue(
            token=self.leading_token[:-1] + [token] + self.trailing_token[1:],
            value=other.value,
        )

    def __len__(self) -> int:
        return len(self.values)

    def __str__(self) -> str:
        return self._to_str(str)

    def __repr__(self) -> str:
        return f"<{self._to_str(repr)}>"

    def _to_str(self, func: Callable[[object], str]) -> str:
        leading = "".join(self.leading_token)
        trailing = "".join(self.trailing_token)
        values = ",".join(func(value) for value in self.values)

        return f"{leading}{values}{trailing}"


@dataclass
class JSONDict:
    leading_token: list[str] = field(default_factory=list)
    trailing_token: list[str] = field(default_factory=list)
    values: dict[JSONValue, JSON] = field(default_factory=dict)

    def update(self, other: JSONValue) -> JSONValue:
        token = "".join(other.token).strip()
        return JSONValue(
            token=self.leading_token[:-1] + [token] + self.trailing_token[1:],
            value=other.value,
        )

    def items(self) -> Iterable[tuple[JSONValue, JSON]]:
        return self.values.items()

    def __getitem__(self, key: JSONValue) -> JSON:
        return self.values[key]

    def __setitem__(self, key: JSONValue, value: JSON) -> None:
        self.values[key] = value

    def __contains__(self, key: JSONValue) -> bool:
        return key in self.values

    def __len__(self) -> int:
        return len(self.values)

    def __str__(self) -> str:
        return self._to_str(str)

    def __repr__(self) -> str:
        return f"<{self._to_str(repr)}>"

    def _to_str(self, func: Callable[[object], str]) -> str:
        leading = "".join(self.leading_token)
        trailing = "".join(self.trailing_token)
        values: list[str] = []
        for idx, (key, value) in enumerate(self.values.items()):
            values += [func(key), ":", str(value)]
            if idx + 1 < len(self.values):
                values.append(",")

        return f"{leading}{''.join(values)}{trailing}"


class Parsing(Enum):
    ExpectingValue = auto()
    ExpectingComma = auto()

    InString = auto()


class JSONParser:
    @classmethod
    def parse(cls, text: str) -> JSON:
        reader = TextIter(text)
        value = cls._parse_next(reader)
        if reader.get() is not None:
            raise JSONError(f"unexpected charcter {reader.get()!r} after root value")

        return value

    @classmethod
    def _parse_next(cls, reader: TextIter) -> JSON:
        token: list[str] = []
        for char in reader:
            token.append(char)

            if char.isspace():
                continue
            elif char == "[":
                return cls._parse_list(reader, token)
            elif char == "{":
                return cls._parse_dict(reader, token)
            elif char == '"':
                return cls._parse_string(reader, token)
            else:
                return cls._parse_scalar(reader, token)

        raise JSONError("no object(s) found")

    @classmethod
    def _parse_dict(cls, reader: TextIter, token: list[str]) -> JSONDict:
        values: dict[JSONValue, JSON] = {}

        expects_value = False
        while reader.get() not in ("}", None):
            key = cls._parse_next(reader)

            if not (isinstance(key, JSONValue) and isinstance(key.value, str)):
                raise JSONError(f"found non-string key {key!r}")
            elif key in values:
                raise JSONError(f"duplicate key {key!r}")
            elif reader.get() != ":":
                raise JSONError(f"unexpected char after key {reader.get()!r}")

            reader.next()
            values[key] = cls._parse_next(reader)

            expects_value = False
            if reader.get() == ",":
                reader.next()
                expects_value = True

        if expects_value:
            raise JSONError("trailing comma in dict")

        char = reader.get()
        if char is None:
            raise JSONError("EOF while parsing dict")

        reader.next()  # consume "}"
        return JSONDict(
            leading_token=token,
            trailing_token=cls._append_trailing_whitespace(reader, [char]),
            values=values,
        )

    @classmethod
    def _parse_list(cls, reader: TextIter, token: list[str]) -> JSONList:
        values: list[JSON] = []

        expects_value = False
        while reader.get() not in ("]", None):
            value = cls._parse_next(reader)
            expects_value = False

            char = reader.get()
            if char is None:
                raise JSONError("EOF while parsing list")
            elif char == ",":
                expects_value = True
                reader.next()
            elif char != "]":
                # trailing whitespace consumed by parse next
                raise JSONError(f"Invalid char in list {char!r}")

            values.append(value)

        if expects_value:
            raise JSONError("trailing comma in list")

        char = reader.get()
        if char is None:
            raise JSONError("EOF while parsing list")

        reader.next()  # consume "]"
        return JSONList(
            leading_token=token,
            trailing_token=cls._append_trailing_whitespace(reader, [char]),
            values=values,
        )

    @classmethod
    def _parse_string(cls, reader: TextIter, token: list[str]) -> JSONValue:
        in_escape = False
        for char in reader:
            token.append(char)
            if char == "\\":
                in_escape = True
            elif char == '"' and not in_escape:
                cls._append_trailing_whitespace(reader, token)
                value = "".join(token).strip()[1:-1]
                return JSONValue(token=token, value=value)
            elif ord(char) not in range(0x20, 0x10FFFF):
                raise JSONError("Invalid character in string")

        raise JSONError("EOF while parsing string")

    @classmethod
    def _parse_scalar(cls, reader: TextIter, token: list[str]) -> JSONValue:
        for char in reader:
            if char in "+-.0123456789Eaeflnrstu":
                token.append(char)
            else:
                reader.rewind()
                break

        # Include only contiguous trailing whitespace
        cls._append_trailing_whitespace(reader, token)

        value = "".join(token).strip()
        if value == "true":
            return JSONValue(token=token, value=True)
        elif value == "false":
            return JSONValue(token=token, value=False)
        elif value == "null":
            return JSONValue(token=token, value=None)
        elif value.isdigit():
            return JSONValue(token=token, value=int(value))

        try:
            return JSONValue(token=token, value=float(value))
        except ValueError:
            raise JSONError(f"invalid scalar value {value!r}") from None

    @classmethod
    def _append_trailing_whitespace(
        cls,
        reader: TextIter,
        token: list[str],
    ) -> list[str]:
        for char in reader:
            if not char.isspace():
                reader.rewind()
                break

            token.append(char)

        return token


def load_json(filepath: Path) -> JSON:
    text = filepath.read_text()
    obj = JSONParser.parse(text)
    if str(obj) != text:
        raise RuntimeError(f"parsing of {filepath} is not lossless")

    return obj


WILDCARDS = {
    "...",
    "...str",
    "...float",
    "...[str]",
    "...[str] | None",
}


def is_wildcard(value: JSON) -> bool:
    return isinstance(value, JSONValue) and value.value in WILDCARDS


def update(source: JSON, destination: JSON) -> JSON:
    if isinstance(source, JSONValue) and source.value in WILDCARDS:
        return destination.update(source)
    elif isinstance(source, JSONDict) and isinstance(destination, JSONDict):
        dict_values: dict[JSONValue, JSON] = dict(destination.values)
        for key, value in dict_values.items():
            if key in source:
                dict_values[key] = update(source[key], value)

        return JSONDict(
            leading_token=destination.leading_token,
            trailing_token=destination.trailing_token,
            values=dict_values,
        )
    elif isinstance(source, JSONList) and isinstance(destination, JSONList):
        list_values: list[JSON] = []
        for src, dst in zip(source.values, destination.values, strict=False):
            list_values.append(update(src, dst))

        return JSONList(
            leading_token=destination.leading_token,
            trailing_token=destination.trailing_token,
            values=list_values + destination.values[len(list_values) :],
        )

    return destination


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=functools.partial(
            argparse.ArgumentDefaultsHelpFormatter,
            width=79,
        ),
        allow_abbrev=False,
    )

    parser.add_argument(
        "reference",
        type=Path,
        help="Original JSON report from regression test",
    )
    parser.add_argument(
        "updated",
        type=Path,
        help="JSON report resulting from (intentional) breaking changes",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save changes to the updated test instead of printing to stdout",
    )

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    reference = load_json(args.reference)
    updated = load_json(args.updated)

    text = update(source=reference, destination=updated)
    if args.save:
        args.updated.write_text(str(text))
    else:
        print(text, end="")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
