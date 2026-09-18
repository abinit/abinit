#!/usr/bin/env python3
"""Translate legacy ABINIT ROBODoc comments for Doxygen without changing source files.

Doxygen invokes this program with a Fortran filename and parses the filtered text
written to standard output. The filter deliberately has no option or code path that
writes to the input file. It also preserves the original number of lines so that
Doxygen source links continue to point to the correct locations.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

_HEADER_START_RE = re.compile(r"^\s*!!\*{4}[a-z]\*", re.IGNORECASE)
_HEADER_END_RE = re.compile(r"^\s*!!\*{3}\s*$")
_SECTION_RE = re.compile(r"^\s*!!\s+(?P<name>[A-Z][A-Z ]+)\s*$")
_PARAMETER_RE = re.compile(
    r"^(?P<prefix>\s*!!)(?P<space>\s+)"
    r"\[?(?P<name>[A-Za-z][A-Za-z0-9_]*)"
    r"(?:\([^)]*\)|<[^>]*>)*\]?\s*=\s*(?P<description>.*)$"
)

_PARAMETER_SECTIONS = {
    "ARGUMENTS": "in",
    "INPUTS": "in",
    "OUTPUT": "out",
    "OUTPUTS": "out",
    "SIDE EFFECTS": "inout",
}

_SECTION_COMMANDS = {
    "FUNCTION": "@brief",
    "NOTES": "@note",
    "TODO": "@todo",
    "WARNINGS": "@warning",
}


def _routine_parameters(lines: list[str], header_start: int) -> set[str]:
    """Return dummy-argument names for the routine following a ROBODoc header."""
    source_index = next(
        (index for index in range(header_start, len(lines)) if lines[index].strip() == "!! SOURCE"),
        None,
    )
    if source_index is None:
        return set()

    declaration = ""
    for line in lines[source_index + 1 : source_index + 101]:
        if _HEADER_END_RE.match(line.rstrip("\r\n")):
            break
        code = line.split("!", maxsplit=1)[0].strip()
        if not declaration and not re.search(r"\b(?:subroutine|function)\s+[A-Za-z][A-Za-z0-9_]*\s*\(", code, re.IGNORECASE):
            continue
        declaration += " " + code.replace("&", " ")
        opening = declaration.find("(")
        if opening >= 0 and ")" in declaration[opening:]:
            arguments = declaration[opening + 1 : declaration.find(")", opening)]
            return {argument.strip().lower() for argument in arguments.split(",") if argument.strip()}

    return set()


def transform_text(text: str) -> str:
    """Return a line-preserving Doxygen view of ROBODoc comments."""
    lines = text.splitlines(keepends=True)
    result: list[str] = []
    routine_parameters: set[str] = set()
    in_header = False
    bridge_to_declaration = False
    section = ""

    for line_index, line in enumerate(lines):
        content = line.rstrip("\r\n")
        ending = line[len(content) :]

        if _HEADER_START_RE.match(content):
            in_header = True
            bridge_to_declaration = False
            section = ""
            routine_parameters = _routine_parameters(lines, line_index)
            content = f"{content[: len(content) - len(content.lstrip())]}!>"
        elif in_header and _HEADER_END_RE.match(content):
            in_header = False
            bridge_to_declaration = False
            section = ""

        if in_header:
            section_match = _SECTION_RE.match(content)
            if section_match:
                section = section_match.group("name").strip()
                if section == "SOURCE":
                    content = "!!"
                    bridge_to_declaration = True
                elif command := _SECTION_COMMANDS.get(section):
                    content = f"!! {command}"
            elif bridge_to_declaration and not content.strip():
                content = "!!"
            elif section in _PARAMETER_SECTIONS:
                parameter_match = _PARAMETER_RE.match(content)
                if parameter_match and parameter_match.group("name").lower() in routine_parameters:
                    direction = _PARAMETER_SECTIONS[section]
                    name = parameter_match.group("name")
                    description = parameter_match.group("description")
                    content = f"{parameter_match.group('prefix')} @param[{direction}] {name} {description}".rstrip()

            if bridge_to_declaration and content.strip() and content != "!!":
                bridge_to_declaration = False

        result.append(content + ending)

    transformed = "".join(result)
    if transformed.count("\n") != text.count("\n"):
        raise RuntimeError("ROBODoc filtering changed the number of source lines")
    return transformed


def main(argv: list[str] | None = None) -> int:
    """Write the filtered view of one source file to standard output."""
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        print(f"Usage: {Path(sys.argv[0]).name} FORTRAN_FILE", file=sys.stderr)
        return 2

    source_path = Path(args[0])
    source_bytes = source_path.read_bytes()
    source_text = source_bytes.decode("utf-8", errors="surrogateescape")
    output_bytes = transform_text(source_text).encode("utf-8", errors="surrogateescape")
    sys.stdout.buffer.write(output_bytes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
