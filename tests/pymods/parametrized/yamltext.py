"""
Text-based comparison of the YAML documents of two output files, with per-path overrides.

The YAML sections (``--- !Tag`` ... ``...``) are extracted as plain text and split into groups:
lines whose path has an ``ignore`` rule are dropped, lines matching an override go to their own group,
everything else goes to the default group. Each group is then compared with the unchanged
``FileToTest.compare`` (line-based fldiff) using the tolerances of that group.

Paths are ``Tag`` (a whole document) or ``Tag.key`` (a top-level key of the document, including its
continuation lines). The most specific rule wins.

This module must depend only on the standard library and ``testsuite``: the test objects that call it
are pickled in the test database, so it must not import pydantic or yaml.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Any, TextIO

from ..testsuite import FileToTest

# Rules applied to every yaml_docs comparison, merged under the template's overrides.
# Only `ignore` rules: these entries describe the run, not the result.
DEFAULT_OVERRIDES: dict[str, str | dict[str, Any]] = {
    # dimensions.mpw is per process with paral_kgb 1; meta.paral_kgb and solver.wfoptalg depend on the sweep.
    "DatasetInfo": "ignore",
    # Repeats the eigensolver (wfoptalg), which paral_kgb 1 switches from 0 to 114.
    "BeginCycle.solver": "ignore",
}

_TOLERANCE_KEYS = ("tolnlines", "tolabs", "tolrel", "fld_options")

_KEY_RE = re.compile(r"^([A-Za-z_][\w ]*?)\s*:")


@dataclass
class YamlLine:
    """A line of a YAML document with the tag of the document and the top-level key it belongs to."""
    tag: str
    key: str | None
    text: str


def extract_yaml_lines(path: str) -> list[YamlLine]:
    """
    Return the lines inside the YAML documents of the file `path`, excluding the document delimiters.
    """
    lines: list[YamlLine] = []
    tag, key = None, None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("--- !"):
                tag, key = line[5:].split("#")[0].strip(), None
                continue
            if tag is None:
                continue
            if line.rstrip() == "...":
                tag, key = None, None
                continue
            m = _KEY_RE.match(line)
            if m:
                key = m.group(1).strip()
            lines.append(YamlLine(tag, key, line))
    return lines


def merge_overrides(overrides: dict[str, Any] | None,
                    defaults: dict[str, Any] | None = None) -> dict[str, Any]:
    """Return the package defaults updated with the template overrides."""
    rules = dict(DEFAULT_OVERRIDES if defaults is None else defaults)
    rules.update(overrides or {})
    return rules


def rule_path(yline: YamlLine, rules: dict[str, Any]) -> str | None:
    """Return the most specific path in `rules` matching `yline`, or None (default group)."""
    if yline.key is not None:
        path = f"{yline.tag}.{yline.key}"
        if path in rules:
            return path
    if yline.tag in rules:
        return yline.tag
    return None


def group_lines(ylines: list[YamlLine], rules: dict[str, Any]) -> dict[str | None, list[str]]:
    """Split the lines into groups (None is the default group). Ignored lines are dropped."""
    groups: dict[str | None, list[str]] = {}
    for yl in ylines:
        path = rule_path(yl, rules)
        if path is not None and rules[path] == "ignore":
            continue
        groups.setdefault(path, []).append(yl.text)
    return groups


def fldiff_text(lines: list[str]) -> list[str]:
    """
    fldiff gives the first column of a line a special meaning (e.g. `-` ignores the line).
    Prepend a blank so that every YAML line is compared as an ordinary numeric line,
    as DataExtractor does for the lines of a document when YAML parsing is disabled.
    """
    return [" " + line if line.endswith("\n") else " " + line + "\n" for line in lines]


#: Files with all the lines compared by yaml_docs (the groups one after the other), in the workdir of the
#: comparison: what the text and HTML diffs of the report show instead of the whole output files.
EXTRACT_OUTPUT, EXTRACT_REFERENCE = "output.txt", "reference.txt"


def _group_order(ref_groups, out_groups, rules) -> list[str | None]:
    """The default group first, then the overrides in the order of the rules."""
    return [None] + [p for p in rules if p in ref_groups or p in out_groups]


def write_extracts(workdir: str, ref_groups, out_groups, rules) -> None:
    """Write EXTRACT_OUTPUT and EXTRACT_REFERENCE in `workdir`, each group preceded by its label."""
    os.makedirs(workdir, exist_ok=True)
    for fname, groups in ((EXTRACT_OUTPUT, out_groups), (EXTRACT_REFERENCE, ref_groups)):
        with open(os.path.join(workdir, fname), "w") as fh:
            for path in _group_order(ref_groups, out_groups, rules):
                if path in ref_groups or path in out_groups:
                    fh.write(f"# yaml_docs[{'default' if path is None else path}]\n")
                    fh.writelines(groups.get(path, []))


@dataclass
class GroupResult:
    path: str | None
    tolerances: dict[str, Any]
    isok: bool
    status: str
    msg: str
    has_line_count_error: bool = False


@dataclass
class YamlTextResult:
    isok: bool
    status: str
    msg: str
    has_line_count_error: bool
    groups: list[GroupResult] = field(default_factory=list)


def _status_of(groups: list[GroupResult]) -> str:
    stats = {g.status for g in groups}
    if "failed" in stats:
        return "failed"
    if "passed" in stats:
        return "passed"
    return "succeeded"


def compare_yaml_text(ref_path: str, out_path: str, base: dict[str, Any], rules: dict[str, Any],
                      workdir: str, fldiff_path: str = "", outf: TextIO | None = None) -> YamlTextResult:
    """
    Compare the YAML documents of `ref_path` and `out_path`.

    Args:
        ref_path: Reference file.
        out_path: File produced by the test.
        base: Tolerances of the default group: tolnlines, tolabs, tolrel, fld_options.
        rules: Path -> "ignore" | dict with a subset of the tolerance keys (see merge_overrides).
        workdir: Directory where the per-group files are written.
        fldiff_path: Passed to FileToTest.compare (unused by the current implementation).
        outf: Stream receiving the fldiff details of every group.
    """
    ref_groups = group_lines(extract_yaml_lines(ref_path), rules)
    out_groups = group_lines(extract_yaml_lines(out_path), rules)
    write_extracts(workdir, ref_groups, out_groups, rules)

    if not ref_groups and not out_groups:
        msg = "yaml_docs: no YAML documents to compare (empty extract)"
        return YamlTextResult(False, "failed", msg, False)

    results: list[GroupResult] = []
    for path in _group_order(ref_groups, out_groups, rules):
        if path not in ref_groups and path not in out_groups:
            continue
        tolerances = dict(base)
        if path is not None:
            tolerances.update({k: v for k, v in rules[path].items() if v is not None})

        gname = "default" if path is None else re.sub(r"\W", "_", path)
        gdir = os.path.join(workdir, gname)
        os.makedirs(os.path.join(gdir, "ref"), exist_ok=True)
        fname = "yaml_docs.txt"
        with open(os.path.join(gdir, "ref", fname), "w") as fh:
            fh.writelines(fldiff_text(ref_groups.get(path, [])))
        with open(os.path.join(gdir, fname), "w") as fh:
            fh.writelines(fldiff_text(out_groups.get(path, [])))

        ftt = FileToTest({"name": fname, **{k: tolerances[k] for k in _TOLERANCE_KEYS if k in tolerances}})
        label = "yaml_docs[default]" if path is None else f"yaml_docs[{path}]"
        if outf is not None:
            outf.write(f"\n==== {label} with {tolerances}\n")
        isok, status, msg = ftt.compare(fldiff_path, os.path.join(gdir, "ref"), gdir, yaml_test={}, outf=outf)
        results.append(GroupResult(path, tolerances, isok, status, f"{label}: {msg}", ftt.has_line_count_error))

    status = _status_of(results)
    isok = all(g.isok for g in results)
    msg = "; ".join(g.msg for g in results)
    return YamlTextResult(isok, status, msg, any(g.has_line_count_error for g in results), results)
