"""
Documentation of the %YAML block of a parametrized template, generated from the pydantic models
(models.py) and from the live TEST_INFO options (testsuite.TESTCNF_KEYWORDS), so that it cannot drift
from the code. Printed by `run_params.py doc`.
"""
from __future__ import annotations

import inspect
import os
import sys
import textwrap
import types
import typing
from typing import Any, TextIO

from pydantic import BaseModel
from pydantic_core import PydanticUndefined

from ..testsuite import TESTCNF_KEYWORDS, TESTCNF_SECTIONS
from .models import (ALLOWED_SECTIONS, FORBIDDEN_KEYS, REACHABLE_NPROCS, TEMPLATE_KEY_NOTES, UNSWEEPABLE,
                     ExplicitAxis, FileToTestSpec, FixtureBlock, FixtureSpec, GridAxis, TemplateSpec, TestInfoSpec,
                     ToleranceSpec, Value)

#: Models documented by `doc_template_format`, with the title of their section.
DOCUMENTED_MODELS: list[tuple[str, type[BaseModel]]] = [
    ("Top level of the %YAML block", TemplateSpec),
    ("axes: grid axis", GridAxis),
    ("axes: explicit axis", ExplicitAxis),
    ("test_info", TestInfoSpec),
    ("test_info.files.files_to_test: one entry", FileToTestSpec),
    ("files_to_test overrides: tolerances of a group", ToleranceSpec),
    ("Fixture file Input/<name>.fixture.abi", FixtureSpec),
    ("fixture: block of a fixture file", FixtureBlock),
]

#: Rules checked by the validators and by the expansion, which the field types cannot express.
RULES = [
    "Only the sections " + ", ".join(ALLOWED_SECTIONS) + " can appear in test_info: "
    "[yaml_test] is not supported and [NCPU_N] is generated from nproc.",
    "files_to_test lists at least one file, each file once, unless setup.no_check is yes.",
    "The mode 'yaml' of TEST_INFO is obsolete: use yaml_docs or fldiff.",
    f"nproc values must be reachable by testbot: {sorted(REACHABLE_NPROCS)}.",
    "A parameter appears in one axis only, and exclude uses only swept parameters.",
    f"Pseudopotentials are fixed: {', '.join(sorted(UNSWEEPABLE))} cannot be swept.",
    "A swept variable must not appear in the body in any form (plain, with a dataset suffix, or in series syntax).",
    "udtset is not supported. With ndtset > 1, dataset is mandatory; without it, dataset is not allowed.",
    "A prerequisite template (an entry of setup.test_chain) has no axes, exclude, dataset, paral_info.max_nprocs "
    "or setup.test_chain of its own.",
    "A fixture is Input/<name>.fixture.abi next to the template; names use lowercase letters, digits and _. "
    "fixture.needs must not form a cycle, and two fixtures used by a template must not declare the same output.",
    "fixtures and setup.test_chain cannot be used together.",
    "Keywords are registered in tests/known_keywords.json and contain no blanks.",
    "Duplicate keys in the YAML block are errors (plain YAML would keep the last value).",
    "The body must not contain a #%% TEST_INFO block: the expansion generates it.",
]

EXAMPLE = """\
%YAML
test_info:
  setup: {executable: abinit}
  files:
    files_to_test:
      - {name: "{id}.abo", mode: yaml_docs, tolnlines: 8, tolabs: 1.0e-6, tolrel: 1.0e-6, fld_options: -medium,
         overrides: {ResultsGS.convergence: ignore}}   # reference: Refs/tNN.abo
  extra_info: {keywords: [NC], authors: [J. Doe], description: "GS results must not depend on useylm or paral_kgb"}
axes:
  - {type: grid, parameters: {useylm: [0, 1], paral_kgb: [0, 1]}}
%ENDYAML
# ABINIT input (the body) ...
"""


def type_name(tp: Any) -> str:
    """Readable name of a field annotation, e.g. 'int | all | null' or 'list of GridAxis | ExplicitAxis'."""
    if tp == Value:
        return "value (number, string, or list of them)"
    origin, args = typing.get_origin(tp), typing.get_args(tp)
    if origin is typing.Annotated:
        if args[0] is int and any(getattr(m, "gt", None) == 0 for m in tp.__metadata__):
            return "positive integer"
        return type_name(args[0])
    if origin is typing.Literal:
        return " | ".join(str(a) for a in args)
    if origin in (typing.Union, types.UnionType):
        return " | ".join(type_name(a) for a in args)
    if origin is list:
        return f"list of {type_name(args[0])}" if args else "list"
    if origin is dict:
        return f"mapping {type_name(args[0])} -> {type_name(args[1])}" if args else "mapping"
    if tp is type(None):
        return "null"
    if tp is Any:
        return "any"
    if inspect.isclass(tp):
        if issubclass(tp, BaseModel):
            return tp.__name__
        return {"int": "integer", "float": "float", "str": "string", "bool": "boolean"}.get(tp.__name__, tp.__name__)
    return str(tp)


def _default(field) -> str:
    if field.is_required():
        return "required"
    if field.default_factory is not None:
        value = field.default_factory()
        return "default: " + ("empty" if not value else repr(value))
    if field.default is PydanticUndefined:
        return "required"
    return f"default: {field.default!r}" if field.default not in (None, "") else (
        "default: none" if field.default is None else "default: ''")


def _wrap(text: str, indent: str, width: int = 110, first: str | None = None) -> str:
    return textwrap.fill(" ".join(text.split()), width=width, initial_indent=indent if first is None else first,
                         subsequent_indent=indent)


def doc_model(model: type[BaseModel], title: str, fh: TextIO) -> None:
    fh.write(f"\n== {title} ({model.__name__}) ==\n")
    if model.__doc__:
        fh.write(_wrap(inspect.cleandoc(model.__doc__), "") + "\n")
    for name, field in model.model_fields.items():
        if field.exclude:   # e.g. test_info.files_to_test: written under files: in the template
            continue
        tname = type_name(field.annotation)
        if field.annotation is int and any(getattr(m, "gt", None) == 0 for m in field.metadata):
            tname = "positive integer"   # PositiveInt: pydantic keeps the constraint in the metadata
        fh.write(f"\n  {name}: {tname}  ({_default(field)})\n")
        if field.description:
            fh.write(_wrap(field.description, "      ") + "\n")


def doc_test_info_options(fh: TextIO) -> None:
    """The TEST_INFO options, as for runtests.py --doc, with the template-specific notes."""
    fh.write("\n== TEST_INFO options available in test_info (same keys and meaning as the INI TEST_INFO) ==\n")
    for section in TESTCNF_SECTIONS:
        fh.write(f"\n  {section}:" + ("" if section in ALLOWED_SECTIONS else "  NOT SUPPORTED in templates") + "\n")
        if section not in ALLOWED_SECTIONS:
            continue
        for key, (_parser, default, sect, desc) in TESTCNF_KEYWORDS.items():
            if sect != section:
                continue
            if key in FORBIDDEN_KEYS:
                fh.write(f"    {key}: NOT ALLOWED ({FORBIDDEN_KEYS[key]})\n")
                continue
            dflt = "mandatory" if default is None else (f"default: {default}" if default else "")
            fh.write(f"    {key}" + (f" ({dflt})" if dflt else "") + "\n")
            fh.write(_wrap(desc, "        ") + "\n")
            if key in TEMPLATE_KEY_NOTES:
                fh.write(_wrap("In a template: " + TEMPLATE_KEY_NOTES[key], "        ") + "\n")


def doc_template_format(fh: TextIO = sys.stdout) -> None:
    """Write the documentation of the %YAML block of a template."""
    fh.write("Format of the %YAML block of a parametrized template (tests/parametrized/Input/tNN.tmpl.abi).\n"
             "Generated from tests/pymods/parametrized/models.py; design in tests/plans/parametrized_tests_plan.md.\n")
    fh.write("\n== Example ==\n" + textwrap.indent(EXAMPLE, "  "))
    for title, model in DOCUMENTED_MODELS:
        doc_model(model, title, fh)
    doc_test_info_options(fh)
    fh.write("\n== Other rules checked when the template is loaded or expanded ==\n")
    for rule in RULES:
        fh.write(_wrap(rule, "    ", first="  - ") + "\n")


FIXTURE_EXAMPLE = """\
%YAML
test_info:
  setup: {executable: abinit, output_prefix: Si_444k_8b}
  extra_info: {keywords: [NC], authors: [J. Doe], description: "NSCF WFK of Si with 8 bands on the 4x4x4 k-mesh"}
fixture: {outputs: [Si_444k_8b_WFK], nproc: 4}
%ENDYAML
# input of the executable (here an NSCF abinit run) ...
"""

FIXTURE_RULES = [
    "The file is Input/<name>.fixture.abi, next to the templates; templates list the fixtures by name (fixtures: [...]). "
    "Names use lowercase letters, digits and _ and start with a letter: names differing only by case are not portable.",
    "A fixture is never compared: files_to_test is not allowed and no_check is always yes.",
    "The number of processes is fixture.nproc (paral_info.max_nprocs is not allowed); mrgddb and mrgdv need nproc 1.",
    "fixture.needs must not form a cycle; two fixtures used by a template must not declare the same output.",
    "A fixture is built once per workdir root, before the tests, and its outputs are linked read-only into the "
    "workdir of each instance. The cache key is a hash of its input, of its executable, of the commit of its git "
    "submodule, of its nproc and of the fixtures it needs.",
]


def doc_fixture_format(fh: TextIO = sys.stdout) -> None:
    """Write the documentation of the fixture files (Input/<name>.fixture.abi)."""
    fh.write("Format of a fixture file (tests/parametrized/Input/<name>.fixture.abi).\n"
             "Generated from tests/pymods/parametrized/models.py; see also tests/parametrized/README.md.\n")
    fh.write("\n== Example ==\n" + textwrap.indent(FIXTURE_EXAMPLE, "  "))
    for title, model in DOCUMENTED_MODELS:
        if model in (FixtureSpec, FixtureBlock):
            doc_model(model, title, fh)
    fh.write("\n== Rules ==\n")
    for rule in FIXTURE_RULES:
        fh.write(_wrap(rule, "    ", first="  - ") + "\n")


def fixture_catalog(input_dir: str, known_keywords: dict[str, str] | None = None) -> list[dict[str, Any]]:
    """
    The fixtures of `input_dir` (*.fixture.abi) with their metadata, and the templates (*.tmpl.abi) that use them,
    directly or through fixture.needs. A fixture that cannot be loaded is listed with its error.
    """
    import glob
    from .template import FIXTURE_SUFFIX, TEMPLATE_SUFFIX, load_fixture, load_template

    catalog: dict[str, dict[str, Any]] = {}
    for path in sorted(glob.glob(os.path.join(input_dir, "*" + FIXTURE_SUFFIX))):
        name = os.path.basename(path)[: -len(FIXTURE_SUFFIX)]
        entry: dict[str, Any] = {"name": name, "path": path, "error": "", "used_by": []}
        try:
            fx = load_fixture(path, known_keywords=known_keywords)
            ti, block = fx.spec.test_info, fx.spec.fixture
            entry.update(executable=ti.setup.get("executable", ""), nproc=block.nproc, nproc_strict=block.nproc_strict,
                         outputs=list(block.outputs), needs=list(block.needs),
                         description=" ".join(str(ti.extra_info.get("description", "")).split()),
                         submodule=ti.files.get("use_git_submodule", ""))
        except Exception as exc:
            entry["error"] = str(exc)
        catalog[name] = entry

    def closure(names: list[str]) -> set[str]:
        seen: set[str] = set()
        todo = list(names)
        while todo:
            n = todo.pop()
            if n in seen:
                continue
            seen.add(n)
            todo += catalog.get(n, {}).get("needs", [])
        return seen

    for path in sorted(glob.glob(os.path.join(input_dir, "*" + TEMPLATE_SUFFIX))):
        try:
            spec = load_template(path, known_keywords=known_keywords).spec
        except Exception:
            continue
        for n in closure(list(spec.fixtures)):
            if n in catalog:
                catalog[n]["used_by"].append(os.path.basename(path)[: -len(TEMPLATE_SUFFIX)])
    return list(catalog.values())


def describe_fixtures(input_dir: str, fh: TextIO = sys.stdout, known_keywords: dict[str, str] | None = None) -> int:
    """Print the fixtures of input_dir with their documentation. Return the number of fixtures that cannot be loaded."""
    catalog = fixture_catalog(input_dir, known_keywords=known_keywords)
    if not catalog:
        fh.write(f"No fixture (*.fixture.abi) in {input_dir}\n")
        return 0
    fh.write(f"Fixtures in {input_dir}:\n")
    nbad = 0
    for e in catalog:
        fh.write(f"\n{e['name']}\n")
        if e["error"]:
            nbad += 1
            fh.write(_wrap(f"INVALID: {e['error']}", "    ") + "\n")
            continue
        if e["description"]:
            fh.write(_wrap(e["description"], "    ") + "\n")
        nproc = f"{e['nproc']} (strict)" if e["nproc_strict"] else str(e["nproc"])
        fh.write(f"    executable: {e['executable']}, nproc: {nproc}\n")
        fh.write(f"    outputs:    {', '.join(e['outputs'])}\n")
        if e["needs"]:
            fh.write(f"    needs:      {', '.join(e['needs'])}\n")
        if e["submodule"]:
            fh.write(f"    submodule:  {e['submodule']}\n")
        fh.write(f"    used by:    {', '.join(e['used_by']) or '(no template)'}\n")
    return nbad


__all__ = ["DOCUMENTED_MODELS", "describe_fixtures", "doc_fixture_format", "doc_template_format", "fixture_catalog",
           "type_name"]
