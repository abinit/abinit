"""
Expand a parametrized template into test objects.

This is the function called by the hook in testsuite.make_abitests_from_inputs for inputs ending
with .tmpl.abi. Nothing is written in the repository: each instance keeps its generated input in memory
and writes it into its workdir when the test runs.
"""
from __future__ import annotations

import itertools
import json
import os
import random
import re
import tempfile
import zlib
from dataclasses import dataclass, field
from typing import Any

from ..testsuite import AbinitTestInfoParser, exec2class
from .models import MAIN_OUTPUT, REACHABLE_NPROCS, FileToTestSpec, GridAxis, TemplateSpec, ini_value
from .paramtests import (BrokenTemplateTest, Fixture, ParamFileToTest, ParametrizedAbinitTest, ParametrizedChain,
                         ParametrizedTest)
from .template import Template, TemplateError, fixture_path, load_fixture, load_template

#: Keyword added to every instance (registered in tests/known_keywords.json).
PARAMETRIZED_KEYWORD = "PARAMETRIZED"


class ExpansionError(Exception):
    """Raised when a valid template cannot be expanded (e.g. conflict with the body)."""


@dataclass
class Instance:
    index: int
    params: dict[str, Any]
    #: Number of MPI processes (1 for a sequential run).
    nproc: int
    #: Test id, including the _MPI<N> suffix for parallel runs.
    test_id: str
    input_text: str = ""
    files: list[tuple[FileToTestSpec, str, str]] = field(default_factory=list)  # (spec, name, ref_name)
    #: A prerequisite of a chain: runs on 1 process whatever the nproc of its chain.
    serial: bool = False


# ---------------------------------------------------------------------------------------------------
# Combinations
# ---------------------------------------------------------------------------------------------------

def axis_dicts(axis) -> list[dict[str, Any]]:
    if isinstance(axis, GridAxis):
        names = list(axis.parameters)
        lists = [axis.parameters[n] for n in names]
        if axis.mode == "zip":
            return [dict(zip(names, vals)) for vals in zip(*lists)]
        return [dict(zip(names, vals)) for vals in itertools.product(*lists)]
    return [dict(v) for v in axis.values]


def combinations(spec: TemplateSpec) -> list[dict[str, Any]]:
    """Cartesian product across the axes, minus the excluded combinations. Deterministic order."""
    combos: list[dict[str, Any]] = [{}]
    for axis in spec.axes:
        combos = [{**c, **d} for c in combos for d in axis_dicts(axis)]

    def excluded(combo: dict[str, Any]) -> bool:
        return any(all(combo.get(k, object()) == v for k, v in ex.items()) for ex in spec.exclude)

    return [c for c in combos if not excluded(c)]


def draw_nproc(template_name: str, max_nprocs: int) -> int:
    """The max_nprocs shorthand: one reachable value <= max_nprocs, seeded by the template name."""
    choices = sorted(n for n in REACHABLE_NPROCS if n <= max_nprocs)
    if not choices:
        raise ExpansionError(f"max_nprocs={max_nprocs}: no reachable process count")
    return random.Random(zlib.crc32(template_name.encode())).choice(choices)


# ---------------------------------------------------------------------------------------------------
# Body scanning and variable insertion
# ---------------------------------------------------------------------------------------------------

def strip_comments(body: str) -> str:
    return "\n".join(re.split(r"[#!]", line, maxsplit=1)[0] for line in body.splitlines())


def scan_datasets(body: str) -> tuple[int, list[int] | None]:
    """Return (ndtset, jdtset). Only these few tokens are scanned: this is not an ABINIT parser."""
    text = strip_comments(body)
    if re.search(r"(?<![\w])udtset(?![\w])", text, re.IGNORECASE):
        raise ExpansionError("udtset is not supported in templates")
    m = re.search(r"(?<![\w])ndtset\s+(\d+)", text, re.IGNORECASE)
    ndtset = int(m.group(1)) if m else 1
    jdtset = None
    m = re.search(r"(?<![\w])jdtset\s+((?:\d+\s*)+)", text, re.IGNORECASE)
    if m:
        jdtset = [int(t) for t in m.group(1).split()][:ndtset]
    return ndtset, jdtset


def find_conflicts(body: str, names: list[str]) -> list[str]:
    """Swept variables must not appear in the body in any form: plain, suffixed or in series syntax."""
    text = strip_comments(body)
    conflicts = []
    for name in names:
        regex = re.compile(r"(?<![\w])" + re.escape(name) + r"\d*[:+*?]?(?=\s|$)", re.IGNORECASE | re.MULTILINE)
        if regex.search(text):
            conflicts.append(name)
    return conflicts


def format_value(value: Any) -> str:
    if isinstance(value, (list, tuple)):
        return " ".join(str(v) for v in value)
    return str(value)


def insertion_lines(params: dict[str, Any], dataset: int | str | None, ndtset: int,
                    jdtset: list[int] | None) -> list[str]:
    """Lines appended to the body for the swept variables (nproc is never inserted)."""
    if ndtset > 1 and dataset is None:
        raise ExpansionError(f"the body has ndtset {ndtset}: the template must set `dataset` (a dataset index or all)")
    if ndtset <= 1 and dataset is not None:
        raise ExpansionError("`dataset` can be used only with a multi-dataset body (ndtset > 1)")
    if isinstance(dataset, int):
        valid = jdtset if jdtset is not None else list(range(1, ndtset + 1))
        if dataset not in valid:
            raise ExpansionError(f"dataset {dataset} is not one of the datasets of the body: {valid}")
    suffix = str(dataset) if isinstance(dataset, int) else ""
    return [f"{name}{suffix} {format_value(value)}" for name, value in params.items() if name != "nproc"]


# ---------------------------------------------------------------------------------------------------
# TEST_INFO compilation
# ---------------------------------------------------------------------------------------------------

def _ini_line(key: str, value: Any) -> list[str]:
    text = ini_value(key, value)
    first, *rest = text.splitlines() or [""]
    return [f"{key} = {first}"] + [f"  {line}" for line in rest]


def _sentinel(line: str) -> str:
    """Prefix a TEST_INFO line with #%%. The parser requires #%%<BEGIN/END TEST_INFO> without blank."""
    return "#%%" + line if line.startswith("<") else "#%% " + line


def files_to_test_ini(files: list[tuple[FileToTestSpec, str, str]]) -> str:
    entries = []
    for spec, name, _ref in files:
        entry = f"{name}, tolnlines = {spec.tolnlines}, tolabs = {spec.tolabs}, tolrel = {spec.tolrel}"
        if spec.fld_options:
            entry += f", fld_options = {spec.fld_options}"
        if spec.verbose_report == "yes":
            entry += ", verbose_report = yes"
        entries.append(entry)
    return "; ".join(entries)


def compile_test_info(template: Template, inst: Instance) -> list[str]:
    """Return the lines (without the #%% sentinel) of the TEST_INFO block of an instance."""
    ti = template.spec.test_info
    sections: dict[str, list[str]] = {s: [] for s in ("setup", "files", "shell", "paral_info", "extra_info")}

    for key, value in ti.setup.items():
        if key == "test_chain":
            continue  # the chain is built by the expansion (ParametrizedChain), not by the INI parser
        sections["setup"] += _ini_line(key, value)
    for key, value in ti.files.items():
        sections["files"] += _ini_line(key, value)
    for key, value in ti.shell.items():
        sections["shell"] += _ini_line(key, value)

    extra = dict(ti.extra_info)
    keywords = extra.pop("keywords", []) or []
    if isinstance(keywords, str):
        keywords = [k.strip() for k in keywords.split(",") if k.strip()]
    extra["keywords"] = list(dict.fromkeys([*keywords, PARAMETRIZED_KEYWORD]))
    suffix = ", ".join(f"{k}={format_value(v)}" for k, v in inst.params.items())
    description = str(extra.pop("description", "Parametrized test")).rstrip()
    extra["description"] = f"{description} [{suffix}]" if suffix else description
    for key, value in extra.items():
        sections["extra_info"] += _ini_line(key, value)

    ftt = files_to_test_ini(inst.files) if inst.files else ""
    ncpu_section: list[str] = []
    if inst.nproc > 1:
        # generate_testinfo_nprocs requires the [NCPU_N] section even when there is nothing to compare (no_check).
        sections["paral_info"] += [f"nprocs_to_test = {inst.nproc}", f"max_nprocs = {inst.nproc}"]
        ncpu_section = [f"[NCPU_{inst.nproc}]"] + ([f"files_to_test = {ftt}"] if ftt else [])
        if inst.serial:
            # [NCPU_N] values override the ones derived from nprocs_to_test: the prerequisite is selected
            # by the passes with N processes but runs on min(max_nprocs, N) = 1 process in dynamic mode.
            ncpu_section.append("max_nprocs = 1")
    elif ftt:
        sections["files"].append(f"files_to_test = {ftt}")

    lines = ["<BEGIN TEST_INFO>"]
    for section, body in sections.items():
        if body:
            lines += [f"[{section}]", *body]
    lines += ncpu_section
    lines.append("<END TEST_INFO>")
    return lines


# ---------------------------------------------------------------------------------------------------
# Instances
# ---------------------------------------------------------------------------------------------------

def build_instances(template: Template) -> list[Instance]:
    spec = template.spec
    combos = combinations(spec)

    max_nprocs = spec.test_info.paral_info.get("max_nprocs")
    if max_nprocs is not None:
        nproc = draw_nproc(os.path.basename(template.path), int(max_nprocs))
        combos = [{**c, "nproc": nproc} for c in combos]

    swept = sorted({k for c in combos for k in c if k != "nproc"})
    conflicts = find_conflicts(template.body, swept)
    if conflicts:
        raise ExpansionError(f"swept variables also appear in the body of the template: {conflicts}")
    ndtset, jdtset = scan_datasets(template.body)

    instances = []
    for index, params in enumerate(combos):
        nproc = int(params.get("nproc", 1))
        test_id = f"ti{index:04d}_{template.number}" + (f"_MPI{nproc}" if nproc > 1 else "")
        inst = Instance(index=index, params=params, nproc=nproc, test_id=test_id)
        inst.files = _instance_files(template, test_id)

        header = (f"# parametrized: template = {os.path.basename(template.path)}, instance = {index}, "
                  f"params = {json.dumps(params)}")
        inst.input_text = render_input(template, inst, header, insertion_lines(params, spec.dataset, ndtset, jdtset))
        instances.append(inst)

    return instances


def render_input(template: Template, inst: Instance, header: str, inserted: list[str]) -> str:
    """
    The input of an instance: header comment, body, swept variables, TEST_INFO block.
    The header goes at the end for executables other than abinit: mrgddb and mrgdv read their input
    line by line from stdin, and the first line is the name of the output file.
    """
    body = template.body if template.body.endswith("\n") else template.body + "\n"
    at_top = template.spec.test_info.setup.get("executable") == "abinit"
    text = (header + "\n" + body) if at_top else body
    if inserted:
        text += "\n# Swept variables (inserted by the parametrized expansion)\n" + "\n".join(inserted) + "\n"
    if not at_top:
        text += "\n" + header + "\n"
    text += "\n" + "\n".join(_sentinel(line) for line in compile_test_info(template, inst)) + "\n"
    return text


def _instance_files(template: Template, test_id: str) -> list[tuple[FileToTestSpec, str, str]]:
    files = []
    if not template.spec.test_info.no_check:
        for fspec in template.spec.test_info.files_to_test:
            files.append((fspec, fspec.name.replace("{id}", test_id), fspec.name.replace("{id}", template.stem)))
    return files


def load_prerequisites(template: Template, known_keywords: dict[str, str] | None = None) -> list[Template]:
    """The templates listed in setup.test_chain, in execution order (same directory as `template`)."""
    prereqs = []
    for name in template.spec.test_info.test_chain:
        path = os.path.join(template.directory, name)
        if os.path.abspath(path) == template.path:
            raise ExpansionError(f"setup.test_chain: {name} is the template itself")
        if not os.path.isfile(path):
            raise ExpansionError(f"setup.test_chain: prerequisite template {name} not found in {template.directory}")
        prereqs.append(load_template(path, known_keywords=known_keywords, prerequisite=True))
    return prereqs


def load_fixtures(template: Template, known_keywords: dict[str, str] | None = None) -> dict[str, Template]:
    """
    The fixtures used by `template` and the ones they need, as {name: fixture template}, dependencies first.
    Fixtures are Input/<name>.fixture.abi in the directory of the template.
    """
    loaded: dict[str, Template] = {}
    visiting: list[str] = []

    def visit(name: str) -> None:
        if name in loaded:
            return
        if name in visiting:
            raise ExpansionError(f"cycle in fixture.needs: {' -> '.join(visiting[visiting.index(name):] + [name])}")
        visiting.append(name)
        fx = load_fixture(fixture_path(template.directory, name), known_keywords=known_keywords)
        for need in fx.spec.fixture.needs:
            visit(need)
        visiting.pop()
        loaded[name] = fx

    for name in template.spec.fixtures:
        visit(name)

    producers: dict[str, str] = {}
    for name, fx in loaded.items():
        for out in fx.spec.fixture.outputs:
            if out in producers:
                raise ExpansionError(f"fixtures {producers[out]} and {name} both declare the output {out}")
            producers[out] = name
    return loaded


def build_fixtures(template: Template, loaded: dict[str, Template], make) -> list[Fixture]:
    """The Fixture objects of the fixtures listed by the template (with their needs), sharing the dependencies."""
    objects: dict[str, Fixture] = {}
    for name, fx in loaded.items():   # dependencies first
        inst = Instance(index=-1, params={}, nproc=1, test_id=name)
        header = f"# parametrized: fixture {name} ({os.path.basename(fx.path)})"
        inst.input_text = render_input(fx, inst, header, [])
        test = make(fx, inst, True)
        test.max_nprocs = fx.spec.fixture.nproc   # the build runs on up to nproc processes (static mode checks it)
        objects[name] = Fixture(name=name, path=fx.path, outputs=fx.spec.fixture.outputs,
                                needs=[objects[n] for n in fx.spec.fixture.needs], test=test,
                                nproc=fx.spec.fixture.nproc, nproc_strict=fx.spec.fixture.nproc_strict)
    return [objects[name] for name in template.spec.fixtures]


def build_prerequisites(template: Template, prereqs: list[Template], nproc: int) -> list[Instance]:
    """
    The prerequisites of the chain with `nproc` processes. Id: tp<K>_<NN>, plus _MPI<N>, where K is the
    position in setup.test_chain and NN the number of the swept template.
    """
    out = []
    for k, prereq in enumerate(prereqs, start=1):
        test_id = f"tp{k}_{template.number}" + (f"_MPI{nproc}" if nproc > 1 else "")
        inst = Instance(index=-k, params={}, nproc=nproc, test_id=test_id, serial=True)
        inst.files = _instance_files(prereq, test_id)
        header = (f"# parametrized: prerequisite {k} ({os.path.basename(prereq.path)}) "
                  f"of {os.path.basename(template.path)}")
        inst.input_text = render_input(prereq, inst, header, [])
        out.append(inst)
    return out


def _parse_test_info(text: str, nproc: int):
    """Parse the TEST_INFO block of `text` with the existing AbinitTestInfoParser."""
    fd, path = tempfile.mkstemp(suffix=".abi", prefix="parametrized_")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        return AbinitTestInfoParser(path).generate_testinfo_nprocs(nproc)
    finally:
        os.remove(path)


def _known_keywords(abenv) -> dict[str, str] | None:
    try:
        with open(os.path.join(abenv.tests_dir, "known_keywords.json"), encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return None


def expand_template(path: str, abenv, keywords=None, need_cpp_vars=None, refs_dir: str | None = None,
                    reference: str | None = None, no_check: bool = False) -> list:
    """
    Return the test objects of the template `path`: one test per instance or, if the template lists
    prerequisites in setup.test_chain, one ParametrizedChain per distinct nproc (plan 8.2).

    Args:
        path: Template file (*.tmpl.abi).
        abenv: AbinitEnvironment (tests.abenv).
        keywords, need_cpp_vars: Suite-wide values, as for make_abitests_from_inputs.
        refs_dir: Directory with the golden references. Default: ../Refs relative to the template.
        reference: Explicit golden reference of the main output "{id}.abo" (run_params.py --reference).
        no_check: Build the instances without files to compare (no reference yet).

    Errors never propagate: a template that cannot be expanded becomes one failing test.
    """
    path = os.path.abspath(path)
    try:
        known = _known_keywords(abenv)
        template = load_template(path, known_keywords=known)
        prereqs = load_prerequisites(template, known_keywords=known)
        if no_check:
            template.spec.test_info.setup["no_check"] = "yes"
            for prereq in prereqs:
                prereq.spec.test_info.setup["no_check"] = "yes"
        instances = build_instances(template)
        if reference is not None and not template.spec.test_info.no_check and \
                MAIN_OUTPUT not in [f.name for f in template.spec.test_info.files_to_test]:
            raise ExpansionError(f"--reference replaces the reference of {MAIN_OUTPUT}, which is not in files_to_test")
        if refs_dir is None:
            refs_dir = os.path.join(os.path.dirname(template.directory), "Refs")

        def make(tmpl: Template, inst: Instance, is_prerequisite: bool):
            return _make_test(tmpl, inst, abenv, keywords, need_cpp_vars, refs_dir,
                              None if is_prerequisite else reference, is_prerequisite)

        if not prereqs:
            fixtures = build_fixtures(template, load_fixtures(template, known_keywords=known), make)
            tests = [make(template, inst, False) for inst in instances]
            for t in tests:
                t.fixtures = tuple(fixtures)
            return tests

        # One chain per distinct nproc: the prerequisites followed by the instances with that nproc (plan 8.2).
        chains = []
        for nproc in sorted({inst.nproc for inst in instances}):
            members = [make(p, pi, True) for p, pi in zip(prereqs, build_prerequisites(template, prereqs, nproc))]
            members += [make(template, inst, False) for inst in instances if inst.nproc == nproc]
            chain_id = f"tc_{template.number}" + (f"_MPI{nproc}" if nproc > 1 else "")
            chains.append(ParametrizedChain(members, chain_id=chain_id, nproc=nproc))
        return chains

    except Exception as exc:
        return [_broken_test(path, abenv, exc)]


def _make_test(template: Template, inst: Instance, abenv, keywords, need_cpp_vars, refs_dir: str,
               reference: str | None, is_prerequisite: bool):
    """The test object of an instance (or of a prerequisite) of `template`."""
    test_info = _parse_test_info(inst.input_text, inst.nproc)
    test_info.inp_fname = template.path
    test_info.add_cpp_vars(need_cpp_vars)
    test_info.add_keywords(keywords)

    cls = ParametrizedAbinitTest if exec2class(test_info.executable).__name__ == "AbinitTest" else ParametrizedTest
    test = cls(test_info, abenv)
    test.setup_instance(instance_id=inst.test_id, input_text=inst.input_text, template_path=template.path,
                        index=None if is_prerequisite else inst.index, params=inst.params,
                        suite_refs_dir=refs_dir, is_prerequisite=is_prerequisite)

    files = []
    for fspec, name, ref_name in inst.files:
        if fspec.name == MAIN_OUTPUT and reference is not None:
            ref_name = os.path.abspath(reference)   # absolute: os.path.join keeps it as is
        dic = {"name": name, "tolnlines": fspec.tolnlines, "tolabs": fspec.tolabs,
               "tolrel": fspec.tolrel, "fld_options": fspec.fld_options,
               "verbose_report": fspec.verbose_report}
        overrides = {p: (r if r == "ignore" else r.model_dump(exclude_none=True))
                     for p, r in fspec.overrides.items()}
        files.append(ParamFileToTest(dic, mode=fspec.mode, overrides=overrides, ref_name=ref_name))
    test.files_to_test = files
    return test


def _broken_test(path: str, abenv, exc: Exception) -> BrokenTemplateTest:
    kind = "" if isinstance(exc, (TemplateError, ExpansionError)) else f"{type(exc).__name__}: "
    text = "\n".join(_sentinel(line) for line in [
        "<BEGIN TEST_INFO>", "[setup]", "executable = abinit", "no_check = yes",
        "[extra_info]", f"keywords = {PARAMETRIZED_KEYWORD}", "description = Template that could not be expanded",
        "<END TEST_INFO>"]) + "\n"
    test_info = _parse_test_info(text, 1)
    test_info.inp_fname = path
    test = BrokenTemplateTest(test_info, abenv)
    test.set_error(f"{kind}{exc}")
    return test
