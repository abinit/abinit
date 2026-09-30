"""
Read a parametrized template: split the %YAML block from the ABINIT body and validate the block.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass

import yaml

from .models import FixtureSpec, TemplateSpec, check_fixture_name, check_keywords

TEMPLATE_SUFFIX = ".tmpl.abi"
FIXTURE_SUFFIX = ".fixture.abi"


class TemplateError(Exception):
    """Raised when a template cannot be parsed or validated."""


class _StrictLoader(yaml.SafeLoader):
    """SafeLoader that rejects duplicate keys (plain safe_load silently keeps the last value)."""


def _construct_mapping_no_dups(loader, node, deep=False):
    keys = [loader.construct_object(k, deep=deep) for k, _ in node.value]
    dups = sorted({str(k) for k in keys if keys.count(k) > 1})
    if dups:
        raise yaml.constructor.ConstructorError(None, None, f"duplicate keys {dups}", node.start_mark)
    return yaml.SafeLoader.construct_mapping(loader, node, deep)


_StrictLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping_no_dups)


def load_yaml_strict(text: str):
    return yaml.load(text, Loader=_StrictLoader)


@dataclass
class Template:
    path: str
    stem: str
    #: Number of the template in the suite (01 in t01.tmpl.abi), or the stem if the name is not tNN.
    #: For a fixture: its name.
    number: str
    spec: TemplateSpec | FixtureSpec
    #: ABINIT input without the %YAML block.
    body: str

    @property
    def directory(self) -> str:
        return os.path.dirname(self.path)


def split_template(text: str) -> tuple[str, str]:
    """Return (yaml_source, body). The %YAML ... %ENDYAML block is removed from the body."""
    lines = text.splitlines(keepends=True)
    starts = [i for i, l in enumerate(lines) if l.strip() == "%YAML"]
    ends = [i for i, l in enumerate(lines) if l.strip() == "%ENDYAML"]
    if len(starts) != 1 or len(ends) != 1 or ends[0] < starts[0]:
        raise TemplateError("a template must contain exactly one %YAML ... %ENDYAML block")
    i, j = starts[0], ends[0]
    body = "".join(lines[:i] + lines[j + 1:])
    if any(l.startswith("#%%") for l in body.splitlines()):
        raise TemplateError("a template must not contain a #%% TEST_INFO block: use test_info in the %YAML block")
    return "".join(lines[i + 1:j]), body


def template_stem(path: str) -> str:
    base = os.path.basename(path)
    if not base.endswith(TEMPLATE_SUFFIX):
        raise TemplateError(f"{path}: a template name must end with {TEMPLATE_SUFFIX}")
    return base[: -len(TEMPLATE_SUFFIX)]


def load_template(path: str, known_keywords: dict[str, str] | None = None, prerequisite: bool = False) -> Template:
    """
    Read and validate a template. Raises TemplateError with a readable message.

    prerequisite: the template is an entry of the setup.test_chain of another template (plan 8.2).
        Its files are not compared unless it lists files_to_test: no_check defaults to yes.
    """
    path = os.path.abspath(path)
    stem = template_stem(path)
    m = re.match(r"^t(\d+)$", stem)
    number = m.group(1) if m else stem

    with open(path, encoding="utf-8") as fh:
        yaml_src, body = split_template(fh.read())

    try:
        data = load_yaml_strict(yaml_src)
    except yaml.YAMLError as exc:
        raise TemplateError(f"{path}: invalid YAML in the %YAML block:\n{exc}") from exc
    if not isinstance(data, dict):
        raise TemplateError(f"{path}: the %YAML block must be a mapping")

    if prerequisite:
        ti = data.get("test_info")
        if isinstance(ti, dict):
            setup, files = ti.get("setup") or {}, ti.get("files") or {}
            if isinstance(setup, dict) and isinstance(files, dict) and "no_check" not in setup \
                    and not files.get("files_to_test"):
                data["test_info"] = {**ti, "setup": {**setup, "no_check": "yes"}}

    try:
        spec = TemplateSpec.model_validate(data)
        if prerequisite:
            spec.check_prerequisite()
        check_keywords(spec, known_keywords)
    except Exception as exc:
        raise TemplateError(f"{path}: invalid template:\n{exc}") from exc

    return Template(path=path, stem=stem, number=number, spec=spec, body=body)


def fixture_path(directory: str, name: str) -> str:
    return os.path.join(directory, name + FIXTURE_SUFFIX)


def load_fixture(path: str, known_keywords: dict[str, str] | None = None) -> Template:
    """Read and validate a fixture file <name>.fixture.abi. Raises TemplateError with a readable message."""
    path = os.path.abspath(path)
    base = os.path.basename(path)
    if not base.endswith(FIXTURE_SUFFIX):
        raise TemplateError(f"{path}: a fixture name must end with {FIXTURE_SUFFIX}")
    name = base[: -len(FIXTURE_SUFFIX)]
    try:
        check_fixture_name(name, "fixture file")
    except ValueError as exc:
        raise TemplateError(f"{path}: {exc}") from exc
    if not os.path.isfile(path):
        raise TemplateError(f"fixture {name}: file {path} not found")

    with open(path, encoding="utf-8") as fh:
        yaml_src, body = split_template(fh.read())
    try:
        data = load_yaml_strict(yaml_src)
    except yaml.YAMLError as exc:
        raise TemplateError(f"{path}: invalid YAML in the %YAML block:\n{exc}") from exc
    if not isinstance(data, dict):
        raise TemplateError(f"{path}: the %YAML block must be a mapping")

    # A fixture is never compared: no_check is always yes.
    ti = data.get("test_info")
    if isinstance(ti, dict) and isinstance(ti.get("setup") or {}, dict):
        setup = dict(ti.get("setup") or {})
        if str(setup.get("no_check", "yes")).lower() not in ("yes", "true"):
            raise TemplateError(f"{path}: a fixture is not compared: remove setup.no_check")
        setup["no_check"] = "yes"
        data["test_info"] = {**ti, "setup": setup}

    try:
        spec = FixtureSpec.model_validate(data)
        check_keywords(spec, known_keywords)
    except Exception as exc:
        raise TemplateError(f"{path}: invalid fixture:\n{exc}") from exc

    return Template(path=path, stem=name, number=name, spec=spec, body=body)
