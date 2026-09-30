"""
Pydantic models validating the %YAML block of a parametrized template.

The test_info section is a one-to-one YAML transcription of the INI TEST_INFO block:
the allowed sections and keys are taken live from TESTCNF_SECTIONS/TESTCNF_KEYWORDS.
"""
from __future__ import annotations

import re
from typing import Annotated, Any, Literal, Union

from pydantic import BaseModel, ConfigDict, Field, PositiveInt, field_validator, model_validator

from ..testsuite import TESTCNF_KEYWORDS, TESTCNF_SECTIONS, TESTBOT_NP_LISTS

Scalar = Union[int, float, str]
Value = Union[Scalar, list[Scalar]]

#: Process counts that testbot.py can reach (see plan section 9.1).
REACHABLE_NPROCS: frozenset[int] = frozenset({1}).union(*(set(v) for v in TESTBOT_NP_LISTS.values()))

#: Sections that can appear in test_info. [yaml_test] and [NCPU_N] are not supported.
ALLOWED_SECTIONS = tuple(s for s in TESTCNF_SECTIONS if s != "yaml_test")

#: Keys computed by the expansion, or not supported.
FORBIDDEN_KEYS = {
    "nprocs_to_test": "computed from nproc",
    "exclude_nprocs": "computed from nproc",
    "psp_files": "pseudopotentials are fixed and declared in the body with pseudos/pp_dirpath",
    "use_files_file": "templates use the command-line invocation",
}

#: Variables that cannot be swept: pseudopotentials are fixed.
UNSWEEPABLE = {"pseudos", "pp_dirpath"}

#: Name of the main output of an instance in files_to_test. Its reference is Refs/tNN.abo.
MAIN_OUTPUT = "{id}.abo"

#: TEST_INFO keys whose meaning or syntax differs in a template (used by docs.py).
TEMPLATE_KEY_NOTES = {
    "files_to_test": "a YAML list of entries (see files_to_test below), not the INI string",
    "test_chain": "the prerequisite templates (*.tmpl.abi in the same directory), in execution order; "
                  "expansion builds one chain per distinct nproc and runs the prerequisites on 1 process",
    "max_nprocs": "shorthand: ONE process count <= max_nprocs, drawn from the template name, for all the instances; "
                  "not an upper bound as in TEST_INFO, and incompatible with an nproc axis (plan 4.4, P18 open)",
    "no_check": "yes: the instances are run but not compared; the default of a prerequisite without files_to_test",
}

#: Executables that abort on more than one MPI process: their fixtures must have nproc 1.
SERIAL_EXECUTABLES = {"mrgddb", "mrgdv"}

#: Name of a fixture: the file Input/<name>.fixture.abi, referenced by `fixtures` and `needs`.
#: Lowercase only: on a case-insensitive filesystem (default on macOS) `MgO_DDB` would find mgo_ddb.fixture.abi,
#: and fail on Linux, and two names differing only by case could not coexist in a checkout.
_FIXTURE_NAME_RE = re.compile(r"^[a-z][a-z0-9_]*$")

_PATH_RE = re.compile(r"^[A-Za-z][\w]*(\.[A-Za-z_][\w ]*)?$")
_FORBID = ConfigDict(extra="forbid")


class ToleranceSpec(BaseModel):
    """Tolerances of an override: a subset of the files_to_test tolerances. Missing values are inherited from the entry."""
    model_config = _FORBID
    tolnlines: int | None = Field(None, description="Maximum number of differing lines in the group.")
    tolabs: float | None = Field(None, description="Maximum absolute difference in the group.")
    tolrel: float | None = Field(None, description="Maximum relative difference in the group.")
    fld_options: str | None = Field(None, description="fldiff options for the group, e.g. -medium or -easy.")


class FileToTestSpec(BaseModel):
    """
    An entry of test_info.files.files_to_test: a file produced by every instance, compared with its golden
    reference in Refs/.
    """
    model_config = _FORBID
    name: str = Field(
        description="File produced by the instance, e.g. \"{id}.abo\" (the main output) or \"{id}o_DS2_DDB\". "
                    "It must contain {id}, the only placeholder: {id} becomes the instance id in the output "
                    "and the template stem in the reference, e.g. Refs/t01.abo and Refs/t01o_DS2_DDB.")
    mode: Literal["fldiff", "yaml_docs"] = Field(
        "fldiff", description="fldiff: the whole file is compared by fldiff. yaml_docs: only the YAML documents "
                              "are extracted as text and compared by fldiff, with per-path overrides.")
    tolnlines: int = Field(0, description="Maximum number of differing lines (text differences count too).")
    tolabs: float = Field(0.0, description="Maximum absolute difference.")
    tolrel: float = Field(0.0, description="Maximum relative difference.")
    fld_options: str = Field(
        "", description="fldiff options, e.g. -medium (lines differing by < 1.01e-8 are equal) or -easy (1.01e-5). "
                        "A line is flagged only if both its absolute and relative differences exceed the threshold.")
    verbose_report: Literal["yes", "no"] = Field("no", description="Print a detailed fldiff report.")
    overrides: dict[str, Literal["ignore"] | ToleranceSpec] = Field(
        default_factory=dict,
        description="yaml_docs only. Maps a document tag (e.g. ResultsGS) or tag.key (e.g. ResultsGS.convergence) "
                    "to ignore or to its own tolerances. DatasetInfo and BeginCycle.solver are ignored by default.")

    @field_validator("mode", mode="before")
    @classmethod
    def _check_mode(cls, v: Any) -> Any:
        if v == "yaml":
            raise ValueError("mode 'yaml' is obsolete and not supported: use 'yaml_docs' or 'fldiff'")
        return v

    @field_validator("fld_options")
    @classmethod
    def _check_fld_options(cls, v: str) -> str:
        for opt in v.split():
            if not opt.startswith("-"):
                raise ValueError(f"wrong fldiff option: {opt}")
        return v

    @field_validator("name")
    @classmethod
    def _check_name(cls, v: str) -> str:
        if re.search(r"\{(?!id\})[^}]*\}", v):
            raise ValueError("the only placeholder allowed in a file name is {id}")
        if "{id}" not in v or "/" in v:
            raise ValueError(f"file name {v!r}: it must contain {{id}} (e.g. \"{{id}}.abo\") and no directory")
        return v

    @field_validator("overrides")
    @classmethod
    def _check_paths(cls, v: dict[str, Any]) -> dict[str, Any]:
        for path in v:
            if not _PATH_RE.match(path):
                raise ValueError(f"invalid override path {path!r}: expected Tag or Tag.key")
        return v

    @model_validator(mode="after")
    def _overrides_need_yaml_docs(self) -> FileToTestSpec:
        if self.overrides and self.mode != "yaml_docs":
            raise ValueError("overrides are allowed only with mode: yaml_docs")
        return self


class TestInfoSpec(BaseModel):
    """
    YAML transcription of the INI TEST_INFO sections: same section names, same keys, same meaning
    (keys and values are checked against the TEST_INFO options, see the list below).
    """
    model_config = _FORBID
    setup: dict[str, Any] = Field(default_factory=dict, description="[setup] options. executable is mandatory.")
    files: dict[str, Any] = Field(default_factory=dict,
                                  description="[files] options, including files_to_test (a list of entries).")
    shell: dict[str, Any] = Field(default_factory=dict, description="[shell] options (pre_commands, post_commands).")
    paral_info: dict[str, Any] = Field(default_factory=dict, description="[paral_info] options.")
    extra_info: dict[str, Any] = Field(default_factory=dict,
                                       description="[extra_info] options. keywords must be registered in "
                                                   "tests/known_keywords.json; PARAMETRIZED is added automatically.")

    # files_to_test is extracted from `files` and validated with its own model.
    files_to_test: list[FileToTestSpec] = Field(default_factory=list, exclude=True,
                                                description="Written under files: in the template (see above).")

    @model_validator(mode="before")
    @classmethod
    def _check_sections(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            raise ValueError("test_info must be a mapping")
        for section in data:
            if section == "yaml_test":
                raise ValueError("[yaml_test] is not supported: use files_to_test with mode/overrides")
            if section.upper().startswith("NCPU_"):
                raise ValueError(f"[{section}] cannot be written in a template: it is generated from nproc")
            if section not in ALLOWED_SECTIONS:
                raise ValueError(f"unknown test_info section {section!r}, expected one of {ALLOWED_SECTIONS}")
        files = dict(data.get("files", {}) or {})
        data = dict(data)
        data["files_to_test"] = files.pop("files_to_test", None) or []
        data["files"] = files
        return data

    @model_validator(mode="after")
    def _check_keys(self) -> TestInfoSpec:
        for section in ALLOWED_SECTIONS:
            for key, value in getattr(self, section).items():
                if key in FORBIDDEN_KEYS:
                    raise ValueError(f"{section}.{key} is not allowed in a template: {FORBIDDEN_KEYS[key]}")
                if key not in TESTCNF_KEYWORDS:
                    raise ValueError(f"unknown TEST_INFO option {section}.{key}")
                if TESTCNF_KEYWORDS[key][2] != section:
                    raise ValueError(f"{key} belongs to section [{TESTCNF_KEYWORDS[key][2]}], not [{section}]")
                if key == "files_to_test":
                    continue
                # The value must be accepted by the same converter used by the INI parser.
                converter = TESTCNF_KEYWORDS[key][0]
                try:
                    converter(ini_value(key, value))
                except Exception as exc:
                    raise ValueError(f"invalid value for {section}.{key}: {value!r} ({exc})") from exc

        if "executable" not in self.setup:
            raise ValueError("setup.executable is mandatory")

        chain = self.setup.get("test_chain")
        if chain is not None:
            names = chain if isinstance(chain, list) else [chain]
            if not names:
                raise ValueError("setup.test_chain must list at least one prerequisite template")
            for name in names:
                if not isinstance(name, str) or "/" in name or not name.endswith(".tmpl.abi"):
                    raise ValueError(f"setup.test_chain: {name!r} must be the file name of a template "
                                     "(*.tmpl.abi) in the same directory")
            if len(set(names)) != len(names):
                raise ValueError(f"setup.test_chain contains repeated templates: {names}")

        names = [f.name for f in self.files_to_test]
        if len(set(names)) != len(names):
            raise ValueError(f"files_to_test lists the same file more than once: {names}")
        if not self.no_check and not names:
            raise ValueError("files_to_test must list at least one file (or set setup.no_check: yes)")
        return self

    @property
    def test_chain(self) -> list[str]:
        """The prerequisite templates (setup.test_chain), in execution order."""
        chain = self.setup.get("test_chain") or []
        return list(chain) if isinstance(chain, list) else [chain]

    @property
    def no_check(self) -> bool:
        """no_check interpreted with the same converter as the INI parser."""
        value = self.setup.get("no_check", "no")
        return bool(TESTCNF_KEYWORDS["no_check"][0](ini_value("no_check", value)))


class GridAxis(BaseModel):
    """An axis built from lists of values, e.g. {type: grid, parameters: {useylm: [0, 1]}}."""
    model_config = _FORBID
    type: Literal["grid"] = Field(description="Selects this kind of axis.")
    mode: Literal["product", "zip"] = Field(
        "product", description="product: all the combinations of the lists. zip: the i-th values of the lists together "
                               "(the lists must have the same length).")
    parameters: dict[str, list[Value]] = Field(
        description="ABINIT variable (or nproc) -> list of values. A value can be a list, e.g. for eph_np_pqbks.")

    @model_validator(mode="after")
    def _check(self) -> GridAxis:
        if not self.parameters:
            raise ValueError("grid axis without parameters")
        for name, values in self.parameters.items():
            if not values:
                raise ValueError(f"grid parameter {name} has no values")
        if self.mode == "zip" and len({len(v) for v in self.parameters.values()}) != 1:
            raise ValueError("zip grid: all the lists must have the same length")
        return self

    def keys(self) -> set[str]:
        return set(self.parameters)


class ExplicitAxis(BaseModel):
    """An axis given as a list of combinations, for variables that must change together (e.g. nproc and npband)."""
    model_config = _FORBID
    type: Literal["explicit"] = Field(description="Selects this kind of axis.")
    values: list[dict[str, Value]] = Field(
        description="One mapping per point of the axis, e.g. {nproc: 4, npkpt: 2, npband: 2, npfft: 1}.")

    @model_validator(mode="after")
    def _check(self) -> ExplicitAxis:
        if not self.values:
            raise ValueError("explicit axis without values")
        return self

    def keys(self) -> set[str]:
        return set().union(*(set(v) for v in self.values))


Axis = Annotated[Union[GridAxis, ExplicitAxis], Field(discriminator="type")]


class TemplateSpec(BaseModel):
    """
    The %YAML ... %ENDYAML block at the top of a template tNN.tmpl.abi. The rest of the file is the ABINIT input
    (the body). The instances are the combinations of the axes (Cartesian product across axes, minus exclude):
    the swept variables are appended to the body, and each instance is compared with the same reference.
    """
    model_config = _FORBID
    test_info: TestInfoSpec = Field(description="The TEST_INFO of every instance, written in YAML.")
    dataset: PositiveInt | Literal["all"] | None = Field(
        None, description="Mandatory when the body has ndtset > 1: the dataset k whose variables are swept "
                          "(var<k> is appended), or all (var without suffix, for all the datasets).")
    axes: list[Axis] = Field(default_factory=list,
                             description="The sweep. nproc is a reserved name: the number of MPI processes of the "
                                         "instance, not an ABINIT variable. No axes: a single instance.")
    exclude: list[dict[str, Value]] = Field(
        default_factory=list, description="Combinations to drop: an instance is dropped if it matches all the "
                                          "values of an entry, e.g. {useylm: 1, nproc: 4}.")
    fixtures: list[str] = Field(
        default_factory=list,
        description="Names of the fixtures (Input/<name>.fixture.abi) whose declared outputs every instance needs, "
                    "e.g. [mgo_ddb_444q, mgo_dvdb_444q, mgo_wfk_444k_70b]. A fixture runs once per workdir root, before the tests, on fixture.nproc processes, and "
                    "its outputs are linked read-only into the workdir of each instance. The fixtures they need "
                    "are added automatically. Cannot be combined with setup.test_chain.")

    @model_validator(mode="after")
    def _check(self) -> TemplateSpec:
        seen: dict[str, int] = {}
        for i, axis in enumerate(self.axes):
            for key in axis.keys():
                if key in seen:
                    raise ValueError(f"parameter {key!r} appears in axes {seen[key]} and {i}")
                seen[key] = i
        swept = set(seen)

        bad = swept & UNSWEEPABLE
        if bad:
            raise ValueError(f"pseudopotential variables cannot be swept: {sorted(bad)}")

        for axis in self.axes:
            values = (axis.parameters.get("nproc", []) if isinstance(axis, GridAxis)
                      else [v["nproc"] for v in axis.values if "nproc" in v])
            for n in values:
                if not isinstance(n, int) or isinstance(n, bool) or n not in REACHABLE_NPROCS:
                    raise ValueError(f"nproc={n!r}: must be one of {sorted(REACHABLE_NPROCS)} (reachable by testbot)")

        if "max_nprocs" in self.test_info.paral_info and "nproc" in swept:
            raise ValueError("paral_info.max_nprocs and an nproc axis cannot be used together")

        for entry in self.exclude:
            unknown = set(entry) - swept
            if unknown:
                raise ValueError(f"exclude uses parameters not defined in any axis: {sorted(unknown)}")

        _check_fixture_names(self.fixtures, "fixtures")
        if self.fixtures and self.test_info.test_chain:
            raise ValueError("fixtures and setup.test_chain cannot be used together")
        return self

    def check_prerequisite(self) -> None:
        """
        A prerequisite (an entry of setup.test_chain of another template) runs once per chain, always on 1 process:
        no axes, no process counts, no chain of its own (plan 8.2).
        """
        if self.axes or self.exclude or self.dataset is not None:
            raise ValueError("a prerequisite template cannot have axes, exclude or dataset")
        if self.test_info.test_chain:
            raise ValueError("a prerequisite template cannot have its own setup.test_chain")
        if "max_nprocs" in self.test_info.paral_info:
            raise ValueError("a prerequisite template always runs on 1 process: remove paral_info.max_nprocs")
        if self.fixtures:
            raise ValueError("a prerequisite template cannot use fixtures")


def check_fixture_name(name: str, what: str = "fixture") -> None:
    """A fixture name has lowercase letters, digits and _, and starts with a letter."""
    if not _FIXTURE_NAME_RE.match(name):
        raise ValueError(f"{what}: invalid fixture name {name!r}: use lowercase letters, digits and _, starting with "
                         "a letter (the file is Input/<name>.fixture.abi; names differing only by case are not portable)")


def _check_fixture_names(names: list[str], what: str) -> None:
    for name in names:
        check_fixture_name(name, what)
    if len(set(names)) != len(names):
        raise ValueError(f"{what} contains repeated fixtures: {names}")


class FixtureBlock(BaseModel):
    """The fixture: block of a fixture file (Input/<name>.fixture.abi)."""
    model_config = _FORBID
    outputs: list[str] = Field(
        description="Files produced by the fixture that are linked into the workdir of its users, e.g. [MgO_444q_DDB]. "
                    "Literal names in the workdir of the fixture: set them in the body (mrgddb, mrgdv) or with "
                    "setup.output_prefix (abinit). The other files it writes stay in the cache. The names may "
                    "describe the content for readability (e.g. 444q for the q-mesh): they are never parsed.")
    needs: list[str] = Field(
        default_factory=list,
        description="Fixtures whose outputs this fixture reads, e.g. an NSCF run that needs the DEN of a GS fixture. "
                    "They are built first and their outputs are linked into its build directory.")
    nproc: PositiveInt = Field(
        1, description="Number of MPI processes of the fixture, fixed by the fixture and not by the tests using it. "
                       "Fixtures are built before the tests, within the CPUs of the run. If the run cannot use nproc "
                       "processes (no MPI, fewer CPUs), the fixture runs on fewer processes, with a warning: its results "
                       "must not depend on the number of processes beyond the tolerances of its users. Prefer inputs "
                       "without distribution variables (npkpt, ...), so that any number of processes works.")
    nproc_strict: bool = Field(
        False, description="yes: never run on fewer than nproc processes; if the run cannot provide them, the fixture "
                           "and the tests using it are skipped.")

    @model_validator(mode="after")
    def _check(self) -> FixtureBlock:
        if not self.outputs:
            raise ValueError("fixture.outputs must list at least one file")
        for name in self.outputs:
            if "/" in name or "{" in name or not name.strip():
                raise ValueError(f"fixture.outputs: {name!r} must be a plain file name (no directory, no placeholder)")
        if len(set(self.outputs)) != len(self.outputs):
            raise ValueError(f"fixture.outputs contains repeated files: {self.outputs}")
        _check_fixture_names(self.needs, "fixture.needs")
        return self


class FixtureSpec(BaseModel):
    """
    The %YAML block of a fixture file Input/<name>.fixture.abi: a step that prepares files for parametrized
    templates (e.g. mrgddb, mrgdv, an NSCF run). The rest of the file is the input of its executable.
    A fixture is built at most once per workdir root and cache key (hash of its input, of its executable, of the
    commit of its git submodule, of its nproc and of the fixtures it needs), even by concurrent processes,
    before the tests and on fixture.nproc processes.
    """
    model_config = _FORBID
    test_info: TestInfoSpec = Field(
        description="As for a template, without files_to_test (a fixture is not compared: no_check is always yes), "
                    "paral_info.max_nprocs and setup.test_chain.")
    fixture: FixtureBlock = Field(description="Outputs and dependencies of the fixture.")

    @model_validator(mode="after")
    def _check(self) -> FixtureSpec:
        if self.test_info.files_to_test:
            raise ValueError("a fixture is not compared: remove files_to_test")
        if "max_nprocs" in self.test_info.paral_info:
            raise ValueError("the number of processes of a fixture is fixture.nproc: remove paral_info.max_nprocs")
        if self.test_info.test_chain:
            raise ValueError("a fixture cannot have setup.test_chain: use fixture.needs")
        exe = self.test_info.setup.get("executable")
        if exe in SERIAL_EXECUTABLES and self.fixture.nproc != 1:
            raise ValueError(f"{exe} runs on 1 process only: fixture.nproc must be 1")
        return self

    @property
    def dataset(self) -> None:
        """Fixtures sweep nothing (duck typing with TemplateSpec for the input rendering)."""
        return None


def ini_value(key: str, value: Any) -> str:
    """Convert a YAML value to the string used in the INI TEST_INFO block."""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (list, tuple)):
        sep = "; " if key in ("pre_commands", "post_commands") else ", "
        return sep.join(str(v) for v in value)
    return str(value)


def check_keywords(spec: TemplateSpec | FixtureSpec, known_keywords: dict[str, str] | None) -> None:
    """Keywords must be registered in known_keywords.json and must not contain blanks."""
    keywords = spec.test_info.extra_info.get("keywords", []) or []
    if isinstance(keywords, str):
        keywords = [k.strip() for k in keywords.split(",") if k.strip()]
    for kw in keywords:
        if " " in kw:
            raise ValueError(f"keyword {kw!r} contains blanks")
        if known_keywords is not None and kw not in known_keywords:
            raise ValueError(f"keyword {kw!r} is not registered in tests/known_keywords.json")
