"""
Test objects for the instances of a parametrized template.

The instances are ordinary BaseTest/AbinitTest objects with a few overrides:

- the generated input is stored on the object and written into the workdir when the test runs;
- the golden reference of the template is linked into ``<workdir>/refs/`` under the name expected by fldiff;
- each file to test has its own comparison mode (``fldiff`` or ``yaml_docs``, see ParamFileToTest).

This module is imported when a pickled test database is loaded, so it must depend only on the
standard library and on the existing testsuite modules: never on pydantic or yaml.
"""
from __future__ import annotations

import copy
import hashlib
import os
import re
import shutil
import subprocess
import sys
from typing import Any

from ..jobrunner import JobRunner
from ..testsuite import AbinitTest, BaseTest, ChainOfTests, FileToTest
from ..tools import unzip
from . import fixture_cache
from .yamltext import EXTRACT_OUTPUT, EXTRACT_REFERENCE, compare_yaml_text, merge_overrides


class ParamFileToTest(FileToTest):
    """
    A file to test with its own comparison mode.

    mode == "fldiff": the standard whole-file comparison of FileToTest.
    mode == "yaml_docs": the YAML documents are compared as text, with per-path overrides (yamltext).
    """

    def __init__(self, dic: dict[str, Any], mode: str = "fldiff", overrides: dict[str, Any] | None = None,
                 ref_name: str | None = None) -> None:
        super().__init__(dic)
        if mode not in ("fldiff", "yaml_docs"):
            raise ValueError(f"Invalid parametrized mode: {mode}")
        self.param_mode = mode
        self.overrides = dict(overrides or {})
        #: yaml_docs: result of each group of lines (path, tolerances, status, msg), for the report.
        self.yaml_groups: list[dict[str, Any]] = []
        # Name of the golden reference in the Refs directory of the suite.
        self.ref_name = ref_name or self.name

    def compare(self, fldiff_path: str, ref_dir: str, workdir: str, yaml_test: dict[str, Any],
                timebomb: Any = None, outf: Any = sys.stdout, simplified_yaml_test: bool = False,
                forced_tolerance: str = "default") -> tuple[bool, str, str]:
        if self.param_mode == "fldiff":
            return super().compare(fldiff_path, ref_dir, workdir, yaml_test, timebomb=timebomb, outf=outf,
                                   simplified_yaml_test=False, forced_tolerance=forced_tolerance)

        ref_path = os.path.abspath(os.path.join(ref_dir, self.name))
        out_path = os.path.abspath(os.path.join(workdir, self.name))
        base = {"tolnlines": self.tolnlines, "tolabs": self.tolabs, "tolrel": self.tolrel,
                "fld_options": " ".join(self.fld_options)}
        try:
            res = compare_yaml_text(ref_path, out_path, base, merge_overrides(self.overrides),
                                    workdir=os.path.join(workdir, self.name + ".yaml_docs"),
                                    fldiff_path=fldiff_path, outf=outf)
            isok, status, msg, line_count_error = res.isok, res.status, res.msg, res.has_line_count_error
            self.yaml_groups = [{"path": g.path or "default", "tolerances": g.tolerances, "status": g.status,
                                 "msg": g.msg} for g in res.groups]
        except Exception as exc:
            isok, status, line_count_error = False, "failed", False
            msg = f"Internal error in yaml_docs comparison:\n{type(exc).__name__}: {exc!s}"

        msg += f" [file={self.name}]"
        self.fld_isok, self.fld_status, self.fld_msg = isok, status, msg
        self.has_line_count_error = line_count_error
        return isok, status, msg


_MPI_NPROC_RE = re.compile(r"\bmpi_nproc:\s*(\d+)")


def mpi_nprocs_from_output(path: str) -> int | None:
    """Number of MPI processes reported by ABINIT in its main output (`mpi_nproc: N`), None if not found."""
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                m = _MPI_NPROC_RE.search(line)
                if m:
                    return int(m.group(1))
    except OSError:
        pass
    return None


def _submodule_revision(tests_dir: str, module: str) -> str:
    """Commit of a git submodule (part of the cache key of the fixtures reading it); empty if unknown."""
    if not module:
        return ""
    path = os.path.join(tests_dir, "modules_with_data", module)
    try:
        out = subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True, text=True, timeout=30)
        return out.stdout.strip() if out.returncode == 0 else ""
    except (OSError, subprocess.SubprocessError):
        return ""


class Fixture:
    """
    A fixture (Input/<name>.fixture.abi, plan 8.4): a step that prepares files for the instances of templates.
    It is pickled with the instances that use it, and built on demand by the first instance that runs,
    through the process-safe cache of fixture_cache, always on 1 process.
    """

    def __init__(self, name: str, path: str, outputs: list[str], needs: list[Fixture], test: BaseTest,
                 nproc: int = 1, nproc_strict: bool = False) -> None:
        self.name = name
        self.path = path
        self.outputs = list(outputs)
        self.needs = list(needs)
        self.test = test
        self.nproc = nproc
        self.nproc_strict = nproc_strict

    def effective_nproc(self, runner, max_cpus: int | None) -> int | None:
        """
        Number of MPI processes of the build: nproc, or fewer if the run cannot provide them
        (no MPI launcher: 1; max_cpus: the CPU budget of the run). None if nproc_strict forbids it.
        """
        if self.nproc == 1:
            return 1
        can_mpi = runner is not None and any(getattr(runner, a, False) for a in ("has_mpirun", "has_srun", "has_poe"))
        n = min(self.nproc, max_cpus) if (can_mpi and max_cpus) else (self.nproc if can_mpi else 1)
        if n != self.nproc and self.nproc_strict:
            return None
        return n

    def all_fixtures(self) -> list[Fixture]:
        """This fixture and the ones it needs, dependencies first."""
        out: list[Fixture] = []
        for n in self.needs:
            out += [f for f in n.all_fixtures() if f not in out]
        return out + [self]

    def cache_key(self, build_env) -> str:
        """
        Hash of the input, of the executable (path, size, mtime), of the submodule commit and of the needs.
        The requested nproc is part of it, the number of processes actually used is not: the first build (normally the one before
        the tests, see prebuild.py) decides it, records it in status.json, and every user shares that entry.
        """
        h = hashlib.sha256()
        h.update(self.test.input_text.encode())
        h.update(f"nproc={self.nproc}".encode())   # the requested number: changing it gives a new entry
        exe = build_env.path_of_bin(self.test.executable)
        try:
            st = os.stat(exe)
            h.update(f"{exe}:{st.st_size}:{st.st_mtime_ns}".encode())
        except OSError:
            h.update(exe.encode())
        h.update(_submodule_revision(self.test.abenv.tests_dir, self.test.use_git_submodule).encode())
        for n in self.needs:
            h.update(n.cache_key(build_env).encode())
        return h.hexdigest()

    def ensure(self, build_env, root: str, runner=None, max_cpus: int | None = None,
               timebomb=None) -> tuple[str, dict[str, Any]]:
        """
        Return (cache entry, status) of the fixture, building it (and what it needs) if necessary.
        runner: JobRunner able to launch MPI jobs (needed for nproc > 1). max_cpus: CPU budget of the run.
        """
        need_dirs = []
        for n in self.needs:
            d, st = n.ensure(build_env, root, runner=runner, max_cpus=max_cpus, timebomb=timebomb)
            if st["status"] != "ok":
                # Not cached under our key: the failure belongs to the needed fixture.
                return d, {**st, "message": f"needed fixture {n.name} {st['status']}: {st['message']}"}
            need_dirs.append((n, d))

        key = self.cache_key(build_env)
        nproc = self.effective_nproc(runner, max_cpus)
        if nproc is None:
            # Not cached: another run (more CPUs, with MPI) can build it.
            entry = os.path.join(root, fixture_cache.entry_name(self.name, key))
            return entry, {"name": self.name, "status": "skipped", "nproc_requested": self.nproc,
                           "message": f"nproc_strict: the run cannot provide {self.nproc} MPI processes"}

        def build(tmp: str) -> tuple[str, str, dict[str, Any]]:
            extra = {"nproc_requested": self.nproc, "nproc": nproc}
            for n, d in need_dirs:
                n.link_outputs(d, tmp)
            test = copy.deepcopy(self.test)
            job_runner = runner if nproc > 1 else JobRunner.sequential(timebomb=timebomb)
            test.run(build_env, job_runner, tmp, mpi_nprocs=nproc, runmode="static", erase_files=0, make_html_diff=0)
            if test._status == "skipped":
                return "skipped", (test.skip_msg or "skipped").strip(), extra
            errors = [str(e) for e in test.exceptions]
            missing = [o for o in self.outputs if not os.path.exists(os.path.join(tmp, o))]
            if missing:
                errors.append(f"declared outputs not produced: {', '.join(missing)} (see {tmp})")
            if test.executable == "abinit":
                ran = mpi_nprocs_from_output(os.path.join(tmp, test.id + ".abo"))
                extra["nproc_run"] = ran
                if ran is not None and ran != nproc:
                    errors.append(f"the output reports mpi_nproc: {ran}, expected {nproc}")
            if errors or getattr(test, "exec_error", False):
                return "failed", "; ".join(errors) or "execution failed", extra
            message = "" if nproc == self.nproc else f"ran on {nproc} of the {self.nproc} requested MPI processes"
            return "ok", message, extra

        return fixture_cache.ensure(root, self.name, key, build)

    def link_outputs(self, entry: str, workdir: str) -> None:
        """Symlink the declared outputs of the cache entry into workdir."""
        for o in self.outputs:
            link = os.path.join(workdir, o)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.join(entry, o), link)


def all_fixtures(fixtures: list[Fixture]) -> list[Fixture]:
    """The fixtures and everything they need, dependencies first, without repetitions."""
    out: list[Fixture] = []
    for fx in fixtures:
        out += [f for f in fx.all_fixtures() if f not in out]
    return out


class ParametrizedMixin:
    """Behaviour shared by the instances of a template (see the module docstring)."""

    #: Fixtures used by the instance (plan 8.4). Set by the expansion.
    fixtures: tuple = ()
    #: Cache directory of the fixtures. None: <parent of the workdir>/fixtures.
    fixtures_root: str | None = None
    #: Number of MPI processes reported by ABINIT in the main output (None if not available).
    mpi_nprocs_run: int | None = None

    def setup_instance(self, *, instance_id: str, input_text: str, template_path: str, index: int | None,
                       params: dict[str, Any], suite_refs_dir: str, is_prerequisite: bool = False) -> None:
        """
        Called by the expansion right after BaseTest.__init__.
        Prerequisites (the members of setup.test_chain) have index None and no params.
        """
        self.id = instance_id
        self.input_text = input_text
        self.template_path = template_path
        self.instance_index = index
        self.params = dict(params)
        self.suite_refs_dir = suite_refs_dir
        self.is_prerequisite = is_prerequisite

    def run(self, build_env, runner, workdir, print_lock=None, mpi_nprocs=1, runmode="static", **kwargs):
        workdir = os.path.abspath(workdir)
        os.makedirs(workdir, exist_ok=True)

        # Write the generated input: prepare_new_cli_invokation reads and patches it in place.
        inp_path = os.path.join(workdir, self.id + ".abi")
        with open(inp_path, "w") as fh:
            fh.write(self.input_text)
        self.inp_fname = inp_path

        self.ref_dir = self.link_references(workdir)
        self.keep_files([inp_path])

        failures = self.prepare_fixtures(build_env, runner, workdir, mpi_nprocs, runmode)
        out = super().run(build_env, runner, workdir, print_lock=print_lock, mpi_nprocs=mpi_nprocs,
                          runmode=runmode, **kwargs)
        for msg in failures:
            self.exceptions.append(self.Error(msg))
            self.cprint(msg=f"{self.full_id}: {msg}", color="red")
        self.check_mpi_nprocs()
        return out

    def check_mpi_nprocs(self) -> None:
        """
        Cross-check with the output: ABINIT must report `mpi_nproc: N` with the N of the instance.
        Otherwise the instance did not test the MPI configuration it claims to test: it fails.
        """
        if self.executable != "abinit" or self._status in ("skipped", "disabled"):
            return
        self.mpi_nprocs_run = mpi_nprocs_from_output(os.path.join(self.workdir, self.id + ".abo"))
        if getattr(self, "is_prerequisite", False):
            # Chain prerequisites and fixtures: selected with the nproc of their chain, or of the fixture,
            # but launched on mpi_nprocs processes (1 for a chain, possibly fewer than nproc for a fixture).
            expected = self.mpi_nprocs
        else:
            expected = self.nprocs_to_test[0] if self.nprocs_to_test else 1
        if self.mpi_nprocs_run is not None and self.mpi_nprocs_run != expected:
            msg = (f"the main output reports mpi_nproc: {self.mpi_nprocs_run}, "
                   f"but the test must run on {expected} MPI processes")
            self.exceptions.append(self.Error(msg))
            self._status = "failed"
            self.cprint(msg=f"{self.full_id}: {msg}", color="red")

    def prepare_fixtures(self, build_env, runner, workdir, mpi_nprocs, runmode) -> list[str]:
        """
        Build (or reuse) the fixtures and link their outputs into workdir, but only if this pass runs the instance.
        The cache is <parent of workdir>/fixtures: Test_suite/fixtures for runtests.py, <workdir-root>/fixtures
        for run_params.py. Return the error messages of the fixtures that failed.
        """
        if not self.fixtures or self._status == "disabled":
            return []
        nprocs, _ = self.compute_nprocs(build_env, mpi_nprocs, runmode=runmode)
        if nprocs == 0:
            return []
        # Normally a cache hit: the fixtures are built before the tests (prebuild.py). Otherwise they are
        # built here, with the runner of the test and no CPU budget (possible oversubscription).
        root = self.fixtures_root or os.path.join(os.path.dirname(os.path.abspath(workdir)), "fixtures")
        failures = []
        for fx in self.fixtures:
            entry, st = fx.ensure(build_env, root, runner=runner, timebomb=getattr(runner, "timebomb", None))
            if st["status"] == "ok":
                fx.link_outputs(entry, workdir)
            elif st["status"] == "skipped":
                self.force_skip = True
                self.cprint(msg=f"{self.full_id}: fixture {fx.name} skipped: {st['message']}", color="cyan")
            else:
                failures.append(f"fixture {fx.name} failed: {st['message']} (cache entry: {entry})")
        return failures

    @property
    def inputs_used(self):
        """The template and the fixture files (for the stale-input check)."""
        inputs = super().inputs_used
        return inputs + [fx.path for fx in all_fixtures(list(self.fixtures)) if fx.path not in inputs]

    def link_references(self, workdir: str) -> str:
        """
        Link the golden reference(s) into <workdir>/refs under the names expected by fldiff.
        Return the path of the refs directory. Used by run() and by `run_params.py compare`.
        """
        refs_dir = os.path.join(os.path.abspath(workdir), "refs")
        os.makedirs(refs_dir, exist_ok=True)
        for f in self.files_to_test:
            link = os.path.join(refs_dir, f.name)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.abspath(os.path.join(self.suite_refs_dir, getattr(f, "ref_name", f.name))), link)
        return refs_dir

    def clean_workdir(self, other_test_files=None):
        """
        BaseTest.clean_workdir refuses directories in the workdir. Remove the directories created by
        the parametrized machinery (refs/ and the yaml_docs groups) if the test is OK and files are erased,
        otherwise keep them for inspection.
        """
        own = [os.path.join(self.workdir, "refs")]
        own += [os.path.join(self.workdir, f.name + ".yaml_docs") for f in self.files_to_test]
        own = [d for d in own if os.path.isdir(d)]
        erasing = (self.erase_files == 1 and self.isok) or self.erase_files == 2
        if erasing and self.isok:
            for d in own:
                shutil.rmtree(d, ignore_errors=True)
        keep = list(other_test_files or []) + [d for d in own if os.path.isdir(d)]
        return super().clean_workdir(other_test_files=keep)

    # Results passed from the worker processes of runtests.py -j to the parent, which writes the template pages.
    def get_results(self, *args, **kwargs):
        d = super().get_results(*args, **kwargs)
        d["parametrized"] = {
            "mpi_nprocs_run": self.mpi_nprocs_run,
            "exceptions": [str(e) for e in getattr(self, "exceptions", [])],
            "files": [{"fld_status": f.fld_status, "fld_msg": f.fld_msg, "fld_isok": f.fld_isok,
                       "yaml_groups": getattr(f, "yaml_groups", [])} for f in self.files_to_test],
        }
        return d

    def load_results(self, d):
        super().load_results(d)
        extra = d.get("parametrized")
        if not extra:
            return
        self.mpi_nprocs_run = extra["mpi_nprocs_run"]
        self.result_errors = extra["exceptions"]
        for f, fd in zip(self.files_to_test, extra["files"], strict=True):
            f.fld_status, f.fld_msg, f.fld_isok = fd["fld_status"], fd["fld_msg"], fd["fld_isok"]
            f.yaml_groups = fd["yaml_groups"]

    # HTML report of the instance (BaseTest.write_html_report): see html.py.
    #: Link from the page of the instance to the page of its template. None: the page written by runtests.py.
    template_page_href: str | None = None

    def html_extra_sections(self) -> str:
        from .htmlreport import instance_section
        return instance_section(self)

    def html_file_options(self, f) -> str:
        from .htmlreport import file_options
        return file_options(f)

    def diff_sources(self, f) -> tuple[str, str]:
        """yaml_docs: the diffs show the compared lines (yamltext.write_extracts), not the whole files."""
        if getattr(f, "param_mode", "fldiff") == "yaml_docs":
            d = os.path.join(self.workdir, f.name + ".yaml_docs")
            return os.path.join(d, EXTRACT_OUTPUT), os.path.join(d, EXTRACT_REFERENCE)
        return super().diff_sources(f)

    def html_diff_allowed(self, f) -> bool:
        return getattr(f, "param_mode", "fldiff") == "yaml_docs" or super().html_diff_allowed(f)

    def get_extra_inputs(self):
        """Copy the extra inputs in the workdir. The main input is already there (see run)."""
        for extra in self.extra_inputs:
            src = os.path.join(self.inp_dir, extra)
            dest = os.path.join(self.workdir, extra)
            if not os.path.isfile(src):
                self.exceptions.append(self.Error(f"{src}: no such file"))
                continue
            shutil.copy(src, dest)
            if dest.endswith(".gz"):
                unzip(dest)


class ParametrizedAbinitTest(ParametrizedMixin, AbinitTest):
    """Instance of a template whose executable is abinit."""


class ParametrizedTest(ParametrizedMixin, BaseTest):
    """Instance of a template for the other executables."""


class ParametrizedChain(ChainOfTests):
    """
    The prerequisites of a template followed by its instances with the same number of MPI processes
    (plan 8.2). All the members run in the same workdir: the instances read the files produced by
    the prerequisites. The prerequisites always run on 1 process ([NCPU_N] max_nprocs = 1).

    The id is short (tc_<NN>, plus _MPI<N>) because it names the workdir; the default id of a
    ChainOfTests joins the ids of all the members.
    """

    def __init__(self, tests, chain_id: str, nproc: int) -> None:
        super().__init__(tests)
        self._chain_id = chain_id
        self.nproc = nproc

    @property
    def id(self):
        return self._chain_id

    @property
    def prerequisites(self) -> list:
        return [t for t in self if getattr(t, "is_prerequisite", False)]

    @property
    def instances(self) -> list:
        return [t for t in self if not getattr(t, "is_prerequisite", False)]


class BrokenTemplateTest(BaseTest):
    """
    Placeholder for a template that could not be expanded.

    Expansion errors must not break the construction of the test database: the error is reported
    as a failed test with the validation message.
    """

    def set_error(self, message: str) -> None:
        self.error_message = message

    def run(self, build_env, runner, workdir, print_lock=None, mpi_nprocs=1, runmode="static", **kwargs):
        self.workdir = os.path.abspath(workdir)
        os.makedirs(self.workdir, exist_ok=True)
        with open(os.path.join(self.workdir, "template_error.txt"), "w") as fh:
            fh.write(self.error_message + "\n")
        self.keep_files([os.path.join(self.workdir, "template_error.txt")])
        self.exceptions = [self.Error(self.error_message)]
        self._status = "failed"
        self._isok = False
        self.run_etime = 0.0
        self.tot_etime = 0.0
        self._executed = True
        self.cprint(msg=f"{self.full_id}: [failed] cannot expand template:\n{self.error_message}", color="red")
