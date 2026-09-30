"""
Run the instances of a template, or compare the outputs already present in their workdirs,
and summarize each instance as an outcome for the report (plan sections 5 and 6).

Execution and comparison go through the existing machinery: BaseTest.run() and FileToTest.compare().
This module does not import pydantic.
"""
from __future__ import annotations

import os
import shutil
from dataclasses import asdict, dataclass, field
from typing import Any

from ..jobrunner import JobRunner, TimeBomb
from .paramtests import BrokenTemplateTest, ParametrizedChain, mpi_nprocs_from_output

#: Statuses of an instance in the report.
STATUSES = ("ok", "differs", "error", "skipped", "pending")

#: The page of the template written by run_params.py in the workdir root, next to the workdirs of the instances.
TEMPLATE_PAGE = "index.html"


@dataclass
class FileOutcome:
    name: str
    mode: str
    #: fldiff status: succeeded, passed or failed.
    fld_status: str
    fld_msg: str
    #: Path of the file with the fldiff details (may not exist).
    fldiff_path: str = ""


@dataclass
class InstanceOutcome:
    test_id: str
    index: int | None
    params: dict[str, Any]
    nproc: int
    workdir: str
    status: str
    message: str = ""
    files: list[FileOutcome] = field(default_factory=list)
    #: "instance", or "prerequisite" for a member of setup.test_chain (plan 8.2).
    role: str = "instance"
    #: Number of MPI processes reported by ABINIT in the main output (`mpi_nproc: N`), None if not available.
    nproc_run: int | None = None

    @property
    def input_file(self) -> str:
        return self.test_id + ".abi"

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["input_file"] = self.input_file
        return d


def nproc_of(test) -> int:
    return test.nprocs_to_test[0] if getattr(test, "nprocs_to_test", None) else 1


def _file_mode(f) -> str:
    return getattr(f, "param_mode", "fldiff")


def _base(test, workdir: str, nproc: int | None = None) -> dict[str, Any]:
    return dict(test_id=test.id, index=getattr(test, "instance_index", None), params=dict(getattr(test, "params", {})),
                nproc=nproc_of(test) if nproc is None else nproc, workdir=workdir,
                role="prerequisite" if getattr(test, "is_prerequisite", False) else "instance")


def units(tests):
    """
    Yield (unit, members, workdir name, nproc). A unit is what runs in one workdir:
    a single test, or a ParametrizedChain (its prerequisites and instances).
    """
    for test in tests:
        if isinstance(test, ParametrizedChain):
            yield test, list(test), test.id, test.nproc
        else:
            yield test, [test], test.id, nproc_of(test)


def _outcome_from_test(test, workdir: str, compared: bool, nproc: int | None = None) -> InstanceOutcome:
    """Classify a test after BaseTest.run() or after compare_instances()."""
    base = _base(test, workdir, nproc)
    if getattr(test, "executable", "") == "abinit":
        base["nproc_run"] = mpi_nprocs_from_output(os.path.join(workdir, test.id + ".abo"))

    if isinstance(test, BrokenTemplateTest):
        return InstanceOutcome(status="error", message=test.error_message, **base)

    files = [FileOutcome(f.name, _file_mode(f), f.fld_status, f.fld_msg, os.path.join(workdir, f.name + ".fldiff"))
             for f in test.files_to_test] if compared else []

    errors = [str(e) for e in getattr(test, "exceptions", [])]
    if getattr(test, "_status", None) == "skipped":
        return InstanceOutcome(status="skipped", message=getattr(test, "skip_msg", "") or "skipped", **base)
    if errors or getattr(test, "exec_error", False):
        return InstanceOutcome(status="error", message="; ".join(errors) or "execution failed", files=files, **base)

    missing = [f.name for f in test.files_to_test if not os.path.exists(os.path.join(workdir, f.name))]
    if missing:
        return InstanceOutcome(status="error", message=f"output not found: {', '.join(missing)}", files=files, **base)

    expected = 1 if base["role"] == "prerequisite" else base["nproc"]
    if base.get("nproc_run") is not None and base["nproc_run"] != expected:
        return InstanceOutcome(status="error", files=files, **base,
                               message=f"the main output reports mpi_nproc: {base['nproc_run']}, expected {expected}")

    if not compared or not test.files_to_test:
        return InstanceOutcome(status="ok", message="not compared (no reference)", files=files, **base)

    if any(f.fld_status == "failed" for f in files):
        return InstanceOutcome(status="differs", message="; ".join(f.fld_msg for f in files if f.fld_status == "failed"),
                               files=files, **base)
    return InstanceOutcome(status="ok", files=files, **base)


def run_instances(tests, build_env, workdir_root: str, max_cpus: int, timeout: float = 900.0) -> list[InstanceOutcome]:
    """
    Run the tests sequentially, each in <workdir_root>/<id>, and return the outcomes of the instances
    and of the prerequisites. A chain runs in dynamic mode, so that its prerequisites
    ([NCPU_N] max_nprocs = 1) run on 1 process and its instances on N.
    """
    timebomb = TimeBomb(timeout)
    outcomes = []
    for unit, members, name, nproc in units(tests):
        workdir = os.path.join(os.path.abspath(workdir_root), name)
        if nproc > max_cpus and not isinstance(unit, BrokenTemplateTest):
            for m in members:
                m._status, m.skip_msg = "skipped", f"nproc={nproc} > max_cpus={max_cpus}"
                outcomes.append(InstanceOutcome(status="skipped", message=m.skip_msg, **_base(m, workdir, nproc)))
            continue
        runner = JobRunner.generic_mpi(timebomb=timebomb) if nproc > 1 else JobRunner.sequential(timebomb=timebomb)
        runmode = "dynamic" if isinstance(unit, ParametrizedChain) else "static"
        # Start from an empty workdir: files left by a previous run can change the results
        # (e.g. EPH restarts from an existing SIGEPH.nc). The fixture cache is not inside it.
        shutil.rmtree(workdir, ignore_errors=True)
        os.makedirs(workdir, exist_ok=True)   # ChainOfTests.run creates only the last component
        for m in members:
            m.template_page_href = "../" + TEMPLATE_PAGE
        unit.run(build_env, runner, workdir, mpi_nprocs=nproc, runmode=runmode, erase_files=0, make_html_diff=1)
        if not isinstance(unit, BrokenTemplateTest):     # its error is shown by the template page
            unit.write_html_report()
        for m in members:
            outcomes.append(_outcome_from_test(m, workdir, compared=not getattr(m, "no_check", False), nproc=nproc))
    return outcomes


def compare_instances(tests, workdir_root: str, fldiff_path: str = "",
                      max_cpus: int | None = None) -> list[InstanceOutcome]:
    """
    Compare the outputs already present in <workdir_root>/<id> with the golden reference(s), without running
    anything. Works the same way whether the outputs come from `run_params.py run` or from a batch job.
    """
    outcomes = []
    for test, workdir_name, nproc in _members(tests):
        workdir = os.path.join(os.path.abspath(workdir_root), workdir_name)
        base = _base(test, workdir, nproc)

        if isinstance(test, BrokenTemplateTest):
            outcomes.append(InstanceOutcome(status="error", message=test.error_message, **base))
            continue
        if max_cpus is not None and base["nproc"] > max_cpus:
            outcomes.append(InstanceOutcome(status="skipped", message=f"nproc={base['nproc']} > max_cpus={max_cpus}", **base))
            continue
        if not os.path.isdir(workdir):
            outcomes.append(InstanceOutcome(status="error", message=f"workdir not found: {workdir}", **base))
            continue

        missing = [f.name for f in test.files_to_test if not os.path.exists(os.path.join(workdir, f.name))]
        if missing:
            outcomes.append(InstanceOutcome(status="error", message=f"output not found: {', '.join(missing)}", **base))
            continue

        refs_dir = test.link_references(workdir)
        for f in test.files_to_test:
            with open(os.path.join(workdir, f.name + ".fldiff"), "w") as fh:
                f.compare(fldiff_path, refs_dir, workdir, yaml_test={}, outf=fh)
        # Only the comparison is evaluated: reset what BaseTest.run() would have set.
        test.exceptions = []
        test.exec_error = False
        outcomes.append(_outcome_from_test(test, workdir, compared=True, nproc=nproc))
    return outcomes


def _members(tests):
    """(test, workdir name, nproc) of the tests to compare. Prerequisites without files to test are omitted."""
    for unit, members, name, nproc in units(tests):
        for m in members:
            if getattr(m, "is_prerequisite", False) and not m.files_to_test:
                continue
            yield m, name, nproc


def is_failure(outcomes: list[InstanceOutcome]) -> bool:
    return any(o.status in ("differs", "error") for o in outcomes)


def runtests_status(outcome: InstanceOutcome) -> str:
    """
    The runtests.py status of an instance: failed (differs or error), passed (ok, some file differs within
    the tolerances), succeeded (ok, no difference), skipped.
    """
    if outcome.status in ("differs", "error"):
        return "failed"
    if outcome.status == "ok":
        return "passed" if any(f.fld_status == "passed" for f in outcome.files) else "succeeded"
    return outcome.status


def write_html_pages(template: str, tests, outcomes: list[InstanceOutcome], workdir_root: str,
                     command: str = "", timeout: float = 900.0) -> str:
    """
    Write the pages of runtests.py for run_params.py: the page of each instance (<workdir>/index.html,
    BaseTest.write_html_report) and the page of the template (<workdir_root>/index.html, htmlreport.py).
    The instances only compared (`compare`) get a page too, rewritten with the new result.
    Return the path of the template page.
    """
    from .htmlreport import write_template_page
    by_id = {o.test_id: o for o in outcomes}
    timebomb = TimeBomb(timeout)
    for test, workdir_name, nproc in _members(tests):
        o = by_id.get(test.id)
        if o is None or getattr(test, "_executed", False) or isinstance(test, BrokenTemplateTest):
            continue      # run: the page was written by run_instances
        test._status = runtests_status(o)
        if o.status == "skipped" or not os.path.isdir(o.workdir):
            test.skip_msg = o.message
            continue
        _prepare_compared(test, o, nproc, timebomb)
        try:
            test.write_html_report()
        except Exception as exc:        # the page is a convenience: never fail the comparison because of it
            print(f"Cannot write the page of {test.id}: {type(exc).__name__}: {exc}")
    root = os.path.abspath(workdir_root)
    return write_template_page(os.path.join(root, TEMPLATE_PAGE), template, tests,
                               fixtures_root=os.path.join(root, "fixtures"), command=command)


def _prepare_compared(test, o: InstanceOutcome, nproc: int, timebomb) -> None:
    """Set what BaseTest.write_html_report reads for a test compared by compare_instances but not run."""
    wd = o.workdir
    test.workdir, test.ref_dir, test.timebomb = wd, os.path.join(wd, "refs"), timebomb
    test.inp_fname = os.path.join(wd, test.id + ".abi")
    test.stdout_fname, test.stderr_fname = (os.path.join(wd, test.id + ext) for ext in (".stdout", ".stderr"))
    test.run_etime = test.tot_etime = 0.0
    test.mpi_nprocs, test.mpi_nprocs_run = nproc, o.nproc_run
    test.exceptions = [test.Error(o.message)] if o.status == "error" and o.message else []
    test.make_html_diff, test._executed, test._isok = 1, True, o.status == "ok"
    test.fld_isok = all(f.fld_isok for f in test.files_to_test)
    test.template_page_href = "../" + TEMPLATE_PAGE
    for f in test.files_to_test:
        f.fldiff_fname = os.path.join(wd, f.name + ".fldiff")


__all__ = ["STATUSES", "TEMPLATE_PAGE", "FileOutcome", "InstanceOutcome", "compare_instances", "is_failure",
           "run_instances", "runtests_status", "units", "write_html_pages"]
