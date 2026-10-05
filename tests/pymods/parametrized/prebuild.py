"""
Build the fixtures before the tests (plan 8.4).

A fixture can run on several MPI processes (fixture.nproc). Built inside a test, it would use CPUs that the
worker did not reserve (runtests.py reserves the CPUs of the test only). So the fixtures of the tests that
will run are built here first, sequentially, within the CPU budget of the run, with the MPI runner of the run.
The tests then find them in the cache (fixture_cache.py).

Called by AbinitTestSuite.run_tests (runtests.py, testbot.py) and by run_params.py.
Only the standard library and the pickled modules are used.
"""
from __future__ import annotations

import os
from typing import Any, Callable

from ..termcolor import cprint
from .paramtests import Fixture, all_fixtures


def fixtures_of(tests) -> list[Fixture]:
    """The fixtures used by `tests` and what they need, dependencies first, one object per name."""
    seen: dict[str, Fixture] = {}
    for t in tests:
        for fx in all_fixtures(list(getattr(t, "fixtures", ()) or ())):
            seen.setdefault(fx.name, fx)
    return list(seen.values())


def prebuild_fixtures(tests, build_env, runner, root: str, max_cpus: int | None,
                      printer: Callable[..., Any] = cprint) -> dict[str, dict[str, Any]]:
    """
    Build (or find in the cache) the fixtures of `tests`, the tests that will actually run.
    Return {name: status}. Failures are reported here, once; the tests using a failed fixture
    report it again when they run.
    """
    fixtures = fixtures_of(tests)
    if not fixtures:
        return {}
    budget = max_cpus if max_cpus and max_cpus > 0 else (os.cpu_count() or 1)
    printer(f"Building {len(fixtures)} fixture(s) in {root} (CPU budget: {budget})", color="yellow")

    statuses = {}
    for fx in fixtures:
        # ensure() builds the needs first; they are already in the cache when their turn comes.
        _entry, st = fx.ensure(build_env, root, runner=runner, max_cpus=budget)
        statuses[fx.name] = st
        nproc = st.get("nproc", st.get("nproc_requested", fx.nproc))
        line = f"  fixture {fx.name}: {st['status']} (np={nproc}, {st.get('build_time', 0.0):.1f} s)"
        if st.get("message"):
            line += f": {st['message']}"
        color = {"ok": "green", "skipped": "cyan"}.get(st["status"], "red")
        if st["status"] == "ok" and nproc != fx.nproc:
            color = "yellow"      # fallback to fewer processes
        printer(line, color=color)
    return statuses


def prebuild_for_suite(suite, build_env, runner, workdir: str, mpi_nprocs: int, runmode: str,
                       max_cpus: int) -> dict[str, dict[str, Any]]:
    """
    Hook of AbinitTestSuite.run_tests: the fixtures of the tests that this pass runs
    (compute_nprocs > 0), in <workdir>/fixtures, the cache used by the tests (ParametrizedMixin).
    """
    selected = []
    for test in suite:
        if not getattr(test, "fixtures", None) or getattr(test, "_status", None) == "disabled":
            continue
        nprocs, _ = test.compute_nprocs(build_env, mpi_nprocs, runmode=runmode)
        if nprocs > 0:
            selected.append(test)
    return prebuild_fixtures(selected, build_env, runner, os.path.join(workdir, "fixtures"), max_cpus)


__all__ = ["fixtures_of", "prebuild_fixtures", "prebuild_for_suite"]
