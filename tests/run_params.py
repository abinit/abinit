#!/usr/bin/env python
"""
Command line interface for the parametrized ABINIT tests.

    run_params.py doc
    run_params.py fixtures [INPUT_DIR] [--doc]
    run_params.py expand  TEMPLATE.tmpl.abi [-o OUT]
    run_params.py run     TEMPLATE.tmpl.abi [--build-dir BUILD] [--max-cpus N] [--reference REF.abo]
                          [--workdir-root W] [--update-ref] [-o comparison_report.json]
    run_params.py compare TEMPLATE.tmpl.abi [--reference REF.abo] [--workdir-root W] [--max-cpus N]
                          [-o comparison_report.json] [-p STATUS]

`run` and `compare` accept `-p STATUS` with the syntax and statuses of runtests.py -p (failed, passed, succeeded, all,
not_<status>): open the output and the golden reference of the files of those instances in $PATCHER (default: vimdiff).
`doc` prints the format of the %YAML block, generated from the pydantic models.
`fixtures` lists the fixtures with their documentation and the templates that use them (`--doc`: the file format).
`run` executes the instances and compares them with the golden reference.
`compare` only compares the outputs already present in the workdirs (e.g. produced by `run` or by a batch job).
Both write the HTML pages of runtests.py: <workdir-root>/index.html for the template (all the instances,
their parameters, statuses and differences) and <workdir-root>/<id>/index.html for each instance,
plus comparison_report.json, and exit with a non-zero status if an instance differs from the reference or fails.

The same expansion (tests.pymods.parametrized.expand.expand_template) is used by
`runtests.py --with-parametrized`: the instances are ordinary test objects run through BaseTest.run().
See tests/plans/parametrized_tests_plan.md.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

# Same path set-up as runtests.py, so that `tests` resolves to this checkout.
pack_dir, _ = os.path.split(os.path.abspath(__file__))
pack_dir, _ = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)

try:
    import pydantic  # noqa: F401
    import yaml  # noqa: F401
except ImportError as exc:
    raise SystemExit(f"run_params.py requires pyyaml and pydantic ({exc}). Install them with: pip install pyyaml pydantic")

import tests  # noqa: E402
from tests.pymods.parametrized.execute import (compare_instances, is_failure, run_instances,  # noqa: E402
                                               runtests_status, write_html_pages)
from tests.pymods.parametrized import fixture_cache  # noqa: E402
from tests.pymods.parametrized.expand import (build_instances, expand_template, load_fixtures,  # noqa: E402
                                              load_prerequisites)
from tests.pymods.jobrunner import JobRunner  # noqa: E402
from tests.pymods.parametrized.paramtests import BrokenTemplateTest, ParametrizedChain  # noqa: E402
from tests.pymods.parametrized.prebuild import prebuild_fixtures  # noqa: E402
from tests.pymods.parametrized.report import build_report, write_reports  # noqa: E402
from tests.pymods.parametrized.template import TemplateError, load_template  # noqa: E402
from tests.pymods.testsuite import BuildEnvironment, find_top_build_tree  # noqa: E402

abenv = tests.abenv

STATUS_ORDER = ("ok", "differs", "error", "skipped", "pending")


def reference_paths(tests) -> list[str]:
    """The golden references of the files compared by the instances (the same for every instance)."""
    for t in tests:
        members = t.instances if isinstance(t, ParametrizedChain) else [t]
        for m in members:
            if not isinstance(m, BrokenTemplateTest):
                return [os.path.join(m.suite_refs_dir, f.ref_name) for f in m.files_to_test]
    return []


def report_path(args) -> str:
    return os.path.abspath(args.output or os.path.join(args.workdir_root, "comparison_report.json"))


def print_outcomes(outcomes) -> None:
    print()
    for o in outcomes:
        what = "prerequisite (1 process)" if o.role == "prerequisite" else o.params
        ran = f"ran {o.nproc_run}" if o.nproc_run is not None else "ran ?"
        print(f"{o.test_id:20s} {o.status:8s} np={o.nproc:<2d} {ran:7s} {what}")
        if o.status in ("differs", "error"):
            print(f"{'':20s} {o.message}")
    counts = {s: sum(1 for o in outcomes if o.status == s) for s in STATUS_ORDER}
    print("\nSummary: " + ", ".join(f"{k}={v}" for k, v in counts.items() if v or k in ("ok", "differs", "error")))


def finish(args, command: str, reference: str | None, outcomes, tests=None) -> int:
    print_outcomes(outcomes)
    report = build_report(args.template, reference, args.workdir_root, command, outcomes)
    json_path = write_reports(report, report_path(args))
    if tests is not None:
        html_path = write_html_pages(os.path.abspath(args.template), tests, outcomes, args.workdir_root,
                                     command=" ".join(["run_params.py"] + sys.argv[1:]))
        print(f"Report: {html_path}\n        {json_path}")
    else:
        print(f"Report: {json_path}")
    if getattr(args, "patch", "") and tests is not None:
        if reference is None:
            print("No reference: nothing to patch (use --update-ref to create it).")
        else:
            patch_references(tests, outcomes, parse_patch_statuses(args.patch))
    return 1 if is_failure(outcomes) else 0


#: Statuses accepted by -p/--patch: the test statuses of runtests.py -p.
PATCH_STATUSES = ("failed", "passed", "succeeded", "skipped", "disabled")


def parse_patch_statuses(spec: str) -> set[str]:
    """Same syntax as runtests.py -p: `failed`, `failed+passed`, `all` (failed+passed+succeeded), `not_succeeded`."""
    if spec == "all":
        return {"failed", "passed", "succeeded"}
    if spec.startswith("not_"):
        excluded = spec[4:]
        if excluded not in PATCH_STATUSES:
            raise SystemExit(f"-p {spec}: {excluded} is not one of {PATCH_STATUSES}")
        return set(PATCH_STATUSES) - {excluded}
    statuses = {t.strip() for t in spec.split("+")}
    bad = statuses - set(PATCH_STATUSES)
    if bad:
        raise SystemExit(f"-p {spec}: invalid status {sorted(bad)}, expected {PATCH_STATUSES}, all or not_<status>")
    return statuses


#: The runtests.py status of an outcome (-p selects the instances by these statuses).
test_status = runtests_status


def patch_pairs(tests, outcomes, statuses: set[str]) -> list[tuple[str, str]]:
    """
    (output, golden reference) of all the compared files of the instances whose runtests.py status is in
    `statuses`, as runtests.py -p (Results.outref_files). Missing outputs are skipped.
    """
    refs = {}
    for t in tests:
        for m in (t if isinstance(t, ParametrizedChain) else [t]):
            for f in getattr(m, "files_to_test", ()):
                # ref_name is relative to the Refs directory, or absolute with --reference.
                refs[(m.id, f.name)] = os.path.abspath(os.path.join(m.suite_refs_dir, getattr(f, "ref_name", f.name)))
    pairs = []
    for o in outcomes:
        if test_status(o) not in statuses:
            continue
        for f in o.files:
            out, ref = os.path.join(o.workdir, f.name), refs.get((o.test_id, f.name))
            if ref is not None and os.path.exists(out):
                pairs.append((out, ref))
    return pairs


def patch_references(tests, outcomes, statuses: set[str]) -> int:
    """Open each (output, reference) pair in $PATCHER (default vimdiff), as runtests.py -p does."""
    from tests.pymods.tools import Patcher
    pairs = patch_pairs(tests, outcomes, statuses)
    if not pairs:
        print(f"Nothing to patch for status {'+'.join(sorted(statuses))}.")
        return 0
    print(f"Opening {len(pairs)} file(s): output on the left, golden reference on the right.")
    print("Note: all the instances share the golden reference, so patching it changes the comparison of every instance.")
    return Patcher().patch_files([o for o, _ in pairs], [r for _, r in pairs])


def cmd_doc(args) -> int:
    from tests.pymods.parametrized.docs import doc_template_format
    doc_template_format(sys.stdout)
    return 0


def cmd_fixtures(args) -> int:
    from tests.pymods.parametrized.docs import describe_fixtures, doc_fixture_format
    input_dir = os.path.abspath(args.input_dir or os.path.join(abenv.tests_dir, "parametrized", "Input"))
    nbad = describe_fixtures(input_dir, sys.stdout, known_keywords=tests.KNOWN_KEYWORDS)
    if args.doc:
        print()
        doc_fixture_format(sys.stdout)
    return 1 if nbad else 0


def cmd_expand(args) -> int:
    try:
        template = load_template(args.template, known_keywords=tests.KNOWN_KEYWORDS)
        prereqs = load_prerequisites(template, known_keywords=tests.KNOWN_KEYWORDS)
        fixtures = load_fixtures(template, known_keywords=tests.KNOWN_KEYWORDS)
        instances = build_instances(template)
    except Exception as exc:
        print(f"Cannot expand {args.template}:\n{exc}", file=sys.stderr)
        return 1

    print(f"{args.template}: {len(instances)} instance(s)")
    if prereqs:
        names = ", ".join(os.path.basename(p.path) for p in prereqs)
        nprocs = sorted({inst.nproc for inst in instances})
        print(f"  prerequisites (setup.test_chain, on 1 process): {names}")
        print(f"  one chain per nproc {nprocs}: the prerequisites run once in each chain")
    if fixtures:
        print(f"  fixtures (built once per workdir root, before the tests): {', '.join(template.spec.fixtures)}")
        for name, fx in fixtures.items():
            needs = f", needs {', '.join(fx.spec.fixture.needs)}" if fx.spec.fixture.needs else ""
            print(f"    {name}: {fx.spec.test_info.setup['executable']} (np={fx.spec.fixture.nproc}) "
                  f"-> {', '.join(fx.spec.fixture.outputs)}{needs}")
    for inst in instances:
        print(f"  {inst.test_id:20s} nproc={inst.nproc:<3d} {inst.params}")

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        for inst in instances:
            with open(os.path.join(args.output_dir, inst.test_id + ".abi"), "w") as fh:
                fh.write(inst.input_text)
        print(f"Inputs written to {args.output_dir} (for inspection only)")
    return 0


def find_instance(tests, index: int):
    """Return (unit to run, test) for the instance `index`: a chain keeps only its prerequisites and that instance."""
    for t in tests:
        if isinstance(t, ParametrizedChain):
            for inst in t.instances:
                if inst.instance_index == index:
                    return ParametrizedChain(t.prerequisites + [inst], chain_id=t.id, nproc=t.nproc), inst
        elif getattr(t, "instance_index", None) == index:
            return t, t
    raise ValueError(f"instance {index} not found")


def set_fixtures_root(tests, root: str) -> None:
    """All the runs of this command share the fixture cache <workdir-root>/fixtures (also the --update-ref run)."""
    for t in tests:
        for m in (t if isinstance(t, ParametrizedChain) else [t]):
            m.fixtures_root = root


def prebuild(tests, build_env, root: str, max_cpus: int) -> None:
    """Build the fixtures of the instances that will run (nproc <= max_cpus) before running them (plan 8.4)."""
    selected = [t for t in tests if not isinstance(t, (ParametrizedChain, BrokenTemplateTest))
                and (t.nprocs_to_test[0] if t.nprocs_to_test else 1) <= max_cpus]
    prebuild_fixtures(selected, build_env, JobRunner.generic_mpi(), root, max_cpus)


def cmd_run(args) -> int:
    build_env = BuildEnvironment(find_top_build_tree(args.build_dir))
    os.environ.setdefault("ABI_PSPDIR", abenv.psps_dir)

    template_path = os.path.abspath(args.template)
    reference = os.path.abspath(args.reference) if args.reference else None
    workdir_root = os.path.abspath(args.workdir_root)
    os.makedirs(workdir_root, exist_ok=True)
    fixtures_root = os.path.join(workdir_root, "fixtures")
    if args.rebuild_fixtures:
        fixture_cache.remove_all(fixtures_root)

    if args.update_ref:
        # The golden references are the files produced by instance 0000 (plan P9):
        # the main output and every other file listed in files_to_test.
        full = expand_template(template_path, abenv, reference=reference)
        first = expand_template(template_path, abenv, no_check=True)
        broken = [t for t in full + first if isinstance(t, BrokenTemplateTest)]
        if broken:
            print(broken[0].error_message, file=sys.stderr)
            return 1
        _, full_test = find_instance(full, 0)
        set_fixtures_root(first, fixtures_root)
        prebuild([t for t in first if getattr(t, "instance_index", None) == 0], build_env, fixtures_root, args.max_cpus)
        unit, _ = find_instance(first, 0)
        outcomes = run_instances([unit], build_env, os.path.join(workdir_root, "update_ref"), args.max_cpus, args.timeout)
        bad = [o for o in outcomes if o.status != "ok"]
        if bad:
            for o in bad:
                print(f"Cannot create the references from instance 0 ({o.test_id}): {o.status} {o.message}",
                      file=sys.stderr)
            return 1
        outcome = next(o for o in outcomes if o.index == 0)
        refs_dir = full_test.suite_refs_dir
        for f in full_test.files_to_test:
            src = os.path.join(outcome.workdir, f.name)
            dest = os.path.join(refs_dir, f.ref_name)   # ref_name is absolute for --reference
            if not os.path.exists(src):
                print(f"Cannot create the reference {dest}: instance {outcome.test_id} did not produce {f.name}",
                      file=sys.stderr)
                return 1
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copy(src, dest)
            print(f"Reference written: {dest} (from {f.name} of instance {outcome.test_id})")

    refs = reference_paths(expand_template(template_path, abenv, reference=reference))
    missing = [p for p in refs if not os.path.exists(p)]
    no_check = bool(missing)
    if no_check:
        print(f"Missing reference(s) {', '.join(missing)}: the instances are run without comparison "
              "(use --update-ref to create them).")

    test_list = expand_template(template_path, abenv, reference=reference, no_check=no_check)
    set_fixtures_root(test_list, fixtures_root)
    prebuild(test_list, build_env, fixtures_root, args.max_cpus)
    outcomes = run_instances(test_list, build_env, workdir_root, args.max_cpus, args.timeout)

    if no_check:
        for o in outcomes:
            out = os.path.join(o.workdir, o.test_id + ".abo")
            if os.path.exists(out):
                print(f"\n==== YAML documents of {o.test_id}")
                with open(out) as fh:
                    inside = False
                    for line in fh:
                        inside = inside or line.startswith("--- !")
                        if inside:
                            print(line, end="")
                        if line.rstrip() == "...":
                            inside = False

    return finish(args, "run", None if no_check else ", ".join(refs), outcomes, tests=test_list)


def cmd_compare(args) -> int:
    template_path = os.path.abspath(args.template)
    reference = os.path.abspath(args.reference) if args.reference else None
    test_list = expand_template(template_path, abenv, reference=reference)
    refs = reference_paths(test_list)
    missing = [p for p in refs if not os.path.exists(p)]
    if missing:
        print(f"Reference(s) not found: {', '.join(missing)}", file=sys.stderr)
        return 1

    outcomes = compare_instances(test_list, args.workdir_root, fldiff_path=abenv.fldiff_path, max_cpus=args.max_cpus)
    return finish(args, "compare", ", ".join(refs), outcomes, tests=test_list)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("doc", help="Print the format of the %%YAML block of a template (from the pydantic models).")
    p.set_defaults(func=cmd_doc)

    p = sub.add_parser("fixtures", help="List the fixtures (Input/*.fixture.abi) with their documentation and users.")
    p.add_argument("input_dir", nargs="?", default=None,
                   help="Directory with the fixtures (default: tests/parametrized/Input).")
    p.add_argument("--doc", action="store_true", help="Also print the format of the fixture files (from the pydantic models).")
    p.set_defaults(func=cmd_fixtures)

    p = sub.add_parser("expand", help="Validate a template and list (or write) its instances.")
    p.add_argument("template")
    p.add_argument("-o", "--output-dir", default=None, help="Write the instance inputs here (inspection only).")
    p.set_defaults(func=cmd_expand)

    def add_common(p):
        p.add_argument("template")
        p.add_argument("--reference", default=None, help="Golden reference of the main output {id}.abo (default: Refs/tNN.abo of the suite).")
        p.add_argument("--workdir-root", default="run_params_workdir", help="Root of the workdirs (default: %(default)s).")
        p.add_argument("-o", "--output", default=None,
                       help="Path of the JSON report; the HTML report is written next to it "
                            "(default: <workdir-root>/comparison_report.json).")
        p.add_argument("-p", "--patch", default="", metavar="STATUS",
                       help="Patch the golden references of the instances with the status specified by -p, as runtests.py -p: "
                            "the output and the reference of every compared file are opened in $PATCHER (default: vimdiff). "
                            "Examples: `-p failed` (instances that differ or are in error), `-p all`, "
                            "`-p failed+passed` or, equivalently, `-p not_succeeded`. "
                            "All the instances share the reference: patching it changes the comparison of every instance.")

    p = sub.add_parser("run", help="Run the instances of a template and compare them with the golden reference.")
    add_common(p)
    p.add_argument("--build-dir", default=os.curdir, help="Directory inside the ABINIT build tree (default: cwd).")
    p.add_argument("--max-cpus", type=int, default=1, help="Instances with more MPI processes are skipped.")
    p.add_argument("--update-ref", action="store_true", help="Create the reference from instance 0000, then run.")
    p.add_argument("--timeout", type=float, default=900.0, help="Timeout per instance in seconds.")
    p.add_argument("--rebuild-fixtures", action="store_true",
                   help="Remove the fixture cache (<workdir-root>/fixtures) before running.")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("compare", help="Compare the outputs already present in the workdirs (no execution).")
    add_common(p)
    p.add_argument("--max-cpus", type=int, default=None,
                   help="Report instances with more MPI processes as skipped instead of missing.")
    p.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    try:
        return args.func(args)
    except TemplateError as exc:
        print(exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
