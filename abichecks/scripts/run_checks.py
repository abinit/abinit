#!/usr/bin/env python3
"""Pure-Python driver for the abichecks maintenance-check series.

Replaces the legacy Perl chain (run-standard-tests.pl, fldiff.pl,
reportdiff.pl, wrap-standard-tests.sh, Sort.sh). Reads a series' tests.cnf
dispatch table (format unchanged), runs each registered check, and writes a
`report` file plus per-case log files in the same on-disk shape the legacy
runner produced -- this is what abibuildbot's special_processor.py parses,
and that contract is preserved byte-for-byte on purpose (see abichecks/README.md).

Usage (mirrors the legacy positional CLI so the Makefile recipe just swaps
the interpreter/script):

    run_checks.py <machine_name> <series> [<start> [<stop>]]

Only six dispatch verbs are supported -- the only ones any tests.cnf file
actually registers today: warnchk, check_forbidden, check_ascii,
check_inlined_macros, statchk, report. Everything else in the legacy Perl
driver existed to run the unrelated main Fortran numerical test suite.
"""

from __future__ import annotations

import difflib
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
ABICHECKS_DIR = SCRIPTS_DIR.parent

# Verbs dispatched to a same-named checker script under scripts/, diffed
# against Refs/t{id}.out with no stdout filtering.
GOLDEN_FILE_SCRIPTS = {
    "check_forbidden": "check_forbidden.py",
    "check_ascii": "check_ascii.py",
    "check_inlined_macros": "check_inlined_macros.py",
}


def parse_tests_cnf(path: Path) -> list[tuple[str, str, list[str]]]:
    """Parse a tests.cnf dispatch table, returning (test_id, verb, params) entries.

    Mirrors run-standard-tests.pl's two-pass parsing exactly: the first
    non-'#' line in the file is a title, consumed and never treated as a
    case (legacy reads it in a separate loop before reopening nothing --
    the same filehandle position carries over). Every following line is
    dropped as a comment if it starts with '#', and dropped as malformed
    (matching legacy's tolerant `next if (...)`) if it has no test_id, the
    test_id is longer than 10 characters, or there is no verb.
    """
    lines = path.read_text().splitlines()

    idx = 0
    for idx, line in enumerate(lines):
        if not line.startswith("#"):
            idx += 1  # consume the title line itself too
            break
    else:
        idx = len(lines)

    entries: list[tuple[str, str, list[str]]] = []
    for line in lines[idx:]:
        if line.startswith("#"):
            continue
        parts = line.split()
        if not parts:
            continue
        test_id = parts[0]
        verb = parts[1] if len(parts) > 1 else ""
        if not test_id or len(test_id) > 10 or not verb:
            continue
        entries.append((test_id, verb, parts[2:]))
    return entries


def select_range(
    entries: list[tuple[str, str, list[str]]], first: str, last: str
) -> list[tuple[str, str, list[str]]]:
    """Select the inclusive [first, last] range of cases by literal test_id equality.

    Reproduces run-standard-tests.pl's range-matching state machine exactly,
    including the "3" != "03" string-equality behavior (matches Make's
    start=/stop= variables, which are plain text).
    """
    selected: list[tuple[str, str, list[str]]] = []
    cur_test = ""
    for test_id, verb, params in entries:
        if cur_test == "":
            if test_id == first or first == "first":
                cur_test = test_id
            else:
                continue
        elif test_id != cur_test:
            if cur_test == last or (test_id == "end" and last == "end"):
                break
            cur_test = test_id
        selected.append((test_id, verb, params))
    return selected


def parse_report_in(path: Path) -> dict[str, tuple[int, float, float]]:
    """Parse Input/report.in, returning {case_id: (tolnlines, tolabs, tolrel)}.

    Format: "Case_NN  tolnlines=N  tolabs=X  tolrel=X", '#'-prefixed lines
    disabled/ignored. Missing file or missing entry both mean "no recorded
    tolerance" -- callers default to (0, 0.0, 0.0), i.e. exact match.
    """
    tolerances: dict[str, tuple[int, float, float]] = {}
    if not path.is_file():
        return tolerances
    pattern = re.compile(
        r"Case_(\d+)\s+tolnlines=\s*(\d+)\s+tolabs=\s*([0-9.eE+-]+)\s+tolrel=\s*([0-9.eE+-]+)"
    )
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        m = pattern.match(line)
        if m:
            case_id, tolnlines, tolabs, tolrel = m.groups()
            tolerances[case_id] = (int(tolnlines), float(tolabs), float(tolrel))
    return tolerances


def sort_warnchk_body(lines: list[str]) -> list[str]:
    """Reimplement Sort.sh: normalize non-deterministic ordering for warnchk 3/4.

    warningschk.py's fixed output shape is a 3-line header, a variable-length
    body of one entry per warning, and a 3-line footer. Variable-name
    enumeration order in the compiler's own warning output isn't
    deterministic across runs/compiler versions for warning categories 3/4
    ("Unused variable"/"Unused dummy argument"), so the body is sorted while
    the header/footer are kept verbatim -- matching Sort.sh's
    `sed 1,3 / sed 4,/^*/ | head -n-1 | sort / tail -3` pipeline exactly for
    this fixed-shape input.
    """
    if len(lines) < 6:
        return lines
    return lines[:3] + sorted(lines[3:-3]) + lines[-3:]


def classify_golden(
    produced: Path, reference: Path, diff_path: Path, case_id: str, tolerances: dict[str, tuple[int, float, float]]
) -> tuple[str, str]:
    """Diff `produced` against `reference` and classify succeeded/passed/failed.

    Reuses tests/pymods/fldiff.py's Differ/Result -- a tested, existing
    Python port of fldiff.pl/reportdiff.pl's tolerance-aware classification
    -- instead of a new comparison engine. `Input/report.in`'s tolerances
    are 0 for every case today, so in practice this reduces to "exact text
    match, else failed", but reusing the real API stays correct if that
    ever changes.
    """
    from pymods.fldiff import Differ  # local import: needs sys.path set up by caller first

    ref_lines = reference.read_text().splitlines(keepends=True) if reference.is_file() else []
    prod_lines = produced.read_text().splitlines(keepends=True)
    diff_path.write_text(
        "".join(difflib.unified_diff(ref_lines, prod_lines, fromfile=str(reference), tofile=str(produced)))
    )

    tolnlines, tolabs, tolrel = tolerances.get(case_id, (0, 0.0, 0.0))
    result = Differ(ignore=True, ignoreP=True).diff(str(produced), str(reference))
    _is_ok, status, msg = result.passed_within_tols(tolnlines, tolabs, tolrel)
    return status, msg


def _stderr_note_lines(proc: subprocess.CompletedProcess, script_name: str) -> list[str]:
    """Extra lines to append when a checker subprocess wrote to stderr.

    A checker script's normal, successful run never touches stderr -- only a
    real error does (an uncaught exception, a missing-prerequisite
    RuntimeError like warningschk.py's "Cannot find make.log or
    make.stderr"). Without this, such a failure showed up only as an
    unexplained diff mismatch against Refs/ -- the actual reason, printed to
    stderr, was captured by subprocess.run() and then silently discarded.
    Appending it here surfaces the real reason in t{id}.out/.log, which
    special_processor.py aggregates into the Buildbot step's visible
    out.log/log.log for failed cases.
    """
    if not proc.stderr:
        return []
    return ["", f"--- {script_name} stderr ---", *proc.stderr.splitlines()]


def run_warnchk(
    case_id: str,
    warno: str,
    work_dir: Path,
    refs_dir: Path,
    tolerances: dict[str, tuple[int, float, float]],
    abinit_builddir: Path,
) -> tuple[str, str]:
    """Run warningschk.py, apply the legacy SUCCESS-filtering/sort normalization, classify."""
    script = SCRIPTS_DIR / "warningschk.py"
    log_path = work_dir / f"t{case_id}.log"
    out_path = work_dir / f"t{case_id}.out"
    diff_path = work_dir / f"diff.t{case_id}"

    # No home_dir argument (2nd, left blank), matching legacy dochkwarnings:
    # warningschk.py derives it from its own (absolute) argv[0] instead,
    # resolving to the source tree root. build_dir (3rd, new) is passed
    # explicitly so make.log is looked up under the actual build directory
    # instead -- for an out-of-tree build (e.g. configured from a separate
    # _build/ directory) those two differ, and make.log only ever exists
    # under the latter.
    proc = subprocess.run(
        [sys.executable, str(script), warno, "", str(abinit_builddir)],
        cwd=work_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    log_lines = proc.stdout.splitlines()
    if warno in ("3", "4"):
        log_lines = sort_warnchk_body(log_lines)
    log_lines += _stderr_note_lines(proc, "warningschk.py")
    log_path.write_text("\n".join(log_lines) + ("\n" if log_lines else ""))

    out_lines = [ln for ln in log_lines if "SUCCESS" not in ln]
    out_path.write_text("\n".join(out_lines) + ("\n" if out_lines else ""))

    return classify_golden(out_path, refs_dir / f"t{case_id}.out", diff_path, case_id, tolerances)


def run_golden_check(
    case_id: str, verb: str, work_dir: Path, refs_dir: Path, tolerances: dict[str, tuple[int, float, float]]
) -> tuple[str, str]:
    """Run check_forbidden.py/check_ascii.py/check_inlined_macros.py and classify.

    Unlike warnchk, legacy never writes a t{id}.log for these three verbs and
    never filters SUCCESS lines -- raw stdout is the golden-file candidate.
    """
    script = SCRIPTS_DIR / GOLDEN_FILE_SCRIPTS[verb]
    out_path = work_dir / f"t{case_id}.out"
    diff_path = work_dir / f"diff.t{case_id}"

    proc = subprocess.run([sys.executable, str(script)], cwd=work_dir, capture_output=True, text=True, check=False)
    out_lines = proc.stdout.splitlines() + _stderr_note_lines(proc, GOLDEN_FILE_SCRIPTS[verb])
    out_path.write_text("\n".join(out_lines) + ("\n" if out_lines else ""))

    return classify_golden(out_path, refs_dir / f"t{case_id}.out", diff_path, case_id, tolerances)


def run_statchk(case_id: str, script_rel_path: str, work_dir: Path, abinit_srcdir: Path) -> tuple[str, str]:
    """Run an arbitrary script and classify by (exit_code, stderr_size), matching legacy dochkstatus.

    cwd is `abinit_srcdir`, not `work_dir` -- statchk scripts (e.g.
    check-config-h.py) locate the source tree via `os.getcwd()`, so they must
    run from the actual top level, unlike the abirules golden-file checks
    above which tolerate running from deep inside the build tree (their own
    abirules_tools helper walks up parent directories to find it).
    """
    script = abinit_srcdir / script_rel_path
    out_path = work_dir / f"t{case_id}.out"
    err_path = work_dir / f"t{case_id}.err"

    proc = subprocess.run(
        [sys.executable, str(script)], cwd=abinit_srcdir, capture_output=True, text=True, check=False
    )
    out_path.write_text(proc.stdout)
    err_path.write_text(proc.stderr)

    if proc.returncode != 0:
        return "failed", f"exit code {proc.returncode}"
    if proc.stderr:
        return "passed", f"exit 0 with {len(proc.stderr)} byte(s) on stderr"
    return "succeeded", "exit 0, no stderr output"


def main(argv: list[str]) -> int:
    """Parse CLI/env, dispatch every selected tests.cnf case, write the report file."""
    if len(argv) < 3:
        print(f"Usage: {argv[0]} machine_name series [start [stop]]", file=sys.stderr)
        return 16
    machine, series = argv[1], argv[2]

    # A single start= (no stop=) runs exactly that one case -- matches
    # run-standard-tests.pl's own arg parsing and abichecks/README.md's
    # documented "To run one check, specify only start" behavior.
    if len(argv) > 3 and argv[3]:
        first = argv[3]
        last = argv[4] if len(argv) > 4 and argv[4] else first
    else:
        first, last = "first", "end"

    try:
        abinit_srcdir = Path(os.environ["abinit_srcdir"])
    except KeyError:
        print("Error: abinit_srcdir is not set -- source abichecks.env first.", file=sys.stderr)
        return 1
    abinit_inpdir = Path(os.environ.get("abinit_inpdir", abinit_srcdir / "abichecks"))
    abinit_outdir = Path(os.environ.get("abinit_outdir", abinit_srcdir / "abichecks"))
    # Defaults to abinit_srcdir (the historical in-tree-build assumption)
    # when unset, matching warningschk.py's own build_dir default -- see
    # run_warnchk()'s comment for why an out-of-tree build needs this kept
    # separate from abinit_srcdir at all.
    abinit_builddir = Path(os.environ.get("abinit_builddir", abinit_srcdir))

    tests_cnf = abinit_inpdir / series / "tests.cnf"
    if not tests_cnf.is_file():
        print(f"Error: {tests_cnf} not found.", file=sys.stderr)
        return 20

    entries = parse_tests_cnf(tests_cnf)
    selected = select_range(entries, first, last)
    print(f"Testing the Abinit code on the {machine} platform")
    print(f"Following tests will be run: {first} to {last}")

    series_dir = abinit_outdir / series
    series_dir.mkdir(parents=True, exist_ok=True)

    ostype = platform.system()
    work_dir_name = f"tmp-{machine}_{ostype}_{time.strftime('%Y%m%d')}"
    # Remove other same-series tmp-* dirs before starting: special_processor.py
    # globs "{series}/tmp*/report" and uses the first match, which is
    # glob-order-dependent, not "most recent" -- stale sibling dirs from an
    # earlier day are a source of the wrong report being picked up. Legacy
    # never cleaned these up (only `make clean-local` did); doing it here
    # removes that whole class of bug at no cost, since no consumer needs
    # old runs kept around.
    for sibling in series_dir.glob("tmp-*"):
        if sibling.name != work_dir_name and sibling.is_dir():
            shutil.rmtree(sibling)
    work_dir = series_dir / work_dir_name
    work_dir.mkdir(exist_ok=True)

    refs_dir = abinit_inpdir / series / "Refs"
    tolerances = parse_report_in(abinit_inpdir / series / "Input" / "report.in")

    # Needed by classify_golden()'s local `from pymods.fldiff import Differ`.
    # Insert <abinit_srcdir>/tests specifically (not abinit_srcdir itself) so
    # this stays a plain `import pymods.fldiff`, not a dotted `import
    # tests.pymods.fldiff` -- the latter would execute tests/__init__.py
    # (1000+ lines, reads known_keywords.json, imports the full test-suite
    # machinery) as an import-time side effect we don't want here.
    sys.path.insert(0, str(abinit_srcdir / "tests"))

    report_lines: list[str] = []
    for case_id, verb, params in selected:
        if verb == "warnchk":
            status, _msg = run_warnchk(case_id, params[0], work_dir, refs_dir, tolerances, abinit_builddir)
        elif verb in GOLDEN_FILE_SCRIPTS:
            status, _msg = run_golden_check(case_id, verb, work_dir, refs_dir, tolerances)
        elif verb == "statchk":
            status, _msg = run_statchk(case_id, params[0], work_dir, abinit_srcdir)
        elif verb == "report":
            # No-op: the report file below is always (re)written for
            # whatever range was actually processed, regardless of whether
            # this directive was reached -- unlike legacy, where a stop=
            # bound before this trailing tests.cnf line silently skipped
            # report generation entirely for partial runs.
            continue
        else:
            print(f"Unknown verb '{verb}' for case {case_id} in {tests_cnf}", file=sys.stderr)
            continue
        print(f"[{series}][{case_id}] {status}")
        report_lines.append(f"Case_{case_id}   {status}")

    report_path = work_dir / "report"
    report_path.write_text("\n".join(report_lines) + ("\n" if report_lines else ""))

    print(f"End of tests_{series}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
