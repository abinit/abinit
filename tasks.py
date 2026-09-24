"""Developer tasks for building, testing, and maintaining ABINIT.

Install the task dependencies from the repository root with::

    python3 -m pip install --group tasks

Invoke discovers this file when run from the source tree or one of its
subdirectories. Use ``invoke --list`` for the task catalog and
``invoke --help TASK`` for task-specific options.

Common workflows::

    invoke make --jobs=8                 # Build with eight parallel jobs.
    invoke make --clean --binary=abinit  # Clean and rebuild only ABINIT.
    invoke runemall --no-make            # Test an existing build.
    invoke runemall --keywords=fast      # Build and run matching tests.
    invoke abinit --input-name=run.abi --run-make
    invoke config-log                    # Summarize the nearest config.log.
    invoke large-files --size-threshold-mb=10
    invoke doxygen                       # Build the source-code reference.
    invoke io-bench --directories=/tmp,/scratch  # Compare directory I/O throughput.

Git and release tasks can modify branches or remotes. Review their help and
ensure the working tree is clean before using them. ``official-release`` is a
simulation unless explicitly called with ``--no-dry-run``.
"""

from __future__ import annotations

import os
import platform
import subprocess
import sys
import webbrowser
from contextlib import contextmanager
from glob import glob
from pathlib import Path
from shutil import which
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator

try:
    from invoke import Context, task
except ImportError:
    raise ImportError("Cannot import invoke package. Use `pip install invoke` (or `pip install fabric` which includes invoke)")


from tests.pymods.devtools import number_of_cpus
from tests.pymods.termcolor import cprint
from tests.pymods.testsuite import find_top_build_tree

ABINIT_ROOTDIR = os.path.dirname(__file__)
ABINIT_SRCDIR = os.path.join(ABINIT_ROOTDIR, "src")

ALL_BINARIES = [
    "abinit",
    "abitk",
    "aim",
    "anaddb",
    "band2eps",
    "conducti",
    "cut3d",
    "dummy_tests",
    "fftprof",
    "fold2Bloch",
    "ioprof",
    "lapackprof",
    "macroave",
    "mrgddb",
    "mrgdv",
    "mrggkk",
    "mrgscr",
    "multibinit",
    "optic",
    "atdep",
    "testtransposer",
    "lruj",
]


SYSTEM = platform.system()


def which_vim() -> str:
    """
    Find a Vim-compatible editor in the system PATH.

    Returns:
        str: Name of the Vim executable found (mvim, nvim, or vim).

    Raises:
        RuntimeError: If no Vim executable is found.
    """
    if which("mvim") is not None:
        return "mvim"
    if which("nvim") is not None:
        return "nvim"
    if which("vim") is not None:
        return "vim"
    raise RuntimeError("Cannot find vim in $PATH!")


def which_differ() -> str:
    """
    Find a visual diff tool in the system PATH.

    Returns:
        str: Name of the differ executable found (mvimdiff or vimdiff).

    Raises:
        RuntimeError: If no differ executable is found.
    """
    differ = "vimdiff"
    if which("mvimdiff") is not None:
        differ = "mvimdiff"

    if which(differ) is None:
        raise RuntimeError(f"Cannot find {differ=} in $PATH!")

    return differ


def change_output_file(input_file: str | Path, output_file: str) -> None:
    """
    Set or update the `output_file` variable in an ABINIT input file.

    Args:
        input_file (str or Path): Path to the input file to modify.
        output_file (str): New output file name to insert.
    """
    input_path = Path(input_file)
    lines = input_path.read_text(encoding="utf-8").splitlines()
    assignment = f'output_file = "{output_file}"'

    for index, line in enumerate(lines):
        if line.lstrip().startswith("output_file"):
            lines[index] = assignment
            break
    else:
        lines.insert(0, assignment)

    input_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@contextmanager
def cd(path: str | Path) -> Iterator[None]:
    """
    A Fabric-inspired cd context that temporarily changes directory for
    performing some tasks, and returns to the original working directory
    afterwards.

    Example::

        with cd("/my/path/"):
            do_something()

    Args:
        path: Path to cd to.
    """
    # Taken from monty.os
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


def list_from_string(string: str, converter=int) -> list:
    """
    Convert a comma or space-separated string into a list of typed values.

    Args:
        string (str): The input string to parse.
        converter: Callable used to convert each value. Defaults to ``int``.

    Returns:
        list: List of converted values.
    """
    values = string.replace(",", " ").split()
    if not values:
        raise ValueError("Expected at least one value")
    return [converter(value) for value in values]


def _run_make(
    ctx: Context,
    jobs: str | int = "auto",
    touch: bool = False,
    clean: bool = False,
    binary: str = "",
) -> None:
    """Implement the :func:`make` task."""
    if touch:
        with cd(ABINIT_ROOTDIR):
            cmd = "./abisrc.py touch"
            cprint(f"Executing: {cmd}", color="yellow")
            ctx.run(cmd, pty=True)

    top = find_top_build_tree(".", with_abinit=False)
    job_count = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)

    with cd(top):
        if clean:
            ctx.run("make -C src clean", pty=True)
            ctx.run("make -C shared clean", pty=True)

        command = f"make -j{job_count}"
        if binary:
            command += f" {binary}"
        cprint(f"Executing: {command}", color="yellow")
        ctx.run(command, pty=True)


@task
def make(
    ctx: Context,
    jobs: str | int = "auto",
    touch: bool = False,
    clean: bool = False,
    binary: str = "",
) -> None:
    """Recompile the ABINIT source code.

    Examples:
        ``invoke make`` uses half of the available CPUs.
        ``invoke make --jobs=8 --binary=abinit`` builds only ``abinit``.
        ``invoke make --clean`` cleans the build tree first.

    Args:
        ctx: Invoke context.
        jobs (str or int, optional): Number of parallel threads for make.
            Use "auto" to use half of available CPUs. Defaults to "auto".
        touch (bool, optional): If True, touch all modified source files
            before recompilation. Defaults to False.
        clean (bool, optional): If True, execute `make clean` before recompiling.
            Defaults to False.
        binary (str, optional): Specific binary to recompile (e.g., "abinit").
            Defaults to all binaries.
    """
    _run_make(ctx, jobs=jobs, touch=touch, clean=clean, binary=binary)

    # TODO Check for errors in make.stderr
    # cprint("Exit code: %s" % retcode, "green" if retcode == 0 else "red")

    # if SYSTEM == "Darwin":
    #    for binary in ALL_BINARIES:
    #        cmd = f"codesign -v --force --deep src/98_main/{binary}"
    #        cprint("Executing: %s" % cmd, "yellow")
    #        ctx.run(cmd, pty=True)


@task
def clean(ctx: Context) -> None:
    """
    Remove object files in the `src` and `shared` directories.

    Note:
        Does not affect object files in the `fallbacks` directory.

    Args:
        ctx: Invoke context.
    """
    top = find_top_build_tree(".", with_abinit=False)
    with cd(top):
        ctx.run("cd src && make clean && cd ..", pty=True)
        ctx.run("cd shared && make clean && cd ..", pty=True)


@task
def runemall(
    ctx: Context,
    make: bool = True,
    jobs: str | int = "auto",
    touch: bool = False,
    clean: bool = False,
    keywords: str | None = None,
) -> None:
    """Run the sequential and parallel test suites.

    The task first ensures the code is compiled, and then executes the
    test suite. It exits immediately if any critical error occurs.

    Examples:
        ``invoke runemall --no-make`` tests the current binaries.
        ``invoke runemall --jobs=8 --keywords=fast`` selects fast tests.

    Args:
        ctx: Invoke context.
        make (bool, optional): If True, compile the code before running tests.
            Defaults to True.
        jobs (str or int, optional): Parallel threads for compilation. Defaults to "auto".
        touch (bool, optional): If True, touch modified files before recompilation.
            Defaults to False.
        clean (bool, optional): If True, perform a clean build. Defaults to False.
        keywords (str, optional): Keyword filter for selecting specific tests.
            Defaults to None.
    """
    if make:
        _run_make(ctx, jobs=jobs, touch=touch, clean=clean)

    top = find_top_build_tree(".", with_abinit=True)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)
    kws = "" if keywords is None else f"-k {keywords}"

    with cd(os.path.join(top, "tests")):
        cmd = f"./runtests.py -j{jobs} {kws}"
        cprint(f"Executing: {cmd}", color="yellow")
        ctx.run(cmd, pty=True)
        # Now run the parallel tests.
        for n in [2, 4, 10]:
            j = jobs // n
            if j == 0:
                continue
            cmd = f"./runtests.py paral mpiio -j{j} -n{n} {kws}"
            cprint(f"Executing: {cmd}", color="yellow")
            ctx.run(cmd, pty=True)


@task
def bisect(ctx: Context, start: str, end: str, runtests_args: str) -> None:
    """
    Perform a bisection search between two commits to find where a bug was introduced.

    Args:
        ctx: Invoke context.
        start: The last known GOOD commit hash.
        end: The first known BAD commit hash.
        runtests_args: String containing arguments for `runtests.py`.
    """
    from tests.pymods.termcolor import cprint

    cprint(f"Bisection started: start (good)={start}, end (bad)={end}", color="yellow")

    # Get the list of commits between start and end.
    # git rev-list end ^start gives commits from start to end.
    res = ctx.run(f"git rev-list --reverse {end} ^{start}", hide=True)
    commits = res.stdout.strip().split("\n")

    if not commits:
        cprint("No commits found between the specified points.", color="red")
        return

    cprint(f"Testing {len(commits)} commits...", color="yellow")

    low = 0
    high = len(commits) - 1
    first_bad = end

    with open("_bisection.txt", "w") as log:
        log.write(f"Bisection started: start (good)={start}, end (bad)={end}\n")
        log.write(f"Number of commits to check: {len(commits)}\n\n")

        try:
            while low <= high:
                mid = (low + high) // 2
                current_commit = commits[mid]

                msg = f"\n--- Checking commit {mid + 1}/{len(commits)}: {current_commit} ---"
                cprint(msg, color="cyan")
                log.write(msg + "\n")

                # 1. Checkout
                ctx.run(f"git checkout {current_commit}", hide=True)

                # 2. Compile (required for runtests.py)
                cprint("Compiling...", color="yellow")
                make_res = ctx.run("invoke makedeep", warn=True)
                if not make_res.ok:
                    skip_msg = f"Compilation failed for {current_commit}. Skipping this commit."
                    cprint(skip_msg, color="red")
                    log.write(skip_msg + "\n")
                    low = mid + 1
                    continue

                # 3. Run tests
                top = find_top_build_tree(".", with_abinit=True)
                with cd(os.path.join(top, "tests")):
                    cmd = f"./runtests.py {runtests_args}"
                    cprint(f"Running tests: {cmd}", color="yellow")
                    test_res = ctx.run(cmd, warn=True)

                    if test_res.ok:
                        cprint(f"Commit {current_commit} is GOOD", color="green")
                        log.write(f"{current_commit}: GOOD\n")
                        low = mid + 1
                    else:
                        cprint(f"Commit {current_commit} is BAD", color="red")
                        log.write(f"{current_commit}: BAD\n")
                        first_bad = current_commit
                        high = mid - 1

            res_msg = f"\nResult: The bug was introduced in commit {first_bad}"
            cprint(res_msg, color="magenta")
            log.write(res_msg + "\n")

            show_res = ctx.run(f"git show --summary {first_bad}", hide=True)
            log.write("\nCommit Details:\n")
            log.write(show_res.stdout)
            ctx.run(f"git show --summary {first_bad}")

        except Exception as e:
            err_msg = f"An error occurred during bisection: {e}"
            cprint(err_msg, color="red")
            log.write(err_msg + "\n")

        finally:
            cprint("\nBisection finished.", color="yellow")


@task
def makemake(ctx: Context) -> None:
    """
    Invoke the `makemake` script to rebuild the build system.

    Args:
        ctx: Invoke context.
    """
    with cd(ABINIT_ROOTDIR):
        ctx.run("./config/scripts/makemake", pty=True)


@task
def makedeep(ctx: Context, jobs: str | int = "auto") -> None:
    """
    Perform a complete rebuild cycle: makemake, clean, and build.

    Args:
        ctx: Invoke context.
        jobs (str or int, optional): Parallel threads for compilation. Defaults to "auto".
    """
    makemake(ctx)
    make(ctx, jobs=jobs, clean=True)


@task
def abichecks(ctx: Context) -> int:
    """
    Execute the Abinit sanity check scripts (abichecks).

    Returns:
        int: The number of failed check scripts.

    Args:
        ctx: Invoke context.
    """
    import time

    retcode = 0
    with cd(ABINIT_ROOTDIR):
        script_dir = os.path.join("abichecks", "scripts")
        exclude = [
            "check-libpaw.py",
            "warningschk.py",
            "abirules_tools.py",
            "__init__.py",
            # Not a standalone pass/fail checker like the others in this
            # directory -- it's the dispatcher that runs them (the pure-Python
            # replacement for the old Perl run-standard-tests.pl), and exits
            # 16 with a usage message when run with no arguments, which this
            # loop would otherwise count as a failed check.
            "run_checks.py",
        ]
        for py_script in [f for f in os.listdir(script_dir) if f.endswith(".py")]:
            if py_script in exclude:
                continue
            py_script = os.path.join(script_dir, py_script)
            print("Running", py_script, "... ")
            start = time.time()
            result = ctx.run(py_script, warn=True, pty=True)
            # print(result.ok)
            msg, color = ("[OK]", "green") if result.ok else ("[FAILED]", "red")
            cprint(f"{msg} ({time.time() - start:.2f} s)", color=color)
            if not result.ok:
                retcode += 1

    if retcode != 0:
        cprint(f"{retcode} FAILED TESTS", color="red")
    else:
        cprint("ALL TESTS OK", color="green")

    return retcode


@task
def robodoc(ctx: Context) -> bool | None:
    """
    Build the Robodoc documentation and open the index in the browser.

    Args:
        ctx: Invoke context.
    """
    with cd(ABINIT_ROOTDIR):
        result = ctx.run("./mkrobodoc.sh", pty=True)

        if result.ok:
            cprint("ROBODOC BUILD OK", color="green")
            # https://stackoverflow.com/questions/44447469/cannot-open-an-html-file-from-python-in-a-web-browser-notepad-opens-instead
            html_path = os.path.join(ABINIT_ROOTDIR, "./tmp-robodoc/www/robodoc/masterindex.html")
            print(f"Trying to open {html_path} in browser ...")
            return webbrowser.open_new_tab(html_path)
        cprint("ROBODOC BUILD FAILED", color="red")

        return result.ok


@task
def doxygen(ctx: Context, open_browser: bool = True, warnings_as_errors: bool = False) -> bool:
    """Build the Doxygen source-code reference.

    Args:
        ctx: Invoke context.
        open_browser: Open the generated index after a successful build.
        warnings_as_errors: Fail when Doxygen reports any warnings.
    """
    env = {"DOXYGEN_WARNINGS_AS_ERRORS": "1" if warnings_as_errors else "0"}
    with cd(ABINIT_ROOTDIR):
        result = ctx.run("./mkdoxygen.sh", env=env, pty=True, warn=True)

    if not result.ok:
        cprint("DOXYGEN BUILD FAILED", color="red")
        return False

    cprint("DOXYGEN BUILD OK", color="green")
    if open_browser:
        html_path = Path(ABINIT_ROOTDIR, "doxygen_docs", "html", "index.html")
        print(f"Opening {html_path} in the browser ...")
        webbrowser.open_new_tab(html_path.as_uri())
    return True


@task
def mksite(ctx: Context) -> None:
    """
    Build the Abinit documentation site and serve it locally.

    Args:
        ctx: Invoke context.
    """
    with cd(ABINIT_ROOTDIR):
        webbrowser.open_new_tab("http://127.0.0.1:8000")
        ctx.run("./mksite.py serve --dirtyreload", pty=True)


@task
def links(ctx: Context) -> None:
    """
    Create symbolic links to all Abinit executables in the current directory.

    Args:
        ctx: Invoke context.
    """
    top = find_top_build_tree(".", with_abinit=True)
    main98 = os.path.join(top, "src", "98_main")
    for dest in ALL_BINARIES:
        if os.path.islink(os.path.join(os.getcwd(), dest)):
            continue
        source = os.path.join(main98, dest)
        if os.path.isfile(source):
            os.symlink(source, dest)
        else:
            cprint(f"Cannot find `{source}` in dir `{main98}", color="yellow")


@task
def ctags(ctx: Context) -> None:
    """
    Regenerate the ctags file for the Abinit source tree.

    Args:
        ctx: Invoke context.
    """
    with cd(ABINIT_ROOTDIR):
        cmd = "ctags -R --langmap=fortran:+.finc.f90.F90,c:.c.cpp.cu shared/ src/"
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def fgrep(ctx: Context, pattern: str) -> None:
    """
    Case-insensitive search for a pattern in all Fortran and C/C++ files.

    Args:
        ctx: Invoke context.
        pattern (str): The pattern to search for.
    """
    # grep -r -i --include \*.h
    # Syntax notes:
    #    -r - search recursively
    #    -i - case-insensitive search
    #    --include=\*.${file_extension} - search files that match the extension(s) or file pattern only
    with cd(ABINIT_ROOTDIR):
        cmd = f'grep -r -i --color --include "*[.F90,.f90,.finc,.c,.cu,.cpp,.h]" "{pattern}" src shared'
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def cgrep(ctx: Context, pattern: str) -> None:
    """
    Case-insensitive search for a pattern specifically in C files.

    Args:
        ctx: Invoke context.
        pattern (str): The pattern to search for.
    """
    with cd(ABINIT_ROOTDIR):
        cmd = f'grep -r -i --color --include "*.c" "{pattern}" src shared'
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def tgrep(ctx: Context, pattern: str) -> None:
    """
    Search for a pattern within all test input files.

    Args:
        ctx: Invoke context.
        pattern (str): The pattern to search for.
    """
    with cd(ABINIT_ROOTDIR):
        cmd = f'grep -r -i --color "{pattern}" tests/*/Input/*'
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def vimt(ctx: Context, tagname: str) -> None:
    """
    Open the file defining a ctags tag in Vim.

    Args:
        ctx: Invoke context.
        tagname (str): The tag to jump to.
    """
    vim = which_vim()
    with cd(ABINIT_ROOTDIR):
        cmd = f"{vim} -f {tagname}"
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def env(ctx: Context) -> None:
    """
    Generate shell commands to configure the environment for the current build.

    Prints the `export` commands for $PATH and $ABI_PSPDIR.

    Args:
        ctx: Invoke context.
    """
    cprint("\nExecute the following lines in the shell to set the env:\n", color="green")
    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main")
    print(f"export ABI_PSPDIR={ABINIT_ROOTDIR}/tests/Pspdir")
    print(f"export PATH={binpath}:$PATH")


@task
def diff2(ctx: Context, filename: str = "run.abo") -> None:
    """
    Compare the specified output file with the most recent backup file.

    Args:
        ctx: Invoke context.
        filename (str, optional): The base output filename. Defaults to "run.abo".
    """
    vimdiff = which_differ()
    files = sorted([f for f in os.listdir(".") if f.startswith(filename)])
    if not files:
        return
    cmd = f"{vimdiff} {filename} {files[-1]}"
    cprint(f"Executing {cmd}", color="green")
    ctx.run(cmd, pty=True)


@task
def diff3(ctx: Context, filename: str = "run.abo") -> None:
    """
    Compare the current output file with the two most recent backup files.

    Args:
        ctx: Invoke context.
        filename (str, optional): The base output filename. Defaults to "run.abo".
    """
    differ = which_differ()

    files = sorted([f for f in os.listdir(".") if f.startswith(filename)])
    if not files:
        return

    cmd = f"{differ} {filename} {files[-2]} {files[-1]}" if len(files) > 2 else f"{differ} {filename} {files[-1]}"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def add_trunk(ctx: Context) -> None:
    """
    Add the main Abinit GitLab repository as a git remote named "trunk".

    Args:
        ctx: Invoke context.
    """
    cmd = "git remote add trunk https://gitlab.uliege.be/abinit/developers/trunk/abinit.git"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)
    cmd = "git fetch trunk"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def remote_add(ctx: Context, remote: str) -> None:
    """
    Register a developer's fork as a git remote and fetch it.

    Args:
        ctx: Invoke context.
        remote (str): GitLab username of the developer.
    """
    cmd = f"git remote add {remote} https://gitlab.uliege.be/abinit/developers/{remote}/abinit.git"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)
    cmd = f"git fetch {remote}"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def gdb(ctx: Context, input_name: str, exec_name: str = "abinit", run_make: bool = False) -> None:
    """
    Launch the GDB debugger for a specific executable and input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ABINIT input file.
        exec_name (str, optional): Name of the executable. Defaults to "abinit".
        run_make (bool, optional): If True, build before debugging. Defaults to False.
    """
    if run_make:
        make(ctx)

    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main", exec_name)
    cprint(f"Using binpath: {binpath}", "green")
    cmd = f"gdb {binpath} --eval-command 'set args {input_name}'"
    cprint(f"Executing gdb command: {cmd}", color="green")
    # mpirun -np 2 xterm -e gdb fftprof --command=dbg_file
    # cprint("Type run to start lldb debugger", color="green")
    # cprint("Then use `bt` to get the backtrace\n\n", color="green")
    ctx.run(cmd, pty=True)


@task
def lldb(ctx: Context, input_name: str, exec_name: str = "abinit", run_make: bool = False) -> None:
    """
    Launch the LLDB debugger for a specific executable and input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ABINIT input file.
        exec_name (str, optional): Name of the executable. Defaults to "abinit".
        run_make (bool, optional): If True, build before debugging. Defaults to False.
    """
    if run_make:
        make(ctx)

    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main", exec_name)
    cprint(f"Using binpath: {binpath}", color="green")
    cmd = f"lldb {binpath} --one-line 'settings set target.run-args {input_name}'"
    cprint(f"Executing lldb command: {cmd}", color="green")
    cprint("Type run to start lldb debugger", color="green")
    cprint("Then use `bt` to get the backtrace\n\n", color="green")
    ctx.run(cmd, pty=True)


@task
def mpi_check(
    ctx: Context,
    np_list: str = "1, 2",
    abinit_input_file: str = "run.abi",
    mpi_runner: str = "mpiexec",
    run_make: bool = False,
) -> None:
    """
    Run an ABINIT input file with various MPI process counts and compare results.

    Args:
        ctx: Invoke context.
        np_list (str): List of process counts (e.g., "1, 2, 4"). Defaults to "1, 2".
        abinit_input_file (str, optional): Path to the input file. Defaults to "run.abi".
        mpi_runner (str, optional): Command used to launch MPI. Defaults to "mpiexec".
        run_make (bool, optional): If True, build before running. Defaults to False.
    """
    if run_make:
        make(ctx)

    cprint(f"Will run {abinit_input_file=} with MPI nprocs in {np_list=}", color="yellow")

    differ = which_differ()
    np_list = list_from_string(np_list)

    for np in np_list:
        change_output_file(abinit_input_file, f"run_mpi{np}.abo")
        cmd = f"{mpi_runner} -n {np} abinit {abinit_input_file} | tee run_mpi{np}.log"
        cprint(f"About to execute {cmd=}", color="yellow")
        ctx.run(cmd)

    np_ref = np_list[0]
    for np in np_list[1:]:
        cmd = f"{differ} run_mpi{np_ref}.abo run_mpi{np}.abo"
        cprint(f"About to execute {cmd=}", color="yellow")
        ctx.run(cmd, pty=True)
        cmd = f"{differ} run_mpi{np_ref}.log run_mpi{np}.log"
        cprint(f"About to execute {cmd=}", color="yellow")
        ctx.run(cmd, pty=True)


@task
def omp_check(
    ctx: Context,
    omp_threads: str = "1, 2",
    np: int = 1,
    abinit_input_file: str = "run.abi",
    mpi_runner: str = "mpiexec",
    run_make: bool = False,
) -> None:
    """
    Run an ABINIT input file with various OpenMP thread counts and compare results.

    Args:
        ctx: Invoke context.
        omp_threads (str): List of thread counts (e.g., "1, 2, 4"). Defaults to "1, 2".
        np (int, optional): Number of MPI processes to use. Defaults to 1.
        abinit_input_file (str, optional): Path to the input file. Defaults to "run.abi".
        mpi_runner (str, optional): Command used to launch MPI. Defaults to "mpiexec".
        run_make (bool, optional): If True, build before running. Defaults to False.
    """
    if run_make:
        make(ctx)

    differ = which_differ()
    omp_threads = list_from_string(omp_threads)

    cprint(
        f"Will run {abinit_input_file=} with OMP threads={omp_threads!r} and MPI nprocs={np}",
        color="yellow",
    )
    for nth in omp_threads:
        change_output_file(abinit_input_file, f"run_omp{nth}_mpi{np}.abo")
        cmd = f"OMP_NUM_THREADS={nth} {mpi_runner} -n {np} abinit {abinit_input_file} | tee run_omp{nth}_mpi{np}.log"
        cprint(f"About to execute {cmd=}", color="yellow")
        ctx.run(cmd)

    omp_ref = omp_threads[0]
    for nth in omp_threads[1:]:
        cmd = f"{differ} run_omp{omp_ref}_mpi{np}.abo run_omp{nth}_mpi{np}.abo"
        cprint(f"About to execute {cmd=}", color="yellow")
        ctx.run(cmd, pty=True)
        cmd = f"{differ} run_omp{omp_ref}_mpi{np}.log run_omp{nth}_mpi{np}.log"
        cprint(f"About to execute {cmd=}", color="yellow")
        ctx.run(cmd, pty=True)


@task
def pyenv_clean(ctx: Context) -> None:
    """
    Purge the Conda and Pip caches.

    Args:
        ctx: Invoke context.
    """
    if which("conda") is not None:
        cmd = "conda clean --all --yes"
        cprint(f"About to execute {cmd=}")
        ctx.run(cmd)

    cmd = "pip cache purge"
    cprint(f"About to execute {cmd=}")
    ctx.run(cmd)


@task
def abinit(ctx: Context, input_name: str, run_make: bool = False) -> None:
    """
    Run the `abinit` executable with the specified input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ABINIT input file.
        run_make (bool, optional): If True, build before running. Defaults to False.
    """
    _run(ctx, input_name, exec_name="abinit", run_make=run_make)


@task
def anaddb(ctx: Context, input_name: str, run_make: bool = False) -> None:
    """
    Run the `anaddb` executable with the specified input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ANADDB input file.
        run_make (bool, optional): If True, build before running. Defaults to False.
    """
    _run(ctx, input_name, exec_name="anaddb", run_make=run_make)


def _run(ctx: Context, input_name: str, exec_name: str, run_make: bool):
    """
    Internal helper to execute an Abinit binary with an input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the input file.
        exec_name (str): Name of the executable.
        run_make (bool): If True, build before running.
    """
    if run_make:
        make(ctx)
    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main", exec_name)
    cprint(f"Using binpath: {binpath}", color="green")
    cmd = f"{binpath} {input_name}"
    cprint(f"Executing {cmd}", color="green")
    ctx.run(cmd, pty=True)


@task
def pull_trunk(ctx: Context) -> None:
    """
    Update the current branch from the trunk's develop branch, maintaining
    local changes via stash.

    Args:
        ctx: Invoke context.
    """
    ctx.run("git stash")
    ctx.run("git pull trunk develop")
    ctx.run("git pull trunk develop --tags")
    ctx.run("git commit")
    ctx.run("git push")
    ctx.run("git push --tags")
    ctx.run("git stash apply")


@task
def pull(ctx: Context) -> None:
    """
    Update the current repository and its submodules, maintaining
    local changes via stash.

    Args:
        ctx: Invoke context.
    """
    ctx.run("git stash")
    ctx.run("git pull --recurse-submodules")
    ctx.run("git stash apply")
    makemake(ctx)


@task
def push(ctx: Context) -> None:
    """
    Commit and push current changes including tags to the origin repository.

    Args:
        ctx: Invoke context.
    """
    ctx.run("git commit")
    ctx.run("git push")
    ctx.run("git push --tags")


@task
def submodules(ctx: Context) -> None:
    """
    Update all Abinit submodules to their latest remote versions.

    Args:
        ctx: Invoke context.
    """
    with cd(ABINIT_ROOTDIR):
        # https://stackoverflow.com/questions/1030169/easy-way-to-pull-latest-of-all-git-submodules
        ctx.run("git submodule update --remote --init", pty=True)
        ctx.run("git submodule update --recursive --remote", pty=True)


@task
def branchoff(ctx: Context, start_point: str) -> None:
    """
    Create a new local branch starting from a remote branch (e.g., trunk/develop).

    Automatically sets the upstream to origin for the new branch.

    Args:
        ctx: Invoke context.
        start_point (str): The remote branch to start from (e.g., "trunk/release-9.0").
    """
    if "/" in start_point:
        remote, _branch = start_point.split("/", maxsplit=1)
        remote_ref = start_point
    else:
        remote = "trunk"
        remote_ref = f"{remote}/{start_point}"

    def run(cmd: str):
        """Execute a shell command via context."""
        cprint(f"Executing: `{cmd}`", color="green")
        ctx.run(cmd)

    run(f"git fetch {remote}")
    # Create new branch `test_v9.0` using trunk/release-9.0 as start_point:
    # git checkout [-q] [-f] [-m] [[-b|-B|--orphan] <new_branch>] [<start_point>]
    my_branch = "my_" + remote_ref.replace("/", "-")
    run(f"git checkout -b {my_branch} {remote_ref}")
    # Change default upstream. If you forget this step, you may push to trunk.
    run("git push --set-upstream origin HEAD")


@task
def dryrun_merge(ctx: Context, start_point: str) -> None:
    """
    Perform a dry-run merge of a remote branch into the current branch.

    Args:
        ctx: Invoke context.
        start_point (str): The remote branch to merge from.
    """

    def run(cmd):
        """Execute a shell command via context."""
        cprint(f"Executing: `{cmd}`", color="green")
        ctx.run(cmd)

    run(f"git merge --no-commit --no-ff {start_point}")

    print("""
To examine the staged changes:

    $ git diff --cached

And you can undo the merge, even if it is a fast-forward merge:

$ git merge --abort
""")


@task
def prune_branches(ctx: Context, base: str = "develop", dry_run: bool = True) -> None:
    """
    List (and optionally delete) local branches already merged into `base`.

    Never touches the currently checked-out branch or `develop`/`master`/`main`,
    and only ever runs `git branch -d` (safe delete, refuses an unmerged branch),
    never `-D`.

    Examples:
        ``invoke prune-branches`` lists merged branches without deleting them.
        ``invoke prune-branches --no-dry-run`` deletes them.
        ``invoke prune-branches --base=my_trunk-release-9.0`` checks against another branch.

    Args:
        ctx: Invoke context.
        base (str, optional): Branch to check "merged into". Defaults to "develop".
        dry_run (bool, optional): If True, only list candidates. Defaults to True.
    """
    protected = {"develop", "master", "main", get_current_branch()}

    result = ctx.run(f'git branch --merged {base} --format="%(refname:short)"', hide=True, warn=True)
    if not result.ok:
        cprint(f"Could not list branches merged into {base!r}: {result.stderr.strip()}", color="red")
        return

    candidates = [b.strip() for b in result.stdout.splitlines() if b.strip() and b.strip() not in protected]

    if not candidates:
        cprint(f"No local branches to prune (already merged into {base!r}).", color="green")
        return

    cprint(f"Local branches already merged into {base!r}:", color="yellow")
    for branch in candidates:
        print(f"  {branch}")

    if dry_run:
        cprint("\nDry run: nothing deleted. Re-run with --no-dry-run to delete these branches.", color="yellow")
        return

    for branch in candidates:
        cmd = f"git branch -d {branch}"
        cprint(f"Executing: {cmd}", color="green")
        ctx.run(cmd)


@task
def watchdog(ctx: Context, jobs: str | int = "auto", sleep_time: int = 5) -> None:
    """
    Monitor the source directory for changes and trigger recompilation automatically.

    Args:
        ctx: Invoke context.
        jobs (str or int, optional): Parallel threads for make. Defaults to "auto".
        sleep_time (int, optional): Sleep time in seconds between checks. Defaults to 5.
    """
    cprint(
        "Starting watchdog service to watch F90 files and execute `make` when changes are detected",
        color="green",
    )
    cprint("Enter <CTRL + C> in the terminal to kill the service.", color="green")

    cprint(f"Start watching F90 files with sleep_time {sleep_time} s ....", color="green")
    top = find_top_build_tree(".", with_abinit=True)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)

    # http://thepythoncorner.com/dev/how-to-create-a-watchdog-in-python-to-look-for-filesystem-changes/
    # https://stackoverflow.com/questions/19991033/generating-multiple-observers-with-python-watchdog
    import time

    from watchdog.events import PatternMatchingEventHandler
    from watchdog.observers import Observer

    event_handler = PatternMatchingEventHandler(
        patterns="*.F90",
        ignore_patterns="",
        ignore_directories=False,
        case_sensitive=True,
    )

    def on_created(event):
        """Handle file creation events."""
        print(f"hey, {event.src_path} has been created!")

    def on_deleted(event):
        """Handle file deletion events."""
        print(f"File deleted: {event.src_path}")

    def on_modified(event):
        """Trigger parallel make when a watched file is modified."""
        print(f"hey buddy, {event.src_path} has been modified")
        cmd = f"make -j{jobs} > >(tee -a make.log) 2> >(tee -a make.stderr >&2)"
        cprint(f"Executing: {cmd}", color="yellow")
        with cd(top):
            try:
                result = ctx.run(cmd, pty=True)
                if result.ok:
                    cprint("Make completed successfully", color="green")
                    cprint("Watching for changes ...", color="green")
            except Exception:
                cprint("Make returned non-zero exit status", color="red")
                cprint(
                    "Keep on watching for changes hoping you get it right ...",
                    color="red",
                )

    def on_moved(event):
        """Handle file rename or move events."""
        print(f"ok ok ok, someone moved {event.src_path} to {event.dest_path}")

    event_handler.on_created = on_created
    event_handler.on_deleted = on_deleted
    event_handler.on_modified = on_modified
    event_handler.on_moved = on_moved

    observer = Observer()
    path = ABINIT_SRCDIR
    observer.schedule(event_handler, path, recursive=True)
    observer.start()

    try:
        while True:
            time.sleep(sleep_time)
    except KeyboardInterrupt:
        observer.stop()
        observer.join()


def get_current_branch() -> str:
    """
    Get the name of the currently active git branch.

    Returns:
        str: The branch name.

    Raises:
        RuntimeError: If not in a git repository or git command fails.
    """
    try:
        return subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True).strip()
    except subprocess.CalledProcessError as exc:
        raise RuntimeError("Not inside a Git repository or git command failed") from exc


def get_git_tags() -> list[str]:
    """
    Retrieve a list of all git tags in the repository.

    Returns:
        list[str]: List of tag names.

    Raises:
        RuntimeError: If not in a git repository or git command fails.
    """
    try:
        return subprocess.check_output(["git", "tag"], text=True).splitlines()
    except subprocess.CalledProcessError as exc:
        raise RuntimeError("Not inside a Git repository or git command failed") from exc


@task
def official_release(ctx: Context, new_version: str, dry_run: bool = True) -> None:
    """
    Automate the process of creating a new official Abinit release.

    This involves merging develop into master, tagging, and pushing to
    remote repositories.

    Examples:
        ``invoke official-release --new-version=10.0.0`` prints the plan.
        Add ``--no-dry-run`` to perform the release.

    Args:
        ctx: Invoke context.
        new_version (str): The version string for the new release.
        dry_run (bool, optional): If True, only simulate the steps. Defaults to True.
    """
    # Set variables
    github_user = "gonzex"
    github_repo = "abinit"
    github_url = f"git@github.com:{github_user}/{github_repo}.git"

    _run_kwargs = dict(pty=True, echo=True)

    def _run(command: str):
        """Run a release command, or display it when simulating."""
        if dry_run:
            cprint(f"Would execute: {command}", color="yellow")
            return None
        return ctx.run(command, **_run_kwargs)

    current_branch = get_current_branch()
    if current_branch != "develop":
        raise RuntimeError(f"You are on the '{current_branch}' branch, not 'develop'.")

    old_tags = get_git_tags()
    if new_version in old_tags:
        raise RuntimeError(f"{new_version=} is already in {old_tags=}")

    # List of files that should be added to master and then removed in develop
    configure_paths = [
        "configure",
        "config/gnu/compile",
        "config/gnu/config.guess",
        "config/gnu/config.sub",
        "config/gnu/install-sh",
        "config/gnu/missing",
        "config/gnu/depcomp",
    ]

    with cd(ABINIT_ROOTDIR):
        # The version in .current_version is updated manually.
        # Here we check that the value stored in the file is equal to the command line argument.
        with open(".current_version") as fh:
            old_version = fh.read().strip()

        if old_version != new_version:
            raise ValueError(f"{old_version=} != {new_version=}")

        # Step 1: Checkout master, merge changes and run makemake
        _run("git checkout master")
        _run("git merge develop")
        _run("./config/scripts/makemake")

        # Add files required by configure.
        for path in configure_paths:
            _run(f"git add -f {path}")

        _run(f"git commit -a -m 'v{new_version}'")
        _run(f"git tag -a {new_version} -m 'v{new_version}'")
        _run("git push origin master")

        # Step 2: Push to GitHub
        _run(f"git remote add abinit {github_url} || echo 'Remote already exists but this is not critical'")
        _run("git push -u abinit master --tags")

        _run("git checkout develop")
        _run("git merge master")
        _run("git push --tags")

        # Step 3: Ensure 'configure_paths' are ignored in develop branch and commit changes.
        for path in configure_paths:
            _run(f"git rm --cached {path}")
        _run("git commit -a -m 'Remove configure files from tracking in develop'")
        _run("git push origin develop")


@task
def git_info(ctx: Context, top_n: int = 20) -> None:
    """
    Analyze git history to find the largest files ever committed.

    Args:
        ctx: Invoke context.
        top_n (int, optional): Number of top files to display. Defaults to 20.
    """

    def get_git_objects():
        """Return list of all Git objects (hash, path)."""
        result = subprocess.run(
            ["git", "rev-list", "--objects", "--all"],
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )
        objects = []
        for line in result.stdout.splitlines():
            parts = line.split(" ", 1)
            if len(parts) == 2:
                objects.append((parts[0], parts[1]))
        return objects

    def get_blob_sizes(hashes):
        """Return a dict of {hash: (size_in_bytes, path)} for blobs."""
        input_text = "\n".join(hashes)
        result = subprocess.run(
            [
                "git",
                "cat-file",
                "--batch-check=%(objectname) %(objecttype) %(objectsize)",
            ],
            input=input_text,
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )

        sizes = {}
        for line in result.stdout.splitlines():
            obj_hash, obj_type, obj_size = line.split()
            if obj_type == "blob":
                sizes[obj_hash] = int(obj_size)
        return sizes

    print("Scanning Git history for largest files...")

    objects = get_git_objects()
    hashes = [obj[0] for obj in objects]
    paths = {obj[0]: obj[1] for obj in objects}

    sizes = get_blob_sizes(hashes)

    sorted_blobs = sorted(
        ((size, paths[_hash], _hash) for _hash, size in sizes.items() if _hash in paths),
        reverse=True,
    )

    print(f"\nTop {top_n} largest files ever committed:")
    for size, path, _obj_hash in sorted_blobs[:top_n]:
        print(f"{size / (1024 * 1024):7.2f} MB\t{path}")

    ctx.run("git count-objects -vH", pty=True)


@task
def large_files(ctx: Context, top_dir: str | Path | None = None, size_threshold_mb: int = 5) -> None:
    """
    Find and list files larger than `size_threshold_mb` megabytes under `top_dir`.

    Args:
        ctx: Invoke context.
        top_dir: Root directory to start searching from.
        size_threshold_mb: Minimum file size in MB to report.
    """
    large_files = []
    threshold_bytes = size_threshold_mb * 1024 * 1024

    if top_dir is None:
        top_dir = ABINIT_ROOTDIR

    exclude_dirs = {".git", "__pycache__", ".ruff_cache", "modules_with_data", "site"}
    print(f"Scanning files starting from {top_dir=}")

    for root, dirs, files in os.walk(top_dir):
        dirs[:] = [d for d in dirs if d not in exclude_dirs]
        dirs[:] = [d for d in dirs if not d.startswith("_build")]

        for name in files:
            path = os.path.join(root, name)
            try:
                size = os.path.getsize(path)
                if size > threshold_bytes:
                    large_files.append((size / (1024 * 1024), Path(path)))
            except OSError:
                # skip unreadable files
                continue

    large_files.sort(reverse=True)

    for size_mb, path in large_files:
        print(f"{size_mb:.2f} MB\t{path}")


@task
def disk_usage(ctx: Context, top_dir: str | Path | None = None, depth: int = 1, top_n: int = 20) -> None:
    """
    Summarize disk usage per subdirectory, similar to `du -d depth | sort -rh`.

    Complements `large_files`, which only reports individual oversized files
    and can miss a directory that's large because of many small files
    (e.g. an accumulated pile of old build trees).

    Args:
        ctx: Invoke context.
        top_dir: Root directory to scan. Defaults to the Abinit source root.
        depth (int, optional): Number of path components below `top_dir` to
            group by. Defaults to 1 (immediate subdirectories).
        top_n (int, optional): Number of largest entries to display. Defaults to 20.
    """
    if top_dir is None:
        top_dir = ABINIT_ROOTDIR
    top_dir = Path(top_dir)
    depth = max(1, depth)

    exclude_dirs = {".git", "__pycache__", ".ruff_cache"}

    sizes: dict[Path, int] = {}
    for root, dirs, files in os.walk(top_dir):
        dirs[:] = [d for d in dirs if d not in exclude_dirs]
        root_path = Path(root)
        rel_parts = root_path.relative_to(top_dir).parts
        group_key = top_dir if len(rel_parts) < depth else top_dir.joinpath(*rel_parts[:depth])

        total = 0
        for name in files:
            try:
                total += (root_path / name).stat().st_size
            except OSError:
                continue
        sizes[group_key] = sizes.get(group_key, 0) + total

    ranked = sorted(sizes.items(), key=lambda kv: kv[1], reverse=True)[:top_n]

    print(f"Disk usage under {top_dir} (grouped by {depth} path component(s)):\n")
    for path, size in ranked:
        print(f"{size / (1024**2):10.2f} MB\t{path}")


def _fsync_full(fh) -> None:
    """Flush a file to physical disk, bypassing any OS write-back cache.

    Plain ``os.fsync`` does not guarantee data reaches the physical disk on
    macOS (the kernel may just hand it to the drive's own volatile cache), so
    write timings would otherwise mostly measure that cache instead of the
    directory's real I/O performance. ``F_FULLFSYNC`` is Apple's documented
    way to force a real flush; other platforms are fine with ``os.fsync``.
    """
    fh.flush()
    if SYSTEM == "Darwin":
        import fcntl

        fcntl.fcntl(fh.fileno(), fcntl.F_FULLFSYNC)
    else:
        os.fsync(fh.fileno())


@task
def io_bench(ctx: Context, directories: str, size_gb: float = 1.0, chunk_mb: int = 4, iterations: int = 1) -> None:
    """Compare the I/O performance of one or more directories.

    For each directory, writes `iterations` temporary file(s) totalling
    `size_gb` gigabytes, measures the write throughput, rereads them to
    measure the read throughput, and removes them afterwards. Use
    `iterations` > 1 to approximate a many-small-files access pattern
    (each file gets its own open/close and fsync) instead of one large
    sequential file.

    Note:
        The reread may be served (partly) from the OS page cache rather than
        the physical disk, so read numbers can look better than a cold read
        would be. There is no portable way to drop the cache without root.

    Examples:
        ``invoke io-bench --directories=/tmp``
        ``invoke io-bench --directories=/tmp,/scratch --size-gb=2``
        ``invoke io-bench --directories=/tmp --size-gb=1 --iterations=1000``  # 1000 x ~1MB files

    Args:
        ctx: Invoke context.
        directories (str): Comma-separated list of directories to benchmark.
        size_gb (float, optional): Total size in GB written per directory. Defaults to 1.0.
        chunk_mb (int, optional): Chunk size in MB used for writing/reading each file. Defaults to 4.
        iterations (int, optional): Number of separate files `size_gb` is split across.
            Defaults to 1 (one big file).
    """
    import time
    import uuid

    from tabulate import tabulate

    dirs = [d.strip() for d in directories.split(",") if d.strip()]
    if not dirs:
        raise ValueError("Expected at least one directory")
    if iterations < 1:
        raise ValueError(f"{iterations=} must be >= 1")

    size_bytes = int(size_gb * 1024**3)
    chunk_bytes = int(chunk_mb * 1024**2)
    if size_bytes <= 0:
        raise ValueError(f"{size_gb=} must be > 0")
    if chunk_bytes <= 0:
        raise ValueError(f"{chunk_mb=} must be > 0")

    file_size_bytes = size_bytes // iterations
    if file_size_bytes <= 0:
        raise ValueError(f"{size_gb=} split across {iterations=} gives a 0-byte file; lower iterations or raise size_gb")

    cprint(
        f"Benchmarking {len(dirs)} directory(ies) with {iterations} file(s) of "
        f"{file_size_bytes / 1024**2:.2f} MB each ({size_gb} GB total) ...",
        color="yellow",
    )

    rows = []
    for directory in dirs:
        dir_path = Path(directory)
        if not dir_path.is_dir():
            cprint(f"Skipping {directory!r}: not a directory", color="red")
            continue

        tmp_paths = [dir_path / f".io_bench_{uuid.uuid4().hex}.tmp" for _ in range(iterations)]
        chunk = os.urandom(min(chunk_bytes, file_size_bytes))
        try:
            written = 0
            start = time.perf_counter()
            for tmp_path in tmp_paths:
                remaining = file_size_bytes
                with open(tmp_path, "wb") as fh:
                    while remaining > 0:
                        to_write = chunk if remaining >= len(chunk) else chunk[:remaining]
                        fh.write(to_write)
                        remaining -= len(to_write)
                        written += len(to_write)
                    _fsync_full(fh)
            write_time = time.perf_counter() - start

            read_bytes = 0
            start = time.perf_counter()
            for tmp_path in tmp_paths:
                with open(tmp_path, "rb") as fh:
                    while True:
                        data = fh.read(chunk_bytes)
                        if not data:
                            break
                        read_bytes += len(data)
            read_time = time.perf_counter() - start
        finally:
            for tmp_path in tmp_paths:
                tmp_path.unlink(missing_ok=True)

        write_mb_s = (written / 1024**2) / write_time if write_time > 0 else float("inf")
        read_mb_s = (read_bytes / 1024**2) / read_time if read_time > 0 else float("inf")
        rows.append(
            [directory, iterations, f"{write_mb_s:.1f} MB/s", f"{write_time:.2f} s", f"{read_mb_s:.1f} MB/s", f"{read_time:.2f} s"]
        )

    if not rows:
        cprint("No valid directory to benchmark.", color="red")
        return

    print()
    print(tabulate(rows, headers=["Directory", "Files", "Write speed", "Write time", "Read speed", "Read time"], tablefmt="grid"))


@task
def doctor(ctx: Context) -> bool:
    """
    Check that common tools needed to build, test, and debug Abinit are on PATH.

    Examples:
        ``invoke doctor``

    Args:
        ctx: Invoke context.

    Returns:
        bool: True if every required tool was found, False otherwise.
    """
    from tabulate import tabulate

    required = [
        ("Make", ["make"]),
        ("Git", ["git"]),
        ("Fortran compiler", [os.environ.get("FC", ""), "gfortran", "ifort", "ifx"]),
        ("MPI launcher", ["mpirun", "mpiexec", "srun"]),
        ("Python 3", ["python3"]),
    ]
    optional = [
        ("CMake", ["cmake"]),
        ("Doxygen", ["doxygen"]),
        ("Ctags", ["ctags"]),
        ("GDB", ["gdb"]),
        ("LLDB", ["lldb"]),
        ("Vim-compatible editor", ["mvim", "nvim", "vim"]),
    ]

    def _first_found(names: list[str]) -> str | None:
        """Return the path of the first name in `names` found on PATH, if any."""
        for name in names:
            if name and (path := which(name)):
                return path
        return None

    rows = []
    all_required_ok = True
    for label, names in required:
        path = _first_found(names)
        if not path:
            all_required_ok = False
        rows.append([label, "required", path or "NOT FOUND"])

    for label, names in optional:
        path = _first_found(names)
        rows.append([label, "optional", path or "not found"])

    print(tabulate(rows, headers=["Tool", "Kind", "Location"], tablefmt="grid"))

    if all_required_ok:
        cprint("\nAll required tools found.", color="green")
    else:
        cprint("\nSome required tools are missing -- builds are likely to fail.", color="red")

    return all_required_ok


@task
def system(ctx: Context) -> None:
    """
    Display comprehensive system information as a formatted table.

    Includes OS, Kernel, Architecture, Processor, CPU Cores, Memory,
    and Cache details.

    Args:
        ctx: Invoke context.
    """
    import platform

    import psutil
    from tabulate import tabulate

    info = []
    info.append(["OS", f"{platform.system()} {platform.release()}"])
    info.append(["Kernel", platform.version()])
    info.append(["Architecture", platform.machine()])
    info.append(["Processor", platform.processor()])
    info.append(["CPU Cores (Physical)", psutil.cpu_count(logical=False)])
    info.append(["CPU Cores (Logical)", psutil.cpu_count()])
    info.append(["Memory (Total)", f"{psutil.virtual_memory().total / (1024**3):.2f} GB"])
    for level, size in get_cache_info().items():
        info.append([level, size])

    print(tabulate(info, headers=["Item", "Value"], tablefmt="grid"))


@task
def pid(ctx: Context, pid: int | str) -> None:
    """
    Display detailed information for a specific process ID (PID).

    Args:
        ctx: Invoke context.
        pid (int or str): The PID of the process to inspect.
    """
    import psutil
    from tabulate import tabulate

    pid = int(pid)
    try:
        p = psutil.Process(pid)
        info = []
        info.append(["PID", p.pid])
        info.append(["Name", p.name()])
        info.append(["Executable", p.exe()])
        info.append(["Command Line", " ".join(p.cmdline())])
        info.append(["Status", p.status()])
        info.append(["User", p.username()])
        info.append(["CPU %", f"{p.cpu_percent(interval=0.1):.1f} %"])
        info.append(["Memory %", f"{p.memory_percent():.2f} %"])
        info.append(["Memory RSS", f"{p.memory_info().rss / (1024**2):.2f} MB"])
        info.append(["Threads", p.num_threads()])
        info.append(["CWD", p.cwd()])
        info.append(["Parent PID", p.ppid()])
        info.append(["Start Time (Epoch)", int(p.create_time())])
        print(tabulate(info, headers=["Item", "Value"], tablefmt="grid"))

    except psutil.NoSuchProcess:
        print(f"❌ Process with PID {pid} does not exist.")
        sys.exit(1)


def get_cache_info() -> dict[str, str]:
    """
    Retrieve CPU cache information for the current platform.

    Returns:
        dict: Mapping of cache levels to their sizes.

    Raises:
        RuntimeError: If the platform is not supported.
    """
    system = platform.system()
    if system == "Linux":
        return get_cache_info_linux()
    if system == "Darwin":
        return get_cache_info_mac()
    if system == "Windows":
        return get_cache_info_windows()
    raise RuntimeError(f"Unsupported platform {system}")


def get_cache_info_linux() -> dict[str, str]:
    """
    Retrieve CPU cache information on Linux systems using sysfs.

    Returns:
        dict: Mapping of cache levels to their sizes.
    """
    caches = {}
    base_path = "/sys/devices/system/cpu/cpu0/cache"
    for index_dir in glob(f"{base_path}/index*"):
        try:
            with open(os.path.join(index_dir, "level")) as f:
                level = f.read().strip()
            with open(os.path.join(index_dir, "type")) as f:
                cache_type = f.read().strip()
            with open(os.path.join(index_dir, "size")) as f:
                size = f.read().strip()
            caches[f"L{level} {cache_type}"] = size
        except Exception:
            continue
    return caches


def get_cache_info_mac() -> dict[str, str]:
    """
    Retrieve CPU cache information on macOS using sysctl.

    Returns:
        dict: Mapping of cache levels to their sizes.
    """
    caches = {}
    keys = {
        "hw.l1dcachesize": "L1d",
        "hw.l1icachesize": "L1i",
        "hw.l2cachesize": "L2",
        "hw.l3cachesize": "L3",
    }
    for key, label in keys.items():
        try:
            out = subprocess.check_output(["sysctl", "-n", key]).decode().strip()
            caches[label] = f"{int(out) // 1024} KB"
        except Exception:
            continue
    return caches


def get_cache_info_windows() -> dict[str, str]:
    """
    Retrieve CPU cache information on Windows using WMIC.

    Returns:
        dict: Mapping of cache levels to their sizes.
    """
    caches = {}
    try:
        out = subprocess.check_output(["wmic", "cpu", "get", "L2CacheSize,L3CacheSize"], stderr=subprocess.DEVNULL).decode()
        lines = out.strip().split("\n")
        if len(lines) >= 2:
            _, l2, l3 = lines[1].split()
            if l2:
                caches["L2"] = f"{l2} KB"
            if l3:
                caches["L3"] = f"{l3} KB"
    except Exception:
        pass
    return caches


def _extract_errors(logfile: str | Path, context_lines: int = 5) -> list[str]:
    """
    Parse a log file to extract error messages with surrounding context.

    Args:
        logfile (str or Path): The path to the log file to analyze.
        context_lines (int, optional): Number of context lines to capture
            around each error. Defaults to 5.

    Returns:
        list[str]: A list of formatted error blocks.
    """
    print(f"Extracting error lines and some context from {logfile}...")

    with open(logfile, errors="ignore") as f:
        lines = f.readlines()

    # Common patterns indicating critical problems
    ERROR_PATTERNS = [
        r"error",  # generic errors
        r"fail",  # tests failing
        r"cannot\s+find",  # missing library or header
        r"no\s+such\s+file",  # missing file
        r"not\s+found",  # program not found
        r"undefined\s+reference",  # linking errors
    ]

    import re

    regex = re.compile("|".join(ERROR_PATTERNS), re.IGNORECASE)
    n = len(lines)
    errors = []
    for i, line in enumerate(lines):
        if regex.search(line):
            # Ignore matches such as `sd_yakl_options='optional fail'
            if line.startswith("sd_") and "fail" in line:
                continue
            # Capture context
            start = max(0, i - context_lines)
            end = min(n, i + context_lines + 1)
            context = "".join(lines[start:end])
            errors.append(context.strip())

    return errors


def find_filename(filename: str, start_dir: Path | None = None) -> Path:
    """
    Search for a file by walking upwards from a starting directory.

    Args:
        filename (str): The filename to search for.
        start_dir (Path, optional): Directory to start the search from.
            Defaults to the current working directory.

    Returns:
        Path: The absolute path to the found file.

    Raises:
        FileNotFoundError: If the file is not found before reaching the root.
    """
    current = Path.cwd() if start_dir is None else start_dir.resolve()
    while True:
        candidate = current / filename
        if candidate.exists():
            return candidate
        if current.parent == current:  # reached filesystem root
            raise FileNotFoundError(f"File '{filename}' not found.")
        current = current.parent


@task
def config_log(ctx: Context, log_path: str = "config.log") -> None:
    """
    Analyze the `config.log` file generated by `configure` to find critical errors.

    Args:
        ctx: Invoke context.
        log_path (str, optional): Path to the log file. Defaults to "config.log".
    """
    log_path = find_filename(log_path)
    results = _extract_errors(log_path)
    if not results:
        print("✅ No critical errors detected.")
        return

    print("❌ Critical errors found:\n")
    for idx, block in enumerate(results, start=1):
        print(f"--- Error block #{idx} ---")
        print(block)
        print("-" * 40)


@task
def nvidia_prof(ctx, sh_path="nv_prof.sh"):
    """Generate a shell script to run an Abinit input file with GPU and profile it with nsys."""
    sh_template = r"""\
#!/bin/bash
set -x
set -e

invoke make -b abinit

export OMP_TARGET_OFFLOAD=MANDATORY
#export LIBOMPTARGET_INFO=4     # LLVM
#export NVCOMPILER_OMP_DEBUG=1  # NVHPC

# Set OpenMP environment
export nt=1
echo "Running ABINIT with OMP_NUM_THREADS=${nt}"
export OMP_NUM_THREADS=${nt}
#export OMP_PLACES=cores
#export OMP_PROC_BIND=close

# GPU version
rm profile_*
nsys profile \
  --trace=cuda,openmp,nvtx \
  --cuda-memory-usage=true \
  -o profile_run \
  mpirun -n 1 abinit run_gpu.abi | tee run_gpu.log

nsys stats profile_run.nsys-rep | tee prof.out
#nsys-ui profile_run.nsys-rep

#ncu \
#  --set roofline \
#  --kernel-name your_kernel_name \
#  mpirun -n 1 abinit run_gpu.abi | tee run_gpu.log

#vimdiff run_gpu.abo ref_cpu.abo
#vimdiff run_gpu.abo ref_cpu.log
"""

    with open(sh_path, "w") as fh:
        fh.write(sh_template)
