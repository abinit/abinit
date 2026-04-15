"""
Pyinvoke file for automating build and configuration tasks within the Abinit repository.

This file can be executed from any location within the Abinit directory structure,
including build directories.

General usage:
    To list available commands:
        invoke --list

    To run specific tasks (e.g., abichecks):
        invoke abichecks
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

try:
    from invoke import task
except ImportError:
    raise ImportError("Cannot import invoke package. Use `pip install invoke`")

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
    if which("mvim") is not None: return "mvim"
    if which("nvim") is not None: return "nvim"
    if which("vim") is not None: return "vim"
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
    if which("mvimdiff") is not None: differ = "mvimdiff"

    if which(differ) is None:
        raise RuntimeError(f"Cannot find {differ=} in $PATH!")

    return differ


def change_output_file(input_file, output_file):
    """
    Set or update the `output_file` variable in an ABINIT input file.

    Args:
        input_file (str or Path): Path to the input file to modify.
        output_file (str): New output file name to insert.
    """
    with open(input_file) as fh:
        remove_iline = None
        lines = [l.lstrip() for l in fh]
        for i, l in enumerate(lines):
            if l.lstrip().startswith("output_file"):
                remove_iline = i
                break

    if remove_iline is not None:
        lines.pop(remove_iline)

    lines.insert(0, f'output_file = "{output_file}"')
    with open(input_file, "w") as fh:
        fh.write("\n".join(lines))


@contextmanager
def cd(path):
    """
    A Fabric-inspired cd context that temporarily changes directory for
    performing some tasks, and returns to the original working directory
    afterwards. E.g.,

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


def list_from_string(string, type=int) -> list[str]:
    """
    Convert a comma or space-separated string into a list of typed values.

    Args:
        string (str): The input string to parse.
        type (type, optional): The type to convert elements to. Defaults to int.

    Returns:
        list: List of converted values.
    """
    if "," in string:
        return [type(s) for s in string.split(",")]
    return [type(s) for s in string.split(" ")]


@task
def make(ctx, jobs="auto", touch=False, clean=False, binary=""):
    """
    Recompile the Abinit source code.

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
    if touch:
        with cd(ABINIT_ROOTDIR):
            cmd = "./abisrc.py touch"
            cprint(f"Executing: {cmd}", color="yellow")
            result = ctx.run(cmd, pty=True)
            if not result.ok:
                cprint(f"{cmd=} failed. Aborting now!", color="red")
                return 1

    top = find_top_build_tree(".", with_abinit=False)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)

    with cd(top):
        if clean:
            ctx.run("cd src && make clean && cd ..", pty=True)
            ctx.run("cd shared && make clean && cd ..", pty=True)

        #cmd = f"make -j{jobs} {binary} | tee make.log 2> make.stderr"
        cmd = f"make -j{jobs} {binary}"
        cprint(f"Executing: {cmd}", color="yellow")
        result = ctx.run(cmd, pty=True)
        if not result.ok:
            cprint(f"{cmd=} failed. Aborting now!", color="red")
            sys.exit(1)

        # TODO Check for errors in make.stderr
        #cprint("Exit code: %s" % retcode, "green" if retcode == 0 else "red")

        #if SYSTEM == "Darwin":
        #    for binary in ALL_BINARIES:
        #        cmd = f"codesign -v --force --deep src/98_main/{binary}"
        #        cprint("Executing: %s" % cmd, "yellow")
        #        ctx.run(cmd, pty=True)

@task
def clean(ctx):
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
def runemall(ctx, make=True, jobs="auto", touch=False, clean=False, keywords=None):
    """
    Run all sequential and parallel tests.

    The task first ensures the code is compiled, and then executes the
    test suite. It exits immediately if any critical error occurs.

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
    make(ctx, jobs=jobs, touch=touch, clean=clean)

    top = find_top_build_tree(".", with_abinit=True)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)
    kws = "" if keywords is None else "-k %s" % keywords

    with cd(os.path.join(top, "tests")):
        cmd = "./runtests.py -j%d %s" % (jobs, kws)
        cprint(f"Executing: {cmd}", color="yellow")
        ctx.run(cmd, pty=True)
        # Now run the parallel tests.
        for n in [2, 4, 10]:
            j = jobs // n
            if j == 0: continue
            cmd = "./runtests.py paral mpiio -j%d -n%d %s" % (j, n, kws)
            cprint(f"Executing: {cmd}", color="yellow")
            ctx.run(cmd, pty=True)


@task
def makemake(ctx):
    """
    Invoke the `makemake` script to rebuild the build system.

    Args:
        ctx: Invoke context.
    """
    with cd(ABINIT_ROOTDIR):
        ctx.run("./config/scripts/makemake", pty=True)


@task
def makedeep(ctx, jobs="auto"):
    """
    Perform a complete rebuild cycle: makemake, clean, and build.

    Args:
        ctx: Invoke context.
        jobs (str or int, optional): Parallel threads for compilation. Defaults to "auto".
    """
    makemake(ctx)
    make(ctx, jobs=jobs, clean=True)


@task
def abichecks(ctx):
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
        exclude = ["check-libpaw.py", "warningschk.py", "abirules_tools.py", "__init__.py"]
        for py_script in [f for f in os.listdir(script_dir) if f.endswith(".py")]:
            if py_script in exclude: continue
            py_script = os.path.join(script_dir, py_script)
            print("Running", py_script, "... ")
            start = time.time()
            result = ctx.run(py_script, warn=True, pty=True)
            #print(result.ok)
            msg, color = ("[OK]", "green") if result.ok else ("[FAILED]", "red")
            cprint("%s (%.2f s)" % (msg, time.time() - start), color=color)
            if not result.ok: retcode += 1

    if retcode != 0:
        cprint("%d FAILED TESTS" % retcode, color="red")
    else:
        cprint("ALL TESTS OK", color="green")

    return retcode


@task
def robodoc(ctx):
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
            print("Trying to open %s in browser ..." % html_path)
            return webbrowser.open_new_tab(html_path)
        cprint("ROBODOC BUILD FAILED", color="red")

        return result.ok


@task
def mksite(ctx):
    """
    Build the Abinit documentation site and serve it locally.

    Args:
        ctx: Invoke context.
    """
    with cd(ABINIT_ROOTDIR):
        webbrowser.open_new_tab("http://127.0.0.1:8000")
        ctx.run("./mksite.py serve --dirtyreload", pty=True)


@task
def links(ctx):
    """
    Create symbolic links to all Abinit executables in the current directory.

    Args:
        ctx: Invoke context.
    """
    top = find_top_build_tree(".", with_abinit=True)
    main98 = os.path.join(top, "src", "98_main")
    for dest in ALL_BINARIES:
        if os.path.islink(os.path.join(os.getcwd(), dest)): continue
        source = os.path.join(main98, dest)
        if os.path.isfile(source):
            os.symlink(source, dest)
        else:
            cprint("Cannot find `%s` in dir `%s" % (source, main98), color="yellow")


@task
def ctags(ctx):
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
def fgrep(ctx, pattern):
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
        cmd  = 'grep -r -i --color --include "*[.F90,.f90,.finc,.c,.cu,.cpp,.h]" "%s" src shared' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def cgrep(ctx, pattern):
    """
    Case-insensitive search for a pattern specifically in C files.

    Args:
        ctx: Invoke context.
        pattern (str): The pattern to search for.
    """
    with cd(ABINIT_ROOTDIR):
        cmd  = 'grep -r -i --color --include "*.c" "%s" src shared' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def tgrep(ctx, pattern):
    """
    Search for a pattern within all test input files.

    Args:
        ctx: Invoke context.
        pattern (str): The pattern to search for.
    """
    with cd(ABINIT_ROOTDIR):
        cmd  = 'grep -r -i --color "%s" tests/*/Input/*' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def vimt(ctx, tagname):
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
def env(ctx):
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
def diff2(ctx, filename="run.abo"):
    """
    Compare the specified output file with the most recent backup file.

    Args:
        ctx: Invoke context.
        filename (str, optional): The base output filename. Defaults to "run.abo".
    """
    vimdiff = which_differ()
    files = sorted([f for f in os.listdir(".") if f.startswith(filename)])
    if not files: return
    cmd = f"{vimdiff} {filename} {files[-1]}"
    cprint(f"Executing {cmd}", color="green")
    ctx.run(cmd, pty=True)


@task
def diff3(ctx, filename="run.abo"):
    """
    Compare the current output file with the two most recent backup files.

    Args:
        ctx: Invoke context.
        filename (str, optional): The base output filename. Defaults to "run.abo".
    """
    differ = which_differ()

    files = sorted([f for f in os.listdir(".") if f.startswith(filename)])
    if not files: return

    if len(files) > 2:
        cmd = "%s %s %s %s" % (differ, filename, files[-2], files[-1])
    else:
        cmd = "%s %s %s" % (differ, filename, files[-1])
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def add_trunk(ctx):
    """
    Add the main Abinit GitLab repository as a git remote named "trunk".

    Args:
        ctx: Invoke context.
    """
    cmd = "git remote add trunk git@gitlab.abinit.org:trunk/abinit.git"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)
    cmd = "git fetch trunk"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def remote_add(ctx, remote):
    """
    Register a developer's fork as a git remote and fetch it.

    Args:
        ctx: Invoke context.
        remote (str): GitLab username of the developer.
    """
    cmd = f"git remote add {remote} git@gitlab.abinit.org:{remote}/abinit.git"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)
    cmd = f"git fetch {remote}"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def gdb(ctx, input_name, exec_name="abinit", run_make=False):
    """
    Launch the GDB debugger for a specific executable and input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ABINIT input file.
        exec_name (str, optional): Name of the executable. Defaults to "abinit".
        run_make (bool, optional): If True, build before debugging. Defaults to False.
    """
    if run_make: make(ctx)

    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main", exec_name)
    cprint(f"Using binpath: {binpath}", "green")
    cmd = f"gdb {binpath} --one-line 'settings set target.run-args {input_name}'"
    cprint(f"Executing gdb command: {cmd}", color="green")
    # mpirun -np 2 xterm -e gdb fftprof --command=dbg_file
    #cprint("Type run to start lldb debugger", color="green")
    #cprint("Then use `bt` to get the backtrace\n\n", color="green")
    ctx.run(cmd, pty=True)


@task
def lldb(ctx, input_name, exec_name="abinit", run_make=False):
    """
    Launch the LLDB debugger for a specific executable and input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ABINIT input file.
        exec_name (str, optional): Name of the executable. Defaults to "abinit".
        run_make (bool, optional): If True, build before debugging. Defaults to False.
    """
    if run_make: make(ctx)

    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main", exec_name)
    cprint(f"Using binpath: {binpath}", color="green")
    cmd = f"lldb {binpath} --one-line 'settings set target.run-args {input_name}'"
    cprint(f"Executing lldb command: {cmd}", color="green")
    cprint("Type run to start lldb debugger", color="green")
    cprint("Then use `bt` to get the backtrace\n\n", color="green")
    ctx.run(cmd, pty=True)


@task
def mpi_check(ctx, np_list="1, 2", abinit_input_file="run.abi", mpi_runner="mpiexec", run_make=False):
    """
    Run an ABINIT input file with various MPI process counts and compare results.

    Args:
        ctx: Invoke context.
        np_list (str): List of process counts (e.g., "1, 2, 4"). Defaults to "1, 2".
        abinit_input_file (str, optional): Path to the input file. Defaults to "run.abi".
        mpi_runner (str, optional): Command used to launch MPI. Defaults to "mpiexec".
        run_make (bool, optional): If True, build before running. Defaults to False.
    """
    if run_make: make(ctx)

    cprint(f"Will run {abinit_input_file=} with MPI nprocs in {np_list=}", color="yellow")

    differ = which_differ()
    np_list = list_from_string(np_list)

    for np in np_list:
        change_output_file(abinit_input_file, f"run_mpi{np}.abo")
        cmd = f"{mpi_runner} -n {np} abinit {abinit_input_file} | tee run_mpi{np}.log"
        cprint(f"About to execute {cmd=}", color="yellow" )
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
def omp_check(ctx, omp_threads="1, 2", np=1, abinit_input_file="run.abi", mpi_runner="mpiexec", run_make=False):
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
    if run_make: make(ctx)

    differ = which_differ()

    omp_threads = list_from_string(omp_threads)

    cprint("Will run {abinit_input_file=} with OMP threads={omp_threads=} and MPI nprocs={np}", color="yellow")
    for nth in omp_threads:
        change_output_file(abinit_input_file, f"run_omp{nth}_mpi{np}.abo")
        cmd = f"{mpi_runner} -n {np} abinit -o {nth} {abinit_input_file} | tee run_omp{nth}_mpi{np}.log"
        cprint(f"About to execute {cmd=}", color="yellow" )
        ctx.run(cmd)

    omp_ref = omp_threads[0]
    for nth in omp_threads[1:]:
        cmd = f"{differ} run_omp{omp_ref}.abo run_omp{nth}_mpi{np}.abo"
        cprint(f"About to execute {cmd=}", color="yellow" )
        ctx.run(cmd, pty=True)
        cmd = f"{differ} run_omp{omp_ref}.log run_omp{nth}_mpi{np}.log"
        cprint(f"About to execute {cmd=}", color="yellow" )
        ctx.run(cmd, pty=True)

@task
def pyenv_clean(ctx):
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
    cprint("About to execute {cmd=}")
    ctx.run(cmd)


@task
def abinit(ctx, input_name, run_make=False):
    """
    Run the `abinit` executable with the specified input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ABINIT input file.
        run_make (bool, optional): If True, build before running. Defaults to False.
    """
    _run(ctx, input_name, exec_name="abinit", run_make=run_make)


@task
def anaddb(ctx, input_name, run_make=False):
    """
    Run the `anaddb` executable with the specified input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the ANADDB input file.
        run_make (bool, optional): If True, build before running. Defaults to False.
    """
    _run(ctx, input_name, exec_name="anaddb", run_make=run_make)


def _run(ctx, input_name, exec_name, run_make):
    """
    Internal helper to execute an Abinit binary with an input file.

    Args:
        ctx: Invoke context.
        input_name (str): Path to the input file.
        exec_name (str): Name of the executable.
        run_make (bool): If True, build before running.
    """
    if run_make: make(ctx)
    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main", exec_name)
    cprint(f"Using binpath: {binpath}", color="green")
    cmd = f"{binpath} {input_name}"
    cprint(f"Executing {cmd}", color="green")
    ctx.run(cmd, pty=True)


@task
def pull_trunk(ctx):
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
def pull(ctx):
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
def push(ctx):
    """
    Commit and push current changes including tags to the origin repository.

    Args:
        ctx: Invoke context.
    """
    ctx.run("git commit")
    ctx.run("git push")
    ctx.run("git push --tags")


@task
def submodules(ctx):
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
def branchoff(ctx, start_point):
    """
    Create a new local branch starting from a remote branch (e.g., trunk/develop).

    Automatically sets the upstream to origin for the new branch.

    Args:
        ctx: Invoke context.
        start_point (str): The remote branch to start from (e.g., "trunk/release-9.0").
    """
    try:
        remote, branch = start_point.split("/")
    except:
        remote = "trunk"

    def run(cmd):
        cprint(f"Executing: `{cmd}`", color="green")
        ctx.run(cmd)

    run(f"git fetch {remote}")
    # Create new branch `test_v9.0` using trunk/release-9.0 as start_point:
    # git checkout [-q] [-f] [-m] [[-b|-B|--orphan] <new_branch>] [<start_point>]
    my_branch = "my_" + start_point.replace("/", "-")
    run(f"git checkout -b {my_branch} {start_point}")
    # Change default upstream. If you forget this step, you will be pushing to trunk
    run("git branch --set-upstream-to origin")
    run("git push origin HEAD")


@task
def dryrun_merge(ctx, start_point):
    """
    Perform a dry-run merge of a remote branch into the current branch.

    Args:
        ctx: Invoke context.
        start_point (str): The remote branch to merge from.
    """

    def run(cmd):
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
def watchdog(ctx, jobs="auto", sleep_time=5):
    """
    Monitor the source directory for changes and trigger recompilation automatically.

    Args:
        ctx: Invoke context.
        jobs (str or int, optional): Parallel threads for make. Defaults to "auto".
        sleep_time (int, optional): Sleep time in seconds between checks. Defaults to 5.
    """
    cprint("Starting watchdog service to watch F90 files and execute `make` when changes are detected", color="green")
    cprint("Enter <CTRL + C> in the terminal to kill the service.", color="green")

    cprint(f"Start watching F90 files with sleep_time {sleep_time} s ....", color="green")
    top = find_top_build_tree(".", with_abinit=True)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)

    # http://thepythoncorner.com/dev/how-to-create-a-watchdog-in-python-to-look-for-filesystem-changes/
    # https://stackoverflow.com/questions/19991033/generating-multiple-observers-with-python-watchdog
    import time

    from watchdog.events import PatternMatchingEventHandler
    from watchdog.observers import Observer
    event_handler = PatternMatchingEventHandler(patterns="*.F90", ignore_patterns="",
                                                ignore_directories=False, case_sensitive=True)

    def on_created(event):
        print(f"hey, {event.src_path} has been created!")

    def on_deleted(event):
        print(f"what the f**k! Someone deleted {event.src_path}!")

    def on_modified(event):
        print(f"hey buddy, {event.src_path} has been modified")
        cmd = "make -j%d  > >(tee -a make.log) 2> >(tee -a make.stderr >&2)" % jobs
        cprint("Executing: %s" % cmd, color="yellow")
        with cd(top):
            try:
                result = ctx.run(cmd, pty=True)
                if result.ok:
                    cprint("Make completed successfully", color="green")
                    cprint("Watching for changes ...", color="green")
            except Exception:
                cprint("Make returned non-zero exit status", color="red")
                cprint("Keep on watching for changes hoping you get it right ...", color="red")

    def on_moved(event):
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
        return subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).strip().decode("utf-8")
    except subprocess.CalledProcessError:
        raise RuntimeError("Not inside a git repository or an error occurred")


def get_git_tags() -> list[str]:
    """
    Retrieve a list of all git tags in the repository.

    Returns:
        list[str]: List of tag names.

    Raises:
        RuntimeError: If not in a git repository or git command fails.
    """
    try:
        tags = subprocess.check_output(["git", "tag"]).decode("utf-8").split("\n")
        # Remove empty strings from the list
        return [tag for tag in tags if tag]
    except subprocess.CalledProcessError:
        raise RuntimeError("Not inside a git repository or an error occurred")


@task
def official_release(ctx, new_version, dry_run=True):
    """
    Automate the process of creating a new official Abinit release.

    This involves merging develop into master, tagging, and pushing to
    remote repositories.

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
        return ctx.run(command, **_run_kwargs)

    current_branch = get_current_branch()
    if current_branch != "develop":
        raise RuntimeError(f"You are on the '{current_branch}' branch, not 'develop'.")

    old_tags = get_git_tags()
    if new_version in old_tags and not dry_run:
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

        if not dry_run:
            _run(f"git commit -a -m 'v{version}'")
            _run(f"git tag -a {version} -m 'v{version}'")
            _run("git push origin master")

        # Step 2: Push to GitHub
        _run(f"git remote add abinit {github_url} || echo 'Remote already exists but this is not critical'")
        if not dry_run:
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
def git_info(ctx, top_n=20):
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
            check=True
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
            ["git", "cat-file", "--batch-check=%(objectname) %(objecttype) %(objectsize)"],
            input=input_text,
            stdout=subprocess.PIPE,
            text=True,
            check=True
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
        reverse=True
    )

    print(f"\nTop {top_n} largest files ever committed:")
    for size, path, obj_hash in sorted_blobs[:top_n]:
        print(f"{size / (1024*1024):7.2f} MB\t{path}")

    ctx.run("git count-objects -vH", pty=True)


@task
def large_files(ctx, top_dir=None, size_threshold_mb=5):
    """
    Find and list files larger than `size_threshold_mb` megabytes under `top_dir`.

    Args:
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
def system(ctx):
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
    info.append(["Memory (Total)", f"{psutil.virtual_memory().total / (1024 ** 3):.2f} GB"])
    for level, size in get_cache_info().items():
        info.append([level, size])

    print(tabulate(info, headers=["Item", "Value"], tablefmt="grid"))


@task
def pid(ctx, pid):
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
        info.append(["Memory RSS", f"{p.memory_info().rss / (1024 ** 2):.2f} MB"])
        info.append(["Threads", p.num_threads()])
        info.append(["CWD", p.cwd()])
        info.append(["Parent PID", p.ppid()])
        info.append(["Start Time (Epoch)", int(p.create_time())])
        print(tabulate(info, headers=["Item", "Value"], tablefmt="grid"))

    except psutil.NoSuchProcess:
        print(f"❌ Process with PID {pid} does not exist.")
        sys.exit(1)


def get_cache_info() -> dict:
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


def get_cache_info_linux() -> dict:
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


def get_cache_info_mac() -> dict:
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
        "hw.l3cachesize": "L3"
    }
    for key, label in keys.items():
        try:
            out = subprocess.check_output(["sysctl", "-n", key]).decode().strip()
            caches[label] = f"{int(out) // 1024} KB"
        except Exception:
            continue
    return caches


def get_cache_info_windows() -> dict:
    """
    Retrieve CPU cache information on Windows using WMIC.

    Returns:
        dict: Mapping of cache levels to their sizes.
    """
    caches = {}
    try:
        out = subprocess.check_output(["wmic", "cpu", "get", "L2CacheSize,L3CacheSize"],
                                       stderr=subprocess.DEVNULL).decode()
        lines = out.strip().split("\n")
        if len(lines) >= 2:
            _, l2, l3 = lines[1].split()
            if l2: caches["L2"] = f"{l2} KB"
            if l3: caches["L3"] = f"{l3} KB"
    except Exception:
        pass
    return caches


def _extract_errors(logfile, context_lines=5) -> list[str]:
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
        r"error",                # generic errors
        r"fail",                 # tests failing
        r"cannot\s+find",        # missing library or header
        r"no\s+such\s+file",     # missing file
        r"not\s+found",          # program not found
        r"undefined\s+reference" # linking errors
    ]

    import re
    regex = re.compile("|".join(ERROR_PATTERNS), re.IGNORECASE)
    n = len(lines)
    errors = []
    for i, line in enumerate(lines):
        if regex.search(line):
            # Ignore maches such as `sd_yakl_options='optional fail'
            if line.startswith("sd_") and "fail" in line: continue
            # Capture context
            start = max(0, i - context_lines)
            end = min(n, i + context_lines + 1)
            context = "".join(lines[start:end])
            errors.append(context.strip())

    return errors


def find_filename(filename, start_dir=None) -> Path:
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
    current = start_dir.resolve()
    while True:
        candidate = current / filename
        if candidate.exists():
            return candidate
        if current.parent == current:  # reached filesystem root
            raise FileNotFoundError(f"File '{filename}' not found.")
        current = current.parent


@task
def config_log(ctx, log_path="config.log"):
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
