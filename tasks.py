"""
Pyinvoke file for automating build/config stuff.

Example:

    invoke abichecks
    invoke --list

Can be executed everywhere inside the Abinit directory, including build directories.
"""
from __future__ import annotations

import os
import sys
import webbrowser
import subprocess

from contextlib import contextmanager
try:
    from invoke import task
except ImportError:
    raise ImportError("Cannot import invoke package. Use `pip install invoke`")

from tests.pymods.testsuite import find_top_build_tree
from tests.pymods.devtools import number_of_cpus
from tests.pymods.termcolor import cprint

ABINIT_ROOTDIR = os.path.dirname(__file__)
ABINIT_SRCDIR = os.path.join(ABINIT_ROOTDIR, "src")

# Set ABI_PSPDIR env variable to point to the absolute path of Pspdir
#os.environ["ABI_PSPDIR"] = os.path.abspath(os.path.join(ABINIT_ROOTDIR, "Pspdir"))
#print("ABI_PSPDIR:", os.environ["ABI_PSPDIR"])

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

import platform
SYSTEM = platform.system()


def which_vim() -> str:
    """
    Find vim in $PATH
    """
    if which("mvim") is not None: return "mvim"
    if which("nvim") is not None: return "nvim"
    if which("vim") is not None: return "vim"
    raise RuntimeError("Cannot find vim in $PATH!")


def which_differ() -> str:
    """
    Find differ in $PATH
    """
    differ = "vimdiff"
    if which("mvimdiff") is not None: differ = "mvimdiff"

    if which(differ) is None:
        raise RuntimeError(f"Cannot find {differ=} in $PATH!")

    return differ


def change_output_file(input_file, output_file):
    """
    Change the name of the main output file in the `input_file` using `output_file`
    """
    with open(input_file, "rt") as fh:
        remove_iline = None
        lines = [l.lstrip() for l in fh.readlines()]
        for i, l in enumerate(lines):
            if l.lstrip().startswith("output_file"):
                remove_iline = i
                break

    if remove_iline is not None:
        lines.pop(remove_iline)

    lines.insert(0, f'output_file = "{output_file}"')
    with open(input_file, "wt") as fh:
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
    if "," in string:
        return [type(s) for s in string.split(",")]
    return [type(s) for s in string.split(" ")]


@task
def make(ctx, jobs="auto", touch=False, clean=False, binary=""):
    """
    Touch all modified files and recompile the code

    Args:
        jobs: Use `jobs` threads for make -jNUM
        touch: Touch all changed files
        clean: Issue `make clean` before `make`.
        binary: Binary to recompile, default: all
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
    """Remove object files in src and shared. Do not remove object files in fallbacks"""
    top = find_top_build_tree(".", with_abinit=False)
    with cd(top):
        ctx.run("cd src && make clean && cd ..", pty=True)
        ctx.run("cd shared && make clean && cd ..", pty=True)


@task
def runemall(ctx, make=True, jobs="auto", touch=False, clean=False, keywords=None):
    """Run all tests (sequential and parallel). Exit immediately if errors"""
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
def makemake(ctx, without_chmod=True):
    """Invoke makemake"""
    with cd(ABINIT_ROOTDIR):
        opt = "--without-chmod" if without_chmod else ""
        ctx.run(f"./config/scripts/makemake {opt}", pty=True)


@task
def makedeep(ctx, jobs="auto"):
    """Execute `makemake && make clean && make`"""
    makemake(ctx)
    make(ctx, jobs=jobs, clean=True)


@task
def abichecks(ctx):
    """Execute (some of the) abichecks scripts."""
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
    """Build robodoc documentation."""
    with cd(ABINIT_ROOTDIR):
        result = ctx.run("./mkrobodoc.sh", pty=True)

        if result.ok:
            cprint("ROBODOC BUILD OK", color="green")
            # https://stackoverflow.com/questions/44447469/cannot-open-an-html-file-from-python-in-a-web-browser-notepad-opens-instead
            html_path = os.path.join(ABINIT_ROOTDIR, "./tmp-robodoc/www/robodoc/masterindex.html")
            print("Trying to open %s in browser ..." % html_path)
            return webbrowser.open_new_tab(html_path)
        else:
            cprint("ROBODOC BUILD FAILED", color="red")

        return result.ok


@task
def mksite(ctx):
    """
    Build the Abinit documentation by running the mksite.py script and open the main page in the browser.
    """
    with cd(ABINIT_ROOTDIR):
        webbrowser.open_new_tab("http://127.0.0.1:8000")
        ctx.run("./mksite.py serve --dirtyreload", pty=True)


@task
def links(ctx):
    """
    Create symbolic links to Abinit executables in current working directory.
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
    Update ctags file.
    """
    with cd(ABINIT_ROOTDIR):
        cmd = "ctags -R --langmap=fortran:+.finc.f90.F90,c:.c.cpp shared/ src/"
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)
        #ctx.run('ctags -R --exclude="_*"', pty=True)

@task
def fgrep(ctx, pattern):
    """
    Grep for `pattern` in all F90 files contained in `src` and `shared` directories.
    """
    # grep -r -i --include \*.h
    # Syntax notes:
    #    -r - search recursively
    #    -i - case-insensitive search
    #    --include=\*.${file_extension} - search files that match the extension(s) or file pattern only
    with cd(ABINIT_ROOTDIR):
        cmd  = 'grep -r -i --color --include "*[.F90,.f90,.finc,.c,.cu,.cpp]" "%s" src shared' % pattern
        #cmd  = 'grep -r -i --color --include "*.F90" "%s" src shared' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def cgrep(ctx, pattern):
    """
    Grep for `pattern` in all C files contained in `src` and `shared` directories.
    """
    with cd(ABINIT_ROOTDIR):
        cmd  = 'grep -r -i --color --include "*.c" "%s" src shared' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def tgrep(ctx, pattern):
    """
    Grep for `pattern` in all input files contained in the `tests` directory.
    """
    with cd(ABINIT_ROOTDIR):
        cmd  = 'grep -r -i --color "%s" tests/*/Input/*' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def vimt(ctx, tagname):
    """
    Execute `vim -t tagname` with tagname as ctags tag.
    """
    vim = which_vim()
    with cd(ABINIT_ROOTDIR):
        cmd = f"{vim} -f {tagname}"
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def env(ctx):
    """Print sh code to set $PATH and $ABI_PSPDIR in order to work with build directory."""
    cprint("\nExecute the following lines in the shell to set the env:\n", color="green")
    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main")
    print(f"export ABI_PSPDIR={ABINIT_ROOTDIR}/tests/Pspdir")
    print(f"export PATH={binpath}:$PATH")


@task
def diff2(ctx, filename="run.abo"):
    """
    Execute `vimdiff` to compare run.abo with the last run.abo0001 found in the cwd.
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
    Execute `vimdiff` to compare run.abo with the last run.abo0001 found in the cwd.
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
    """Register trunk as remote."""
    cmd = "git remote add trunk git@gitlab.abinit.org:trunk/abinit.git"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def remote_add(ctx, remote):
    """Register `remote` as remote branch and fetch it"""
    cmd = f"git remote add {remote} git@gitlab.abinit.org:{remote}/abinit.git"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)
    cmd = f"git fetch {remote}"
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)


@task
def gdb(ctx, input_name, exec_name="abinit", run_make=False):
    """
    Execute `gdb` debugger with the given `input_name`.
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
    Execute `lldb` debugger with the given `input_name`.
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
    Args:
        np_list: List of MPI procs
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
    if which("conda") is not None:
        cmd = f"conda clean --all --yes"
        cprint(f"About to execute {cmd=}")
        ctx.run(cmd)

    cmd = f"pip cache purge"
    cprint("About to execute {cmd=}")
    ctx.run(cmd)

@task
def abinit(ctx, input_name, run_make=False):
    """
    Execute `abinit` with the given `input_name`.
    """
    _run(ctx, input_name, exec_name="abinit", run_make=run_make)


@task
def anaddb(ctx, input_name, run_make=False):
    """"execute `anaddb` with the given `input_name`."""
    _run(ctx, input_name, exec_name="anaddb", run_make=run_make)


def _run(ctx, input_name, exec_name, run_make):
    """"Execute `exec_name input_name`"""
    if run_make: make(ctx)
    top = find_top_build_tree(".", with_abinit=True)
    binpath = os.path.join(top, "src", "98_main", exec_name)
    cprint(f"Using binpath: {binpath}", color="green")
    cmd = f"{binpath} {input_name}"
    cprint(f"Executing {cmd}", color="green")
    ctx.run(cmd, pty=True)


@task
def pull_trunk(ctx):
    """"Execute `git stash && git pull trunk develop && git stash apply`"""
    ctx.run("git stash")
    ctx.run("git pull trunk develop")
    ctx.run("git pull trunk develop --tags")
    ctx.run("git commit")
    ctx.run("git push")
    ctx.run("git push --tags")
    ctx.run("git stash apply")

@task
def pull(ctx):
    """"Execute `git stash && git pull --recurse-submodules && git stash apply && makemake`"""
    ctx.run("git stash")
    ctx.run("git pull --recurse-submodules")
    ctx.run("git stash apply")
    makemake(ctx)


@task
def push(ctx):
    """"Execute `git commit && git push && git push --tags`"""
    ctx.run("git commit")
    ctx.run("git push")
    ctx.run("git push --tags")


@task
def submodules(ctx):
    """Update submodules."""
    with cd(ABINIT_ROOTDIR):
        # https://stackoverflow.com/questions/1030169/easy-way-to-pull-latest-of-all-git-submodules
        ctx.run("git submodule update --remote --init", pty=True)
        ctx.run("git submodule update --recursive --remote", pty=True)

@task
def branchoff(ctx, start_point):
    """"Checkout new branch from start_point e.g. `trunk/release-9.0` and set default upstream to origin."""
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
    """"Merge `remote/branch` in dry-run mode."""

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
    Start watchdog service to watch F90 files and execute `make` when changes are detected.
    """
    cprint("Starting watchdog service to watch F90 files and execute `make` when changes are detected", color="green")
    cprint("Enter <CTRL + C> in the terminal to kill the service.", color="green")

    cprint(f"Start watching F90 files with sleep_time {sleep_time} s ....", color="green")
    top = find_top_build_tree(".", with_abinit=True)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)

    # http://thepythoncorner.com/dev/how-to-create-a-watchdog-in-python-to-look-for-filesystem-changes/
    # https://stackoverflow.com/questions/19991033/generating-multiple-observers-with-python-watchdog
    import time
    from watchdog.observers import Observer
    from watchdog.events import PatternMatchingEventHandler
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
                cprint(f"Make returned non-zero exit status", color="red")
                cprint(f"Keep on watching for changes hoping you get it right ...", color="red")

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


def which(cmd):
    """
    Returns full path to a executable.

    Args:
        cmd (str): Executable command to search for.

    Returns:
        (str) Full path to command. None if it is not found.

    Example::

        full_path_to_python = which("python")

    Take from monty.path. See https://github.com/materialsvirtuallab/monty
    """

    def is_exe(fp):
        return os.path.isfile(fp) and os.access(fp, os.X_OK)

    fpath, fname = os.path.split(cmd)
    if fpath:
        if is_exe(cmd):
            return cmd
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, cmd)
            if is_exe(exe_file):
                return exe_file
    return None


def get_current_branch() -> str:
    """Run git command to get the current branch"""
    try:
        return subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).strip().decode("utf-8")
    except subprocess.CalledProcessError:
        raise RuntimeError("Not inside a git repository or an error occurred")


def get_git_tags() -> list[str]:
    """Run git command to list tags"""
    try:
        tags = subprocess.check_output(["git", "tag"]).decode("utf-8").split("\n")
        # Remove empty strings from the list
        return [tag for tag in tags if tag]
    except subprocess.CalledProcessError:
        raise RuntimeError("Not inside a git repository or an error occurred")


@task
def official_release(ctx: Context, new_version: str, dry_run: bool = True) -> None:
    """
    Build new officiale release ...

    Example usage:

        invoke official-release 10.2.4
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
        'config/gnu/compile',
        'config/gnu/config.guess',
        'config/gnu/config.sub',
        'config/gnu/install-sh',
        'config/gnu/missing',
        'config/gnu/depcomp',
    ]

    with cd(ABINIT_ROOTDIR):
        # The version in .current_version is updated manually.
        # Here we check that the value stored in the file is equal to the command line argument.
        with open(".current_version", "rt") as fh:
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
def git_info(ctx: Context, top_n=20) -> None:
    """Scan git history for largest top_n files"""

    def get_git_objects():
        """Return list of all Git objects (hash, path)."""
        result = subprocess.run(
            ['git', 'rev-list', '--objects', '--all'],
            stdout=subprocess.PIPE,
            text=True,
            check=True
        )
        objects = []
        for line in result.stdout.splitlines():
            parts = line.split(' ', 1)
            if len(parts) == 2:
                objects.append((parts[0], parts[1]))
        return objects

    def get_blob_sizes(hashes):
        """Return a dict of {hash: (size_in_bytes, path)} for blobs."""
        input_text = '\n'.join(hashes)
        result = subprocess.run(
            ['git', 'cat-file', '--batch-check=%(objectname) %(objecttype) %(objectsize)'],
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
