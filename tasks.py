"""
Pyinvoke file for automating build/config stuff.

Example:

    invoke abichecks
    invoke --list

Can be executed everywhere inside the Abinit directory, including build directories.

Use: `pip install invoke --user` to install invoke package.
"""
import os
import sys
import webbrowser

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
    "tdep",
    "testtransposer",
    "ujdet",
]


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


@task
def make(ctx, jobs="auto", touch=False, clean=False):
    """Touch all modified files and recompile the code with -jNUM."""
    if touch:
        with cd(ABINIT_ROOTDIR):
            cmd = "./abisrc.py touch"
            cprint("Executing: %s" % cmd, "yellow")
            result = ctx.run(cmd, pty=True)
            if not result.ok:
                cprint("`%s` failed. Aborting now!" % cmd, "red")
                return 1

    top = find_top_build_tree(".", with_abinit=False)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)

    with cd(top):
        if clean:
            ctx.run("cd src && make clean && cd ..", pty=True)
            ctx.run("cd shared && make clean && cd ..", pty=True)
        cmd = "make -j%d  > >(tee -a make.log) 2> >(tee -a make.stderr >&2)" % jobs
        cprint("Executing: %s" % cmd, "yellow")
        results = ctx.run(cmd, pty=True)
        # TODO Check for errors in make.stderr
        #cprint("Exit code: %s" % retcode, "green" if retcode == 0 else "red")


@task
def clean(ctx):
    """Remove object files in src and shared. Do not object files in fallbacks"""
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
        cprint("Executing: %s" % cmd, "yellow")
        ctx.run(cmd, pty=True)
        # Now run the parallel tests.
        for n in [2, 4, 10]:
            j = jobs // n
            if j == 0: continue
            cmd = "./runtests.py paral mpiio -j%d -n%d %s" % (j, n, kws)
            cprint("Executing: %s" % cmd, "yellow")
            ctx.run(cmd, pty=True)


@task
def makemake(ctx):
    """Invoke makemake"""
    with cd(ABINIT_ROOTDIR):
        ctx.run("./config/scripts/makemake", pty=True)


@task
def makedeep(ctx, jobs="auto"):
    """Execute `makemake && make clean && make`"""
    makemake(ctx)
    make(ctk, jobs=jobs, clean=True)


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
        cprint("%d FAILED TESTS" % retcode, "red")
    else:
        cprint("ALL TESTS OK", "green")

    return retcode


@task
def robodoc(ctx):
    with cd(ABINIT_ROOTDIR):
        result = ctx.run("./mkrobodoc.sh", pty=True)

        if result.ok:
            cprint("ROBODOC BUILD OK", "green")
            # https://stackoverflow.com/questions/44447469/cannot-open-an-html-file-from-python-in-a-web-browser-notepad-opens-instead
            html_path = os.path.join(ABINIT_ROOTDIR, "./tmp-robodoc/www/robodoc/masterindex.html")
            print("Trying to open %s in browser ..." % html_path)
            return webbrowser.open_new_tab(html_path)
        else:
            cprint("ROBODOC BUILD FAILED", "red")

        return result.ok


@task
def mksite(ctx):
    """
    Build the Abinit documentation by running the mksite.py script and open the main page in the browser.
    """
    with cd(ABINIT_ROOTDIR):
        ctx.run("./mksite.py serve --dirtyreload", pty=True)
        return webbrowser.open_new_tab("http://127.0.0.1:8000")


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
            cprint("Cannot find `%s` in dir `%s" % (source, main98), "yellow")


@task
def ctags(ctx):
    """
    Update ctags file.
    """
    with cd(ABINIT_ROOTDIR):
        cmd = "ctags -R shared/ src/"
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
        cmd  = 'grep -r -i --color --include "*.F90" "%s" src shared' % pattern
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
    Execute `vim -t tagname` with tagname a ctags tag.
    """
    with cd(ABINIT_ROOTDIR):
        if which("mvim") is not None:
            cmd  = "mvim -t %s" % tagname
        else:
            cmd  = "vim -t %s" % tagname
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def pull_trunk(ctx):
    """"Execute `git stash && git pull trunk develop && git stash apply`"""
    ctx.run("git stash")
    ctx.run("git pull trunk develop")
    ctx.run("git stash apply")


@task
def branchoff(ctx, start_point):
    """"Checkout new branch from start_point e.g. `trunk/release-9.0` and set default upstream to origin."""
    remote, branch = start_point.split("/")
    def run(cmd):
        cprint(f"Executing: `{cmd}`", "green")
        ctx.run(cmd)

    run(f"git fetch {remote}")
    # Create new branch `test_v9.0` using trunk/release-9.0 as start_point:
    # git checkout [-q] [-f] [-m] [[-b|-B|--orphan] <new_branch>] [<start_point>]
    my_branch = "my_" + branch
    run(f"git checkout -b {my_branch} {start_point}")
    # Change default upstream. If you forget this step, you will be pushing to trunk
    run("git branch --set-upstream-to origin")
    run("git push origin HEAD")


@task
def watchdog(ctx, jobs="auto", sleep_time = 5):
    """
    Start watchdog service to watch F90 files and execute `make` when changes are detected.
    """
    cprint("Starting watchdog service to watch F90 files and execute `make` when changes are detected", "green")
    cprint("Enter <CTRL + C> in the terminal to kill the service.", "green")

    cprint(f"Start watching F90 files with sleep_time {sleep_time} s ....", "green")
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
        cprint("Executing: %s" % cmd, "yellow")
        with cd(top):
            try:
                result = ctx.run(cmd, pty=True)
                if result.ok:
                    cprint("Make completed successfully", "green")
                    cprint("Watching for changes ...", "green")
            except Exception:
                cprint(f"Make returned non-zero exit status", "red")
                cprint(f"Keep on watching for changes hoping you get it right ...", "red")

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
