"""
Pyinvoke tasks.py file for automating build/config stuff.

Example:

    invoke abichecks

    invoke --list

Can be executed everywhere inside the Abinit directory, including build directories.
"""
import os
import sys
PY2 = sys.version_info[0] <= 2

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


#@task
#def vim(ctx, tag):
#    """Invoke """
#    with cd(ABINIT_SRCDIR):
#        cmd = "mvim -t %s" % tag
#        #cmd = "vim -t %s" % tag
#        ctx.run(cmd)


@task
def make(ctx, jobs="auto", touch=False, clean=False):
    """
    Touch all modified files and recompile the code with -jNUM.
    """
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
        retcode = ctx.run(cmd, pty=True)
        # TODO Check for errors in make.stderr
        #cprint("Exit code: %s" % retcode, "green" if retcode == 0 else "red")


@task
def clean(ctx):
    """Remove object files in src and shared."""
    top = find_top_build_tree(".", with_abinit=False)
    with cd(top):
        ctx.run("cd src && make clean && cd ..", pty=True)
        ctx.run("cd shared && make clean && cd ..", pty=True)


@task
def runemall(ctx, make=True, jobs="auto", touch=False, clean=False, keywords=None):
    """
    Run all tests (sequential and parallel).
    Exit immediately if errors
    """
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
    """makemake && make clean && make"""
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
            #if PY2:
            print("Running", py_script, "... ")
            #else:
            #    print("Running", py_script, "... ", end="")
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
        if os.path.isfile(dest): continue
        source = os.path.join(main98, dest)
        if os.path.isfile(source):
            os.symlink(source, dest)
        else:
            cprint("Cannot find `%s` in dir `%s" % (source, main98), "yellow")

#def pulltrunk(ctx):
#    ctx.run("git stash")
#    ctx.run("git pull trunk develop")
#    ctx.run("git stash apply")
#    ctx.run("git push")

#def move_to_master(ctx):
#    ctx.run("git tag -a v%s -m \"v%s release\"" % (NEW_VER, NEW_VER))
#    ctx.run("git push --tags")
#    ctx.run("git checkout master")
#    ctx.run("git pull")
#    ctx.run("git merge develop")
#    ctx.run("git push")
#    ctx.run("git checkout develop")
