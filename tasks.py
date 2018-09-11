"""
Pyinvoke tasks.py file for automating build/config stuff.

Example:

    invoke abichecks

    invoke --list

Can be executed everywhere inside the Abinit directory, including build directories.
"""
import os
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
def make(ctx, jobs="auto", clean=False):
    """
    Touch all modified files and recompile the code with -jNUM.
    """
    with cd(ABINIT_SRCDIR):
        cmd = "./abisrc.py touch"
        cprint("Executing: %s" % cmd, "yellow")
        result = ctx.run(cmd, pty=True)
        if not result.ok:
            cprint("`%s` failed. Aborting now!" % cmd, "red")
            return 1

    top = find_top_build_tree(".", with_abinit=False)
    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)

    with cd(top):
        if clean: ctx.run("make clean", pty=True)
        #cmd = "make -j%d > make.log 2> make.stderr" % jobs
        cmd = "make -j%d  > >(tee -a make.log) 2> >(tee -a make.stderr >&2)" % jobs
        cprint("Executing: %s" % cmd, "yellow")
        ctx.run(cmd, pty=True)
        # TODO Check for errors in make.stderr


@task
def makemake(ctx):
    """Invoke makemake"""
    with cd(ABINIT_ROOTDIR):
        ctx.run("./config/scripts/makemake", pty=true)


@task
def makedeep(ctx, jobs="auto"):
    """makemake && make clean && make"""
    makemake(ctx)
    make(ctk, jobs=jobs, clean=True)


@task
def abilint(ctx):
    """Invoke abilint"""
    with cd(ABINIT_ROOTDIR):
        ctx.run("./config/scripts/abilint . .", pty=True)

@task
def abichecks(ctx):
    """Execute (some of the) abichecks scripts."""
    import time
    retcode = 0
    with cd(ABINIT_ROOTDIR):
        script_dir = os.path.join("abichecks", "scripts")
        exclude = ["check-input-vars.py", "check-libpaw.py", "warningschk.py"]
        for py_script in [f for f in os.listdir(script_dir) if f.endswith(".py")]:
            if py_script in exclude: continue
            py_script = os.path.join(script_dir, py_script)
            print("Running", py_script, "... ", end="")
            start = time.time()
            result = ctx.run(py_script, pty=True)
            msg, color = ("[OK]", "green") if result.ok else ("[FAILED]", "red")
            cprint("%s (%.2f s)" % (msg, time.time() - start), color=color)
            if not result.ok: retcode += 1

    if retcode != 0:
        cprint("%d FAILED TESTS", "red")
    else:
        cprint("ALL TESTS OK", "green")

    return retcode

#@task
#def robodoc(ctx):
#    with cd(ABINIT_ROOTDIR):
#        cmd = "./abisrc.py robodoc"
#        result = ctx.run(cmd, pty=True)
#
#    if result.ok
#        return 0
#        return webbrowser.open_new_tab("http://127.0.0.1:8000")
#    else:
#        print("robodoc build FAILED")

#@task
#def mksite(ctx):
#    """
#    Build the Abinit documentation by invoking mksite and open the main page in the browser.
#    """
#    with cd(ABINIT_ROOTDIR):
#        ctx.run("./mksite.py serve --dirtyreload", pty=True)
#        return webbrowser.open_new_tab("http://127.0.0.1:8000")


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
