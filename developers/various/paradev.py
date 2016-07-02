#!/usr/bin/env python
"""
This script runs the given abinit input file with different number of processors 
and compare the results
"""
from __future__ import print_function, division

import sys
import os
import argparse 
import warnings
import tempfile

from collections import namedtuple


def browse_localpath(filepath):
    """
    Open filepath in the default browser.

    .. warning: This code is not portable since we should pass a url.
    """
    import webbrowser
    if not filepath.startswith("file://"): filepath = "file://" + filepath

    try:
        webbrowser.open(filepath)
    except webbrowser.Error as exc:
        # Warn the user and ignore the exception.
        warnings.warn(str(exc))


def make_html_diff(fromfile, tofile):
    import difflib

    with open(fromfile) as fh: fromlines = fh.readlines()
    with open(tofile) as fh: tolines = fh.readlines()

    diff = difflib.HtmlDiff().make_file(fromlines, tolines, fromfile, tofile, context=True, numlines=3)

    # writelines because diff is a generator
    _, tmp_path = tempfile.mkstemp(text=True, suffix=".html")
    with open(tmp_path, "w") as fh: fh.writelines(diff)

    return tmp_path


def find_nextout(root="."):
    """
    Returns the path the next output file. Assumes Abinit conventions + gmatteo's convention 
    for the extension i.e. run.abo, run.aboA ...
    This function is the most complicated of the entire module!
    """
    from string import ascii_letters
    root = os.path.abspath(root)
                                                                                                               
    # Find all files with extensions `.abo?`
    chars, found = [], False
    for f in os.listdir(root):
        head, ext = os.path.splitext(f)
        if ext.startswith(".abo"):
            found = True
            if len(ext) > 4: chars.append(ext[-1])

    if not found:
        # First run
        fname = "run.abo"
    else:
        if not chars: 
            fname = "run.aboA"
        else:
            last_ch = sorted(chars)[-1]
            if last_ch == "Z": raise RunTimeError("Got Z as last character, clean your directory!")
            print("last",last_ch)
            i = ascii_letters.index(last_ch)
            next_ch =  ascii_letters[i+1]
            fname = "run.abo" + next_ch
                                                                                                               
    return os.path.join(root, fname)

def main():
    def str_examples():
        examples = (
          "\n"
          "Usage example:\n\n" 
          "paradev.py 1 2 3    ==> Run mpirun -n# abinit < files > log 2> err for n in [1,2,3] and compare the output files\n"
        )
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: 
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                     help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument("-b", "--binary", type=str, default="abinit", help="Executable to be launched")

    parser.add_argument("-f", "--file-type", type=str, default="out", help="File to compared (out from main output or log for log file)")

    #parser.add_argument("-s", "--strict", type=bool, default=True, action="store_if_true, 
    #                    help="Exit immediately if the subprocess returns nonzero exit status")

    parser.add_argument("nprocs", nargs="*", help="List of MPI nodes to be tested")

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    if options.nprocs: 
        options.nprocs = map(int, options.nprocs)
    else:
        options.nprocs = [1,2]

    class Run(namedtuple("Run", "out log err")):
        """Stores the paths to the stderr, the log file and the output file of the run."""

    run_list = []
    for np in options.nprocs:
        out = find_nextout()
        print("Expecting next output %s: " % os.path.basename(out))
        run = Run(out=out, log="log_np%d" % np, err="err_np%d" % np)
        run_list.append(run)

        files_file = "files"

        cmd = "mpirun -n %d %s < %s > %s 2> %s" % (np, options.binary, files_file, run.log, run.err)
        print("Executing: %s" % cmd)

        retcode = os.system(cmd)

        if retcode: 
            print("Return code: %d" % retcode)
            with open(run.log, "r") as fh: print(fh.read())
            with open(run.err, "r") as fh: print(fh.read())
            #if options.strict
            #    return retcode
            #else

    # Generate HTML diffs and open the (temporary) files in the default browser
    # TODO: Add call to fldiff!

    # This function will select the file we want to analyze.
    gat = {
        "log": lambda x: getattr(x, "log"),
        "out": lambda x: getattr(x, "out"),
    }[options.file_type]

    run0, html_files = run_list[0], []
    for run in run_list[1:]:
        hfile = make_html_diff(gat(run0), gat(run))
        html_files.append(hfile)

    for html_file in html_files:
        browse_localpath(html_file)

    return 0


if __name__ == "__main__":
    sys.exit(main())

