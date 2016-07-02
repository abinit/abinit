#!/usr/bin/env python
"""
This script automates the steps needed to debug a MPI executable with the GNU debugger
Note that GNU gdb is not a parallel debugger hence this crude approach may not work properly.
"""
from __future__ import print_function, division

import sys
import os
import argparse 
import tempfile

from subprocess import Popen, PIPE

def str_examples():
    examples = (
      "\n"
      "Usage example:\n\n" 
         "pmdbg.py -n 2 abinit files_file  ==> (Try) to run Abinit in parallel under the control of gdb\n"
         "                                      Use 2 MPI processes, read input from files_file "
    )
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def main():
    """mpirun -n# xterm gdb binary -command=file"""
    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                     help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('-n', default=1, type=int, help='Number of MPI processes')

    parser.add_argument("args", nargs="+", help='executable stdin_file')

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    binary, stdin_fname= options.args[0], options.args[1]

    _, dbg_fname = tempfile.mkstemp()

    with open(dbg_fname, "w") as fh:
        fh.write("run < %s" % stdin_fname)

    # mpirun -n# xterm gdb binary -command=file.
    cmd = "mpirun -n %i xterm -e gdb %s -command=%s" % (options.n, binary, dbg_fname)
    print(cmd)

    p = Popen(cmd, shell=True) #, cwd=cwd, env=env)
    retcode = p.wait()

    try:
        os.remove(dbg_fname)
    except IOError:
        pass

    return retcode


if __name__ == "__main__":
    sys.exit(main())
