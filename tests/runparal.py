#!/usr/bin/env python
"""
Script for running the parallel tests of the test suite with different number of MPI nodes.
Note: This script is just a small wrapper around runtests.py. It receives a list of paral
tests, register them and run all of them in parallel using distinct subprocesses. 
The user is responsible for providing a list of tests that do not cause an excessive overload of the OS.
"""
from __future__ import print_function, division, absolute_import #, unicode_literals

import sys
import os
import subprocess

script = os.path.join(os.path.dirname(__file__), "runtests.py")

MPI_NPROCS_LIST = [1, 2, 4, 10]

def str_examples():
    examples = """
      Usage example (assuming the script is executed within a build tree):
      \n
      runparal.py 1        => Run paral[1] with 1,2,4,10 MPI nodes
      runparal.py 1:3      => To select both test 1 and 2 
      runparal.py 1:3  4   => Same as above but add test number 4
    """
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    "Display the usage of the script."
    sys.stderr.write(str_examples())
    if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)

def main():
    try:
        args = sys.argv[1:]

        tests = []
        for arg in args:
            if ":" in arg:
                # String defining a python range
                tokens = map(int, arg.split(":"))
                tests.extend(list(range(*tokens)))
            else:
                # Assume integer
                num = int(arg)
                tests.append(num)

        # Build the list of strings to pass to runtest.py
        test_strings = []
        for t in tests:
            test_strings.append("paral[" + str(t) + "]")

    except Exception, exc:
        show_examples_and_exit(str(exc))

    # Lauch the tests in a subprocess (non blocking algorithm)
    processes = []
    for s in test_strings:
        for mpi_nprocs in MPI_NPROCS_LIST:
            workdir = os.path.join(os.getcwd(), s + "_MPI_" + str(mpi_nprocs))
            cmd_args = [script, s, "-n", str(mpi_nprocs), "-w", workdir]
            #print "About to execute", args
            p = subprocess.Popen(cmd_args)
            processes.append(p)

    # Parent waits here.
    for p in processes:
        p.wait()

    retcode = max([p.returncode for p in processes])

    print("Final return code %s" % retcode)
    return retcode 


if __name__ == "__main__":
    sys.exit(main())
