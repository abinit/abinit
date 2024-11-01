#!/usr/bin/env python
"""Search for inlined CPP macros in ABINIT src files"""
from __future__ import unicode_literals, division, print_function, absolute_import

__author__  = "M. Giantomassi"

import os
import re
import sys

from abirules_tools import find_src_dirs

# Files that will be checked.
re_srcfile = re.compile("\.([Ff]|[Ff]90|finc)$")

def is_srcfile(dirpath, fname):
    return re_srcfile.search(fname)

# List of macros that should not be inlined with if statements or other instructions.
# These macros are defined in shared/common/src/incs/abi_common.h
#
MACRO_NAMES = [
  "ABI_MALLOC",
  "ABI_MALLOC_SCALAR",
  "ABI_FREE_SCALAR",
  "ABI_FREE",
  "ABI_STAT_MALLOC",
  "ABI_MALLOC_OR_DIE",
  "ABI_MOVE_ALLOC",
  "ABI_MALLOC_TYPE_SCALAR",
  "ABI_CALLOC",
  "ABI_ICALLOC",
  "ABI_CALLOC_OR_DIE",
  "ABI_FREE_NOCOUNT",
  "ABI_MALLOC_IFNOT",
  "ABI_SFREE",
  "ABI_REMALLOC",
  "ABI_RECALLOC",
  "ABI_MALLOC_RAND",
]

# Regular expressions for each macro.
#regexps = dict()
#
#for name in MACRO_NAMES:
#    regexps[name] = re.compile(name + " ?\(.*\)")
#    #regexps[name] = re.compile(name + "\(.*\)") blanks between macro name and () are not permitted


def wrong_string(string):
    """Return 0 if input string does not contain inlined macros else 0"""
    s = string.strip().upper()
    if s.startswith("!"): return 0
    for macro_name in MACRO_NAMES:
        if macro_name not in s: continue
        #pattern = regexps[macro_name]
        if s.startswith("IF") and "THEN " not in s:
            return 1

    return 0


def test_wrong_string():

    strings_and_retcode = [
       ("if (.True.) ABI_MALLOC(dummy_atifc(3))", 1),
       ("!if (.True.) ABI_MALLOC(dummy_atifc(3))", 0),
       ("if (allocated(dummy_atifc)) ABI_FREE(dummy_atifc)", 1)
    ]

    failed = 0
    for (string, expected_retcode) in strings_and_retcode:
        retcode = wrong_string(string)
        if retcode != expected_retcode:
            failed += 1
            print("Failed test for string:", string)

    if failed:
        raise RuntimeError("test_wrong_string failed. Number of failed tests: %d" % failed)


def main():
    print("-------------------------------------------------------")
    print(" Searching for inlined CPP macros in ABINIT src files  ")
    print("-------------------------------------------------------")
    exit_status = 0
    for top in find_src_dirs():
        assert os.path.isdir(top)
        for dirpath, dirnames, files in os.walk(top):
            for src in files:
                if not is_srcfile(dirpath, src): continue
                fpath = os.path.join(dirpath,src)
                with open(fpath, "rt") as fh:
                    for lno, line in enumerate(fh):
                        retcode = wrong_string(line)
                        if retcode:
                            print("Found INLINED MACRO at %s:%d:  %s " % (src, lno + 1, line))
                            exit_status += 1

    if exit_status > 0:
        err_msg = """
Please, avoid instructions like:

    if (npts > 0) ABI_FREE(arr)

When the code is compiled in profiling mode, indeed,
ABI_FREE expands to the set of Fortran instructions:

    deallocate(arr)
    call memocc_abi()

These instructions MUST be placed inside an "if then" "end if" block else you get:

    if (npts > 0) deallocate(arr)
    call memocc_abi()

that is cleary wrong!

For the time being, one has to use the more verbose form:

  if (npts > 0) then
    ABI_FREE(arr)
  end if

To deallocate an array only if it's already allocated, consider using the Safe FREE Macro:

    ABI_SFREE(arr)

that expands to:

    if (allocated(arr)) ABI_FREE(arr)

and is compatible with memory profiling.

This is the list of macros that cannot be inlined:

%(MACRO_NAMES)s
""" % globals()
        print(err_msg)

    return exit_status


if __name__ == "__main__":
    #test_wrong_string()
    sys.exit(main())
