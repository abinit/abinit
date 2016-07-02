#!/usr/bin/env python
"""Search for non-ASCII characters in the ABINIT src files"""
from __future__ import unicode_literals, division, print_function, absolute_import

import os
import re
import sys

__author__ = "M. Giantomassi"


def isascii(string, verbose=False):
    """False if string contains non-ASCII characters."""
    return all(ord(c) < 128 for c in string)
    #try:
    #    string.decode('ascii')
    #    return True
    #except UnicodeDecodeError:
    #    if verbose: print("not a ascii-encoded unicode string")
    #    return False
    #else:
    #    if verbose: print("It may have been an ascii-encoded unicode string")
    #    return False


re_srcfile = re.compile("\.([Ff]|[Ff]90|finc|h|c|cu)$")


def is_srcfile(dirpath, fname):
    return re_srcfile.search(fname)


def abinit_test_generator():
    def test_func(abenv):
        "Search for non-ASCII characters in the ABINIT src files"
        top = abenv.apath_of("src")
        try:
            return main(top)
        except Exception:
            import sys
            raise sys.exc_info()[1]  # Reraise current exception (py2.4 compliant)

    return {"test_func": test_func}


def main(top):
    exit_status = 0
    for dirpath, dirnames, files in os.walk(top):
        for src in files:
            if is_srcfile(dirpath, src):
                fpath = os.path.join(dirpath, src)
                with open(fpath, "rt") as fh:
                    lines = fh.readlines()
                for lno, line in enumerate(lines):
                    if not isascii(line):
                        exit_status += 1
                        print(">>> Non-ASCII character at: ", fpath, ":", lno + 1)
                        print(line)
    return exit_status


if __name__ == "__main__":

    if len(sys.argv) == 1:
        top = "../../../src"
        print("--------------------------------------------------------")
        print(" Searching for non-ASCII characters in ABINIT src files ")
        print("--------------------------------------------------------")
    else:
        top = sys.argv[1]

    exit_status = main(top)
    sys.exit(exit_status)
