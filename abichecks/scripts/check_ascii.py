#!/usr/bin/env python
"""Search for non-ASCII characters in the ABINIT src files"""
from __future__ import unicode_literals, division, print_function, absolute_import

import os
import re
import sys

from abirules_tools import find_src_dirs

__author__ = "M. Giantomassi"


def isascii(string, verbose=False):
    """False if string contains non-ASCII characters."""
    return all(ord(c) < 128 for c in string)

re_srcfile = re.compile("\.([Ff]|[Ff]90|finc|h|c|cu)$")


def is_srcfile(dirpath, fname):
    return re_srcfile.search(fname)


def main():
    print("--------------------------------------------------------")
    print(" Searching for non-ASCII characters in ABINIT src files ")
    print("--------------------------------------------------------")

    exit_status = 0
    for top in find_src_dirs():
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
    sys.exit(main())
