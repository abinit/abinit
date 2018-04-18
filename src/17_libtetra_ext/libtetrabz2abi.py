#!/usr/bin/env python
"""
Change libtetrabz source file so that they are compatible with the Abinit build system.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import os

here = os.path.dirname(__file__)
f90_files = [f for f in os.listdir(here) if f.startswith("libtetrabz") and f.endswith(".F90")]
paths = [os.path.join(here, f) for f in f90_files]

def abi_sanitize(line):
    line = line.replace("__MPI", "HAVE_MPI")
    # Bind(C) with optional args requires F2015
    line = line.replace("BIND(C)", "")
    return line

print("Replacing __MPI with HAVE_MPI in the following routines:")
for path in paths:
    print("\t", path)
    lines = ["""
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
"""]
    lines = []

    with open(path, "rt") as fh:
        lines += [abi_sanitize(l) for l in fh]

    print(lines)
    with open(path, "wt") as fh:
        fh.writelines(lines)
