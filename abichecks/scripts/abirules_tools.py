# coding: utf-8
from __future__ import print_function, division, absolute_import, unicode_literals

import os

def find_abinit_src_directory(start_path=None, ntrials=10):
    top = find_abinit_toplevel_directory(start_path=start_path, ntrials=ntrials)
    return os.path.join(top, "src")


def find_src_dirs(start_path=None, ntrials=10):
    """
    Return list of directories containing source files
    taking into account the new division in shared and ~abinit/src.
    """
    top = find_abinit_toplevel_directory(start_path=start_path, ntrials=ntrials)
    dlist= [os.path.join(top, "shared", "common", "src"),
            os.path.join(top, "shared", "libpaw", "src"),
            os.path.join(top, "src")]
    if not all(os.path.isdir(d) for d in dlist):
        raise RuntimeError("Non-existent directory in list:" % str(dlist))
    return dlist


def find_abinit_toplevel_directory(start_path=None, ntrials=10):
    """
    Returns the absolute path of the ABINIT top level directory.
    Use current working directory is start_path is None

    Raises:
        `RuntimeError` if build tree is not found after ntrials attempts.
    """
    if start_path is None: start_path = os.getcwd()
    abs_path = os.path.abspath(start_path)

    trial = 0
    while trial <= ntrials:
        config_ac = os.path.join(abs_path, "configure.ac")
        abinit_f90 = os.path.join(abs_path, "src", "98_main", "abinit.F90")
        # Check if we are in the top of the ABINIT source tree
        found = os.path.isfile(config_ac) and os.path.isfile(abinit_f90)

        if found:
            return abs_path
        else:
            abs_path, tail = os.path.split(abs_path)
            trial += 1

    raise RuntimeError("Cannot find the ABINIT top level directory after %s trials" % ntrials)
