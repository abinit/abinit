# coding: utf-8
from __future__ import print_function, division, absolute_import, unicode_literals

import os


def find_abinit_src_directory(start_path=None, ntrials=10):
    top = find_abinit_toplevel_directory()
    return os.path.join(top, "src")


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
