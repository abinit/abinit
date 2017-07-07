from __future__ import print_function, division, absolute_import #, unicode_literals

import sys
import os
import imp
import platform
import re

from tests.pymods.termcolor import cprint

# Handle py2, py3k differences.
py2 = sys.version_info[0] <= 2
if py2:
    import cPickle as pickle
    from cStringIO import StringIO
else:
    import pickle
    from io import StringIO

from socket import gethostname
from pprint import pprint

#from tests.pymods.testsuite import ChainOfTests, AbinitTestSuite
#from tests.pymods.devtools import FileLock

import logging
logger = logging.getLogger(__name__)


__all__ = [
]


class AbinitEnvironment(object):
    """
    Container with information on the abinit source tree.
    Provide helper functions to construct the absolute path of directories.
    """
    def __init__(self):
        self.uname = platform.uname()
        self.hostname = gethostname()
        try:
            self.username = os.getlogin()
        except:
            self.username = "No_username"

        filepath = os.path.abspath(__file__)
        # Follows symlink (if any)
        filepath = os.path.realpath(filepath)

        # Paths of important dirs.
        self.tests_dir, tail = os.path.split(filepath)
        self.home_dir, tail = os.path.split(self.tests_dir)

        self.src_dir = os.path.join(self.home_dir, "src")
        self.psps_dir = os.path.join(self.tests_dir, "Psps_for_tests")
        self.fldiff_path = os.path.join(self.tests_dir, "Scripts", "fldiff.pl")

    def __str__(self):
        return "\n".join([str(k) + " : " + str(v) for (k, v) in self.__dict__.items()])

    def apath_of(self, *p):
        """
        Return the absolute path of p where p is one or more pathname components
        relative to the top level directory of the package.
        """
        return os.path.join(self.home_dir, *p)

    def isbuild(self):
        """True if we have built the code in the top level directory."""
        configh_path = os.path.join(self.home_dir, "config.h")
        abinit_path = os.path.join(self.home_dir, "src", "98_main", "abinit")
        return os.path.isfile(configh_path) and os.path.isfile(abinit_path)

    def touch_srcfiles(self, patterns):
        """
        Touch all Abinit source files containing the specified list of `patterns`.
        """
        def touch(fname):
            """
            Python touch
            See also http://stackoverflow.com/questions/1158076/implement-touch-using-python
            for a race-condition free version based on py3k
            """
            try:
                os.utime(fname, None)
            except:
                open(fname, 'a').close()

        top = os.path.join(abenv.home_dir, "src")
        print("Walking through directory:", top)
        print("Touching all files containing:", str(patterns))

        for root, dirs, files in os.walk(self.src_dir):
            for path in files:
                path = os.path.join(root, path)
                with open(path, "rt") as fh:
                    for line in fh:
                        if any(p in line for p in patterns):
                            print("Touching %s" % path)
                            touch(path)
                            break

    def start_watching_sources(self):
        def is_source(fname):
            # NB: `.finc` include files are ignored on purpose because
            # make is not tracking them. one should change the Fortran file
            # that includes .finc to force recompilation.
            _, ext = os.path.splitext(fname)
            return ext.lower() in (".f", ".f90", ".c") # ".finc"

        self.srcpath_stat = {}
        for root, dirs, files in os.walk(self.src_dir):
            for fname in files:
                if not is_source(fname): continue
                path = os.path.join(root, fname)
                self.srcpath_stat[path] = os.stat(path)

        return self.srcpath_stat

    def changed_sources(self):
        """Return list of files that have been modified."""
        changed = []
        for path, old_stat in self.srcpath_stat.items():
            now_stat = os.stat(path)
            if now_stat.st_mtime != old_stat.st_mtime:
               changed.append(path)
               self.srcpath_stat[path] = now_stat
        return changed


def exclude_path(p):
    p = os.path.basename(p)
    if p.startswith(".") or p.endswith("~"): return True
    return False


def path2str(path):
    head, fname = os.path.split(path)
    head, x = os.path.split(head)
    _, dirname = os.path.split(head)

    return "["+dirname+"]["+fname+"]"
