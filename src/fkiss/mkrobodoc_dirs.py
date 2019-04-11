#!/usr/bin/env python
"""
This script generates the ROBODOC headers located in the Abinit directories (e.g src/70_gw/_70_gw_)
Usage: mkrobodoc_dirs.py abinit/src/
"""
from __future__ import print_function

import sys
import os
import fnmatch

def is_string(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        dummy = s + " "
        return True

    except TypeError:
        return False


def list_strings(arg):
    """
    Always return a list of strings, given a string or list of strings as
    input.

    :Examples:

    >>> list_strings('A single string')
    ['A single string']

    >>> list_strings(['A single string in a list'])
    ['A single string in a list']

    >>> list_strings(['A','list','of','strings'])
    ['A', 'list', 'of', 'strings']
    """
    if is_string(arg):
        return [arg]
    else:
        return arg


class WildCard(object):
    """
    This object provides an easy-to-use interface for
    filename matching with shell patterns (fnmatch).

    .. example:

    >>> w = WildCard("*.nc|*.pdf")
    >>> w.filter(["foo.nc", "bar.pdf", "hello.txt"])
    ['foo.nc', 'bar.pdf']

    >>> w.filter("foo.nc")
    ['foo.nc']
    """
    def __init__(self, wildcard, sep="|"):
        """
        Args:
            wildcard:
                String of tokens separated by sep.
                Each token represents a pattern.
            sep:
                Separator for shell patterns.
        """
        self.pats = ["*"]
        if wildcard:
            self.pats = wildcard.split(sep)

    def __str__(self):
        return "<%s, patterns = %s>" % (self.__class__.__name__, self.pats)

    def filter(self, names):
        """
        Returns a list with the names matching the pattern.
        """
        names = list_strings(names)

        fnames = []
        for f in names:
            for pat in self.pats:
                if fnmatch.fnmatch(f, pat):
                    fnames.append(f)

        return fnames

    def match(self, name):
        """
        Returns True if name matches one of the patterns.
        """
        for pat in self.pats:
            if fnmatch.fnmatch(name, pat):
                return True

        return False


def robodoc_dheader(dirname):
    """Return a string with the ROBODOC header for the specified directory."""
    dirname = os.path.basename(dirname)

    return """\
!!****d* ABINIT/%(dirname)s
!! NAME
!! %(dirname)s
!!
!! DESCRIPTION
!!  FIXME: Description is missing
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT Group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! CHILDREN
""" % locals()


def mkrobodoc_files(top):
    """
    Generate the ROBODOC files in all the ABINIT directories
    located within the top level directory top.

    Returns:
        Exit status.
    """
    top = os.path.abspath(top)
    # Select files with these extensions.
    wildcard = WildCard("*.F90|*.finc")

    # Walk the directory tree starting from top
    # Find the source files in the Abinit directories
    # and add their name to the _dirname_ file used by robodoc
    wrong_dirpaths = []

    for dirpath, dirnames, filenames in os.walk(top):
        dirname = os.path.basename(dirpath)
        if "dirname" in ("__pycache__", ): continue
        robo_dfile = "_" + dirname  + "_"

        if robo_dfile not in filenames:
            # Not an abinit dir.
            wrong_dirpaths.append(os.path.abspath(dirpath))
            continue

        robo_dfile = os.path.abspath(os.path.join(dirpath, robo_dfile))
        with open(robo_dfile, "r") as f:
            robotext = []
            for line in f:
                robotext.append(line.strip())
                if line.startswith("!! CHILDREN"): break
            else:
                raise ValueError("File %s does not have a valid '!! CHILDREN' section'" % robo_dfile)

        # Get source files, sort them in alphabetical order and generate new robodoc file.
        src_files = wildcard.filter(os.listdir(dirpath))
        if not src_files: continue
        src_files = sorted(src_files)

        robotext +=  ["!! " + src_file for src_file in src_files]
        robotext += ["!!", "!!***"]
        #for line in robotext: print(line)

        # Write new robodoc file.
        with open(robo_dfile, "w") as f:
            f.write("\n".join(robotext) + "\n")

    if wrong_dirpaths:
        # Some dirs should not have a ROBODOC file.
        # This code does not work if we have the same basename but oh well.
        EXCLUDE_BASEDIRS = set([
            "src",
            "replacements",
            "01_interfaces_ext",
            "ptgroup_data",
            "incs",
            "libs",
            "mods",
        ])

        # Remove dirs in EXCLUDE_BASEDIRS
        wrong_dirpaths = [d for d in wrong_dirpaths if os.path.basename(d) not in EXCLUDE_BASEDIRS]

        if wrong_dirpaths:
            print("dirs_without_files")
            for i, dirpath in enumerate(wrong_dirpaths):
                print("%d) %s" % (i, dirpath))
                robodoc_dfile = os.path.join(dirpath, "_" + os.path.basename(dirpath) + "_")
                assert not os.path.isfile(robodoc_dfile)
                with open(robodoc_dfile, "w") as f:
                    f.write(robodoc_dheader(dirpath))

    return len(wrong_dirpaths)


def main():
    try:
        top = os.path.abspath(sys.argv[1])
    except:
        raise ValueError("Top level directory must be specified.\nEx: mkrobodoc_dirs.py ./70_gw")

    return mkrobodoc_files(top)


if __name__ == "__main__":
    sys.exit(main())
