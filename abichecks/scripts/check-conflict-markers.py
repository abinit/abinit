#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import re
import os
import sys

from abirules_tools import find_abinit_toplevel_directory

re_markers = re.compile("^(<<<<<<< TREE|=======|>>>>>>> MERGE-SOURCE)$")
re_fbktop  = re.compile("fallbacks$")
re_fbkdir  = re.compile("(exports|sources|stamps)")
re_tmpdir  = re.compile("^tmp")
re_tmpfile = re.compile("\.(orig|rej)$")
re_rstfile = re.compile("\.rst$")

# TODO: Should look at ..gitignore
exclude_exts = set([
".gz", ".tgz", ".png", ".nc", ".jar", ".xcf", ".pyc",
".gif", ".pct", ".jpg", ".jpeg", ".gif", ".pdf", ".svg",
".a", ".o", ".mod", ".cpkl", ".pickle", ".tar",
".swp", ".swo", ".odt",
])

exclude_bins = set([
  "abinit", "anaddb", "mrgddb", "aim", "fftprof", "mrgdv", "mrgddb", "mrggkk", "ujdet",
  "band2eps", "abitk", "cut3d", "fold2Bloch", "conducti", "ioprof", "lapackprof",
  "macroave", "optic", "vdw_kernelgen", "vdw_kernelgen", "mrgscr",
])

def check_item(item):
  "True if item has to be analyzed."
  if re_tmpfile.search(item): return False
  if re_rstfile.search(item): return False
  if item in exclude_bins: return False

  # check extension
  _, ext = os.path.splitext(item)
  if ext and ext.lower() in exclude_exts: return False
  return True


def main():
  retval = 0
  top = find_abinit_toplevel_directory()
  assert os.path.exists(top)

  for root, dirs, files in os.walk(top):
    # Ignore Makefiles
    if "Makefile.am" in files: files.remove("Makefile.am")
    if "Makefile.in" in files: files.remove("Makefile.in")
    if "Makefile" in files: files.remove("Makefile")

    # Ignore Autotools subdirs
    if "autom4te.cache" in dirs: dirs.remove("autom4te.cache")

    # Ignore temporary dirs
    garb_dirs = [item for item in dirs if re_tmpdir.match(item)]
    for d in garb_dirs: dirs.remove(d)

    # Ignore installed fallbacks
    if re_fbktop.search(root):
      garb_dirs = [item for item in dirs if re_fbkdir.match(item)]
      for d in garb_dirs: dirs.remove(d)

    # Display conflict markers found
    for item in files:
      path = os.path.join(root, item)
      if not check_item(item): continue

      try:
          if sys.version_info >= (3, 0):
            with open(os.path.join(root, item), "rt", encoding="ISO-8859-1") as fh:
              chk_data = fh.readlines()
          else:
            with open(path, "r") as fh:
              chk_data = fh.readlines()

          chk_stat = False
          for line in chk_data:
            if re_markers.match(line):
              chk_stat = True
              retval += 1
              break

          if chk_stat:
            sys.stderr.write("Found conflict markers in:\n" % path)

      except Exception as exc:
        retval += 1
        sys.stderr.write("Exception while testing: %s\n%s\n" % (path, str(exc)))

  return retval


if __name__ == "__main__":
  sys.exit(main())
