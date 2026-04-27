#!/usr/bin/env python
"""
Check for unresolved git conflict markers.

This script recursively scans the ABINIT source tree to identify any files
that still contain unresolved git conflict markers (e.g., `<<<<<<< TREE`,
`=======`, `>>>>>>> MERGE-SOURCE`).
"""

import os
import re
import sys

from abirules_tools import find_abinit_toplevel_directory

re_markers = re.compile("^(<<<<<<< TREE|=======|>>>>>>> MERGE-SOURCE)$")
re_fbktop  = re.compile("fallbacks$")
re_fbkdir  = re.compile("(exports|sources|stamps)")
re_tmpdir  = re.compile("^tmp")
re_tmpfile = re.compile(r"\.(orig|rej)$")
re_rstfile = re.compile(r"\.rst$")

# TODO: Should look at ..gitignore
exclude_exts = set([
".gz", ".tgz", ".png", ".nc", ".jar", ".xcf", ".pyc",
".gif", ".pct", ".jpg", ".jpeg", ".gif", ".pdf", ".svg",
".a", ".o", ".mod", ".cpkl", ".pickle", ".tar",
".swp", ".swo", ".odt",
])

exclude_bins = set([
  "abinit", "anaddb", "mrgddb", "aim", "fftprof", "mrgdv", "mrgddb", "mrggkk", "lruj",
  "band2eps", "abitk", "cut3d", "fold2Bloch", "conducti", "ioprof", "lapackprof",
  "macroave", "optic", "vdw_kernelgen", "vdw_kernelgen", "mrgscr", "multibinit",
])

def check_item(item: str) -> bool:
  """
  Determine whether a file should be analyzed for conflict markers.

  Args:
      item: The name of the file to check.

  Returns:
      True if the file should be analyzed, False if it matches an ignore pattern.
  """
  if re_tmpfile.search(item): return False
  if re_rstfile.search(item): return False
  if item in exclude_bins: return False

  # check extension
  _, ext = os.path.splitext(item)
  if ext and ext.lower() in exclude_exts: return False
  return True


def main() -> int:
  """
  Main logic for validating conflict markers.

  Returns:
      Number of files found containing conflict markers (0 if OK).
  """
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

    # Ignore hidden directories
    hidden_dirs = [d for d in dirs if d.startswith('.')]
    for d in hidden_dirs: dirs.remove(d)

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
      #print("Checking path:", path)

      try:
          if sys.version_info >= (3, 0):
            with open(path, encoding="ISO-8859-1") as fh:
              chk_data = fh.readlines()
          else:
            with open(path) as fh:
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
