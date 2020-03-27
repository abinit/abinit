#!/usr/bin/env python
"check build refs"
#
# Copyright (C) 2011-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#
from __future__ import unicode_literals, division, print_function, absolute_import

from abirules_tools import find_abinit_toplevel_directory

from time import gmtime,strftime

import os
import re
import sys

def getstatusoutput(cmd):
    """
    Return (status, output) of executing cmd in a shell.

    Execute the string 'cmd' in a shell with 'check_output' and
    return a 2-tuple (status, output). Universal newlines mode is used,
    meaning that the result with be decoded to a string.

    A trailing newline is stripped from the output.
    The exit status for the command can be interpreted
    according to the rules for the function 'wait'. Example:
    """
    from subprocess import check_output, STDOUT, CalledProcessError
    try:
        data = check_output(cmd, shell=True, universal_newlines=True, stderr=STDOUT)
        status = 0
    except CalledProcessError as ex:
        data = ex.output
        status = ex.returncode
    if data[-1:] == '\n':
        data = data[:-1]
    return status, data


def main():
  home_dir = find_abinit_toplevel_directory()
  # Init
  nerr = 0
  bex_diffs = list()
  bex_missing = list()

  bex_dir = os.path.join(home_dir,"doc/build/config-examples")
  ref_dir = os.path.join(home_dir,"abichecks/buildsys/Refs")
  assert os.path.exists(bex_dir) and os.path.exists(ref_dir)

  # Check files
  ref_list = os.listdir(ref_dir)
  ref_list.sort()
  for ref_file in ref_list:
    if os.path.exists("%s/%s" % (bex_dir,ref_file)):
      ret, tmp = getstatusoutput("diff  %s/%s %s/%s" % (ref_dir,ref_file,bex_dir,ref_file))
      if ret != 0:
        bex_diffs.append(ref_file)
        sys.stdout.write(tmp)
    else:
      bex_missing.append(ref_file)

  nerr = len(bex_diffs) + len(bex_missing)

  # Report any mismatch
  if nerr > 0:
    sys.stderr.write("%s: reporting wrongly generated build examples\n\n" % (os.path.basename(sys.argv[0])))
    sys.stderr.write("Reference files are in ~abinit/abichecks/buildsys")
    sys.stderr.write("X: D=Difference detected / M=Missing File\n\n")
    sys.stderr.write("%s  %-64s\n" % ("X","File"))
    sys.stderr.write("%s  %s\n" % ("-","-" * 64))

    for bex in bex_diffs:
      sys.stderr.write("%s  %-64s\n" % ("D",bex))
    for bex in bex_missing:
      sys.stderr.write("%s  %-64s\n" % ("M",bex))

    sys.stderr.write("\n")

  return nerr

if __name__ == "__main__":
  sys.exit(main())
