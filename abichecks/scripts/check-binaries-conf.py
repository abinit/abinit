#!/usr/bin/env python
"Check binaries configuration"
#
# Copyright (C) 2012-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#
from __future__ import unicode_literals, division, print_function, absolute_import

from abirules_tools import find_abinit_toplevel_directory

try:
    from ConfigParser import ConfigParser,NoOptionError
except ImportError:
    from configparser import ConfigParser, NoOptionError

import os
import re
import sys

class MyConfigParser(ConfigParser):

  def optionxform(self,option):
    return str(option)

# ---------------------------------------------------------------------------- #

#
# Auxiliary data
#

dep_levels = {
  "bigdft":10,
  "fft":8,
  "levmar":9,
  "libpsml":9,
  "libxc":9,
  "libxml2":0,
  "linalg":7,
  "mpi":1,
  "gpu":2,
  "hdf5":4,
  "netcdf":5,
  "netcdf_fortran":6,
  "papi":3,
  "triqs":3, 
  "wannier90":9,
  "xmlf90":3,
}

def main():
  home_dir = find_abinit_toplevel_directory()

  # Init
  cnf_bin = MyConfigParser()
  cnf_fname = os.path.join(home_dir, "config/specs/binaries.conf")
  assert os.path.exists(cnf_fname)
  cnf_bin.read(cnf_fname)
  bin_list = cnf_bin.sections()
  bin_list.sort()
  dep_order = {}
  lib_order = {}
  re_num = re.compile("[0-9][0-9]_")

  # Check order of dependencies and libraries
  for prg in bin_list:
    if cnf_bin.has_option(prg,"dependencies"):
      bin_deps = cnf_bin.get(prg,"dependencies").split()
    else:
      bin_deps = []
    dep_old = 100
    dep_new = 100
    for dep in bin_deps:
      if dep in dep_levels:
        dep_new = dep_levels[dep]
      else:
        sys.stderr.write("Error: unregistered dependency '%s'\n" %  dep)
        sys.exit(10)
      if dep_new > dep_old:
        if prg not in dep_order:
          dep_order[prg] = list()
        dep_order[prg].append(dep)
      dep_old = dep_new

    if cnf_bin.has_option(prg,"libraries"):
      bin_libs = cnf_bin.get(prg,"libraries").split()
    else:
      bin_libs = list()
    lib_old = 100
    lib_new = 100
    for lib in bin_libs:
      if re_num.match(lib):
        lib_new = int(re.sub("_.*","",lib))
        if lib_new > lib_old:
          if prg not in lib_order:
            lib_order[prg] = list()
          lib_order[prg].append(lib)
      lib_old = lib_new

  # Report any disorder
  nerr = len(dep_order) + len(lib_order)

  if nerr > 0 :
    sys.stderr.write("%s: reporting disordered libraries\n\n" % (os.path.basename(sys.argv[0])))
    sys.stderr.write("X: D=Dependency / L=Library\n\n")
    sys.stderr.write("%s  %-24s  %-24s\n" % ("X","Program","Dependency/Library"))
    sys.stderr.write("%s  %s  %s\n" % ("-","-" * 24,"-" * 24))

    dep_keys = list(dep_order.keys())
    dep_keys.sort()
    for prg in dep_keys:
      for dep in dep_order[prg]:
        sys.stderr.write("%s  %-24s  %-24s\n" % ("D",prg,dep))
    lib_keys = list(lib_order.keys())
    lib_keys.sort()
    for prg in lib_keys:
      for lib in lib_order[prg]:
        sys.stderr.write("%s  %-24s  %-24s\n" % ("L",prg,lib))
    sys.stderr.write("\n")

  return nerr


if __name__ == "__main__":
  sys.exit(main())
