#!/usr/bin/env python
"Check forbidden CPP flags"
#
# Copyright (C) 2010-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#
from __future__ import unicode_literals, division, print_function, absolute_import

from abirules_tools import find_abinit_src_directory

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser
from time import gmtime,strftime

import os
import re
import sys

class MyConfigParser(ConfigParser):

  def optionxform(self,option):
    return str(option)

env_ignore = ["CFLAGS","CXXFLAGS","FCFLAGS","NVCC_*","fcflags_opt_*"]

def is_ignored(keyword):
  for env in env_ignore:
    if "*" in env:
      if re.match(env,keyword):
        return True
    elif env == keyword:
        return True
  return False


def main():
  home_dir = find_abinit_src_directory()
  # Init
  re_dbgflags = re.compile("(^-g|[^0-9A-Za-z]-g)")
  re_optflags = re.compile("(-O[0-9]|-xO[0-9])")

  # Extract environment variables from config file
  cnf_env = MyConfigParser()
  cnf_env.read(os.path.join(home_dir,"config/specs/environment.conf"))
  env_config = list()
  for env in cnf_env.sections():
    if cnf_env.get(env,"reset") == "no":
      if not is_ignored(env):
        env_config.append(env)
  env_config.sort()

  # Extract information from build example config file
  bex_ok = True
  cnf_bex = MyConfigParser()
  cnf_bex.read( os.path.join(home_dir,"config/specs/testfarm.conf") )
  env_forbidden = dict()
  for bot in cnf_bex.sections():
    env_forbidden[bot] = list()
    for var in cnf_bex.options(bot):
      if var in env_config:
        val = cnf_bex.get(bot,var)
        if re_dbgflags.search(val):
          env_forbidden[bot].append(("D",var))
          bex_ok = False
        if re_optflags.search(val):
          env_forbidden[bot].append(("O",var))
          bex_ok = False
  env_fkeys = list(env_forbidden.keys())
  env_fkeys.sort()

  # Report any match
  my_exitcode = 0
  if not bex_ok:
    sys.stderr.write("%s: reporting use of forbidden flags\n\n" % (os.path.basename(sys.argv[0])))
    sys.stderr.write("X: D=debug / O=Optimization\n\n")
    sys.stderr.write("%s  %-24s  %-48s\n" % ("X","Variable","Bot"))
    sys.stderr.write("%s  %s  %s\n" % ("-","-" * 24,"-" * 48))

    for bot in env_fkeys:
      if len(env_forbidden[bot]) > 0:
        my_exitcode = 1
        for (tag,var) in env_forbidden[bot]:
          sys.stderr.write("%s  %-24s  %-48s\n" % (tag,var,bot))

    sys.stderr.write("\n")

  return my_exitcode

if __name__ == "__main__":
  sys.exit(main())
