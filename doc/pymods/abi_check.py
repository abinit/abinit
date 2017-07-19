#!/usr/bin/env python
from __future__ import print_function

import sys

from variables import *

try:
  import yaml
except ImportError:
  raise ImportError("yaml package is not available. Install it with `pip install pyyaml`")


def main():
  path = '../input_variables/origin_files/'
  with open(path+'abinit_vars.yml', 'r') as f:
    variables = yaml.load(f)
  with open(path+'characteristics.yml', 'r') as f:
    list_chars = yaml.load(f)
  with open(path+'varfiles.yml', 'r') as f:
    list_varfiles = yaml.load(f)
    varfile_names=[]
    for varfile in list_varfiles:
      varfile_names.append(varfile.name)

  retcode = 0
  for var in variables:

    if var.abivarname is None:
      print('Variable %s has no name. This is forbidden' % str(var))
      retcode += 1
      continue

    abivarname = var.abivarname

    if var.vartype is None:
      print('FAIL: ', abivarname, ' does not have a vartype')
      retcode += 1

    if var.characteristics is not None:
      if not isinstance(var.characteristics, list):
        print('FAIL: the field characteristics of ', abivarname, ' is not a list')
        retcode += 1
      else:
        for cat in var.characteristics:
          if cat.replace("[[", "").replace("]]", "") not in list_chars:
            print('FAIL: the characteristics ', cat, ' of ', abivarname, ' is not valid')
            retcode += 1

    if var.dimensions is not None:
      if var.dimensions != "scalar":
        if not isinstance(var.dimensions, list) \
              and not isinstance(var.dimensions, ValueWithConditions):
          print('FAIL: the field dimensions of ', abivarname, ' is not a list neither a valuewithconditions')
          retcode += 1

    if var.varfile is None:
      print('FAIL: ', abivarname, ' does not have a varfile')
      retcode += 1
    else:
      if not isinstance(var.varfile, str) or var.varfile not in varfile_names:
        print('FAIL: the field varfile of ', abivarname, ' should be one of the valid varfiles')
        retcode += 1

  return retcode


if __name__ == "__main__":
  sys.exit(main())
