#!/usr/bin/env python
from __future__ import print_function

import sys

from variables import *

try:
    import yaml
except ImportError:
    raise ImportError("yaml package is not available. Install it with `pip install pyyaml`")


def main():
    path = 'abinit_vars.yml'
    with open(path, 'r') as f:
        variables = yaml.load(f)

    retcode = 0
    for var in variables:

        if var.varname is None:
            print('Variable %s has no name. This is forbidden' % str(var))
            retcode += 1
            continue

        varname = var.varname

        if var.vartype is None:
            print('FAIL: ', varname, ' does not have a vartype')
            retcode += 1

        if var.characteristics is not None:
            if not isinstance(var.characteristics, list):
                print('FAIL: the field characteristics of ', varname, ' is not a list')
                retcode += 1
            else:
                for cat in var.characteristics:
                    if cat.replace("[[", "").replace("]]", "") not in list_chars:
                        print('FAIL: the characteristics ', cat, ' of ', varname, ' is not valid')
                        retcode += 1

        if var.dimensions is not None:
            if var.dimensions != "scalar":
                if not isinstance(var.dimensions, list) \
                        and not isinstance(var.dimensions, ValueWithConditions):
                    print('FAIL: the field dimensions of ', varname, ' is not a list neither a valuewithconditions')
                    retcode += 1

        if var.section is None:
            print('FAIL: ', varname, ' does not have a section')
            retcode += 1
        else:
            if not isinstance(var.section, str) or var.section not in list_sections:
                print('FAIL: the field section of ', varname, ' should be one of the valid sections')
                retcode += 1

    return retcode


if __name__ == "__main__":
    sys.exit(main())
