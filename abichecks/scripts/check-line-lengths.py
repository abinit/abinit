#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import re
import os
import sys

from abirules_tools import find_src_dirs

# Init
re_srcfile = re.compile("\.([Ff]|[Ff]90|finc)$")
len_limit = 132


def main():
    retval = 0
    for top in find_src_dirs():
        for root, dirs, files in os.walk(top):
            # Check line lengths in Fortran source files
            for item in files:
                if re_srcfile.search(item) and item != "m_build_info.F90":
                    lineno = 1
                    path = os.path.join(root, item)
                    with open(path, "rt") as fh:
                        for line in fh:
                            line = re.sub("!.*", "", line)
                            line = re.sub("\n", "", line)
                            if len(line) > len_limit:
                                sys.stderr.write(
                                    "%s: line %d has more than %d characters\n" % (path, lineno, len_limit))
                                sys.stdout.write(
                                    "%s: line %d has more than %d characters\n" % (path, lineno, len_limit))
                                retval = 1
                            lineno += 1
    return retval


if __name__ == "__main__":
    sys.exit(main())
