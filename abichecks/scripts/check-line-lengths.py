#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import re
import os
import sys

# Init
re_srcfile = re.compile("\.([Ff]|[Ff]90|finc)$")
len_limit = 132


def main(top):
    retval = 0
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

    if len(sys.argv) == 1:
        top = "src"
    else:
        top = sys.argv[1]

    if not os.path.exists(top):
        raise ValueError("Path %s does not exist" % top)

    sys.exit(main(top))
