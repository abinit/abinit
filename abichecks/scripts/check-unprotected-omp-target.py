#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import re
import os
import sys
from subprocess import PIPE, run

from abirules_tools import find_src_dirs

# Init
re_srcfile = re.compile(r"\.([Ff]|[Ff]90|finc)$")
len_limit = 132

black_list = {
"m_build_info.F90",
#"minimax_omega.F90",
#"minimax_tau.F90",
}


def main():
    retval = 0
    out = "out.F90"
    for top in find_src_dirs():
        for root, dirs, files in os.walk(top):
            # Check line lengths in Fortran source files
            for item in files:
                if re_srcfile.search(item) and item not in black_list:
                    lineno = 1
                    path = os.path.join(root, item)
                    cmd_list = ["gfortran", "-E", path]
                    cmd_list.append(f"-I{root}")
                    cmd_list.append(f"-I{top}/incs")
                    cmd_list.append(f"-I{top}/../shared/common/src/incs")
                    cmd_list.append(f"-I{top}/../shared/libpaw/incs")
                    p = run(cmd_list, stdin=None, timeout=5, capture_output=True, encoding="utf-8")
                    output = p.stdout.split('\n')
                    for line in output:
                        line = re.sub("!!.*", "", line)
                        line = re.sub("\n", "", line)
                        if re.search(r"!\$OMP TARGET", line):
                            sys.stderr.write(
                                "%s: line %d has an unprotected OMP TARGET directive:\n\n%s\n" % (path, lineno, line))
                            sys.stdout.write(
                                "%s: line %d has an unprotected OMP TARGET directive:\n\n%s\n" % (path, lineno, line))
                            retval = 1
                        lineno += 1
    return retval


if __name__ == "__main__":
    sys.exit(main())
