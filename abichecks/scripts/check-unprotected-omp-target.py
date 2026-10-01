#!/usr/bin/env python

import os
import re
import shutil
import sys
from subprocess import run

from abirules_tools import find_src_dirs

# Init
re_srcfile = re.compile(r"\.([Ff]|[Ff]90|finc)$")
re_omp_target = re.compile(r"!\$OMP TARGET")
re_fortran_comment = re.compile(r"!!.*")

black_list = {
"m_build_info.F90",
#"minimax_omega.F90",
#"minimax_tau.F90",
}


def main():
    if shutil.which("gfortran") is None:
        sys.stderr.write(
            "check-unprotected-omp-target.py: 'gfortran' not found in PATH -- "
            "it is used here only to preprocess Fortran sources (-E), not to compile. "
            "Load a GNU toolchain module (or otherwise put gfortran on PATH) before running this check.\n"
        )
        return 1

    retval = 0
    for top in find_src_dirs():
        for root, dirs, files in os.walk(top):
            # Check line lengths in Fortran source files
            for item in files:
                if re_srcfile.search(item) and item not in black_list:
                    lineno = 1
                    path = os.path.join(root, item)

                    # Most source files cannot produce a finding. Avoid the
                    # relatively expensive compiler startup for those files.
                    with open(path, encoding="utf-8") as fh:
                        if not any(re_omp_target.search(re_fortran_comment.sub("", line)) for line in fh):
                            continue

                    cmd_list = ["gfortran", "-cpp", "-E", path]
                    cmd_list.append(f"-I{root}")
                    cmd_list.append(f"-I{top}/incs")
                    cmd_list.append(f"-I{top}/../shared/common/src/incs")
                    cmd_list.append(f"-I{top}/../shared/libpaw/incs")
                    p = run(cmd_list, stdin=None, timeout=5, capture_output=True, encoding="utf-8", check=False)
                    if p.returncode != 0:
                        sys.stderr.write(
                            "%s: gfortran preprocessing failed with exit status %d:\n%s\n"
                            % (path, p.returncode, p.stderr)
                        )
                        retval = 1
                        continue

                    for line in p.stdout.splitlines():
                        line = re_fortran_comment.sub("", line)
                        if re_omp_target.search(line):
                            sys.stderr.write(
                                "%s: line %d has an unprotected OMP TARGET directive:\n\n%s\n" % (path, lineno, line))
                            sys.stdout.write(
                                "%s: line %d has an unprotected OMP TARGET directive:\n\n%s\n" % (path, lineno, line))
                            retval = 1
                        lineno += 1
    return retval


if __name__ == "__main__":
    sys.exit(main())
