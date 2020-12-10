#!/usr/bin/env python
# coding: utf-8
"""
This script replaces strings in the source files. See replace_string function.
"""
import sys
import os

TOPDIR = os.path.dirname(os.path.realpath(__file__))


def source_paths_from_abinit_src(top=None):
    """
    Return list with absolute paths of source files extracted from abinit.src
    """
    import imp
    source_paths = []
    black_list = ["m_build_info.F90", "m_optim_dumper.F90", os.path.basename(__file__)]
    if top is None: top = TOPDIR
    for root, dirs, files in os.walk(top):
        files = [f for f in files if f == "abinit.src"]
        if not files: continue
        assert len(files) == 1
        abinit_src = os.path.join(root, "abinit.src")
        mod = imp.load_source(abinit_src, abinit_src)
        if hasattr(mod, "sources"):
            source_paths.extend((os.path.join(root, s) for s in mod.sources if s not in black_list))

    return source_paths


def all_source_files(top=None, types=("fortran", "c", "h")):
    """
    Return list with the absolute paths of all the files in the project, exclude binary files
    or files that should not be modified.
    """
    if top is None: top = TOPDIR
    ext_list = []
    if "fortran" in types: ext_list += [".F90", ".f90", "finc"]
    if "c" in types: ext_list += [".c"]
    if "h" in types: ext_list += [ ".h"]
    def select_basename(f):
        if any(f.endswith(b) for b in ext_list): return True
        if f ==  os.path.basename(__file__): return False
        return False

    all_files = []
    for root, dirs, files in os.walk(top):
        all_files.extend(os.path.join(root, f) for f in files if select_basename(f))
    return all_files


def replace_string(s):
    """
    Main entry point for users.
    Change old2new dict to replace old_expression with new one
    using replace method of python string.
    """
    old2new = {
         # Malloc/free macros
         "ABI_ALLOCATE(": "ABI_MALLOC(",
         "ABI_DEALLOCATE(": "ABI_FREE(",
         "ABI_DATATYPE_ALLOCATE(": "ABI_MALLOC(",
         "ABI_DATATYPE_DEALLOCATE(": "ABI_FREE(",
         "ABI_STAT_ALLOCATE(": "ABI_STAT_MALLOC(",
         "ABI_DATATYPE_ALLOCATE_SCALAR(": "ABI_MALLOC_TYPE_SCALAR(",
         "ABI_DATATYPE_DEALLOCATE_SCALAR(": "ABI_FREE(",
         # Logging macros
         "MSG_COMMENT(": "ABI_COMMENT(",
         "MSG_WARNING(": "ABI_WARNING(",
         "MSG_COMMENT_UNIT(": "ABI_COMMENT_UNIT(",
         "MSG_WARNING_UNIT(": "ABI_WARNING_UNIT(",
         "MSG_ERROR(": "ABI_ERROR(",
         "MSG_ERROR_CLASS(": "ABI_ERROR_CLASS(",
         "MSG_BUG(": "ABI_BUG(",
         "MSG_STOP(": "ABI_STOP(",
         "MSG_ERROR_NODUMP(": "ABI_ERROR_NODUMP(",
         "MSG_ERROR_NOSTOP(": "ABI_ERROR_NOSTOP(",
         "MSG_WARNING_IF(": "ABI_WARNING_IF(",
    }

    #old2new = {
    #     "LIBPAW_COMMENT(": "ABI_COMMENT(",
    #     "LIBPAW_WARNING(": "ABI_WARNING(",
    #     "LIBPAW_COMMENT_UNIT(": "ABI_COMMENT_UNIT(",
    #     "LIBPAW_WARNING_UNIT(": "ABI_WARNING_UNIT(",
    #     "LIBPAW_ERROR(": "ABI_ERROR(",
    #     "LIBPAW_ERROR_CLASS(": "ABI_ERROR_CLASS(",
    #     "LIBPAW_BUG(": "ABI_BUG(",
    #     "LIBPAW_STOP(": "ABI_STOP(",
    #     "LIBPAW_ERROR_NODUMP(": "ABI_ERROR_NODUMP(",
    #     "LIBPAW_ERROR_NOSTOP(": "ABI_ERROR_NOSTOP(",
    #     "LIBPAW_WARNING_IF(": "ABI_WARNING_IF(",
    #}
    ## invert mapping
    #old2new = {v:k for k, v in old2new.items()}

    # This is just to print the table in github md format.
    #from tabulate import tabulate
    #print(str(tabulate(old2new.items(), headers=["OLD", "NEW"], tablefmt="github")))
    #raise RuntimeError("Just to print table.")

    for old, new in old2new.items():
        s = s.replace(old, new)

    return s


def main():

    if "-h" in sys.argv or "--help" in sys.argv:
        print(" Usage: change_source.py [DIR]")
        print(" to replace strings in all source files inside DIR")
        print(" If DIR not provided, uses current directory.")
        return 0

    #for path in all_files():
    #for path in source_paths_from_abinit_src():
    #for path in all_source_files(types=("fortran", "c", "h")):
    top = None
    if len(sys.argv) > 1:
        top = sys.argv[1]
        print("Replacing strings in files inside directory:", top)

    for path in all_source_files(top=top, types=("fortran", "h", "c")):
        print("Replacing strings in:", path)
        with open(path, "rt") as fh:
            s = replace_string(fh.read())
        with open(path, "wt") as fh:
            fh.write(s)

    return 0


if __name__ == "__main__":
    main()
