#!/usr/bin/env python
# coding: utf-8
"""
This script replaces strings in the source files. See replace_string function.
"""
import sys
import os

TOPDIR = os.path.dirname(os.path.realpath(__file__))


def find_source_paths():
    """
    Return list with absolute paths of source files extracted from abinit.src
    """
    import imp
    source_paths = []
    black_list = ["m_build_info.F90", "m_optim_dumper.F90", os.path.basename(__file__)]
    for root, dirs, files in os.walk(TOPDIR):
        files = [f for f in files if f == "abinit.src"]
        if not files: continue
        assert len(files) == 1
        abinit_src = os.path.join(root, "abinit.src")
        mod = imp.load_source(abinit_src, abinit_src)
        if hasattr(mod, "sources"):
            source_paths.extend((os.path.join(root, s) for s in mod.sources if s not in black_list))

    return source_paths


def all_files():
    """
    Return list with the absolute paths of all the files in the project, exclude binary files
    or files that should not be modified.
    """
    def ignore_basename(f):
        black_list = [".pdf", ".pickle", ".swp", ".cpkl", ".gz", ".nc", "_HIST",
                      ".pyc", ".png", "_DEN", ".pdf", ".nc", "_WFK"]
        if any(f.endswith(b) for b in black_list): return True
        if f ==  os.path.basename(__file__): return True
        return False

    all_files = []
    for root, dirs, files in os.walk(TOPDIR):
        #if os.path.basename(root).startswith("_"): continue
        all_files.extend(os.path.join(root, f) for f in files if not ignore_basename(f))
    return all_files


def replace_string(s):
    """
    Main entry point for users.
    Change old2new dict to replace old_expression with new one
    using replace method of python string.
    """
    old2new = {
         "ABI_ALLOCATE(": "ABI_MALLOC(",
         "ABI_DEALLOCATE(": "ABI_FREE(",
         "ABI_DATATYPE_ALLOCATE(": "ABI_MALLOC(",
         "ABI_DATATYPE_DEALLOCATE(": "ABI_FREE(",
         "ABI_STAT_ALLOCATE(": "ABI_STAT_MALLOC(",
         "ABI_DATATYPE_ALLOCATE_SCALAR(": "ABI_MALLOC_SCALAR(",
         "ABI_DATATYPE_DEALLOCATE_SCALAR(": "ABI_FREE(",
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
    old2new = {"informations": "information"}

    for old, new in old2new.items():
        s = s.replace(old, new)

    return s


def main():
    #for path in all_files():
    for path in find_source_paths():
        print("Replacing strings in:", path)
        with open(path, "rt") as fh:
            s = replace_string(fh.read())
        with open(path, "wt") as fh:
            fh.write(s)

    return 0


if __name__ == "__main__":
    main()
