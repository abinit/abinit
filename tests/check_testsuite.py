#!/usr/bin/env python
"""Script for checking the ABINIT automatic tests."""
from __future__ import print_function, division, absolute_import #, unicode_literals

import sys
import os

from pprint import pprint
from optparse import OptionParser

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

#try:
#    import tests
#except ImportError:
# Add the directory [...]/abinit/tests to $PYTHONPATH
pack_dir, tail = os.path.split(os.path.abspath(__file__))
pack_dir, tail = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)

import tests
from tests import abitests
abenv = tests.abenv

__version__ = "0.1"
__author__ = "Matteo Giantomassi"


def check_authors(suite):
    def first_second_name(string):
        idx = string.rfind(".")
        if idx == -1:
            first, second = "", string
        else:
            first, second = string[:idx], string[idx+1:]

        return first.strip(), second.strip()

    second_names = []
    for test in suite:
        if not hasattr(test, "authors"):
            authors = []
            for t in test:
                authors.extend(t.authors)
            authors = set(authors)
        else:
            authors = test.authors

        for string in authors:
            f, s = first_second_name(string)
            if not f and s and s != "Unknown":
                print("author(s) first name is missing in file %s, string = %s " %(test.full_id, s))
            second_names.append(s)
        #print(test.id, first_second_name(test.authors[0])[1])

    return set(second_names)

def get_allowed_cpp_vars():
    """
    Inspect the libpaw header file, the autoconf macros and config.ac
    Extract and return the set of allowed CPP options, used to check
    the exclude_cpp_vars TEST_INFO section for possible typos.

    Based on ~abinit/abichecks/scripts/check-cpp-options.
    """
    import re
    re_m4file = re.compile("\.m4$")
    re_hdrfile = re.compile("\.h$")
    re_acdef = re.compile("AC_DEFINE\\(")
    re_cppdef = re.compile("^([ ]?)+#([ ]?)+define [0-9A-Z_]*")
    
    abidir = os.path.abspath("../")

    # Extract CPP options from the libPAW header files
    cpp_libpaw = set()
    for root, dirs, files in os.walk(os.path.join(abidir, "src/39_libpaw")):
        for src in files:
            if not re_hdrfile.search(src): continue
            with open(os.path.join(root, src), "rt") as fh:
                for line in fh:
                    if not re_cppdef.search(line): continue
                    tmp_def = re.sub("^[# ]*define[ ]*([0-9A-Z_]*).*","\\1", line).strip()
                    cpp_libpaw.add(tmp_def)

    # Extract CPP options from the build system
    cpp_buildsys = set()
    for root, dirs, files in os.walk(os.path.join(abidir, "config/m4")):
        for src in files:
            if not re_m4file.search(src): continue
            with open(os.path.join(root, src), "rt") as fh:
                for line in fh:
                    if not re_acdef.search(line): continue
                    tmp_def = re.sub(".*AC_DEFINE\\([\\[]?([^\\],]*).*","\\1",line).strip()
                    cpp_buildsys.add(tmp_def)

    with open(os.path.join(abidir, "configure.ac"), "rt") as fh:
        for line in fh:
            if not re_acdef.search(line): continue
            tmp_def = re.sub(".*AC_DEFINE\\([\\[]?([^\\],]*).*","\\1",line).strip()
            cpp_buildsys.add(tmp_def)

    return cpp_buildsys.union(cpp_libpaw)

def main():
    usage = "usage: %prog [suite_name] [options] [-h|--help] for help)"
    version = "%prog "+ str(__version__)
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="verbose mode")

    options, args = parser.parse_args()

    # Get the full database.
    # TODO should use with_disabled=True
    full_database = abitests.build_database(with_disabled=False)

    retcode = 0

    print("Stale or lost reference files...  ", end="")
    err = full_database.find_stale_or_lost_refs()

    if err:
        retcode += len(err)
        print("FAILED")
        sys.stderr.write(err)
    else:
        print("OK")

    print("Stale or lost inputs...  ", end="")
    err = full_database.find_stale_or_lost_inputs()

    if err:
        retcode += len(err)
        print("FAILED")
        sys.stderr.write(err)
    else:
        print("OK")

    unknowns, wrong = full_database.find_unknown_wrong_keywords()

    print("Unknown keywords...  ", end="")
    if unknowns:
        retcode += len(unknowns)
        print("FAILED")
        print("The following keywords are not documented:\n\t%s" % unknowns)
        print("ACTION: Add the corresponding documentation to the KNOWN_KEYWORDS dictionary defined in tests/__init__.py")
    else:
        print("OK")

    print("Wrong keywords ...  ", end="")
    if wrong:
        retcode += len(wrong)
        print("FAILED")
        print("The following keywords contain blank spaces:\n\t%s" % wrong)
        print("ACTION: Replace blank spaces with underscores")
    else:
        print("OK")

    print("Testing whether important TEST_INFO entries are present...  ", end="")
    errstr = full_database.check_testinfo_options()
    if errstr:
        retcode += 1
        print("FAILED")
        print(errstr)
    else:
        print("OK")

    # Check authors.
    #print("Testing whether authors are defined...  ", end="")
    #second_names = set()
    #for suite_name, suite in full_database.items():
    #  second_names = second_names.union(check_authors(suite))

    #if second_names:
    #    retcode += len(second_names)
    #    print("FAILED")
    #    pprint(second_names)
    #else:
    #    print("OK")
    
    # Add test on CPP options
    allowed_cpp_vars = get_allowed_cpp_vars()
    #print(allowed_cpp_vars)
    for suite_name, suite in full_database.items():
        for test in suite:
            # Remove ! from string e.g. !HAVE_MPI
            tvars = set(v[1:] if v.startswith("!") else v for v in test.need_cpp_vars)
            diff = tvars.difference(allowed_cpp_vars)
            if diff:
                print("in test: ", test)
                print("diff", diff)


    return retcode


if __name__ == "__main__":
    sys.exit(main())
