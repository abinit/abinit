#!/usr/bin/env python
"""
This script analyzes the `abimem_rank.mocc` files produced by Abinit when
the code is executed in memory-profiling mode.
"""
from __future__ import print_function, division, unicode_literals

__version__ = "0.1.0"
__author__ = "Matteo Giantomassi"

import sys
import os
import re
try:
    import argparse
except ImportError:
    print("abimem.py requires argparse module and python >= 2.7")
    raise

import logging
logger = logging.getLogger(__name__)

# We don't install with setup.py hence we have to add the directory [...]/abinit/tests to $PYTHONPATH
# TODO: Use Absolute imports and rename tests --> abitests to
# avoid possible conflicts with the packages in PYTHONPATH
# monty installs the subpackage paths and this breaks the import below
pack_dir, x = os.path.split(os.path.abspath(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)

try:
    from tests.pymods.memprof import AbimemParser
except ImportError:
    print("Cannot find tests.pymods.memprof module in python path")
    print("Very likely, you have copied the scripts from ~abinit/test/Scripts to another directory")
    print("In this case, you have to add the abinit directory to the python path by executing:\n")
    print("     export PYTHONPATH=~/abinit_directory\n")
    print("before running the script.")
    print("Note that there's no need to change PYTHONPATH if you invoke the script as ~abinit/tests/Scripts/abimem.py")
    raise


def leaks(options, parsers):
    """Find possible memory leaks"""
    retcode = 0
    for parser in parsers:
        parser.find_memleaks()
    return retcode


def small(options, parsers):
    """Find small allocations."""
    retcode = 0
    for parser in parsers:
        smalls = parser.find_small_allocs() #nbytes=options.nbytes)
    return retcode


def intens(options, parsers):
    """Find routines with intensive allocations."""
    retcode = 0
    for parser in parsers:
        intensive = parser.find_intensive() #threshold=options.threshold)
    return retcode


def peaks(options, parsers):
    """Find memory peaks."""
    retcode = 0
    for parser in parsers:
        parser.find_peaks()
    return retcode


def weird(options, parsers):
    """Find weird allocations. i.e. pointers <= 0"""
    retcode = 0
    for parser in parsers:
        parser.find_weird_ptrs()
    return retcode


def zerosized(options, parsers):
    """Find zero-sized allocations."""
    retcode = 0
    for parser in parsers:
        parser.find_zerosized()
    return retcode


def plot(options, parsers):
    """Plot data with matplotlib."""
    for parser in parsers:
        parser.plot_memory_usage()
    return 0


def get_parsers(options):
    if not hasattr(options, "paths"):
        raise ValueError("paths argument is missing")

    if os.path.isfile(options.paths[0]):
        return [AbimemParser(path) for path in options.paths]

    else:
        # Walk directory tree and find abimem files.
        top = options.paths[0]
        if len(options.paths) != 1:
           raise ValueError("Expecting one argument with dirname, got %s" % len(options.paths))
        if not os.path.isdir(top):
           raise ValueError("Expecting existenting directory, got %s" % top)

        re_abimem = re.compile("^abimem_rank(\d+)\.mocc$")

        paths = []
        for dirpath, dirnames, filenames in os.walk(top):
            for f in filenames:
                if not re_abimem.match(f): continue
                paths.append(os.path.join(dirpath, f))

        options.paths = paths
        print("Will analyze %s abimem file(s)" % len(paths))
        if not paths: sys.exit(1)
        return [AbimemParser(path) for path in paths]


def main():
    def str_examples():
        return """\
usage example:
    abimem.py leaks [FILES]       => Find possible memory leaks in FILE(s)
    abimem.py small [FILES]       => Find small memory allocations in FILE(s)
    abimem.py intens [FILES]      => Find periods of intense memory allocation in FILE(s)
    abimem.py peaks [FILES]       => Find peaks in memory allocation in FILE(s)
    abimem.py plot [FILES]        => Plot memory allocations in FILE(s) with matplotlib

    FILES could be either a list of files or a single directory containing abimem_ran.mocc files.

    TIP: To profile the python code add `prof` before the command e.g.
         abimem.py prof leaks [FILES]
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    paths_selector_parser = argparse.ArgumentParser(add_help=False)
    paths_selector_parser.add_argument('paths', nargs="+", help="List of files or directory containing abimem files.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Build Subparsers for commands
    p_leaks = subparsers.add_parser('leaks', parents=[paths_selector_parser], help=leaks.__doc__)
    # Subparser for small
    p_small = subparsers.add_parser('small', parents=[paths_selector_parser], help=small.__doc__)
    # Subparser for intens
    p_intensive = subparsers.add_parser('intens', parents=[paths_selector_parser], help=intens.__doc__)
    # Subparser for peaks command.
    p_peaks = subparsers.add_parser('peaks', parents=[paths_selector_parser], help=peaks.__doc__)
    # Subparser for weird command.
    p_weird = subparsers.add_parser('weird', parents=[paths_selector_parser], help=weird.__doc__)
    # Subparser for zerosized command.
    p_zerosized = subparsers.add_parser('zerosized', parents=[paths_selector_parser], help=zerosized.__doc__)
    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[paths_selector_parser], help=plot.__doc__)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    parsers = get_parsers(options)

    # Dispatch
    return globals()[options.command](options, parsers)


if __name__ == "__main__":
    # Check whether we are in profiling mode
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof: sys.argv.pop(1)
    except:
        do_prof = False

    if not do_prof:
        sys.exit(main())
    else:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
