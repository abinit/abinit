#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals

__version__ = "0.1.0"
__author__ = "Matteo Giantomassi"

import sys
import os
import argparse

import logging
logger = logging.getLogger(__name__)

# We don't install with setup.py hence we have to add the directory [...]/abinit/tests to $PYTHONPATH
# TODO: Use Absolute imports and rename tests --> abitests to 
# avoid possible conflicts with the packages in PYTHONPATH
pack_dir, x = os.path.split(os.path.abspath(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)

import tests.pymods.etsf_specs as etsf

def validate(options, paths):
    """Validate etsf file(s)."""
    retcode = 0
    for path in paths:
        errors = etsf.validate_ncfile(path)
        retcode += len(errors)
        #title = "=== GROUPS IN FILE %s ===" % os.path.relpath(path)
        #bar =  "=" * len(title)
        #print(bar); print(title); print(bar)
    print("OK" if retcode == 0 else "FAILED")
    return retcode


def groups(options, paths):
    """Show etsf groups present in file(s)."""
    retcode = 0
    for path in paths:
        groups = etsf.find_groups(path)
        title = "=== GROUPS IN FILE %s ===" % os.path.relpath(path)
        bar =  "=" * len(title)
        print(bar); print(title); print(bar)
        #for g in groups: print(g)

    return 0


def ncvars(options, paths):
    """
    Validate etsf variables present in file(s).
    """
    all_errors = []
    for path in paths:
        errors = etsf.validate_vars(path)
        all_errors.append(errors)
        #title = "=== VARS IN FILE %s ===" % os.path.relpath(path)
        #bar =  "=" * len(title)
        #print(bar); print(title); print(bar)

    return len(all_errors)


def get_paths(options):
    if os.path.isfile(options.paths[0]):
        return [path for path in options.paths if path.endswith(".nc")]

    else:
        # Walk directory tree and find netcdf files.
        top = options.paths[0]
        if len(options.paths) != 1:
           raise ValueError("Expecting one argument with dirname, got %s" % len(options.paths))
        if not os.path.isdir(top):
           raise ValueError("Expecting existenting directory, got %s" % top)

        paths = []
        for dirpath, dirnames, filenames in os.walk(top):
            for f in filenames:
                if not f.endswith(".nc"): continue
                paths.append(os.path.join(dirpath, f))

        options.paths = paths
        print("Will analyze %s netcdf file(s)" % len(paths))
        if not paths: sys.exit(1)
        return paths


def main():
    def str_examples():
        return """\
Usage example:

    ncheck.py validate [FILES]     => Validate list of nc files.
    ncheck.py groups [FILES]       => Show groups in list of nc files.

    FILES could be either a single file or a list of files.
 
    TIP: To profile the python code add `prof` before the command.
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
    p_validate = subparsers.add_parser('validate', parents=[paths_selector_parser], help=validate.__doc__)

    p_vars = subparsers.add_parser('ncvars', parents=[paths_selector_parser], help=ncvars.__doc__)

    p_groups = subparsers.add_parser('groups', parents=[paths_selector_parser], help=groups.__doc__)

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

    paths = get_paths(options)

    # Dispatch
    return globals()[options.command](options, paths)


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
