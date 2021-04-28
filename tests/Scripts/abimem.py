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


from pprint import pprint

# We don't install with setup.py hence we have to add the directory [...]/abinit/tests to $PYTHONPATH
# TODO: Use Absolute imports and rename tests --> abitests to
# avoid possible conflicts with the packages in PYTHONPATH
# monty installs the subpackage paths and this breaks the import below
pack_dir, x = os.path.split(os.path.abspath(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)

try:
    from tests.pymods.memprof import AbimemFile
except ImportError:
    print("Cannot find tests.pymods.memprof module in PYTHONPATH.")
    print("Very likely, you have copied the scripts from ~abinit/test/Scripts to another directory")
    print("In this case, you have to add the abinit directory to the python path by executing:\n")
    print("     export PYTHONPATH=~/abinit_directory\n")
    print("before running the script.")
    print("Note that there's no need to change PYTHONPATH if you invoke the script as ~abinit/tests/Scripts/abimem.py")
    raise


from tests.pymods.termcolor import cprint


def summarize(options):
    """Print basic info to terminal."""
    for memfile in options.memfiles:
        print(memfile)
    return 0


def leaks(options):
    """Find possible memory leaks."""
    retcode = 0
    for memfile in options.memfiles:
        retcode += memfile.find_memleaks()
    return retcode


def small(options):
    """Find small allocations."""
    retcode = 0
    for memfile in options.memfiles:
        smalls = memfile.find_small_allocs() #nbytes=options.nbytes)
    return retcode


def large(options):
    """Find large allocations."""
    retcode = 0
    for memfile in options.memfiles:
        larges = memfile.find_large_allocs() #nbytes=options.nbytes)
    return retcode


def intense(options):
    """Find routines with intensive allocations."""
    retcode = 0
    for memfile in options.memfiles:
        intens = memfile.get_intense_dataframe()
        print(intens)
        #df = memfile.get_hotspots_dataframe()
        #print(df)
    return retcode


def peaks(options):
    """Find memory peaks."""
    retcode = 0
    maxlen = 20
    for memfile in options.memfiles:
        #print(memfile.get_peaks(maxlen=maxlen, as_dataframe=True))
        for i, peak in enumerate(memfile.get_peaks(maxlen=maxlen)):
            print("[%d] %s" % (i, peak))
        memfile.plot_peaks(maxlen=maxlen, title="Memory peaks")

    return retcode


def weird(options):
    """Find weird allocations. i.e. pointers <= 0."""
    retcode = 0
    for memfile in options.memfiles:
        memfile.find_weird_ptrs()
    return retcode


def zerosized(options):
    """Find zero-sized allocations."""
    retcode = 0
    for memfile in options.memfiles:
        elist = memfile.find_zerosized()

        if elist:
            print("Found %d zero-sized entries:" % len(elist))
            pprint(elist)
        else:
            print("No zero-sized found")


    return retcode


def plot(options):
    """Plot data with matplotlib."""
    for memfile in options.memfiles:
        memfile.expose()
    return 0


def panel(options):
    """
    Open GUI in web browser, requires panel package
    """
    try:
        import panel  # noqa: F401
    except ImportError as exc:
        cprint("Use `conda install panel` or `pip install panel` to install the python package.", "red")
        raise exc

    import matplotlib
    matplotlib.use("Agg")

    for memfile in options.memfiles:
       memfile.get_panel().show()
    return 0


def ipython(options):
    """Open file in ipython shell."""
    # Start ipython shell with namespace
    # Use embed because I don't know how to show a header with start_ipython.
    import IPython
    abifile = options.memfiles[0]
    IPython.embed(header="""
The Abinit file is bound to the `parsers` variable.
Use `abifile.<TAB>` to list available methods.
Use e.g. `abifile.plot?` to access docstring and `abifile.plot??` to visualize source.
Use `print(abifile)` to print the object.
""")
    return 0


def get_parser(with_epilog=False):

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    #copts_parser.add_argument('paths', nargs="+", help="List of ABIMEM files.")
    copts_parser.add_argument('paths', nargs="+", help="List of ABIMEM files or directory containing ABIMEM files.")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='Verbose, can be supplied multiple times to increase verbosity.')
    copts_parser.add_argument('-sns', "--seaborn", const="paper", default=None, action='store', nargs='?', type=str,
        help='Use seaborn settings. Accept value defining context in ("paper", "notebook", "talk", "poster"). Default: paper')
    copts_parser.add_argument('-mpl', "--mpl-backend", default=None,
        help=("Set matplotlib interactive backend. "
              "Possible values: GTKAgg, GTK3Agg, GTK, GTKCairo, GTK3Cairo, WXAgg, WX, TkAgg, Qt4Agg, Qt5Agg, macosx."
              "See also: https://matplotlib.org/faq/usage_faq.html#what-is-a-backend."))
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG.")

    # Parent parser for commands supporting (ipython/jupyter)
    #ipy_parser = argparse.ArgumentParser(add_help=False)
    #ipy_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    #ipy_parser.add_argument('--foreground', action='store_true', default=False,
    #    help="Run jupyter notebook in the foreground.")
    #ipy_parser.add_argument('-ipy', '--ipython', default=False, action="store_true", help='Invoke ipython terminal.')

    # Parent parser for commands supporting (jupyter notebooks)
    #nb_parser = argparse.ArgumentParser(add_help=False)
    #nb_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    #nb_parser.add_argument('--foreground', action='store_true', default=False,
    #    help="Run jupyter notebook in the foreground.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Build Subparsers for commands
    p_summarize = subparsers.add_parser('summarize', parents=[copts_parser], help=summarize.__doc__)

    # Build Subparsers for commands
    p_leaks = subparsers.add_parser('leaks', parents=[copts_parser], help=leaks.__doc__)

    # Subparser for small
    p_small = subparsers.add_parser('small', parents=[copts_parser], help=small.__doc__)

    # Subparser for large
    p_large = subparsers.add_parser('large', parents=[copts_parser], help=large.__doc__)

    # Subparser for intense
    p_intense = subparsers.add_parser('intense', parents=[copts_parser], help=intense.__doc__)

    # Subparser for peaks command.
    p_peaks = subparsers.add_parser('peaks', parents=[copts_parser], help=peaks.__doc__)

    # Subparser for weird command.
    p_weird = subparsers.add_parser('weird', parents=[copts_parser], help=weird.__doc__)

    # Subparser for zerosized command.
    p_zerosized = subparsers.add_parser('zerosized', parents=[copts_parser], help=zerosized.__doc__)

    # Subparser for plot command.
    p_plot = subparsers.add_parser('plot', parents=[copts_parser], help=plot.__doc__)

    # Subparser for panel command.
    p_panel = subparsers.add_parser('panel', parents=[copts_parser], help=panel.__doc__)

    # Subparser for ipython command.
    p_ipython = subparsers.add_parser('ipython', parents=[copts_parser], help=ipython.__doc__)

    return parser


def get_epilog():
        return """\
Usage example:

    abimem.py leaks [FILES]       => Find possible memory leaks in FILE(s)
    abimem.py small [FILES]       => Find small memory allocations in FILE(s)
    abimem.py large [FILES]       => Find large memory allocations in FILE(s)
    abimem.py intense [FILES]     => Find periods of intense memory allocation in FILE(s)
    abimem.py peaks [FILES]       => Find peaks in memory allocation in FILE(s)
    abimem.py plot [FILES]        => Plot memory allocations in FILE(s) with matplotlib
    abiopen.py FILE               => Open file in ipython shell.

    FILES could be either a list of files or a single directory containing abimem_ran.mocc files.

    TIP: To profile the python code add `prof` before the command e.g.
         abimem.py prof leaks [FILES]
"""


def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

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

    if options.seaborn:
        # Use seaborn settings.
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)

    if os.path.isfile(options.paths[0]):
        options.memfiles = [AbimemFile(path) for path in options.paths]

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
        options.memfiles = [AbimemFile(path) for path in paths]

    # Dispatch
    return globals()[options.command](options)


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
