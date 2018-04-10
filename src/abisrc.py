#!/usr/bin/env python
# coding: utf-8
"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import argparse

from fkiss.project import FortranFile, AbinitProject


def get_epilog():
    return """\

Usage example:

  abisrc.py parse 41_geometry/m_crystal.F90   ==> Parse file, print results

  abisrc.py print 41_geometry/m_crystal.F90   ==> Print info about file
  abisrc.py print 41_geometry                 ==> Print info about directory
  abisrc.py print crystal_init                ==> Print info about file

########
# Graphs
########

  abisrc.py graph 41_geometry/m_crystal.F90   => Plot dependency graph for module.
  abisrc.py graph fourdp                      => Plot dependecy graph for function.

#############
# Developers
#############

  abisrc.py makemake              => Generate files required by build system.
"""

def get_parser():
    """Build and return parser object."""

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity.')
    copts_parser.add_argument('-r', '--regenerate', default=False, action="store_true",
        help='Parse files, generate new pickle file')
    #copts_parser.add_argument('--no-colors', default=False, action="store_true", help='Disable ASCII colors.')
    #copts_parser.add_argument('--no-logo', default=False, action="store_true", help='Disable AbiPy logo.')
    #copts_parser.add_argument('--loglevel', default="ERROR", type=str,
    #    help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument('flowdir', nargs="?", help=("File or directory containing the ABINIT flow/work/task. "
    #                                                "If not given, the flow in the current workdir is selected."))
    #parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for print.
    p_parse = subparsers.add_parser('parse', parents=[copts_parser],
        help="Parse file.")
    p_parse.add_argument("what", type=str, help="File to parse.")

    # Subparser for print.
    p_print = subparsers.add_parser('print', parents=[copts_parser],
        help="Show children of module/procedure.")
    p_print.add_argument("what", nargs="?", default=None, help="File or procedure name")

    # Subparser for graph.
    p_graph = subparsers.add_parser('graph', parents=[copts_parser],
        help=("Draw flow and node dependencies with graphviz package. Accept (FLOWDIR|WORKDIR|TASKDIR). "
             "See https://graphviz.readthedocs.io/."))
    p_graph.add_argument("-e", "--engine", type=str, default="automatic",
        help=("graphviz engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']. "
            "Default: automatic i.e. the engine is automatically selected. See http://www.graphviz.org/pdf/dot.1.pdf "
            "Use `conda install python-graphviz` or `pip install graphviz` to install the python package"))
    #p_graph.add_argument("-d", '--dirtree', default=False, action="store_true",
    #    help='Visualize files and directories in workdir instead of tasks/works.')
    p_graph.add_argument("what", nargs="?", default=None, help="File of directory to visualize")

    return parser


def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # This to avoid RecursionError in pickle as we have a highly recursive datastructure.
    sys.setrecursionlimit(sys.getrecursionlimit() * 3)

    parser = get_parser()

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    if not options.command:
        show_examples_and_exit(error_code=1)

    if options.command == "parse":
        fort_file = FortranFile.from_path(options.what, options.verbose)
        print(fort_file)
        print(fort_file.stree())
        return 0

    needs_reload = True
    if not options.regenerate and os.path.exists(AbinitProject.DEFAULT_PICKLE_FILE):
        proj = AbinitProject.pickle_load()
        needs_reload = proj.needs_reload()
        if needs_reload:
            print("Source tree changed. Will parse source files to rebuild dependency graph...")

    if needs_reload:
        proj = AbinitProject(".", verbose=options.verbose)
        proj.pickle_dump()

    # FIXME
    #node = proj.fort_files["abinit"]
    #print("abinit parents", node.parents)

    if options.command == "makemake":
        proj.write_buildsys_files()

    elif options.command == "print":

        if options.what is None:
            print(proj.to_string(verbose=options.verbose))

        elif os.path.isdir(options.what):
            proj.print_dir(options.what, verbose=options.verbose)

        elif os.path.isfile(options.what):
            fort_file = proj.fort_files[os.path.basename(options.what)]
            print(fort_file.to_string(verbose=options.verbose))
            #print(fort_file.stree())
            #print(fort_file.modules[0].contains[getng]

        else:
            obj = proj.find_public_entity(options.what)
            if obj is not None:
                print(obj.to_string(verbose=options.verbose))
            else:
                print("Cannot find public entity `%s`" % str(options.what))
                return 1

    elif options.command == "graph":

        if options.what is None:
            raise NotImplementedError()
            #print(proj.to_string(verbose=options.verbose))
            graph = proj.get_graphviz(engine=options.engine)

        elif os.path.isdir(options.what):
            graph = proj.get_graphviz_dir(options.what, engine=options.engine)

        elif os.path.isfile(options.what):
            fort_file = proj.fort_files[os.path.basename(options.what)]
            print(fort_file.to_string(verbose=options.verbose))
            graph = fort_file.get_graphviz(engine=options.engine)

        else:
            obj = proj.find_public_entity(options.what)
            if obj is not None:
                #print(obj.to_string(verbose=options.verbose))
                graph = obj.get_graphviz(engine=options.engine)
            else:
                print("Cannot find public entity `%s`" % str(options.what))
                return 1

        import tempfile
        directory = tempfile.mkdtemp()
        print("Producing source files in:", directory)
        graph.view(directory=directory, cleanup=False)

    else:
        raise ValueError("Don't know how to handle command: %s" % options.command)

    #elif options.command == "plot":
    #    if os.path.isdir(options.what):
    #        proj.plot_dir(options.what)
    #    else:
    #        key = os.path.basename(options.what)
    #        proj.plot_file(key)

    #elif options.command == "canimove":
    #   return proj.canimove(src, dest)

    #elif options.command in ("edit_parents", "edit_children"):
    #   obj = proj.find_public_entity(name)
    #   if obj is None
    #       print("Cannot find public entity `%s`" % str(name))
    #       return 1

    #   # Find files with procedures.
    #   relation = options.command.split("_")[1]
    #   plist = dict(parents=obj.parents, children=obj.children)[relation]
    #   paths = sorted(set(p.path for p in plist))
    #   from pymods.tools import Editor
    #   return Editor().edit_files(paths, ask_for_exit=True)

    #elif options.command == "abirules":
    #   proj.abirules()

    # Dispatch
    #return globals()["abisrc_" + options.command](options)

if __name__ == "__main__":
    sys.exit(main())
