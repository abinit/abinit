#!/usr/bin/env python
# coding: utf-8
"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import argparse

from fkiss.project import FortranFile, AbinitProject

master_foo = """\

Master Foo and the Hardware Designer

On one occasion, as Master Foo was traveling to a conference
with a few of his senior disciples, he was accosted by a hardware designer.

The hardware designer said:
“It is rumored that you are a great programmer. How many lines of code do you write per year?”

Master Foo replied with a question:
“How many square inches of silicon do you lay out per year?”

“Why...we hardware designers never measure our work in that way,” the man said.

“And why not?” Master Foo inquired.

“If we did so,” the hardware designer replied, “we would be tempted to design chips
so large that they cannot be fabricated - and, if they were fabricated,
their overwhelming complexity would make it be impossible to generate proper test vectors for them.”

Master Foo smiled, and bowed to the hardware designer.

In that moment, the hardware designer achieved enlightenment.

From http://www.catb.org/esr/writings/unix-koans/
"""

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
  abisrc.py graph fourdp                      => Plot dependency graph for function.

#############
# Developers
#############

  abisrc.py makemake              => Generate files required by build system.
  abisrc.py abirules              =>
  abisrc.py master                => Master the Abinit source tree.
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

    # Subparser for parse.
    p_parse = subparsers.add_parser('parse', parents=[copts_parser], help="Parse file.")
    p_parse.add_argument("what", type=str, help="File to parse.")

    # Subparser for makemake.
    p_makemake = subparsers.add_parser('makemake', parents=[copts_parser],
        help="Generate configuration files required by the build system.")

    # Subparser for touch.
    p_touch = subparsers.add_parser('touch', parents=[copts_parser],
        help="Change timestamp of all files that.")
    p_touch.add_argument("what_list", nargs="*", default=None, help="List of files or empty for auto.")

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

    # Subparser for validate.
    p_validate = subparsers.add_parser('validate', parents=[copts_parser],
        help="Validate source tree.")

    p_master = subparsers.add_parser('master', parents=[copts_parser], help="Master.")

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
        print(fort_file.to_string(verbose=options.verbose))
        return 0

    # After this point I operate an AbinitProject instance.
    # Load the object from pickle first and then check if we need to parse the source again.
    needs_reload = True
    if not options.regenerate and os.path.exists(AbinitProject.DEFAULT_PICKLE_FILE):
        proj = AbinitProject.pickle_load()
        needs_reload = proj.needs_reload()
        if needs_reload:
            print("Source tree changed. Will parse source files to rebuild dependency graph...")

    if needs_reload:
	# Parse the source and save new object.
        proj = AbinitProject(".", verbose=options.verbose)
        proj.pickle_dump()

    #assert "abinit.F90" in proj.fort_files

    if options.command == "makemake":
        proj.write_buildsys_files()

    elif options.command == "touch":
        proj.touch_alldeps(what_list=options.what_list, verbose=options.verbose)

    elif options.command == "print":
        if options.what is None:
            print(proj.to_string(verbose=options.verbose))
        elif os.path.isdir(options.what):
            proj.print_dir(options.what, verbose=options.verbose)
        elif os.path.isfile(options.what):
            fort_file = proj.fort_files[os.path.basename(options.what)]
            print(fort_file.to_string(verbose=options.verbose))
        else:
            obj = proj.find_public_entity(options.what)
            if obj is not None:
                print(obj.to_string(verbose=options.verbose))
            else:
                print("Cannot find public entity `%s`" % str(options.what))
                return 1

    elif options.command == "graph":
        if options.what is None:
            graph = proj.get_graphviz(engine=options.engine)
        elif os.path.isdir(options.what):
            graph = proj.get_graphviz_dir(options.what, engine=options.engine)
        elif os.path.isfile(options.what):
            fort_file = proj.fort_files[os.path.basename(options.what)]
            print(fort_file.to_string(verbose=options.verbose))
            graph = fort_file.get_graphviz(engine=options.engine)
        else:
            graph = proj.get_graphviz_dir(options.what, engine=options.engine)
            if graph is None:
                return 1

        # Visualize graph
        import tempfile
        directory = tempfile.mkdtemp()
        print("Producing source files in:", directory)
        graph.view(directory=directory, cleanup=False)

    elif options.command == "validate":
       proj.validate(verbose=options.verbose)

    elif options.command == "master":
       print(master_foo)

    #elif options.command == "plot":
    #    if os.path.isdir(options.what):
    #        proj.plot_dir(options.what)
    #    else:
    #        key = os.path.basename(options.what)
    #        proj.plot_file(key)

    #elif options.command == "canimove":
    #   return proj.canimove(src, dest)

    elif options.command in ("edit_parents", "edit_children"):
        relation = options.command.split("_")[1]
        proj.edit_connections(options.what, relation)

    elif options.command == "stats":
        if options.what is None:
            proj.stats(verbose=options.verbose)
        elif os.path.isdir(options.what):
            proj.stats_dir(options.what, verbose=options.verbose)
        elif os.path.isfile(options.what):
            proj.stats_file(options.what, verbose=options.verbose)
        else:
            raise TypeError("Don't know how to produce stats for %s" % str(options.what))

    else:
        raise ValueError("Don't know how to handle command: %s" % options.command)

    return 0


if __name__ == "__main__":
    sys.exit(main())