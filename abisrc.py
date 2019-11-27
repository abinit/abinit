#!/usr/bin/env python
# coding: utf-8
"""
This script analyzes the Abinit source tree and generates the dependency graph
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import argparse

from fkiss import termcolor
from fkiss.tools import print_dataframe
from fkiss.project import FortranFile, AbinitProject
from fkiss.termcolor import cprint


def get_epilog():
    return """\

Usage example:

################
# Documentation
################

  ./abisrc.py print src/41_geometry/m_crystal.F90   ==> Print info about file.
  ./abisrc.py print src/41_geometry                 ==> Print info about directory.
  ./abisrc.py print crystal_init                    ==> Print info about public procedure.
  ./abisrc.py print m_crystal                       ==> Print info about module.
  ./abisrc.py print crystal_t                       ==> Print info about datatype.
  ./abisrc.py print xmpi_sum                        ==> Print info about interface.
  ./abisrc.py interface xmpi_sum                    ==> Print Fortran interfaces.
  ./abisrc.py interface                             ==> Print ALL Fortran interfaces.

  ./abisrc.py parse src/41_geometry/m_crystal.F90   ==> Parse file, print results (debugging tool)

#################
# Graphviz graphs
#################

  ./abisrc.py graph src/41_geometry/m_crystal.F90   ==> Plot dependency graph for module.
  ./abisrc.py graph src/41_geometry                 ==> Plot dependency graph for directory.
  ./abisrc.py graph fourdp                          ==> Plot dependency graph for public procedure.

Requires graphviz and python graphviz: https://graphviz.readthedocs.io/en/stable/manual.html

#############
# Developers
#############

  ./abisrc.py makemake           ==> Generate files required by the build system.
  ./abisrc.py validate           ==> Validate files in project.
  ./abisrc.py touch              ==> Touch all files that have been changed + parents.
                                     so that make can recompile all the relevant files.
                                     Useful when changing API/ABI.
  ./abisrc.py pedit fourdp       ==> Call $EDITOR to edit all the parents of the fourdp routine.
  ./abisrc.py stats              ==> Print pandas Dataframe with stats about project/dir/file.
                                     Use -c to copy to system clipboard
  ./abisrc.py orphans            ==> Show orphans.
  ./abisrc.py ipython            ==> Open project in ipython terminal.
  ./abisrc.py cpp                ==> List CPP options.
  ./abisrc.py panel              ==> Generate dashboard in browser.


NB: graph and panel commands require external libraries. Use:

    pip install -r abisrc_requirements.txt
or

    conda install --file abisrc_requirements.txt
"""
  #./abisrc.py robodoc            ==> Generate robodoc files.

def get_parser():
    """Build and return parser object."""

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity.')
    copts_parser.add_argument('-r', '--regenerate', default=False, action="store_true",
        help='Parse files, generate new pickle file.')
    #copts_parser.add_argument('--loglevel', default="ERROR", type=str,
    #    help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG.")
    copts_parser.add_argument('--no-colors', default=False, action="store_true", help='Disable ASCII colors')
    copts_parser.add_argument('-j', "--jobs", type=int, default=2,  help="Number of python processes to use.")

    # Parent parser for commands that operating on pandas dataframes
    pandas_parser = argparse.ArgumentParser(add_help=False)
    pandas_parser.add_argument("-c", '--clipboard', default=False, action="store_true",
            help="Copy dataframe to the system clipboard. This can be pasted into Excel, for example")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help.', description="Valid subcommands")

    # Subparser for parse.
    p_parse = subparsers.add_parser('parse', parents=[copts_parser], help="Parse file.")
    p_parse.add_argument("what", type=str, help="File to parse.")

    # Subparser for makemake.
    p_makemake = subparsers.add_parser('makemake', parents=[copts_parser],
        help="Generate configuration files required by the build system.")

    # Subparser for touch.
    p_touch = subparsers.add_parser('touch', parents=[copts_parser],
        help="Change timestamp of all files that.")
    #p_touch.add_argument("what_list", nargs="*", default=None, help="List of files or empty for auto.")

    # Subparser for print.
    p_print = subparsers.add_parser('print', parents=[copts_parser],
        help="Show children of module/procedure.")
    p_print.add_argument("what", nargs="?", default=None, help="File or procedure name")

    # Subparser for interface.
    p_interfaces = subparsers.add_parser('interface', parents=[copts_parser], help="Print interface.")
    p_interfaces.add_argument("what", nargs="?", default=None, help="Name of the interface.")

    # Subparser for pedit.
    p_edit = subparsers.add_parser("pedit", parents=[copts_parser],
        help="Edit parents of public procedure or module.")
    p_edit.add_argument("what", help="File or procedure name.")

    # Subparser for canimove.
    #p_canimove = subparsers.add_parser("canimove", parents=[copts_parser],
    #    help="Check whether file or directory can be moved ")
    #p_canimove.add_argument("what", help="File or procedure name.")
    #p_canimove.add_argument("dest_level", type=int, help="Destination level")

    # notebook option
    p_notebook = subparsers.add_parser("notebook", parents=[copts_parser], help="Analyze project in jupyter notebook")
    p_notebook.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")

    # panel option
    p_panel = subparsers.add_parser("panel", parents=[copts_parser],
                                    help="Analyze project in dashboard. Requires panel package.")

    # Subparser for stats.
    p_stats = subparsers.add_parser("stats", parents=[copts_parser, pandas_parser],
        help="Print statistics about file or directory or project (if no argument is given). Use -c to copy to system clipboard")
    p_stats.add_argument("what", nargs="?", default=None, help="File, directory or empty for project.")

    # Subparser for graph.
    p_graph = subparsers.add_parser('graph', parents=[copts_parser],
        help="Draw dependency graph flow with graphviz package. See https://graphviz.readthedocs.io/.")
    p_graph.add_argument("-e", "--engine", type=str, default="automatic",
        help=("graphviz engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']. "
              "Default: automatic i.e. the engine is automatically selected. See http://www.graphviz.org/pdf/dot.1.pdf "
              "Use `conda install python-graphviz` or `pip install graphviz` to install the python package."))
    #p_graph.add_argument("-d", '--dirtree', default=False, action="store_true",
    #    help='Visualize files and directories in workdir instead of tasks/works.')
    p_graph.add_argument("what", nargs="?", default=None, help="File of directory to visualize.")

    # Subparser for validate.
    p_validate = subparsers.add_parser('validate', parents=[copts_parser], help="Validate source tree.")

    p_abirules = subparsers.add_parser('abirules', parents=[copts_parser], help="Check abirules (still under development).")
    p_abirules.add_argument("what", nargs="?", default=None, help="File, directory or empty for project.")

    p_orphans = subparsers.add_parser('orphans', parents=[copts_parser], help="Print orphans.")
    p_ipython = subparsers.add_parser('ipython', parents=[copts_parser], help="Open project in ipython terminal.")

    p_dtype = subparsers.add_parser('dtype', parents=[copts_parser], help="Generate code for datatype.")
    p_dtype.add_argument("what", nargs="+", default=None,
            help="Name of datatype followed by the kind of routine that should be generated.")

    p_master = subparsers.add_parser('master', parents=[copts_parser], help="How to become a great programmer.")
    #p_robodoc = subparsers.add_parser('robodoc', parents=[copts_parser], help="Generate robodoc files.")
    p_cpp = subparsers.add_parser('cpp', parents=[copts_parser], help="List CPP options.")

    return parser


def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # This to avoid RecursionError in pickle as we have a highly recursive datastructure.
    #sys.setrecursionlimit(sys.getrecursionlimit() * 3)
    parser = get_parser()

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    if not options.command:
        show_examples_and_exit(error_code=1)

    if options.no_colors:
        # Disable colors
        termcolor.enable(False)

    os.chdir(os.path.dirname(__file__))

    if sys.version_info[0] < 3 and options.jobs > 1:
        cprint("py2.x CANNOT USE jobs > 1. Setting jobs to 1. Use py3k", "yellow")
        options.jobs = 1

    #if options.command == "robodoc":
    #    from fkiss.mkrobodoc_dirs import mkrobodoc_files
    #    return mkrobodoc_files(".")

    if options.command == "cpp":
        from fkiss.list_cpp_options import list_cpp_options
        return list_cpp_options(".")

    if options.command == "parse":
        if options.what == ".":
            proj = AbinitProject(".", processes=options.jobs, verbose=options.verbose)
            print(proj)
        else:
            fort_file = FortranFile.from_path(options.what, macros="abinit",  verbose=options.verbose)
            print(fort_file.to_string(verbose=options.verbose))
            #for mod in fort_file.modules:
            #    for dtype in mod.types:
            #        print("Analyzing", repr(dtype))
            #        dtype.analyze()
        return 0

    elif options.command == "touch":
        # Load old project and touch files that have been changed.
        old_proj = AbinitProject.pickle_load()
        ntouch = old_proj.touch_alldeps(verbose=options.verbose)

        if ntouch:
            print("\nTouched %d files. Need to parse source files again and dump new pickle file\n" % ntouch)
            new_proj = AbinitProject(".", processes=options.jobs, verbose=options.verbose)
            new_proj.pickle_dump()
        else:
            print("\nNo change detected. No need to touch files.\n")

        return 0

    # After this point I operate an AbinitProject instance.
    # Load the object from pickle first and then check if we need to parse the source files again.
    needs_reload = True
    if not options.regenerate and os.path.exists(AbinitProject.get_default_pickle_file()):
        proj = AbinitProject.pickle_load()
        needs_reload = proj.needs_reload()
        if needs_reload:
            cprint("Source tree changed. Need to parse source files again to rebuild dependency graph...", "yellow")

    if needs_reload:
	# Parse the source and save new object.
        proj = AbinitProject(".", processes=options.jobs, verbose=options.verbose)
        proj.pickle_dump()

    #assert "abinit.F90" in proj.fort_files

    if options.command == "makemake":
        retcode = proj.validate(verbose=options.verbose)
        if retcode != 0:
            cprint("validate returned retcode: %s. Aborting now" % retcode, "red")
            return retcode

        #all_mods = proj.find_allmods("dummy_tests.F90")
        #for mod in all_mods: print(mod.basename)

        proj.write_binaries_conf(verbose=options.verbose, dryrun=False)
        proj.write_buildsys_files(verbose=options.verbose, dryrun=False)

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
                return 1

    elif options.command == "interface":
        print(proj.print_interfaces(what=options.what, verbose=options.verbose))

    elif options.command == "graph":
        if options.what is None:
            raise NotImplementedError("graph action requires argument!")
            #graph = proj.get_graphviz(engine=options.engine)

        elif os.path.isdir(options.what):
            graph = proj.get_graphviz_dir(options.what, engine=options.engine)

        elif os.path.isfile(options.what):
            fort_file = proj.fort_files[os.path.basename(options.what)]
            print(fort_file.to_string(verbose=options.verbose))
            graph = fort_file.get_graphviz(engine=options.engine)

        else:
            graph = proj.get_graphviz_pubname(options.what, engine=options.engine)
            if graph is None:
                return 1

        # Visualize graph
        import tempfile
        directory = tempfile.mkdtemp()
        print("Producing source files in:", directory)
        graph.view(directory=directory, cleanup=False)

    elif options.command == "validate":
       return proj.validate(verbose=options.verbose)

    #elif options.command == "canimove":
    #    if os.path.isdir(options.what):
    #        return proj.canimove_dir(options.what, options.dest_level)
    #    elif os.path.isfile(options.what):
    #        return proj.canimove_file(options.what, options.dest_level)
    #    else:
    #        raise TypeError("Requiring directory or file but received %s" % str(options.what))

    elif options.command == "notebook":
        return proj.make_and_open_notebook(foreground=options.foreground)

    elif options.command == "panel":
        proj.get_panel().show() #threaded=True)

    elif options.command == "pedit":
        return proj.pedit(options.what, verbose=options.verbose)

    elif options.command == "stats":
        if options.what is None:
            df = proj.get_stats()
        elif os.path.isdir(options.what):
            df = proj.get_stats_dir(options.what)
        elif os.path.isfile(options.what):
            df = proj.get_stats_file(options.what)
        else:
            raise TypeError("Requiring directory or file or None but received %s" % str(options.what))

        print_dataframe(df)

        if options.clipboard:
            cprint("Copying dataframe to the system clipboard.", "green")
            cprint("This can be pasted into Excel, for example", "green")
            df.to_clipboard()

    elif options.command == "orphans":
        proj.print_orphans(verbose=options.verbose)

    elif options.command == "ipython":
        import IPython
        IPython.embed(header="The Abinit project is bound to the `proj` variable.\n")

    elif options.command == "dtype":
        retcode = 0

        items = proj.get_all_datatypes_and_fortfiles().values()
        for dtype, fort_file in items:
            print("Analyzing:", repr(dtype))
            try:
                dtype.analyze(verbose=0)
                print(dtype)
            except Exception as exc:
                import traceback
                cprint(traceback.format_exc(), "red")

        print("%d errors out of %d" % (retcode, len(items)))

        #dtype = proj.find_datatype(options.what[0])
        #if dtype is None: return 1
        #dtype.analyze()
        #if "free" in options.what[1:]: print(dtype.generate_free_routine())
        #if "copy" in options.what[1:] print(dtype.generate_copy_routine())

    elif options.command == "abirules":
        if not options.what:
            return proj.check_abirules(verbose=options.verbose)

        elif os.path.isdir(options.what):
            retcode = 0
            for fort_file in proj.fort_files_indir(options.what):
                retcode += fort_file.check_abirules(verbose=options.verbose)
            return retcode

        elif os.path.isfile(options.what):
            fort_file = proj.fort_files[os.path.basename(options.what)]
            return fort_file.check_abirules(verbose=options.verbose)

    elif options.command == "master":
        print(proj.master())

    else:
        raise ValueError("Don't know how to handle command: %s" % options.command)

    return 0


if __name__ == "__main__":
    # Check whether we are in profiling mode
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof: sys.argv.pop(1)
    except Exception:
        do_prof = False

    if not do_prof:
        sys.exit(main())
    else:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
