# coding: utf-8
"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import sys
import io
import re
import time
import shutil
import pickle
import difflib

from collections import OrderedDict, defaultdict
from textwrap import TextWrapper
from pprint import pprint, pformat
from .parser import FortranKissParser
from .tools import lazy_property, NotebookWriter
from .termcolor import cprint


EXTERNAL_MODS = {
    # Intrinsics
    "iso_fortran_env",
    "iso_c_binding",
    "ieee_arithmetic",
    "ieee_exceptions",
    "ieee_features",
    # Modules provided by compilers.
    "f90_unix_proc",
    "f90_unix_dir",
    "ifcore",
    # MPI modules.
    "mpi",
    "mpi_f08",
    # External libraries.
    "openacc",
    "omp_lib",
    "mkl_dfti",
    "netcdf",
    "etsf_io_low_level",
    "etsf_io",
    "plasma",
    "elpa",
    "elpa1",
    "fox_sax",
    "m_libpaw_libxc_funcs",
    "m_psml",
    "m_psml_api",
    # Bigdft modules.
    "yaml_output",
    "bigdft_api",
    "module_base",
    "module_types",
    "module_xc",
    "poisson_solver",
    "dynamic_memory",
    # Abinit-specific modules.
    #"m_build_info",
    #"m_optim_dumper",
    "libxc_functionals",
}


class FortranFile(object):
    """
    Base class for files containing Fortran source code.
    A FortranFile can have modules, programs, subroutines and functions.

    .. attributes:

        modules:
        programs:
        subroutines
        functions
    """

    @classmethod
    def from_path(cls, path, macros=None, verbose=0):
        if macros == "abinit": macros = AbinitProject.MACROS
        p = FortranKissParser(macros=macros, verbose=verbose).parse_file(path)

        new = cls(path)
        new.all_includes = p.all_includes
        new.programs = p.programs
        new.modules = p.modules
        new.subroutines = p.subroutines
        new.functions = p.functions

        return new

    def __init__(self, path):
        # A file can contain multiples modules but not procedures outside modules
        # A file with program cannot contain other procedures/modules outside program
        # module/program must be followed by end [module|program] to facilitate parsing.
        self.path = os.path.abspath(path)
        # TODO name --> basename and REMOVE self.name
        self.name = os.path.basename(self.path)
        self.basename = self.name
        self.dirname = os.path.dirname(self.path)
        if self.dirname == "src" and os.path.join("shared", "libpaw") in self.path:
            self.dirname = "39_libpaw"

        # Save initial stat values, used to understand if reload is needed.
        self.stat = os.stat(self.path)

        self.programs, self.modules, self.subroutines, self.functions = [], [], [], []

        #self.uses, self.includes = [], []
        self.all_used_mods = []
        self.all_usedby_mods = []
        #self.num_f90lines, self.num_doclines = 0, 0

    def __eq__(self, other):
        """Use path to compare for equality and compute hash value."""
        if other is None: return False
        if not isinstance(other, self.__class__): return False
        return self.path == other.path

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.path)

    def stree(self):
        """
        Return string with textual representation of the tree.
        """
        lines = [repr(self)]; app = lines.append
        for a in ["programs", "modules", "subroutines", "functions"]:
            for p in getattr(self, a):
                app(p.stree())

        return "\n".join(lines)

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.path)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, with_stats=True, width=90):
        """
        String representation with verbosity level `verbose`.
        Text is wrapped at `width` columns.
        """
        w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)
        lines = []; app = lines.append
        app("%s:\n\t%s\n" % (self.__class__.__name__, os.path.relpath(self.path)))
        app("Use modules:\n%s\n" % w.fill(", ".join(mod.name for mod in self.all_used_mods)))
        app("Use dir_levels:\n%s\n" % w.fill(", ".join(map(str, sorted(set(mod.dirlevel for mod in self.all_used_mods))))))
        app("Includes:\n%s\n" % w.fill(", ".join(self.all_includes)))
        app("Used by modules:\n%s\n" % w.fill(", ".join(mod.name for mod in self.all_usedby_mods)))
        app("Used by dir_levels:\n%s\n" % w.fill(", ".join(map(str, sorted(set(mod.dirlevel for mod in self.all_usedby_mods))))))

        if verbose:
            for a in ["programs", "modules", "subroutines", "functions"]:
                plist = getattr(self, a)
                if not plist: continue
                app("Fortran file contains %d %s(s)\n" % (len(plist), a[:-1]))
                for p in plist:
                    app(p.to_string(verbose=verbose, width=width))

        if with_stats:
            df = self.get_stats()
            app(df.to_string())

        if verbose > 1:
            app(self.stree())

        return "\n".join(lines)

    #@lazy_property
    #def all_used_mods(self):
    #    all_mods = []
    #    for a in ["programs", "modules", "subroutines", "functions"]:
    #        for p in getattr(self, a):
    #            all_mods.extend(p.local_uses)
    #    return sorted(set(all_mods))

    #@lazy_property
    #def all_usedby_mods(self):
    #    all_mods = []
    #    for a in ["programs", "modules", "subroutines", "functions"]:
    #        for p in getattr(self, a):
    #            all_mods.extend(p.usedby_mods)
    #    return sorted(set(all_mods))

    @lazy_property
    def all_uses(self):
        """
        Full list of modules (string) used by the Fortran file
        Includes modules used by the procedures declared in the file.
        """
        all_uses = []
        for a in ["programs", "modules", "subroutines", "functions"]:
            for p in getattr(self, a):
                all_uses.extend(p.uses)
        return sorted(set(all_uses))

    @lazy_property
    def dirlevel(self):
        """Integer given the directory level in the Abinit hierarchy."""
        # 72_response --> 72
        dname = os.path.dirname(self.path)
        try:
            return int(os.path.basename(dname).split("_")[0])
        except Exception:
            head, tail = os.path.split(dname)
            if os.path.join("shared", "libpaw") in head:
                #print("Setting dirlevel to 30 for shared/libpaw/src since head:", head)
                return 39
            cprint("Cannot detect dirlevel in: %s" % self.path, "red")
            raise

    @lazy_property
    def min_dirlevel(self):
        """Minimum directory level used by this file."""
        return max(mod.dirlevel for mod in self.all_used_mods) if self.all_used_mods else 999

    @lazy_property
    def max_dirlevel(self):
        """Maximum directory level used by this file."""
        return min(mod.dirlevel for mod in self.all_usedby_mods) if self.all_usedby_mods else 0

    @lazy_property
    def all_public_procedures(self):
        """Dictionary name --> public_procedure."""
        pubs = OrderedDict()
        for a in ["modules", "programs", "subroutines", "functions"]:
            for container in getattr(self, a):
                for p in container.public_procedures:
                    pubs[p.name] = p
        return pubs

    def find_public_entity(self, name, all_names=None):
        """
        Find and return the public procedure or datatype with `name`.
        Return None if not found.
        """
        for a in ["modules", "programs", "subroutines", "functions"]:
            for p in getattr(self, a):
                for e in p.public_procedures:
                    if e.name == name: return e
                    if all_names is not None: all_names.append(e.name)
                if a == "modules":
                    # Try also in datatypes and interfaces.
                    for e in p.types:
                        if e.name == name: return e
                        if all_names is not None: all_names.append(e.name)
                    for e in p.interfaces:
                        if e.name == name: return e
                        if all_names is not None: all_names.append(e.name)

        return None

    def find_datatype(self, name, all_names=None):
        """
        Find and return the public datatype with `name`. Return None if not found.
        """
        for p in getattr(self, "modules"):
            for e in p.types:
                if e.name == name: return e
                if all_names is not None: all_names.append(e.name)

        return None

    def get_stats(self, as_dict=False):
        """
        Return dictionary with FortFile stats or pandas dataframe if as_dict.
        """
        d = OrderedDict([
            ("use", len(self.all_used_mods)),
            ("usedby", len(self.all_usedby_mods)),
            ("min_dirlevel", self.min_dirlevel),
            ("this_dirlevel", self.dirlevel),
            ("max_dirlevel", self.max_dirlevel),
            ("includes", len(self.all_includes)),
            #("subroutines", len(self.subroutines)),
            #("functions", len(self.functions)),
            #("modules", len(self.modules)),
            #("programs", len(self.programs)),
            #("code_lines", self.num_f90lines),
            #("doc_lines", self.num_doclines),
        ])

        if as_dict: return d
        import pandas as pd
        return pd.DataFrame([d], index=[self.basename], columns=d.keys() if d else None)

    def check_abirules(self, verbose=0):
        retcode = 0
        for a in ["modules", "programs", "subroutines", "functions"]:
            for p in getattr(self, a):
                retcode += p.check_abirules(verbose=verbose)

        return retcode

    def get_graphviz(self, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate dependency graph for this Fortran file in the DOT language
        (only all_used_mods and children of this file).

        Args:
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph
        fg = Digraph(name="Fortran file: %s" % self.name,
            graph_attr=graph_attr, node_attr=node_attr, edge_attr=edge_attr,
            engine="dot" if engine == "automatic" else engine,
        )

        # Set graph attributes.
        fg.attr(label=self.name, rankdir="LR", pagedir="BL") # constraint="false", pack="true", packMode="clust")
        fg.node_attr.update(color='lightblue2', style='filled')

        # Build clusters representing directories.
        all_nodes = sorted([self] + self.all_used_mods + self.all_usedby_mods, key=lambda o: o.dirname)
        dir2nodes = defaultdict(list)
        for node in all_nodes:
            dir2nodes[node.dirname].append(node)

        for dirname, nodes in dir2nodes.items():
            # Create cluster.
            cluster_name = "cluster%s" % dirname
            with fg.subgraph(name=cluster_name) as wg:
                wg.attr(label=os.path.basename(dirname))
                for node in nodes:
                    wg.node(node.name)

        for mod in self.all_used_mods:
            fg.edge(self.name, mod.name)

        for child in self.all_usedby_mods:
            fg.edge(child.name, self.name)

        return fg


class AbinitProject(NotebookWriter):
    """
    This object defines the main entry point for client code.
    It contains a dictionary (fort_files) mapping the basename of the Fortran files to FortranFile instances.
    Provides methods to traverse the DAG, print information about the source code
    and generate the configuration files required by the build system.
    """

    # m_optim_dumper and m_build_info are generated by makemake so its source is not available in the src tree
    # if one uses build directories
    # Note from YP (20190514): The .F90.in are however available and the
    # rest of the build system knows about them through the abinit.src files
    # and the ABI_SRC_BLT variable declared there.
    BUILDSYS_POSTPROCESS_FILES = {"m_build_info.F90", "m_optim_dumper.F90"}

    # Macros used to import modules in abinit libraries.
    # Must be consistent with CPP version. See incs/abi_common.h
    MACROS = {
        # Abinit MACROS.
        "ABI_ASYNC": ",asynchronous",
        "ABI_PRIVATE": ",private",
        "ABI_PROTECTED": ",protected",
        "ABI_CONTIGUOUS": "contigous,",
        # Libpaw.
        "USE_DEFS": "use defs_basis",
        "USE_MPI_WRAPPERS": "use m_xmpi",
        "USE_MSG_HANDLING": "use m_errors, only : msg_hndl, netcdf_check; use m_abicore",
        "USE_MEMORY_PROFILING": "use m_profiling_abi",
        # Libtetra.
    }

    @classmethod
    def get_default_pickle_file(cls):
        """
        Return string with the default name of the pickle file used to save the object to disk.
        The string contains the python major version to avoid possible incompatibilites
        in the pickle protocol that may occur when reading a pickle file produced by another python interpreter.
        """
        return "_project_py%s.pickle" % sys.version_info[0]

    @classmethod
    def pickle_load(cls, filepath=None):
        """
        Reconstruct object from pickle file. Default is used if filepath is None.
        """
        if filepath is None: filepath = cls.get_default_pickle_file()
        with open(filepath, "rb") as fh:
            return pickle.load(fh)

    def __init__(self, top, processes=1, verbose=0):
        # Find directories with abinit.src files inside srcdir
        # and get list of files treated by the build system.
        self.top = os.path.abspath(top)
        #self.srcdir = os.path.abspath(os.path.join(self.top, "src"))
        self.dirpaths = self.get_dirpaths()
        self.verbose = verbose

        # Get source files from abinit.src and build mapping basename --> FortranFile
        start = time.time()
        def filter_fortran(files):
            return [f for f in files if f.endswith(".f") or f.endswith(".F90")]

        import imp
        name2path = OrderedDict()
        for d in self.dirpaths:
            if os.path.basename(d) == "98_main":
                # Special treatment for programs (no abinit.src here)
                for basename in filter_fortran(os.listdir(d)):
                    if basename in name2path:
                        raise RuntimeError("Found two Fortran files with same basename `%s`" % basename)
                    name2path[basename] = os.path.join(d, basename)
            else:
                # Get source files from abinit.src.
                abinit_src = os.path.join(d, "abinit.src")
                mod = imp.load_source(abinit_src, abinit_src)
                for basename in filter_fortran(mod.sources):
                    if basename in name2path:
                        raise RuntimeError("Found two Fortran files with same basename `%s` and other in dir `%s`\n"
                          % (name2path[basename], abinit_src))

                    if basename in self.BUILDSYS_POSTPROCESS_FILES:
                        # Here I need a file that provides the public API and looks like valid Fortran
                        # so that abisrc can parse part of it.
                        name2path[basename] = os.path.join(d, basename + ".in")
                    else:
                        name2path[basename] = os.path.join(d, basename)

        print("Using %d processes to analyze %d directories with %d files" % (processes, len(self.dirpaths), len(name2path)))

        self.fort_files = OrderedDict()
        if processes == 1:
            for basename, path in name2path.items():
                #print(path)
                self.fort_files[basename] = FortranFile.from_path(path, macros=self.MACROS, verbose=self.verbose)
        else:
            from multiprocessing import Pool
            pool = Pool(processes=processes)
            results = pool.map(self._pool_f, list(name2path.items()))
            for basename, fort_file in results:
                self.fort_files[basename] = fort_file

        print("Parsing completed in %.2f [s]" % (time.time() - start))
        self.correlate()

    def _pool_f(self, item):
        #sys.stdout, sys.stderr = open(str(os.getpid()) + ".out", "w"), open(str(os.getpid()) + ".err", "w")
        basename, path = item
        return (basename, FortranFile.from_path(path, macros=self.MACROS, verbose=self.verbose))

    def correlate(self):
        print("Building dependency graph ...")
        start = time.time()

        # Build dependency graph TODO: check for cyclic dependencies?
        for fort_file in self.fort_files.values():
            for use_name in fort_file.all_uses:
                if use_name in EXTERNAL_MODS: continue
                try:
                    used_mod = self.all_modules[use_name]
                except KeyError:
                    #print(self.all_modules.keys())
                    raise RuntimeError((
                        "Cannot find Fortran module `%s`\n used by `%s`\nin Abinit project.\n" +
                        "It may be a syntax error or a stale import if you've removed the module.\n" +
                        "If it's an external module (e.g. mpi), add it to the EXTERNAL_MODS list in ~abinit/fkiss/project.py."
                        ) % (use_name, fort_file.path))
                fort_file.all_used_mods.append(used_mod)

                # This trick is needed for F90.in files that will be post-processed by the build system
                key = os.path.basename(used_mod.path)
                if key.endswith(".F90.in"): key = key[:-3]
                self.fort_files[key].all_usedby_mods.append(fort_file)

        pub_procs = self.get_all_public_procedures()
        all_interfaces = self.get_all_interfaces()
        miss = []
        for fort_file in self.fort_files.values():
            for proc in fort_file.all_public_procedures.values():
                for child_name in proc.children:
                    try:
                        pub_procs[child_name].parents.append(proc)
                    except KeyError:
                        # TODO: Could be in subroutine contains
                        #if child_name in all_interfaces:
                        #    print("Found in interfaces:", child_name)
                        #else:
                        miss.append(child_name)

        def is_internal(name):
            return not any(name.startswith(s) for s in
                set(("mpi_", "dfftw_", "mkl_", "papif_", "plasma_", "elpa", "etsf_io_", "blacs_",
                     "libpaw_", "gpu_", "xc_", "bigdft_")))

        miss = filter(is_internal, miss)
        from .check_linalg_calls import blas_routines, lapack_routines
        from .regex import FORTRAN_INTRINSICS
        if miss:
            miss = set(miss)
            # Remove blas and lapack routines.
            miss = miss.difference(blas_routines)
            miss = miss.difference(lapack_routines)
            # Remove pblas and scalalapack routines.
            miss = miss.difference("p" + n for n in blas_routines)
            miss = miss.difference("p" + n for n in lapack_routines)
            # Remove public interfaces.
            miss = miss.difference(self.get_all_interfaces().keys())
            # Remove intrincs.
            miss = miss.difference(FORTRAN_INTRINSICS)
            miss = sorted(miss)
            print("Cannot find %d callees. Use --verbose to show list." % len(miss))
            if self.verbose: pprint(miss)

        print("Graph completed in %.2f [s]" % (time.time() - start))

    def __str__(self):
         return self.to_string()

    def to_string(self, verbose=0, width=90):
        """
        String representation with verbosity level `verbose`.
        Text is wrapped at `width` columns.
        """
        lines = []; app = lines.append
        # TODO: Print project stats
        for fort_file in self.fort_files.values():
            app(fort_file.to_string(verbose=verbose, width=width))

        return "\n".join(lines)

    @lazy_property
    def all_modules(self):
        """
        Mapping modules name --> Module object with all modules in project.
        """
        omods = OrderedDict()
        for fort_file in self.fort_files.values():
            for m in fort_file.modules:
                assert m.name not in omods
                omods[m.name] = m
        return omods

    def pickle_dump(self, filepath=None):
        """
        Save the object in pickle format. Default name is used if filepath is None.
        """
        if filepath is None: filepath = self.get_default_pickle_file()
        with open(filepath, "wb") as fh:
            return pickle.dump(self, fh)

    @lazy_property
    def all_src_dirs(self):
        """List with all top level directories containing subdirectories with F90 file."""
        return [
            os.path.join(self.top, "shared", "common", "src"),
            # Exclude libpaw because shared/common/39_libpaw --> shared/libpaw/src
            # and we don't want duplicated entries.
            #os.path.join(self.top, "shared", "libpaw"),
            os.path.join(self.top, "src"),
        ]

    def get_dirpaths(self):
        """
        Return list of directory names containing source files.
        """
        l = []
        for src_dir in self.all_src_dirs:
            assert os.path.isdir(src_dir)
            #print("src_dir", src_dir) #, os.listdir(src_dir))
            dpaths = [os.path.join(src_dir, d) for d in os.listdir(src_dir)]
            s = [d for d in dpaths if os.path.isdir(d)
                 and os.path.isfile(os.path.join(d, "abinit.src"))
                 #and os.path.basename(d) not in black_list
            ]

            l += s

        #print(sorted(([os.path.basename(f) for f in l])))
        # 98_main does not have abinit.src so we have to add it here.
        return sorted(l + [os.path.join(self.top, "src", "98_main")])

    def needs_reload(self):
        """
        Returns True if source tree must be parsed again because:

            1. new files/directories have been added
            2. source files have been changed
        """
        if set(self.dirpaths) != set(self.get_dirpaths()): return True
        # FIXME: Return immediately if new files have been added...

        # Compare time of most recent content modification.
        try:
            return any(fort_file.stat.st_mtime != os.stat(fort_file.path).st_mtime
                       for fort_file in self.fort_files.values())
        except FileNotFoundError:
            return True

    #def fort_files_indirname(self, dirname):
    #    if dirname.endswith(os.sep): dirname = dirname[:-1]
    #    dir2files = self.groupby_dirname()
    #    return dir2files[dirname]

    def groupby_dirname(self):
        """
        Return dictionary {dirname --> [List of FortranFile in dirname]}

        Note how we use the dirname of the file as key in the dict.
        This is the trick that allows us to access source files in shared/common/src
        using e.g. "02_clib" instead of the shared/common/src/02_clib.
        """
        dir2files = defaultdict(list)
        for fort_file in self.fort_files.values():
            dir2files[fort_file.dirname].append(fort_file)
        dir2files = {k: dir2files[k] for k in sorted(dir2files.keys())}

        return OrderedDict(dir2files.items())

    def iter_dirname_fortfile(self):
        """Iterate over (dirname, fort_file)"""
        dir2files = self.groupby_dirname()
        for dirname, fort_files in dir2files.items():
            for fort_file in fort_files:
                yield dirname, fort_file

    def get_all_public_procedures(self):
        """
        Return Dictionary mapping name --> Procedure
        """
        # FIXME: Check that there's no name collision
        d = {}
        for f in self.fort_files.values():
            for a in ["modules", "programs", "subroutines", "functions"]:
                for p in getattr(f, a):
                    d.update({p.name: p for p in p.public_procedures})

        return {k: d[k] for k in sorted(d.keys())}

    def get_all_datatypes(self):
        """
        Dictionary name -> datatype with all datatypes available in the project.
        """
        dtypes = {}
        for f in self.fort_files.values():
            for p in getattr(f, "modules"):
                for dtype in p.types:
                    assert dtype.name not in dtypes
                    dtypes[dtype.name] = dtype

        return {k: dtypes[k] for k in sorted(dtypes.keys())}

    def get_all_interfaces(self):
        """
        Dictionary mapping name --> Interface
        """
        d = {}
        for f in self.fort_files.values():
            for p in getattr(f, "modules"):
                d.update({i.name: i for i in p.interfaces})

        return {k: d[k] for k in sorted(d.keys())}

    def print_interfaces(self, what=None, verbose=0):
        """
        Print ALL Fortran interfaces defined in project if what is None else interface with name `what`.
        """
        name2interface = self.get_all_interfaces()

        if what is not None:
            if what in name2interface:
                interface = name2interface[what]
                print(interface.to_string(verbose=verbose))
            else:
                cprint("Cannot find interface `%s` in project" % what, "red")
                matches = difflib.get_close_matches(what, name2interface.keys())
                if matches:
                    cprint("Perhaps you meant: {}".format(matches), "red")
        else:
            # Print all interfaces.
            for interface in name2interface.values():
                cprint(repr(interface), "yellow")
                print(interface.to_string(verbose=verbose))
                print("")

    def find_public_entity(self, name):
        """
        Find and return Procedure object with name `name`.
        Assume name is unique in the project.
        """
        all_names = []
        for f in self.fort_files.values():
            obj = f.find_public_entity(name, all_names=all_names)
            if obj is not None: return obj

        # Print closest matches
        cprint("Cannot find public entity `%s`" % str(name), "red")
        matches = difflib.get_close_matches(name, all_names)
        if matches:
            cprint("Perhaps you meant: {}".format(matches), "red")

        return None

    def find_datatype(self, name):
        """
        Find and return Datatype object with name `name`.
        Assume name is unique in the project.
        """
        all_names = []
        for f in self.fort_files.values():
            obj = f.find_datatype(name, all_names=all_names)
            if obj is not None: return obj

        # Print closest matches
        matches = difflib.get_close_matches(name, all_names)
        if matches: cprint("Perhaps you meant: {}".format(matches), "red")

        return None

    def find_module_from_entity(self, name):
        """
        Return the Module object that contains the public entity `name`.
        """
        obj = self.find_public_entity(name)
        if obj.is_module: return obj
        assert obj.ancestor is not None and obj.ancestor.is_module
        return obj.ancestor

    def print_dir(self, dirname, verbose=0):
        if dirname.endswith(os.sep): dirname = dirname[:-1]

        print("Printing info on directory:", dirname)
        dir2files = self.groupby_dirname()

        # Find parent directory
        for src_dir in self.all_src_dirs:
            dirname = os.path.join(src_dir, os.path.basename(dirname))
            if dirname in dir2files: break
        else:
            raise ValueError("Cannot find dirname `%s` in project" % dirname)

        if verbose:
            if verbose > 1:
                for f in dir2files[dirname]:
                    print(f.to_string(verbose=verbose, with_stats=False), "\n\n")

        df = self.get_stats_dir(dirname)
        print(df)

    #def find_parents(self, obj):
    #    """
    #    Find parents of object `obj` where obj can be a file or the name of a procedure.
    #    """
    #    parents = []
    #    print("Finding parents of", obj)
    #    if os.path.isfile(obj):
    #        # Assume module with same name as F90 file
    #        this = FortranFile(obj)
    #        for fort_file in self.fort_files.values():
    #            if this in fort_file.parents:
    #                print(fort_file.basename)
    #                #parents.append(fort_file)
    #    else:
    #        raise NotImplementedError()

    #    return parents

    #def detect_cyclic_deps(self):
    #    cycles = defaultdict(list)
    #    for fort_file in self.fort_files.values():
    #        # TODO: This should be recursive
    #        for child in fort_file.all_usedby_mods:
    #            if child.depends_on(fort_file):
    #                cycles[fort_file].append(child)
    #    return cycles

    def find_allmods(self, head_path, include_files_in_dirs=True):
        """
        Traverse the *entire* graph starting from head_path.
        Return full list of `Module` objects required by head_path.

        Args:
            include_files_in_dirs: True if the dependency graph should include
                the Fortran modules that are located in the same directory
                as the dependencies even if no explicit `use statement` is found.
                This option is needed when generating `binaries.conf` because `make`
                builds all files inside the directory instead of the minimal set
                required by the target. In a nutshell, we need to find the full
                set of dependencies of the directory.
        """
        dir2files = self.groupby_dirname()
        allmods = self._find_allmods(head_path, dir2files=dir2files)
        if not include_files_in_dirs: return allmods
        #print("initial list of modules", "\n".join(mod.basename for mod in allmods), "end initial list"

        # Include **all** modules inside the directories in which dependencies are located
        # so that make can build object files in all dirs.
        # This step is relatively costly, so avoid scanning the same directory twice.
        dirpaths = sorted(set(mod.dirpath for mod in allmods))
        #dirpaths = [os.path.join(self.srcdir, os.path.basename(d)) for d in dirpaths]
        visited = set()

        while dirpaths:
            dirpath = dirpaths.pop(0)
            if dirpath in visited: continue
            visited.add(dirpath)
            #print("In dirpath:", dirpath)
            for fort_file in dir2files[dirpath]:
                #print("Adding", fort_file.basename, "in dir", fort_file.dirpath)
                other_mods = self._find_allmods(fort_file.basename, dir2files=dir2files)
                allmods.update(other_mods)
                #other_dirpaths = set(os.path.join(self.srcdir, os.path.basename(d))
                #        for d in set(mod.dirname for mod in other_mods))
                other_dirpaths = set(mod.dirpath for mod in other_mods)
                dirpaths.extend(other_dirpaths - visited)

        return allmods

    def _find_allmods(self, head_path, dir2files=None, include_files_in_dirs=True):
        dir2files = dir2files if dir2files is not None else self.groupby_dirname()

        # This trick is needed for F90.in files that will be post-processed by the build system
        if head_path.endswith(".F90.in"): head_path = head_path[:-3]
        head = self.fort_files[head_path]

        allmods, queue, visited = set(), set(), set()
        queue.add(head)
        while queue:
            fort_file = queue.pop()
            for mod in fort_file.all_used_mods:
                if mod.basename in visited: continue
                visited.add(mod.basename)
                allmods.add(mod)
                queue.add(self.fort_files[mod.basename])

        return allmods

    def write_binaries_conf(self, dryrun=False, verbose=0):
        """
        Write new binaries.conf file
        """
        # Read binaries.conf and add new list of libraries.
        # NB: To treat `dependencies` in an automatic way, I would need either an
        # explicit "use external_module" or an explicit "include foo.h" so
        # that I can map these names to external libraries.
        # This means that I **cannot generate** the entire file in a programmatic way
        # but only the libraries entries.
        try:
            from ConfigParser import ConfigParser
        except ImportError:
            from configparser import ConfigParser

        binconf_path = os.path.join(self.top, "config", "specs", "binaries.conf")

        # Read INI file.
        config = ConfigParser()
        config.read(binconf_path)

        # Read header with comments to be added to the new conf file.
        # NB: use [DEFAULT] as sentinel.
        header = []
        with io.open(binconf_path, "rt", encoding="utf8") as fh:
            for line in fh:
                line = line.strip()
                if line == "[DEFAULT]":
                    break
                else:
                    header.append(line)
            else:
                raise ValueError("Cannot find `[DEFAULT]` string in %s" % binconf_path)
        header.append("")
        header = "\n".join(header)

        print("Finding all binary dependencies...")
        start = time.time()
        # Find programs
        program_paths = []
        for path, fort_file in self.fort_files.items():
            if fort_file.programs:
                program_paths.append((fort_file, path))

        for prog_file, path in program_paths:
            # Note include_files_in_dirs
            allmods = self.find_allmods(path, include_files_in_dirs=True)
            dirnames = sorted(set(mod.dirname for mod in allmods), reverse=True)
            if verbose:
                print("For program:", prog_file.name)
                pprint(dirnames)

            prog_name = prog_file.programs[0].name
            if prog_name.lower() == "fold2bloch": prog_name = "fold2Bloch"
            config.set(prog_name, "libraries", "\n" + "\n".join(dirnames))
            # py3k
            #config[prog_name]["libraries"] = "\n" + "\n".join(dirnames)

        print("Analysis completed in %.2f [s]" % (time.time() - start))

        try:
            from StringIO import StringIO
        except ImportError:
            from io import StringIO

        fobj = StringIO()
        fobj.write(header)
        config.write(fobj)
        if dryrun:
            print(fobj.getvalue())
        else:
            with open(binconf_path, "wt") as fh:
                fh.write(fobj.getvalue())
        fobj.close()

        return 0

    def write_buildsys_files(self, dryrun=False, verbose=0):
        """
        Write files require by buildsys:

            abinit.dep --> Dependencies inside the directory
            abinit.dir --> Dependencies outside the directory
            abinit.amf --> File with EXTRA_DIST
        """
        # Group Fortfiles by dirname
        dir2files = self.groupby_dirname()


        template = """\
# Dependencies ({kind}) of directory {directory}
#
# This file has been generated by abisrc.py.
# DO NOT edit this file. All changes will be lost.
# Use `abisrc.py makemake` to regenerate the file.

"""
        for dirpath, fort_files in dir2files.items():
            dirname = os.path.basename(dirpath)
            if dirname == "98_main": continue

            # Find dependencies inside this directory (abinit.dep)
            # Remove extensions if ext in (".F90" or "F90.in")
            inside_deps = {}
            for this in fort_files:
                dlist = [f.name for f in fort_files if any(m in this.all_used_mods for m in f.modules)]
                inside_deps[this.name] = sorted(set(d.replace(".F90.in", "").replace(".F90", "") for d in dlist))
            inside_deps = OrderedDict([(k, inside_deps[k]) for k in sorted(inside_deps.keys())])

            lines = []
            cleanfiles = ["CLEANFILES += \\"]
            n = len(inside_deps)
            for i, (k, dlist) in enumerate(inside_deps.items()):
                k = k.replace(".F90", "")
                end = " " if i == n - 1 else " \\"
                cleanfiles.append("\t%s.$(MODEXT)%s" % (k, end))
                if not dlist: continue
                lines.append("%s.$(OBJEXT): %s " % (k, " ".join("%s.$(OBJEXT)" % v for v in dlist)))
            cleanfiles.append("\n")

            # Write abinit.dep
            s = template.format(kind="inside the directory", directory=dirname)
            s += "\n".join(cleanfiles)
            s += "\n\n".join(lines)

            abinitdep_path = os.path.join(dirpath, "abinit.dep")
            if dryrun:
                print("# For dirpath:", dirpath)
                print(s, end=2 * "\n")
            else:
                #shutil.copyfile(abinitdep_path, abinitdep_path + ".bkp")
                with open(abinitdep_path, "wt") as fh:
                    fh.write(s)

            # Find dependencies outside this directory (abinit.dir).
            outside_dir = []
            for fort_file in fort_files:
                outside_dir.extend(m.dirname for m in fort_file.all_used_mods)

	    # Write abinit.dir
            s = template.format(kind="outside the directory", directory=os.path.basename(dirpath))
            s += "include_dirs = \\\n" + pformat(sorted(set(outside_dir)))
            abinitdir_path = os.path.join(dirpath, "abinit.dir")
            if dryrun:
                print(s, end=2 * "\n")
            else:
                #if os.path.exists(abinitdir_path): shutil.copyfile(abinitdir_path, abinitdir_path + ".bkp")
                with open(abinitdir_path, "wt") as fh:
                    fh.write(s)

    def touch_alldeps(self, verbose=0):
        """
        Touch all files that depend on the modules that have been changed.
        Return number of touched files.
        """
        def touch(fname, times=None):
            """Emulate Unix touch."""
            with open(fname, 'a'):
                os.utime(fname, times)

        # TODO: possible problem if new files have been added.
        count = 0
        changed_fort_files = []
        for fort_file in self.fort_files.values():
            if fort_file.stat.st_mtime != os.stat(fort_file.path).st_mtime:
                changed_fort_files.append(fort_file)
                print("Detected changes in:", os.path.relpath(fort_file.path))
                touch(fort_file.path)
                count += 1

        # Get list of modules that have been changed.
        changed_mods = []
        for f in changed_fort_files:
            changed_mods.extend(f.modules)

        for fort_file in self.fort_files.values():
            if any(m in fort_file.all_used_mods for m in changed_mods):
                print("Touching:", os.path.relpath(fort_file.path))
                touch(fort_file.path)
                count += 1

        return count

    def validate(self, verbose=0):
        """
        Validate project. Return exit status.
        """
        white_list = set(["m_psxml2ab.F90",])

        retcode = 0
        for fort_file in self.fort_files.values():
            count = len(fort_file.subroutines) + len(fort_file.functions)
            if count == 0: continue
            if fort_file.name in white_list:
                cprint("WHITE_LIST [%s] Found %d procedure(s) outside modules!" % (fort_file.name, count), "green")
            else:
                cprint("[%s] Found %d procedure(s) outside modules!" % (fort_file.name, count), "red")
                retcode += 1

        # check_dirlevel()
        for fort_file in self.fort_files.values():
            for child in fort_file.all_usedby_mods:
                if child.dirlevel < fort_file.dirlevel:
                    cprint("%s should be below level: %s" % (repr(child), fort_file.dirlevel), "red")
                    retcode += 1

        #retcode += self.check_abirules(verbose=verbose)

        cprint("retcode %d" % retcode, "green" if retcode == 0 else "red")
        return retcode

    def check_abirules(self, verbose=0):
        retcode = 0
        for name, fort_file in self.fort_files.items():
            print("Checking abirules for %s..." % name)
            rv = fort_file.check_abirules(verbose=verbose)
            retcode += rv
            msg = "%s [OK]" if rv == 0 else "%s [FAILED]"
            cprint(msg % name, color="green" if rv == 0 else "red")

        return retcode

    def pedit(self, name, verbose=0):
        """
        Edit all children of a public entity specified by name.
        """
        obj = self.find_public_entity(name)
        if obj is None: return 1
        if verbose: print(obj)

        # Find files with procedures.
        paths = sorted(set(p.path for p in obj.parents))
        if verbose: print(paths)
        # TODO from pymods.tools import Editor
        from fkiss.tools import Editor
        return Editor().edit_files(paths, ask_for_exit=True)

    def get_stats_file(self, filename, as_dict=False):
        return self.fort_files[os.path.basename(filename)].get_stats(as_dict=as_dict)

    def get_stats_dir(self, dirname):
        """
        Return dataframe with statistics about directory.
        """
        dir2files = self.groupby_dirname()
        dirpath = os.path.join(self.top, dirname)
        if dirpath.endswith(os.sep): dirpath = dirpath[:-1]
        index, rows = [], []
        for fort_file in sorted(dir2files[dirpath], key=lambda f: f.basename):
            index.append(fort_file.basename)
            rows.append(fort_file.get_stats(as_dict=True))

        import pandas as pd
        return pd.DataFrame(rows, index=index, columns=list(rows[0].keys() if rows else None))

    def get_stats(self):
        df_list = []
        for dirpath in self.dirpaths:
            dirname = os.path.basename(dirpath)
            if dirname == "incs": continue
            try:
                df = self.get_stats_dir(dirpath)
                df_list.append(df)
            except Exception as exc:
                print("Exception for dirpath: %s\n%s" % (dirpath, exc))
                raise exc

        import pandas as pd
        return pd.concat(df_list)

    def print_orphans(self, verbose=0):
        """Print orphan procedures and modules."""
        # FIXME: this does not work as expected.

        def find_orphans_in_fort_file(fort_file):
            """List with orphan procedures."""
            orphans = []
            for mod in fort_file.modules:
                if all(mod.name not in ffile.all_uses for ffile in self.fort_files.values()):
                    orphans.append(mod)

            for proc in fort_file.all_public_procedures.values():
                # programs are orphan by definition
                # function callers are difficult to detect.
                if proc.is_program or proc.is_function: continue
                if not proc.parents:
                    orphans.append(proc)

            return orphans

        dir2files = self.groupby_dirname()
        for _, fort_file in self.iter_dirname_fortfile():
            orphans = find_orphans_in_fort_file(fort_file)
            if orphans:
                cprint("Found %d orphans in file: %s" % (len(orphans), fort_file.basename), "yellow")
                for o in orphans:
                    print("\t", repr(o))

    def get_parent_dirs_of_dirname(self, dirname):
        """
        {parent_dirname --> [{module_name_in_parent_dirname: [list_of_modules_in_dirname]}]
        """
        if dirname.endswith(os.sep): dirname = dirname[:-1]

        dir2files = self.groupby_dirname()
        for src_dir in self.all_src_dirs:
            dirpath = os.path.join(src_dir, os.path.basename(dirname))
            if dirpath in dir2files: break
        else:
            raise ValueError("Cannot find dirname `%s` in project" % dirname)

        fort_files_indirname = dir2files[dirpath]
        #used_dirs = defaultdict(list)
        #dirlevel = fort_files_indirname[0].dirlevel

        for parent_dirname, parent_fort_files in dir2files.items():
            #parent_dirlevel = parent_fort_files[0].dirlevel
            #if parent_dirlevel < dirlevel: continue
            for parent_fort_file in parent_fort_files:
                for use_name in parent_fort_file.all_uses:
                    if use_name in EXTERNAL_MODS: continue
                    #if use_name in
                    #used_dir[parent_dirname].append(

        # Order keys.
        #return {k: used_dirs[k] for k in sorted(used_dirs.keys())}

    def get_dirs_used_by_dirname(self, dirname):
        """
        {used_dirname --> [{module_name_in_used_dirname: [list_of_modules_in_dirname]}]
        """
        if dirname.endswith(os.sep): dirname = dirname[:-1]

        dir2files = self.groupby_dirname()
        for src_dir in self.all_src_dirs:
            dirpath = os.path.join(src_dir, os.path.basename(dirname))
            if dirpath in dir2files: break
        else:
            raise ValueError("Cannot find dirname `%s` in project" % dirname)

        fort_files_indirname = dir2files[dirpath]

        # {uses_dirname --> [{module_name_in_used_dirname: [list_of_modules_in_parent_dirname]}]
        used_dirs = defaultdict(list)

        for fort_file in fort_files_indirname:
            for use_name in fort_file.all_uses:
                if use_name in EXTERNAL_MODS: continue
                used_mod = self.all_modules[use_name]
                used_dirs[used_mod.dirname].append(used_mod.name)

        for k, v in used_dirs.items():
            d = {k: [] for k in sorted(set(used_dirs[k]))}
            for fort_file in fort_files_indirname:
                for use_name in fort_file.all_uses:
                    if use_name in EXTERNAL_MODS: continue
                    if use_name in d:
                        d[use_name].append(fort_file.name)

            used_dirs[k] = {k: list(set(v)) for k, v in d.items()}

        # Order keys.
        return {k: used_dirs[k] for k in sorted(used_dirs.keys())}

    def get_graphviz_dir(self, dirname, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate dependency graph for directory `dirname` in the DOT language

        Args:
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        if dirname.endswith(os.sep): dirname = dirname[:-1]
        used_dirs = self.get_dirs_used_by_dirname(dirname)

        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph
        fg = Digraph(name="Connections of dir: %s" % dirname,
            graph_attr=graph_attr, node_attr=node_attr, edge_attr=edge_attr,
            engine="dot" if engine == "automatic" else engine,
        )

        # Set graph attributes.
        fg.attr(label="Connections of dir: %s" % dirname,
                rankdir="LR", pagedir="BL") # constraint="false", pack="true", packMode="clust")
        fg.node_attr.update(color='lightblue2', style='filled')

        target = os.path.basename(dirname)
        fg.node(name=target)
        for used_dirname in used_dirs:
            fg.edge(target, used_dirname)
        #for parent_dir in parent_dirs:
        #    fg.edge(target, used_dirname)

        return fg

    def get_graphviz_pubname(self, name, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate dependency graph for public procedure `name` in the DOT language
        (only parents and children modules of this file).

        Args:
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        obj = self.find_public_entity(name)

        if obj is None:
            cprint("Cannot find public entity `%s` in Abinit project." % name, "red")
            return None
        if not obj.is_procedure:
            cprint("Only procedures can be visualized with graphviz. Received class: %s" % obj.__class__, "yellow")
            return None

        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph
        fg = Digraph(name="Connections of %s: %s" % (obj.__class__.__name__, obj.name),
            graph_attr=graph_attr, node_attr=node_attr, edge_attr=edge_attr,
            engine="dot" if engine == "automatic" else engine,
        )

        # Set graph attributes.
        fg.attr(label="Connections of %s: %s" % (obj.__class__.__name__, obj.name),
                rankdir="LR", pagedir="BL") # constraint="false", pack="true", packMode="clust")
        fg.node_attr.update(color='lightblue2', style='filled')

        target = obj.name
        fg.node(name=target)
        for child in obj.children:
            fg.edge(target, child)
        for parent in obj.parents:
            fg.edge(parent.name, target)

        return fg

    def master(self):
        return """\

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

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_markdown_cell("## Abinit Project"),
            nbv.new_code_cell("""
from fkiss.project import AbinitProject
proj = AbinitProject.pickle_load(filepath=None)
""",),
            nbv.new_code_cell('#proj.get_graphviz_dir("41_geometry")'),
            nbv.new_code_cell('#proj.get_graphviz_pubname("crystal_init")'),
            nbv.new_code_cell('#proj.get_graphviz_pubname("fourdp")'),
            nbv.new_code_cell('#proj.get_graphviz_pubname("m_geometry")'),
            nbv.new_code_cell('#proj.get_stats()'),
            nbv.new_code_cell('#proj.get_stats_dir("src/41_geometry")'),
            nbv.new_code_cell("""
try:
    import param
    import panel as pn
    pn.extension()
    proj_viewer = proj.get_viewer()
    proj_viewer.get_panel()
except ImportError as exc:
    print(exc)
    print('Install panel with `pip install panel`\nSee also: https://panel.pyviz.org/index.html#installation')
""")
        ])

        return self._write_nb_nbpath(nb, nbpath)

    def yield_figs(self, **kwargs):  # pragma: no cover
        # TODO: Activate this
        return None

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""
        return ProjectViewer(self).get_panel()

    def get_procedure_viewer(self):
        """Return tabs with widgets to interact with the DDB file."""
        return ProcedureViewer(self).get_tabs()

import param
import panel as pn

class ProjectViewer(param.Parameterized):

    dir_select = pn.widgets.Select(name="Directory")
    file_select = pn.widgets.Select(name="Fortran File")
    pubproc_select = pn.widgets.Select(name="Public procedure")

    def __init__(self, proj, **params):
        super().__init__(**params)
        self.proj = proj

        self.dir2files = proj.groupby_dirname()
        self.dirname2path = {os.path.basename(p): p for p in self.dir2files.keys()}
        self.dir_select.options = list(self.dirname2path.keys())
        #self.all_pubs = list(proj.get_all_public_procedures.keys())

        #self.autocomplete = pn.widgets.AutocompleteInput(
        #        name='Autocomplete Input', options=list(self.all_pubs.keys()),
        #        placeholder='Write something here')

    @param.depends('dir_select.value')
    def view_dirname(self):
        dirpath = self.dirname2path[self.dir_select.value]
        # Change options in file_select.
        self.file_select.options = [f.name for f in self.dir2files[dirpath]]
        #self.file_select.value = self.file_select.options[0]
        self.file_select.param.trigger("value")

        return pn.Row(self.proj.get_stats_dir(dirpath),
                      self.proj.get_graphviz_dir(dirpath))

    @param.depends('file_select.value')
    def view_fort_file(self):
        dirpath = self.dirname2path[self.dir_select.value]
        for fort_file in self.dir2files[dirpath]:
            if fort_file.name == self.file_select.value:
                break
        else:
            raise ValueError("Cannot find fortran file with name: %s" % self.file_select.value)

        # Change options in pubproc_select.
        self.pubproc_select.options = list(fort_file.all_public_procedures.keys())
        #self.pubproc_select.value = self.pubproc_select.options[0]
        self.pubproc_select.param.trigger("value")

        return pn.Row(fort_file.get_stats(), fort_file.get_graphviz())

    @param.depends('pubproc_select.value')
    def view_pubproc(self):
        pubname = self.pubproc_select.value
        if pubname is None: return
        graph = self.proj.get_graphviz_pubname(pubname) #, engine=options.engine)
        return pn.Row(graph)

    def get_panel(self):
        tabs = pn.Tabs()

        tabs.append(("Directory", pn.Column(self.dir_select, pn.Column(self.view_dirname))))

        tabs.append(("File", pn.Column(
            pn.Column(self.dir_select, self.file_select), pn.Column(self.view_fort_file)),
        ))

        tabs.append(("Procedure", pn.Column(
            pn.Column(self.dir_select, self.file_select, self.pubproc_select),
            pn.Column(self.view_pubproc)),
        ))

        return tabs


class ProcedureViewer(param.Parameterized):

    #autocomplete = pn.widgets.AutocompleteInput(
    autocomplete = pn.widgets.TextInput(
                name='Autocomplete Input',
                placeholder='Write something here',
                ) #options=[], #list(self.all_pubs.keys()),

    view_btn = pn.widgets.Button(name="View", button_type='primary')

    def __init__(self, proj, **params):
        super().__init__(**params)
        self.proj = proj

        self.all_pubs = proj.get_all_public_procedures()
        self.autocomplete = pn.widgets.AutocompleteInput(
                name='Autocomplete Input', options=list(self.all_pubs.keys()),
                placeholder='Write something here')

        #self.autocomplete.options = list(self.all_pubs.keys())
        #self.autocomplete.param.watch(self.view_pubname, 'value')
        #self.autocomplete.param.trigger("value")

    @param.depends('view_btn.clicks')
    def view_pubname(self):
        print("in view_pubname with clicks", self.view_btn.clicks)
        if self.view_btn.clicks == 0: return
        pubname = self.autocomplete.value
        print("pubname:", pubname)
        if pubname is None: return # or pubname not in self.all_pubs: return
        #proc = self.all_pubs[pubname]
        graph = self.proj.get_graphviz_pubname(pubname) #, engine=options.engine)
        return graph

    def get_panel(self):
        tabs = pn.Tabs()

        tabs.append(("Procedure", pn.Row(pn.Column(self.autocomplete, self.view_btn), pn.Column(self.view_pubname))))
        return tabs
