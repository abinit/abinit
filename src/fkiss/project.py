"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re
import pickle

from collections import OrderedDict, defaultdict
from textwrap import TextWrapper
from pprint import pprint
from .parser import FortranKissParser
from .tools import lazy_property


EXTERN_MODS = {
    # Intrinsics
    "iso_fortran_env",
    "iso_c_binding",
    "ieee_arithmetic",
    "ieee_exceptions",
    "ieee_features",
    # MPI modules.
    "mpi",
    "mpi_f08",
    # Modules provide by compilers.
    "f90_unix_proc",
    "f90_unix_dir",
    "ifcore",
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
    "m_build_info",
    "m_optim_dumper",
    "libxc_functionals",
}

class FortranFile(object):
    """
    Base class for files containing Fortran source code.

    modules:
    programs:
    subroutine
    functions
    """

    @classmethod
    def from_path(cls, path, verbose=0):
        p = FortranKissParser(verbose=verbose).parse_file(path)

        new = cls(path)
        new.includes = p.includes
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
        self.name = os.path.basename(self.path)
        self.dirname = os.path.dirname(self.path)

        # Save initial stat values, used to understand if reload is needed.
        self.stat = os.stat(self.path)

        self.programs, self.modules, self.subroutines, self.functions = [], [], [], []

        #self.uses, self.includes = [], []
        self.required_mods, self.children_mods = [], []
        #self.num_f90lines, self.num_doclines = 0, 0

    """
    def __eq__(self, other):
        if other is None: return False
        if not isinstance(other, self.__class__): return False
        return self.path == other.path

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.path)
    """

    def stree(self):
        lines = [repr(self)]; app = lines.append
        for a in ["programs", "modules", "subroutines", "functions"]:
            for p in getattr(self, a):
                app(p.stree())

        return "\n".join(lines)

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.path)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, width=90):
        w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)
        lines = []; app = lines.append
        app("%s:\n\t%s" % (self.__class__.__name__, os.path.relpath(self.path)))
        app("Required modules:\n%s\n" % w.fill(", ".join(mod.name for mod in self.required_mods)))
        app("Used by modules:\n%s\n" % w.fill(", ".join(mod.name for mod in self.children_mods)))
        app("Includes:\n%s\n" % w.fill(", ".join(self.includes)))

        for a in ["programs", "modules", "subroutines", "functions"]:
            plist = getattr(self, a)
            if not plist: continue
            app("Fortran file contains %d %s" % (len(plist), a))
            for p in plist:
                app(p.to_string(verbose=verbose, width=width))

        app(self.stree())

        return "\n".join(lines)

    @lazy_property
    def all_uses(self):
        all_uses = []
        for a in ["programs", "modules", "subroutines", "functions"]:
            for p in getattr(self, a):
                all_uses.extend(p.uses)
        return sorted(set(all_uses))

    def iter_procedures(self, visibility="public"):
        for a in ["modules", "programs", "subroutines", "functions"]:
            for p in getattr(self, a):
                for e in p.public_procedures:
                    #if e.name == "gstate": raise RuntimeError("Gstate")
                    yield e

    def find_public_entity(self, name):
        """Find and returh public procedure with `name` or None if not found"""
        for a in ["modules", "programs", "subroutines", "functions"]:
            for p in getattr(self, a):
                for e in p.public_procedures:
                    if e.name == name: return e
        return None

    def get_graphviz(self, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate dependency graph for this file in the DOT language
        (only required_mods and children of this file).

        Args:
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph

        graph_attr = {
            'size':'8.90625,1000.0',
            'rankdir': "LR",
            'concentrate':'true',
            #'id':self.ident
        }
        node_attr = {
            'shape':'box',
            'height':'0.0',
            'margin':'0.08',
            'fontname':'Helvetica',
            'fontsize':'10.5'
        }

        edge_attr = {
            'fontname':'Helvetica',
            'fontsize':'9.5'
        }
        graph_attr = {}
        node_attr =  {}
        edge_attr = {}

        #fg = Digraph("Abinit project", engine="dot" if engine == "automatic" else engine)
        fg = Digraph("Abinit project",
            graph_attr=graph_attr,
            node_attr=node_attr,
            edge_attr=edge_attr,
            #format='svg',
            engine="dot" if engine == "automatic" else engine,
        )

        # Set graph attributes.
        #fg.attr(label="%s@%s" % (self.__class__.__name__, self.relworkdir))
        #fg.attr(label=repr(self))
        fg.attr(label=self.name)
        #fg.attr(fontcolor="white", bgcolor='purple:pink')
        fg.attr(rankdir="LR", pagedir="BL")
        #fg.attr(constraint="false", pack="true", packMode="clust")
        fg.node_attr.update(color='lightblue2', style='filled')

        # Add input attributes.
        #if graph_attr is not None: fg.graph_attr.update(**graph_attr)
        #if node_attr is not None: fg.node_attr.update(**node_attr)
        #if edge_attr is not None: fg.edge_attr.update(**edge_attr)

        def node_kwargs(node):
            return node_attr
            return dict(
                #shape="circle",
                #color=node.color_hex,
                #label=(str(node) if not hasattr(node, "pos_str") else
                #    node.pos_str + "\n" + node.__class__.__name__),
            )

        edge_kwargs = dict(arrowType="vee", style="solid")
        cluster_kwargs = dict(rankdir="LR", pagedir="BL", style="rounded", bgcolor="azure2")

        # Build clusters representing directories
        all_nodes = sorted([self] + self.required_mods + self.children_mods, key=lambda o: o.dirname)
        dir2nodes = defaultdict(list)
        for node in all_nodes:
            dir2nodes[node.dirname].append(node)

        for dirname, nodes in dir2nodes.items():
            #print(dirname, nodes)
            cluster_name = "cluster%s" % dirname
            with fg.subgraph(name=cluster_name) as wg:
                wg.graph_attr.update(**graph_attr)
                wg.node_attr.update(**node_attr)
                wg.edge_attr.update(**edge_attr)
                #wg.attr(**cluster_kwargs)
                wg.attr(label=os.path.basename(dirname))
                for node in nodes:
                    wg.node(node.name) #, **node_attr) #, **node_kwargs(child))

        for mod in self.required_mods:
            fg.edge(mod.name, self.name) #, label=edge_label, # color=self.color_hex,
        for child in self.children_mods:
            fg.edge(self.name, child.name) #, label=edge_label, # color=self.color_hex,

        return fg

class AbinitProject(object):

    DEFAULT_PICKLE_FILE = "_project.pickle"

    IGNORED_FILES = {"m_build_info.F90", "m_optim_dumper.F90"}

    @classmethod
    def pickle_load(cls, filepath=None):
        if filepath is None: filepath = cls.DEFAULT_PICKLE_FILE
        with open(filepath, "rb") as fh:
            return pickle.load(fh)

    def __init__(self, srcdir, verbose=0):
        # Find directories with abinit.src files inside srcdir
        # and get list of files treated by the build system.
        self.srcdir = os.path.abspath(srcdir)

        # FIXME: 98_main does not provide abinit.src so abinit and other exec are not in project
        self.dirpaths = self.get_dirpaths()

        import imp
        def filter_fortran(files):
            return [f for f in files if f.endswith(".f") or f.endswith(".F90")]

        # basename --> FortranFile
        print("Analyzing directories...")
        self.fort_files = OrderedDict()
        for d in self.dirpaths:
            if verbose: print("Analyzing directory:", d)
            # Source files
            abinit_src = os.path.join(d, "abinit.src")
            #abinit_amf = os.path.join(d, "abinit.amf")
            mod = imp.load_source(abinit_src, abinit_src)
            for basename in filter_fortran(mod.sources):
                if basename in self.IGNORED_FILES: continue
                path = os.path.abspath(os.path.join(d, basename))
                fort_file = FortranFile.from_path(path, verbose=verbose)
                if basename in self.fort_files:
                    raise RuntimeError("Found two fortran files with same basename `%s`" % basename)
                self.fort_files[basename] = fort_file

        # Build dependency graph and check for cyclic dependencies.
        for fort_file in self.fort_files.values():
            for use_name in fort_file.all_uses:
                if use_name in EXTERN_MODS: continue
                try:
                    used_mod = self.all_modules[use_name]
                except KeyError:
                    raise RuntimeError(("Fortran module `%s` used by `%s` not found in Abinit project.\n" +
                                        "Add it to EXTERN_MODS set.") % (use_name, fort_file.path))
                fort_file.required_mods.append(used_mod)

                # FIXME
                key = os.path.basename(used_mod.path)
                self.fort_files[key].children_mods.append(fort_file)

        d = self.all_public_procedures()
        assert "gstate" in d
        miss = []
        for fort_file in self.fort_files.values():
            for p in fort_file.iter_procedures():
                for child_name in p.children:
                    try:
                        d[child_name].parents.append(p)
                    except KeyError:
                        miss.append(child_name)

        def is_internal(name):
            return any(name.startswith(s) for s in
                       ["mpi_", "fftw_", "mkl_", "papif_", "plasma_", "etsf_io_"])

        miss = filter(is_internal, miss)
        if miss:
            miss = sorted(set(miss))
            print("Cannot find %d callees" % len(miss))
            #pprint(miss)

    def __str__(self):
         return self.to_string()

    def to_string(self, verbose=0, width=90):
        lines = []; app = lines.append
        # TODO: Print project stats
        for fort_file in self.fort_files.values():
            app(fort_file.to_string(verbose=verbose, width=width))

        return "\n".join(lines)

    @lazy_property
    def all_modules(self):
        """Mapping modules name --> Module object with all modules in project"""
        omods = OrderedDict()
        for fort_file in self.fort_files.values():
            for m in fort_file.modules:
                assert m.name not in omods
                omods[m.name] = m
        return omods

    def pickle_dump(self, filepath=None):
        if filepath is None: filepath = self.DEFAULT_PICKLE_FILE
        with open(filepath, "wb") as fh:
            return pickle.dump(self, fh)

    def get_dirpaths(self):
        return sorted([d for d in os.listdir(self.srcdir) if os.path.isdir(d) and
                       os.path.isfile(os.path.join(d, "abinit.src"))])

    def needs_reload(self):
        if set(self.dirpaths) != set(self.get_dirpaths()): return True

        # TODO
        # Return immediately if new files have been added...

        # Compare time of most recent content modification.
        return any(fort_file.stat.st_mtime != os.stat(fort_file.path).st_mtime
                   for fort_file in self.fort_files.values())

    def groupby_dirname(self):
        """Dictionary dirname --> List of FortranFile."""
        dir2files = defaultdict(list)
        for fort_file in self.fort_files.values():
            dir2files[fort_file.dirname].append(fort_file)
        return dir2files

    def all_public_procedures(self):
        """name --> Procedure"""
        # FIXME: Check that there's no name collision
        d = {}
        for f in self.fort_files.values():
            for a in ["modules", "programs", "subroutines", "functions"]:
                for p in getattr(f, a):
                    d.update({p.name: p for p in p.public_procedures})
        return d

    def find_public_entity(self, name):
        for f in self.fort_files.values():
            obj = f.find_public_entity(name)
            if obj is not None: return obj
        return None

    def print_dir(self, dirname, verbose=0):
        dirname = os.path.basename(dirname).replace(os.sep, "")
        files = [f for f in self.fort_files.values() if os.path.basename(f.dirname) == dirname]
        for f in files:
            print(f)
            print("")

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
    #        for child in fort_file.children_mods:
    #            if child.depends_on(fort_file):
    #                cycles[fort_file].append(child)
    #    return cycles

    def write_buildsys_files(self):
        """
        Write

            abinit.amf --> File with EXTRA_DIST
            abinit.dep --> Dependencies (inside the directory) of the directory ./src/52_fft_mpi_noabirule
            abinit.dir --> Dependencies (outside the directory) of the directory ./src/52_fft_mpi_noabirule
        """
        dir2files = self.groubpy_dirname()

        # Sort...

        header = """\
# Dependencies ({}) of directory {directory}
#
# This file has been generated by abigraph.py.
# Do not edit this file. All changes will be lost.
# Use abigraph.py makemake to regenerate the file.

"""
        for dirname, fort_files in dir2files.items():
            # Dependencies inside the directory.
            for ffile in fort_files:
                parents_in_dir = ffile.filter_parents(fort_files)

            # Dependencies outside the directory
            parents_outside_dir = []
            for ffile in fort_files:
                parents_outside_dir.extends(ffile.required_mods)

    def get_graphviz_dir(self, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate dependency graph for this file in the DOT language
        (only parents and children modules of this file).

        Args:
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph
        fg = Digraph("Abinit project", engine="dot" if engine == "automatic" else engine)

        return fg
