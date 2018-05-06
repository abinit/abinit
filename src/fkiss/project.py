"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re
import time
import pickle

from collections import OrderedDict, defaultdict
from textwrap import TextWrapper
from pprint import pprint, pformat
from .parser import FortranKissParser
from .tools import lazy_property


EXTERNAL_MODS = {
    # Intrinsics
    "iso_fortran_env",
    "iso_c_binding",
    "ieee_arithmetic",
    "ieee_exceptions",
    "ieee_features",
    # MPI modules.
    "mpi",
    "mpi_f08",
    # Modules provided by compilers.
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

    .. attributes:

        modules:
        programs:
        subroutines
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
        # TODO name --> basename
        self.name = os.path.basename(self.path)
        self.dirname = os.path.dirname(self.path)

        # Save initial stat values, used to understand if reload is needed.
        self.stat = os.stat(self.path)

        self.programs, self.modules, self.subroutines, self.functions = [], [], [], []

        #self.uses, self.includes = [], []
        self.used_mods, self.usedby_mods = [], []
        #self.num_f90lines, self.num_doclines = 0, 0

    def __eq__(self, other):
        if other is None: return False
        if not isinstance(other, self.__class__): return False
        return self.path == other.path

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.path)

    def stree(self):
        """Return string with textual representation of the tree."""
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
        """
        String representation with verbosity level `verbose`.
        Text is wrapped at `width` columns.
        """
        w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)
        lines = []; app = lines.append
        app("%s:\n\t%s\n" % (self.__class__.__name__, os.path.relpath(self.path)))
        app("Use modules:\n%s\n" % w.fill(", ".join(mod.name for mod in self.used_mods)))
        app("Includes:\n%s\n" % w.fill(", ".join(self.includes)))
        app("Used by modules:\n%s\n" % w.fill(", ".join(mod.name for mod in self.usedby_mods)))

        for a in ["programs", "modules", "subroutines", "functions"]:
            plist = getattr(self, a)
            if not plist: continue
            app("Fortran file contains %d %s" % (len(plist), a))
            for p in plist:
                app(p.to_string(verbose=verbose, width=width))

        if verbose:
            app(self.stree())

        return "\n".join(lines)

    @lazy_property
    def all_uses(self):
        """
        List with full list of modules (string) used by the Fortran file
        Includes modules used by the procedures declared in the file.
        """
        all_uses = []
        for a in ["programs", "modules", "subroutines", "functions"]:
            for p in getattr(self, a):
                all_uses.extend(p.uses)
        return sorted(set(all_uses))

    def iter_procedures(self, visibility="public"):
        for a in ["modules", "programs", "subroutines", "functions"]:
            for p in getattr(self, a):
                for e in p.public_procedures:
                    yield e

    def find_public_entity(self, name):
        """
        Find and return the public procedure with `name` or None if not found
        """
        for a in ["modules", "programs", "subroutines", "functions"]:
            for p in getattr(self, a):
                for e in p.public_procedures:
                    if e.name == name: return e
        return None

    def get_graphviz(self, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate dependency graph for this Fortran file in the DOT language
        (only used_mods and children of this file).

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
        all_nodes = sorted([self] + self.used_mods + self.usedby_mods, key=lambda o: o.dirname)
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

        for mod in self.used_mods:
            fg.edge(mod.name, self.name) #, label=edge_label, # color=self.color_hex,
        for child in self.usedby_mods:
            fg.edge(self.name, child.name) #, label=edge_label, # color=self.color_hex,

        return fg


class AbinitProject(object):

    DEFAULT_PICKLE_FILE = "_project.pickle"

    IGNORED_FILES = {"m_build_info.F90", "m_optim_dumper.F90"}

    @classmethod
    def pickle_load(cls, filepath=None):
        """
        Reconstruct object from pickle file. Default is used if filepath is None.
        """
        if filepath is None: filepath = cls.DEFAULT_PICKLE_FILE
        with open(filepath, "rb") as fh:
            return pickle.load(fh)

    def __init__(self, srcdir, verbose=0):
        # Find directories with abinit.src files inside srcdir
        # and get list of files treated by the build system.
        self.srcdir = os.path.abspath(srcdir)
        self.dirpaths = self.get_dirpaths()

        import imp
        def filter_fortran(files):
            return [f for f in files if f.endswith(".f") or f.endswith(".F90")]

        # Build dict basename --> FortranFile
        print("Analyzing directories...")
        start = time.time()
        self.fort_files = OrderedDict()
        for d in self.dirpaths:
            if verbose: print("Analyzing directory:", d)
            if os.path.basename(d) == "98_main":
                # Treat executables
                for basename in filter_fortran(os.listdir(d)):
                    path = os.path.join(d, basename)
                    self.fort_files[basename] = FortranFile.from_path(path, verbose=verbose)
            else:
                # Source files
                abinit_src = os.path.join(d, "abinit.src")
                #abinit_amf = os.path.join(d, "abinit.amf")
                mod = imp.load_source(abinit_src, abinit_src)
                for basename in filter_fortran(mod.sources):
                    if basename in self.IGNORED_FILES: continue
                    path = os.path.abspath(os.path.join(d, basename))
                    fort_file = FortranFile.from_path(path, verbose=verbose)
                    if basename in self.fort_files:
                        raise RuntimeError("Found two Fortran files with same basename `%s`" % basename)
                    self.fort_files[basename] = fort_file

        # def correlate()
        # Build dependency graph
        # TODO: check for cyclic dependencies. ?
        for fort_file in self.fort_files.values():
            for use_name in fort_file.all_uses:
                if use_name in EXTERNAL_MODS: continue
                try:
                    used_mod = self.all_modules[use_name]
                except KeyError:
                    raise RuntimeError(("Fortran module `%s` used by `%s` not found in Abinit project.\n" +
                                        "Add it to EXTERNAL_MODS set.") % (use_name, fort_file.path))
                fort_file.used_mods.append(used_mod)

                # FIXME
                key = os.path.basename(used_mod.path)
                self.fort_files[key].usedby_mods.append(fort_file)

        d = self.all_public_procedures()
        #assert "gstate" in d
        miss = []
        for fort_file in self.fort_files.values():
            for p in fort_file.iter_procedures():
                for child_name in p.children:
                    try:
                        d[child_name].parents.append(p)
                    except KeyError:
                        miss.append(child_name)

        def is_internal(name):
            return not any(name.startswith(s) for s in
                          ("mpi_", "dfftw_", "mkl_", "papif_", "plasma_", "etsf_io_"))

        miss = filter(is_internal, miss)
        if miss:
            miss = sorted(set(miss))
            print("Cannot find %d callees. Use --verbose to show list." % len(miss))
            if verbose: pprint(miss)

        print("Analysis completed in %.2f [s]" % (time.time() - start))

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
        if filepath is None: filepath = self.DEFAULT_PICKLE_FILE
        with open(filepath, "wb") as fh:
            return pickle.dump(self, fh)

    def get_dirpaths(self):
        """
        Return list of directory names with source files.
        """
        l = sorted([d for d in os.listdir(self.srcdir) if os.path.isdir(d) and
                    os.path.isfile(os.path.join(d, "abinit.src"))])

        # 98_main does not contain abinit.src so we have to add it explicitely
        return l + [os.path.join(self.srcdir, "98_main")]

    def needs_reload(self):
        """
        Returns True if source tree must be parsed again because:

            1. new files/directories have been added
            2. source files have been changed
        """
        if set(self.dirpaths) != set(self.get_dirpaths()): return True
        # TODO
        # Return immediately if new files have been added...

        # Compare time of most recent content modification.
        try:
            return any(fort_file.stat.st_mtime != os.stat(fort_file.path).st_mtime
                       for fort_file in self.fort_files.values())
        except FileNotFoundError:
            return True

    def groupby_dirname(self):
        """
        Return dictionary mapping dirname --> List of FortranFile.
        """
        dir2files = defaultdict(list)
        for fort_file in self.fort_files.values():
            dir2files[fort_file.dirname].append(fort_file)
        return OrderedDict(dir2files.items())

    def all_public_procedures(self):
        """
        Dictionary mapping name --> Procedure
        """
        # FIXME: Check that there's no name collision
        d = {}
        for f in self.fort_files.values():
            for a in ["modules", "programs", "subroutines", "functions"]:
                for p in getattr(f, a):
                    d.update({p.name: p for p in p.public_procedures})
        return d

    def find_public_entity(self, name):
        """
        Find and return Procedure object with name `name`.
        Assume name is unique in the project.
        """
        for f in self.fort_files.values():
            obj = f.find_public_entity(name)
            if obj is not None: return obj
        return None

    def find_module_from_entity(self, name):
        obj = self.find_public_entity(name)
        if obj.is_module: return obj
        assert obj.ancestor is not None and obj.ancestor.is_module
        return obj.ancestor

    def print_dir(self, dirname, verbose=0):
        #print("Print directory:", dirname)
        #dirname = os.path.basename(dirname).replace(os.sep, "")
        #dirname = dirname.replace(os.sep, "")
        #files = [f for f in self.fort_files.values() if os.path.basename(f.dirname) == dirname]
        print("Print directory:", dirname)
        dir2files = self.groupby_dirname()
        dirname = os.path.join(self.srcdir, os.path.basename(dirname))
        for f in dir2files[dirname]:
            print(f.to_string(verbose=verbose))
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
    #        for child in fort_file.usedby_mods:
    #            if child.depends_on(fort_file):
    #                cycles[fort_file].append(child)
    #    return cycles

    def find_allmods(self, root_path):
        """
        Traverse the *entire* graph starting from root_path.
        and return full list of `Module` required by this file.
        """
        head = self.fort_files[root_path]
        allmods = set()
        queue = set()
        seen = set()
        queue.add(head)
        while queue:
            fort_file = queue.pop()
            for mod in fort_file.used_mods:
                #print(mod.name)
                if mod.basename in seen: continue
                seen.add(mod.basename)
                allmods.add(mod)
                queue.add(self.fort_files[mod.basename])
        return allmods

    def write_binaries_conf(self):
        print("Find all dependenciies of binaries...")
        start = time.time()
        # Find programs
        program_paths = []
        for path, fort_file in self.fort_files.items():
            if fort_file.programs:
                program_paths.append((fort_file, path))

        for prog_file, path in program_paths:
            allmods = self.find_allmods(path)
            dirnames = sorted(set(mod.dirname for mod in allmods), reverse=True)
            print("For program:", prog_file.name)
            pprint(dirnames)

        print("Analysis completed in %.2f [s]" % (time.time() - start))

        import configparser
        config = configparser.ConfigParser()
        config.read(os.path.join(self.srcdir, "..", "config", "specs", "binaries.conf"))
        print(config)
        #for prog_file, path in program_paths:
        #    config[prog_file.name]["libraries"]

        return 0

    def write_buildsys_files(self):
        """
        Write

            abinit.dep --> Dependencies inside the directory
            abinit.dir --> Dependencies outside the directory
            abinit.amf --> File with EXTRA_DIST
        """
        dir2files = self.groupby_dirname()
        # Sort...

        template = """\
# Dependencies ({kind}) of directory {directory}
#
# This file has been generated by abisrc.py.
# DO NOT edit this file. All changes will be lost.
# Use `abisrc.py makemake` to regenerate the file.

"""
        for dirpath, fort_files in dir2files.items():

            # Find dependencies inside this directory (abinit.dep)
            inside_deps = {}
            for this in fort_files:
                dlist = [f.name for f in fort_files if any(m in this.used_mods for m in f.modules)]
                inside_deps[this.name] = sorted(set(d.replace(".F90", "") for d in dlist))

            lines = []
            for k, dlist in inside_deps.items():
                if not dlist: continue
                k = k.replace(".F90", "")
                lines.append("%s.$(OBJEXT): %s" % (k, " ".join("%s.$(OBJEXT)" % v for v in dlist)))

            s = template.format(kind="inside the directory", directory=os.path.basename(dirpath))
            s += "\n\n".join(lines)
            print(s, end=2 * "\n")
            # TODO: Ask Yann whether CLEAN_FILES section is still needed.
            #with open(os.path.join(dirpath, "abinit.dep"), "wt") as fh:
            #    fh.write(s)

            # Find dependencies outside this directory (abinit.dir).
            outside_dir = []
            for fort_file in fort_files:
                outside_dir.extend(m.dirname for m in fort_file.used_mods)

	    # Write abinit.dir
            s = template.format(kind="outside the directory", directory=os.path.basename(dirpath))
            s += "include_dirs = \\\n" + pformat(sorted(set(outside_dir)))
            #print(s, end=2 * "\n")
            #with open(os.path.join(dirpath, "abinit.dir"), "wt") as fh:
            #    fh.write(s)

    def touch_alldeps(self, verbose=0):
        """
        Touch all files that depend on the modules that have been changed.
        """
        def touch(fname, times=None):
            """Emulate Unix touch."""
            import os
            with open(fname, 'a'):
                os.utime(fname, times)

        count = 0

        # TODO: possible problem if new files have been added
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
            if any(m in fort_file.used_mods for m in changed_mods):
                print("Touching:", os.path.relpath(fort_file.path))
                touch(fort_file.path)
                count += 1

        return count

    #def iter_deps_of_fortpath(self, fort_path):
    #    this = self.fort_files[fort_path]
    #    mod_names = [mod.namd for mod in this.modules]
    #    for fort_file in self.fort_files.values():
    #        if any(m in fort_file.all_uses for m in mod_names):
    #            yield fort_file

    def validate(self, verbose=0):
        """
        Validate project. Return exit status.
        """
        retcode = 0
        for fort_file in self.fort_files.values():
            count = len(fort_file.subroutines) + len(fort_file.functions)
            if count:
                print("[%s] Found %d procedure(s) outside module!" % (fort_file.name, count))
                retcode += 1

        print("retcode", retcode)
        return retcode

    def pedit(self, name, verbose=0):
        """
        Edit all children of a public entity specified by name.
        """
        #elif os.path.isfile(options.what):
        obj = self.find_public_entity(name)
        if obj is None:
            print("Cannot find public entity `%s`" % str(name))
            return 1
        if verbose:
            print(obj)

        # Find files with procedures.
        paths = sorted(set(p.path for p in obj.parents))
        if verbose:
           print(paths)
        # TODO
        #from pymods.tools import Editor
        from fkiss.tools import Editor
        return Editor().edit_files(paths, ask_for_exit=True)

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
        dir2files = self.groupby_dirname()
        dirname = os.path.join(self.srcdir, os.path.basename(dirname))

        parent_dirs, child_dirs = [], []
        for fort_file in dir2files[dirname]:
            for use_name in fort_file.all_uses:
                if use_name in EXTERNAL_MODS: continue
                used_mod = self.all_modules[use_name]
                parent_dirs.append(os.path.basename(used_mod.dirname))

        parent_dirs = sorted(set(parent_dirs))
        print(parent_dirs)

        # https://www.graphviz.org/doc/info/
        #from graphviz import Digraph
        #fg = Digraph("Directory: %s" % dirname, engine="dot" if engine == "automatic" else engine)
        #return fg

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
            print("Cannot find public entity `%s` in project" % str(name))
            return None

        obj.parents
        obj.children
        #return obj.get_graphviz(engine=engine, graph_attr=graph_attr, node_attr=node_attr, edge_attr=edge_attr)

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