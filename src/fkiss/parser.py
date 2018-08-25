"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re

from textwrap import TextWrapper
from pprint import pformat
from collections import OrderedDict, defaultdict, deque
from .tools import lazy_property


class Datatype(object):

    def __init__(self, name, lines):
        self.name = name
        self.lines = lines
        self._analyzed = False

    def __getattr__(self, name):
        """Called when an attribute lookup has not found the attribute in the usual places"""
        if not self._analyzed:
            self.analyze()
            return super(Datatype, self).__getattr__(name)
        else:
            raise AttributeError("Canoot find attributed `%s`" % str(name))

    def analyze(self):
        self._analyzed = True


class Interface(object):

    def __init__(self, name, lines):
        self.name = name
        self.lines = lines



class Procedure(object):
    """
    Base class

    contains: List of contained `Procedure`
    local_uses: List of strings with the name of the modules used (explicit) by this procedure
    includes: List of strings with the name of the included file.
    parents: List of `Procedure` calling this one
    children: List of string with the name of the subroutines called by this procedure.
    """

    def __init__(self, name, ancestor, preamble, path=None):
        self.name = name.strip()
        self.ancestor = ancestor
        self.preamble = "\n".join(preamble)
        self.path = path
        self.basename = os.path.basename(self.path)

        self.num_f90lines, self.num_doclines = 0, 0
        #self.num_omp_statements = 0
        self.contains, self.local_uses, self.includes = [], [], []
        self.parents, self.children = [], []
        self.datatypes = []
        self.interfaces = []

        # TODO
        self.visibility = "public"
        #self.has_implicit_none = False

    @lazy_property
    def is_program(self):
        return isinstance(self, Program)

    @lazy_property
    def is_module(self):
        return isinstance(self, Module)

    @lazy_property
    def is_subroutine(self):
        return isinstance(self, Subroutine)

    @lazy_property
    def is_function(self):
        return isinstance(self, Function)

    @lazy_property
    def dirpath(self):
        """Absolute path of the directory in which the procedure is located."""
        return None if self.path is None else os.path.dirname(self.path)

    @lazy_property
    def dirname(self):
        """name of the directory in which the procedure is located."""
        return None if self.path is None else os.path.basename(os.path.dirname(self.path))

    @lazy_property
    def dirlevel(self):
        # 72_response --> 72
        if self.dirname is None:
            return -1
        else:
            return int(self.dirname.split("_")[0])

    @property
    def is_public(self):
        return self.visibility == "public"

    @property
    def is_private(self):
        return self.visibility == "private"

    #@lazy_property
    #def global_uses(self)
    #    """String with all the modules used by this procedure (locals + globals)"""

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.name)

    @lazy_property
    def public_procedures(self):
        """List of public procedures."""
        if self.is_program:
            return [self]
        elif self.is_module:
            return [self] + [p for p in self.contains if p.is_public]
        elif self.is_subroutine or self.is_function:
            return [self]  if self.is_public else []
        raise TypeError("Don't know how to find public entities of type: %s" % type(self))

    def stree(self, level=0):
        lines = [level * "\t" + repr(self)]; app = lines.append
        level += 1
        #if self.is_module:
        for p in self.contains:
            app(p.stree(level=level))

        return "\n".join(lines)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, width=90):
        """
        String representation with verbosity level `verbose`.
        Text is wrapped at `width` columns.
        """
        w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)
        lines = []; app = lines.append

        app("%s: %s\n" % (self.__class__.__name__.upper(), self.name))
        app("Directory: %s" % os.path.basename(self.dirname))

        if self.ancestor is not None:
            app("ANCESTOR:\n\t%s (%s)" % (self.ancestor.name, self.ancestor.ftype))
        if self.uses:
            app("USES:\n%s\n" % w.fill(", ".join(self.uses)))
            diff = sorted(set(self.local_uses) - set(self.uses))
            if diff:
                app("LOCAL USES:\n%s\n" % w.fill(", ".join(diff)))
        if self.includes:
            app("INCLUDES:\n%s\n" % w.fill(", ".join(self.includes)))
        if self.contains:
            app("CONTAINS:\n%s\n" % w.fill(", ".join(c.name for c in self.contains)))

        if self.datatypes:
            app("DATATYPES:\n%s\n" % w.fill(", ".join(d.name for d in self.datatypes)))
        if self.interfaces:
            app("INTERFACES:\n%s\n" % w.fill(", ".join(i.name for i in self.interfaces)))

        app("PARENTS:\n%s\n" % w.fill(", ".join(sorted(p.name for p in self.parents))))
        #if verbose:
        # Add directory of parents
        dirnames = sorted(set(os.path.basename(p.dirname) for p in self.parents))
        app("PARENT_DIRS:\n%s\n" % w.fill(", ".join(dirnames)))

        app("CHILDREN:\n%s\n" % w.fill(", ".join(sorted(c for c in self.children))))
        #if verbose:
        ## Add directory of children
        #dirnames = sorted(set(os.path.basename(p.dirname) for p in self.children))
        #app("CHILDREN_DIRS:\n%s\n" % w.fill(", ".join(dirnames)))

        if verbose:
            app("DOC_STRING:\n%s" % self.preamble)

        return "\n".join(lines)

    @lazy_property
    def uses(self):
        """
        List of strings with the modules used by this procedure.
        The list includes the modules used explicitly inside the procedure as well as
        the modules imported at the module level.
        """
        uses = self.local_uses[:]
        if self.is_module or self.is_program:
            for p in self.contains:
                uses.extend(p.local_uses)
                # TODO: should be recursive
                for c in p.contains:
                    uses.extend(c.local_uses)

        elif self.is_subroutine or self.is_function:
            ancestor = self.ancestor
            while ancestor is not None:
                uses.extend(ancestor.local_uses)
                ancestor = ancestor.ancestor

        else:
            raise TypeError("Don't know how to handle %s" % type(self))

        return sorted(set(uses))

    #def analyze(self):
    #    if self.is_module:
    #        for p in self.contains:
    #            p.visibility = self.visibility
    #    else:
    #        for p in self.contains:
    #            p.visibility = "private"


class Program(Procedure):
    """Fortran program."""
    ftype = "program"


class Function(Procedure):
    """Fortran function."""
    ftype = "function"


class Subroutine(Procedure):
    """Fortran subroutine."""
    ftype = "subroutine"


class Module(Procedure):
    """Fortran module."""
    ftype = "module"

    def to_string(self, verbose=0, width=90):
        lines = []; app = lines.append
        app(super(Module, self).to_string(verbose=verbose, width=width))
        #w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)

        return "\n".join(lines)

# This one is problematic. The firs one is stricter but it does not play well with the other code!!
# I use the old one used by abilint
#RE_FUNC_END = re.compile('^[ \t]*end[ \t]*(function\n)',re.I).match

#Detect call (not after an if or ;)
#re_call = re.compile('(^[ \t]*call[ ]*)(\w+)', re.MULTILINE+re.I)

# Taken from https://github.com/cmacmackin/ford/blob/master/ford/sourceform.py
#VAR_TYPE_STRING = "^integer|real|double\s*precision|character|complex|logical|type|class|procedure|enumerator"
#VARKIND_RE = re.compile("\((.*)\)|\*\s*(\d+|\(.*\))")
#KIND_RE = re.compile("kind\s*=\s*",re.I)
#LEN_RE = re.compile("len\s*=\s*",re.I)
#ATTRIBSPLIT_RE = re.compile(",\s*(\w.*?)::\s*(.*)\s*")
#ATTRIBSPLIT2_RE = re.compile("\s*(::)?\s*(.*)\s*")
#ASSIGN_RE = re.compile("(\w+\s*(?:\([^=]*\)))\s*=(?!>)(?:\s*([^\s]+))?")
#POINT_RE = re.compile("(\w+\s*(?:\([^=>]*\)))\s*=>(?:\s*([^\s]+))?")
#EXTENDS_RE = re.compile("extends\s*\(\s*([^()\s]+)\s*\)")
#DOUBLE_PREC_RE = re.compile("double\s+precision",re.I)
#QUOTES_RE = re.compile("\"([^\"]|\"\")*\"|'([^']|'')*'",re.I)
#PARA_CAPTURE_RE = re.compile("<p>.*?</p>",re.I|re.DOTALL)
#COMMA_RE = re.compile(",(?!\s)")
#NBSP_RE = re.compile(" (?= )|(?<= ) ")
#DIM_RE = re.compile("^\w+\s*(\(.*\))\s*$")

#ATTRIB_RE = re.compile("^(asynchronous|allocatable|bind\s*\(.*\)|data|dimension|external|intent\s*\(\s*\w+\s*\)|optional|parameter|pointer|private|protected|public|save|target|value|volatile)(?:\s+|\s*::\s*)((/|\(|\w).*?)\s*$",re.I)
#END_RE = re.compile("^end\s*(?:(module|submodule|subroutine|function|procedure|program|type|interface|enum|block\sdata|block|associate)(?:\s+(\w.*))?)?$",re.I)
#BLOCK_RE = re.compile("^(\w+\s*:)?\s*block\s*$",re.I)
#BLOCK_DATA_RE = re.compile('^block\s*data\s*(\w+)?\s*$',re.I)
#ASSOCIATE_RE = re.compile("^(\w+\s*:)?\s*associate\s*\((.+)\)\s*$",re.I)
#ENUM_RE = re.compile("^enum\s*,\s*bind\s*\(.*\)\s*$",re.I)
#MODPROC_RE = re.compile("^(module\s+)?procedure\s*(?:::|\s)\s*(\w.*)$",re.I)
#MODULE_RE = re.compile("^module(?:\s+(\w+))?$",re.I)
#SUBMODULE_RE = re.compile("^submodule\s*\(\s*(\w+)\s*(?::\s*(\w+))?\s*\)\s*(?:::|\s)\s*(\w+)$",re.I)
#PROGRAM_RE = re.compile("^program(?:\s+(\w+))?$",re.I)
#SUBROUTINE_RE = re.compile("^\s*(?:(.+?)\s+)?subroutine\s+(\w+)\s*(\([^()]*\))?(?:\s*bind\s*\(\s*(.*)\s*\))?$",re.I)
#TYPE_RE = re.compile("^type(?:\s+|\s*(,.*)?::\s*)((?!(?:is\s*\())\w+)\s*(\([^()]*\))?\s*$",re.I)
#RE_MOD_END = re.compile(r'end(\s*module\s*\w*|)\Z', re.I)
##~ ABS_INTERFACE_RE = re.compile("^abstract\s+interface(?:\s+(\S.+))?$",re.I)
#BOUNDPROC_RE = re.compile("^(generic|procedure)\s*(\([^()]*\))?\s*(.*)\s*::\s*(\w.*)",re.I)
#COMMON_RE = re.compile("^common(?:\s*/\s*(\w+)\s*/\s*|\s+)(\w+.*)",re.I)
#COMMON_SPLIT_RE = re.compile("\s*(/\s*\w+\s*/)\s*",re.I)
#FINAL_RE = re.compile("^final\s*::\s*(\w.*)",re.I)
#USE_RE = re.compile("^use(?:\s*(?:,\s*(?:non_)?intrinsic\s*)?::\s*|\s+)(\w+)\s*($|,.*)",re.I)
#ARITH_GOTO_RE = re.compile("go\s*to\s*\([0-9,\s]+\)",re.I)
#CALL_RE = re.compile("(?:^|[^a-zA-Z0-9_% ]\s*)(\w+)(?=\s*\(\s*(?:.*?)\s*\))",re.I)


class FortranKissParser(object):
    """
    Parse fortran code.
    """
    # PROGRAM [name]
    # END [PROGRAM [name]]
    RE_PROG_START = re.compile(r"^program\s*(?P<name>\w*)\Z", re.I)
    RE_PROG_END = re.compile(r"^end(\s*program\s*(?P<name>\w*)|)\Z", re.I)

    # MODULE <name>
    # END [MODULE [name]]
    RE_MOD_START= re.compile(r"^module\s+(?P<name>\w+)\Z", re.I)
    RE_MOD_END = re.compile(r"^end(\s*module\s*(?P<name>\w*)|)\Z", re.I)

    # [<prefix>] <SUBROUTINE> <name> [(<args>)] [<suffix>]
    # END [SUBROUTINE [name]]
    RE_SUB_START = re.compile(r'(?P<prefix>(recursive|pure|elemental|\s)*)subroutine\s*(?P<name>\w+)', re.I)
    RE_SUB_END = re.compile(r'^end(\s*subroutine\s*(?P<name>\w*)|)\Z', re.I)

    # [<prefix>] FUNCTION <name> ([<dummy-arg-list>]) [<suffix>]
    # END [FUNCTION [name]]
    RE_FUNC_START = re.compile(r"""
(?P<prefix>
(
recursive | pure | elemental |
logical | integer | integer(\s*\(.+\)\s*)? |
double\s+precision | real(\s*\(.+\))? |
complex(\s*\(.+\))? |
character\s*\(len=\w+\s*\) |
type\s*\(\w+\)
)
\s+
)*
\s*function\s+(?P<name>\w+)\s*""",
re.I | re.VERBOSE)
    #RE_FUNC_START = re.compile('^[ \t]*(([^!\'"\n]*?)function)',re.I)
    #RE_FUNC_START = re.compile("^(?:(.+?)\s+)?function\s+(\w+)\s*(\([^()]*\))?(?=(?:.*result\s*\(\s*(\w+)\s*\))?)(?=(?:.*bind\s*\(\s*(.*)\s*\))?).*$", re.I)
    RE_FUNC_END = re.compile(r"^end(\s*function\s*(?P<name>\w*)|)\Z", re.I)

    RE_INTERFACE_START = re.compile("^(abstract\s+)?interface(?:\s+(\S.+))?$", re.I)
    RE_INTERFACE_END = re.compile("^end\s+interface(?:\s+(\S.+))?$", re.I)
    #RE_INTERFACE_START = re.compile(r"(abstract\s+)?interface\s*(?P<name>\w*)", re.I)
    #RE_INTERFACE_END = re.compile(r"end\s+interface\s+(?P<name>\w*)", re.I)

    # [if ()] call <name> [([<dummy-arg-list>])]
    #RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)\s*(?:\(\s*(.*?)\s*\))?\Z", re.I)
    RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)", re.I)

    # TYPE [ [ , access-spec ] :: ] type-name
    # [ PRIVATE ]
    # [ SEQUENCE ]
    # component-declaration
    # [ component-declaration ]...
    # END TYPE [ type-name ]
    #RE_TYPE_START = re.compile(r'type\s*(|.*::)\s*(?P<name>\w+)', re.I)
    RE_TYPE_START = re.compile(r'^type(?:\s+|\s*(,.*)?::\s*)(?P<name>\w+)\Z', re.I)
    #RE_TYPE_END = re.compile(r'^end\s+type', re.I)
    RE_TYPE_END = re.compile(r'^end(\s*type\s*(?P<name>\w*)|)\Z', re.I)

    def __init__(self, macros=None, verbose=0):
        self.verbose = verbose
        self.macros = {} if macros is None else macros

    #@staticmethod
    #def rstrip_comment(s):
    #    for i in range(len(s), -1, -1)
    #        if s[i] == "!"

    def parse_file(self, path):
        with open(path, "rt") as fh:
            string = fh.read()

            # TODO: Include Fortran files?
            #lines = []
            #for line in string.splitlines():
            #    l =  line.strip().replace("'", "").replace('"', "")
            #    if l.startswith("#include") and (l.endswith(".finc") or l.endswith(".F90")):
            #        basename = l.split()[-1]
            #        with open(os.path.join(os.path.dirname(path), basename), "rt") as incfh:
            #            lines.extend(incfh.readlines())
            #    else:
            #        lines.append(line)
            #string = "\n".join(lines)

            return self.parse_string(string, path=path)

    def parse_string(self, string, path=None):
        if path is None: path = "UknownFile"
        self.path = path

        # Replace macros. Need e.g. to treat USE_DEFS macros in libpaw and tetralib.
        for macro, value in self.macros.items():
            string = re.sub(macro, value, string)

        # Perhaps here one should join multiple lines ending with &
        # Get list of lower-case string.
        self.lines = deque(l.strip().lower() for l in string.splitlines())
        #self.warnings = []

        self.num_doclines, self.num_f90lines, self.num_omp_statements = 0, 0, 0
        preamble, self.stack, includes, uses  = [], [], [], []
        ancestor = None

        stack = self.stack

        while self.lines:
            line = self.lines.popleft()
            if not line: continue
            #print(line)

            # Count number of comments and code line
            # Inlined comments are not counted (also because I don't like them)
            if line.startswith("!"):
                if not stack or (stack and stack[-1][1] != "open"):
                    preamble.append(line)
                if line.replace("!", ""): self.num_doclines += 1
                continue
            else:
                self.num_f90lines += 1

            if line.startswith("!$omp"):
                self.num_omp_statements += 1

            # Invokations of Fortran functions are difficult to handle
            # without inspecting locals so we only handle explicit calls to routines.
            # TODO: At this level subname is a string that will be replaced by a Procedure object afterwards
            # should also handle call obj%foo()
            m = self.RE_SUBCALL.match(line)
            if m:
                subname = m.group("name")
                assert subname
                stack[-1][0].children.append(subname)
                continue

            # Interface declaration.
            m = self.RE_INTERFACE_START.match(line)
            if m:
                self.consume_interface(m)
                continue

            # Datatype declaration.
            m = self.RE_TYPE_START.match(line)
            if m and stack[-1][0].is_module and "(" not in line:
            #if m and "(" not in line:
                assert stack[-1][1] == "open"
                self.consume_datatype(m)
                continue

            # Handle include statement (CPP or Fortran version).
            if line.startswith("#include ") or line.startswith("include "):
                what = line.split()[1].replace("'", "").replace('"', "")
                if stack:
                    stack[-1][0].includes.append(what)
                else:
                    includes.append(what)
                continue

            # Find use statements and the corresponding module
            # TODO in principle one could have `use A; use B`
            if line.startswith("use "):
                smod = line.split()[1].split(",")[0].lower()
                #if smod.startswith("interfaces_"): continue
                # Remove comment at the end of the line if present.
                i = smod.find("!")
                if i != -1: smod = smod[:i]
                if self.verbose > 1: print("Found used module:", smod)
                # TODO
                #if stack:
                stack[-1][0].local_uses.append(smod)
                uses.append(smod)
                continue

            if self.simple_match(line, "contains"):
                assert stack
                ancestor = stack[-1][0]
                if self.verbose > 1: print("Setting ancestor to:", repr(ancestor))
                continue

            # Extract name from (module | program).
            m = self.RE_PROG_START.match(line)
            if m:
                name = m.group("name")
                if self.verbose > 1: print("Entering program:", name)
                if not name:
                    raise ValueError("Cannot find program name in line `%s`" % line)
                stack.append([Program(name, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            m = self.RE_MOD_START.match(line)
            if m:
                name = m.group("name")
                if self.verbose > 1: print("Entering module:", name)
                # Ignore module procedure
                #assert name != "procedure"
                # TODO
                #assert ancestor is None
                stack.append([Module(name, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            # Subroutine declaration.
            m = self.RE_SUB_START.match(line)
            if m:
                subname = m.group("name")
                if not subname:
                    raise ValueError("Cannot find procedure name in line `%s`" % line)
                if self.verbose > 1: print("Found subroutine:", subname, "in line:\n\t", line)
                stack.append([Subroutine(subname, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            # Function declaration
            m = self.RE_FUNC_START.match(line)
            if m:
                func_name = m.group("name")
                if not func_name:
                    raise ValueError("Cannot find procedure name in line:\n\t%s" % line)
                stack.append([Function(func_name, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            #isend, ftype, name = self.maybe_end_procedure(line):
            #if isend:
            if line.startswith("end "):
                tokens = line.split()
                end_ftype = None if len(tokens) == 1 else tokens[1]
                if end_ftype not in {"program", "function", "subroutine", "module", None}:
                    #print("end ftype:", end_ftype)
                    continue

                try:
                    end_name = tokens[2]
                except IndexError:
                    # TODO
                    #print("index error in:", line)
                    end_name = None

                #if end_ftype == "module": print("got end", end_name)
                #print(line, end_name)
                #last_name = stack[-1][0].name
                #if last_name != end_name:
                #    print("WARNING: end_name:", end_name, " != last_name:", last_name)
                #    print("line", line, "end_ftype", end_ftype)

                #self.close_stack_entry(end_name, end_ftype)
                if end_name is not None:
                    # Close the last entry in the stack with name == end_name.
                    for item in reversed(stack):
                        if item[0].name == end_name:
                            item[1] = "closed"
                            break
                    else:
                        raise RuntimeError("Cannot find end_name `%s` in stack:\n%s" % (
                            end_name, pformat([s[0].name for s in stack])))
                else:
                    # Close the last entry in the stack with end_ftype.
                    if end_ftype is not None:
                        for item in reversed(stack):
                            if item[0].ftype == end_ftype:
                                item[1] = "closed"
                                break
                        else:
                            raise RuntimeError("Cannot find end_ftype `%s` in stack:\n%s" % (
                                end_ftype, pformat([s[0].ftype for s in stack])))
                    else:
                        stack[-1][1] = "closed"

                #print("Closing procedure", stack[-1][0].name)
                if ancestor is not None and ancestor.name == end_name:
                    ancestor = ancestor.ancestor

        self.includes = sorted(set(includes))
        self.uses = sorted(set(uses))

        # Extract data from stack.
        self.programs, self.modules, self.subroutines, self.functions = [], [], [], []
        while stack:
            p, status = stack.pop(0)
            if status != "closed":
                print("WARNING: unclosed", repr(p), status, "ancestor", repr(p.ancestor))

            # Sort entries here.
            p.local_uses = sorted(p.local_uses)
            p.children = sorted(set(p.children))

            if p.ancestor is not None:
                #print("Adding %s to ancestor %s" % (repr(p), repr(p.ancestor)))
                p.ancestor.contains.append(p)
                #p.ancestor.contains[p.name] = p
            else:
                if p.is_module: self.modules.append(p)
                elif p.is_subroutine: self.subroutines.append(p)
                elif p.is_function: self.functions.append(p)
                elif p.is_program: self.programs.append(p)
                else: raise ValueError("Don't know how to handle type `%s`" % type(p))

        return self

    @staticmethod
    def simple_match(s, token):
        i = s.find("!")
        if i != -1: s = s[:i]
        return s.strip() == token

    def consume_interface(self, start_match):
        buf = [start_match.string]
        #name = start_match.group("name")
        name = "None"
        if self.verbose > 1: print("begin interface", name, "in line", start_match.string)
        while self.lines:
            line = self.lines.popleft()
            buf.append(line)
            end_match = self.RE_INTERFACE_END.match(line)
            if end_match:
                if self.verbose > 1: print("end interface", line)
                #end_name = end_match.group("name")
                #if name != end_name: raise ValueError("%s != %s" % (name, end_name))
                break
        self.stack[-1][0].interfaces.append(Interface(name, buf))

    def consume_datatype(self, start_match):
        buf = [start_match.string]
        name = start_match.group("name")
        if self.verbose > 1: print("begin datatype", name, "in line", start_match.string)
        while self.lines:
            line = self.lines.popleft()
            buf.append(line)
            end_match = self.RE_TYPE_END.match(line)
            if end_match:
                end_name = end_match.group("name")
                if self.verbose > 1: print("end datatype", line)
                if name == end_name:
                    break
        else:
            raise ValueError("Cannot find `end type %s` in %s" % (name, self.path))

        self.stack[-1][0].datatypes.append(Datatype(name, buf))

    #@staticmethod
    #def maybe_end_procedure(s):
    #    if not line.startswith("end"):
    #        return False, False, False
    #    tokens = line.split()
    #    what = tokens[1]
    #    assert what in ("program", "function", "subroutine", "module"):
    #    try:
    #        fname = tokens[2]
    #    except IndexError:
    #        fname = None

    #    return isend, ftype, name
